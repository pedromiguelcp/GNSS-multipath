function [TckResultCT, navSolutionsCT] = trackingCT_POS(file,signal,track,cmn, Acquired,cnslxyz,eph,sbf)
%Purpose:
%   Conventional tracking and positioning using EKF and WLS
%Inputs: 
%	file        - parameters related to the data file to be processed,a structure 
%	signal   	- parameters related to signals,a structure
%	track     	- parameters related to signal tracking,a structure
%	cmn         - parameters commmonly used,a structure
%	Acquired 	- acquisition results
%	cnslxyz 	- initial position in ECEF coordinate
%	eph         - ephemeris
%	sbf         - parameters used for pseudorange estimation
%
%Outputs:
%	TckResultCT         - conventional tracking results
%	navSolutionsCT   	- navigation solutions in conventional tracking mode
%--------------------------------------------------------------------------
%                           SoftXXXGPS v1.0
% 
% Copyright (C) X X
% Written by X X

%% 
global ALPHA BETA

% 2014/07/03 JAPAN
ALPHA = [0.1490e-07   0.2235e-07  -0.1192e-06  -0.1192e-06];
BETA  = [0.1167e+06   0.1802e+06  -0.1311e+06  -0.4588e+06];

sv_clk              = zeros(1,32);
clkBias_kf         	= 0;
usr_clk_wls         = 0;
clkDrift             = 0;
estusr              = [0;0;0];
estusr_wls          = [0;0;0];
estusr_kf           = cnslxyz;%[0;0;0];
estVel              = [0;0;0];
num_state           = 8; 

Spacing = [-track.CorrelatorSpacing 0 track.CorrelatorSpacing];

sv      = Acquired.sv;
f0      = signal.codeFreqBasis;
fL1     = signal.Fc;
fs      = signal.Fs ;
pdi     = 1; % unit: ms
t       = signal.ms;
svlength    = length(sv);
datalength  = 5000;    %track.msToProcessCT;

num_state   = 8;
state = zeros(num_state,1);
Dynamic_Model       = diag([0,0,0,0,0,0,0,0]);
Dynamic_Model(1,4)  = 1;
Dynamic_Model(2,5)  = 1;
Dynamic_Model(3,6)  = 1;
Dynamic_Model(7,8)  = 1;
Transistion_Matrix  = eye(length(state)) + Dynamic_Model * pdi * t;

state_cov = diag([1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e0,1e0]);
% Sv = 1e2; % m^2/s^3,
% process_noise(1:3,1:3) = diag(ones(1,3)*(1/3 * Sv * (pdi*t)^3));
% process_noise(1:3,4:6) = diag(ones(1,3)*(1/2 * Sv * (pdi*t)^2));
% process_noise(4:6,1:3) = diag(ones(1,3)*(1/2 * Sv * (pdi*t)^2));
% process_noise(4:6,4:6) = diag(ones(1,3)*(1/1 * Sv * (pdi*t)^1));
% h0      = 2e-20;
% h2      = 2e-21;
% S_phi   = cmn.cSpeed^2 * h0/2;
% S_f     = cmn.cSpeed^2 * 2 * pi^2 * h2;
% process_noise(7:8,7:8) = [S_phi*pdi*t+S_f*(pdi*t)^3/3    S_f*(pdi*t)^2/2;...
%                                 S_f*(pdi*t)^2/2                 S_f*(pdi*t)];
process_noise = 1e1*diag([0.1 0.1 0.1 0.01 0.01 0.01 10 1]);
mesurement_noise(1:svlength,1:svlength) = eye(svlength)*1e-1; 
mesurement_noise(svlength+1:2*svlength,svlength+1:2*svlength) = eye(svlength)*1e-2; 

% for experimental variance estimation
counterUptR         = 0;

% Calculate filter coefficient values
[tau1code, tau2code] = calcLoopCoef(track.DLLBW,track.DLLDamp,track.DLLGain);

for svindex = 1:length(sv)
    prn                     = sv(svindex);
    codetemp                = generateCAcode(prn);
    Code(svindex,:)         = [codetemp(end) codetemp codetemp(1)];
    AcqCodeDelay(svindex)   = Acquired.codedelay(svindex);
%     file_ptr(svindex)       = Acquired.codedelay(svindex);
    file_ptr(svindex)       = signal.Sample-AcqCodeDelay(svindex) -1;        %%%%%
    carrFreq(svindex)       = Acquired.fineFreq(svindex);
    AcqFreq(svindex)        = Acquired.fineFreq(svindex);
    
    oldcodedelay_pos(svindex) = 0;
    oldabsoluteSample_pos(svindex) = 0;
end

% CNR estimation
WBP             = zeros(1,svlength);
NBP_I           = zeros(1,svlength);
NBP_Q           = zeros(1,svlength);
flag_snr        = zeros(1,svlength);
snrIndex        = zeros(1,svlength);
M               = 20;
int_rem         = zeros(1,svlength);
flag_int_rem    = zeros(1,svlength);
flag_ifint      = zeros(1,svlength);
index_int       = zeros(1,svlength);

eph_idx         = ones(1,svlength);

corrUpdateSec   = .001;
corrUpt         = corrUpdateSec / (pdi * t);
counter_corr    = corrUpt-1 * ones(svlength,1); 
        
% Tracking parameters
carrier_output      = zeros(1,svlength);
carrier_outputLast  = zeros(1,svlength);
PLLdiscriLast       = zeros(1,svlength);
code_output         = zeros(1,svlength);
code_outputLast     = zeros(1,svlength);
DLLdiscriLast       = zeros(1,svlength);
remChip             = zeros(1,svlength);
codePhase           = zeros(1,svlength);
codeFreq            = ones(1,svlength)*f0;
delay               = zeros(1,svlength);
remCarrPhase        = zeros(1,svlength);
PLLdiscri           = zeros(1,svlength);
DLLdiscri           = zeros(1,svlength);
delayValue          = zeros(svlength,datalength);
codeDelay           = zeros(svlength,datalength);

for msIndex = 1: datalength
    if msIndex < 100
        track.PLLBW = 50;
    else
        track.PLLBW = 10;
    end
    [tau1carr, tau2carr] = calcLoopCoef(track.PLLBW,track.PLLDamp,track.PLLGain);
    
    for svindex = 1 :svlength
        prn = sv(svindex);
        remSample = ((signal.codelength-remChip)/(codeFreq/signal.Fs)) ;
        numSample = ceil((signal.codelength-remChip(svindex))/(codeFreq(svindex)/signal.Fs));  
               
        delayValue(svindex,msIndex) = numSample - signal.Sample;
%         delayValue(svindex,msIndex) = (signal.codelength-remChip(svindex))/(codeFreq(svindex)/signal.Fs) - signal.Sample;
        
        if (flag_snr(svindex) == 1)
            index_int(svindex) = index_int(svindex) + 1;
        end
        if (flag_int_rem(svindex) == 1),
            int_rem(svindex) = int_rem(svindex) + 1;
        end
        
        fseek(file.fid, (file_ptr(svindex)+file.skip*signal.Sample)*file.dataPrecision*file.dataType,'bof');  
        if file.dataPrecision == 2
            rawsignal = fread(file.fid,numSample*file.dataType,'int16')';
            sin_rawsignal = rawsignal(1:2:length(rawsignal));
            cos_rawsignal = rawsignal(2:2:length(rawsignal));
            rawsignal = (sin_rawsignal - mean(sin_rawsignal)) + 1i.*(cos_rawsignal-mean(cos_rawsignal));
        else
            rawsignal = fread(file.fid,numSample*file.dataType,'int8')'; %
        end

        file_ptr(svindex)   = file_ptr(svindex) + numSample; %%%%%%
        
        t_CodeEarly 	= (0 + Spacing(1) + remChip(svindex)) : codeFreq(svindex)/signal.Fs : ((numSample -1) * codeFreq(svindex)/signal.Fs + Spacing(1) + remChip(svindex));
        t_CodePrompt 	= (0 + Spacing(2) + remChip(svindex)) : codeFreq(svindex)/signal.Fs : ((numSample -1) * codeFreq(svindex)/signal.Fs + Spacing(2) + remChip(svindex));
        t_CodeLate   	= (0 + Spacing(3) + remChip(svindex)) : codeFreq(svindex)/signal.Fs : ((numSample -1) * codeFreq(svindex)/signal.Fs + Spacing(3) + remChip(svindex));
        CodeEarly 	= Code(svindex,(ceil(t_CodeEarly) + 1));
        CodePrompt	= Code(svindex,(ceil(t_CodePrompt) + 1));
        CodeLate	= Code(svindex,(ceil(t_CodeLate) + 1));
        remChip(svindex) = t_CodePrompt(numSample) + codeFreq(svindex)/signal.Fs - signal.codeFreqBasis*signal.ms;
        
        CarrTime = (0:numSample) ./ signal.Fs;
        Wave = (2*pi*((carrFreq(svindex)) .* CarrTime)) + remCarrPhase(svindex);
        remCarrPhase(svindex) = rem(Wave(numSample+1),2*pi);
        carrsig = exp(1i.* Wave(1:numSample));
        InphaseSignal    = imag(rawsignal .* carrsig);
        QuadratureSignal = real(rawsignal .* carrsig);
        
        E_i = sum(CodeEarly    .*InphaseSignal);  
        E_q = sum(CodeEarly    .*QuadratureSignal);
        P_i	= sum(CodePrompt   .*InphaseSignal);  
        P_q = sum(CodePrompt   .*QuadratureSignal);
        L_i	= sum(CodeLate     .*InphaseSignal);  
        L_q = sum(CodeLate     .*QuadratureSignal);
        
        if (flag_snr(svindex) == 1)
            WBP(svindex) = WBP(svindex) + (P_i^2+P_q^2);
            NBP_I(svindex) = NBP_I(svindex) + P_i;
            NBP_Q(svindex) = NBP_Q(svindex) + P_q;
            if mod(index_int(svindex),20) == 0
                snrIndex(svindex) = snrIndex(svindex) + 1;
                NBP(svindex) = (NBP_I(svindex)^2 + NBP_Q(svindex)^2);
                NP(svindex) = NBP(svindex) / WBP(svindex);
                CN0_CT(snrIndex(svindex),svindex) =  norm(10*log10(1/t * (NP(svindex)-1)/(M-NP(svindex))));
                WBP(svindex) = 0;NBP_I(svindex)=0;NBP_Q(svindex)=0; NP(svindex)=0; NBP(svindex)=0;
            end
        end
        
        % Implement code loop filter and generate NCO command
        E = sqrt(E_i^2+E_q^2);
        L = sqrt(L_i^2+L_q^2);
        DLLdiscri(svindex) = 0.5*(E-L)/(E+L);  % DLL discriminator
        code_output(svindex) = code_outputLast(svindex) + (tau2code/tau1code)*(DLLdiscri(svindex) - DLLdiscriLast(svindex)) + DLLdiscri(svindex)* (0.001/tau1code);
        DLLdiscriLast(svindex) = DLLdiscri(svindex);
        code_outputLast(svindex) = code_output(svindex);
        codeFreq(svindex) = signal.codeFreqBasis - code_output(svindex);
        
        % PLL discriminator
        PLLdiscri(svindex) = atan(P_q/P_i)/(2*pi);
        carrier_output(svindex) = carrier_outputLast(svindex) + (tau2carr/tau1carr)*(PLLdiscri(svindex) - PLLdiscriLast(svindex)) + PLLdiscri(svindex) * (0.001/tau1carr);
        carrier_outputLast(svindex) = carrier_output(svindex); % for Laplace Transform
        PLLdiscriLast(svindex) = PLLdiscri(svindex);
        carrFreq(svindex)  = AcqFreq(svindex) + carrier_output(svindex);  % Modify carrier freq based on NCO command
        
        % Data Record
        TckResultCT(prn).P_i(msIndex)             = P_i;
        TckResultCT(prn).P_q(msIndex)             = P_q;
        TckResultCT(prn).E_i(msIndex)             = E_i;
        TckResultCT(prn).E_q(msIndex)             = E_q;
        TckResultCT(prn).L_i(msIndex)             = L_i;
        TckResultCT(prn).L_q(msIndex)             = L_q;
        TckResultCT(prn).PLLdiscri(msIndex)       = PLLdiscri(svindex);
        TckResultCT(prn).DLLdiscri(msIndex)       = DLLdiscri(svindex);
        TckResultCT(prn).codeFreq(msIndex)        = codeFreq(svindex);
        TckResultCT(prn).carrFreq(msIndex)        = carrFreq(svindex);
        TckResultCT(prn).numSample(msIndex)        = numSample;
        TckResultCT(prn).remChip(msIndex)         = remChip(svindex);        
        TckResultCT(prn).remCarrPhase(msIndex)    = remCarrPhase(svindex);
        TckResultCT(prn).absoluteSample(msIndex)  = (ftell(file.fid));%- codePhase(svindex)/codePhaseStep(svindex);
        TckResultCT(prn).absoluteSampleCodedelay(msIndex)  = mod(TckResultCT(prn).absoluteSample(msIndex)/4,fs*t);
        TckResultCT(prn).codedelay(msIndex)       = mod(TckResultCT(prn).absoluteSample(msIndex)/4,fs*t);%signal.Sample - Acquired.codedelay(svindex)+ sum(delayValue(svindex,(1:msIndex)));   % mod(TckResultCT(prn).absoluteSample(msIndex)/4,fs*t);%
        TckResultCT(prn).delayValue(msIndex)      = delayValue(svindex,msIndex);
        
        if (P_i > 0)
            TckResultCT(prn).Navi(msIndex) = 1;
        else
            TckResultCT(prn).Navi(msIndex) = -1;
        end
        if ( msIndex > 20 && flag_snr(svindex) == 0 && TckResultCT(prn).Navi(msIndex)~=TckResultCT(prn).Navi(msIndex-1))
            if (abs(sum(TckResultCT(prn).Navi(msIndex-20:msIndex-1)))/20 == 1)
                flag_int_rem(svindex) =1;
            end
        end
        if (int_rem(svindex) == (20-1) && flag_ifint(svindex)==0)
            flag_snr(svindex) =1;
            flag_ifint(svindex)=1;
            WBP(svindex) = 0;
        end
    end % end for svindex in Tracking
    
    
	%% Interpolate
    for svindex=1:svlength
    	prn = sv(svindex);
%         x1 = [file.skip*fs*t+(msIndex-1)*pdi*fs*t  TckResultCT(prn).absoluteSample(msIndex)/4];
        x1 = [oldabsoluteSample_pos(svindex)/4 TckResultCT(prn).absoluteSample(msIndex)/4];
      	y2 = [oldcodedelay_pos(svindex) TckResultCT(prn).codedelay(msIndex)];
      	x_s = file.skip*fs*t + msIndex*pdi*fs*t;
        codedelay_pos(msIndex, svindex) = interp1(x1,y2,x_s);
    end
    oldcodedelay_pos = codedelay_pos(msIndex, :);
    oldabsoluteSample_pos(svindex) = TckResultCT(prn).absoluteSample(msIndex);
    
   
    %% Pseudorange calculation and correction   
    [pseudorange relative_time] = pr_est_interp(eph, Acquired, sbf, signal, codedelay_pos(msIndex,:), pdi, msIndex);
 
%   usr_clk = usr_clk_wls ;%%%
% 	estusr = estusr_wls;%%%
    usr_clk = clkBias_kf;%%%
	estusr = estusr_kf;%%%
    
    for svindex = 1 : svlength
        prn = sv(svindex);
       	fs_new = (f0/codeFreq(svindex))*fs;
        
        %
        time = eph(prn).TOW1(eph_idx(svindex)) - sbf.sfb1(prn) * 1/50 ...
              	- sbf.nav1(prn) * 1/1000 + msIndex * 1/1000;
            
        % difference for all satellites in unit of samples
        diff_of_dat_pos(svindex) = sbf.sfb1(prn)*20*fs*t ...     
                                        + sbf.nav1(prn)*fs*t ...
                                        + (codedelay_pos(msIndex,svindex)-1);                                    
        % time of transmition                          
        tot_est_pos(svindex) = eph(prn).TOW1(eph_idx(svindex)) ...   % TOW of current subframe 1
                                    - diff_of_dat_pos(svindex)/fs ...% 
                                    + pdi*msIndex*t ...
                                    + (1/cmn.cSpeed)*(sv_clk(prn)); 
                                    %- ((fs*pdi*t)-fs_new*pdi*t)/fs;% ...

        tot_est(svindex) = tot_est_pos(svindex);
        [svxyz(svindex,:) sv_vel(svindex,:) sv_clk(prn) sv_clk_vel(prn) grpdel] = ...
                svPosVel(prn,eph,tot_est(svindex),eph_idx(svindex));
        prvec(svindex)      = pseudorange(svindex) + sv_clk(prn);% + grpdel*cmn.cSpeed;   % C/A-code pseudorange corrected for satellite clock and Tgd 
        svxyzr(svindex,:)   = erotcorr(svxyz(svindex,:),prvec(svindex) - usr_clk);   % Adjust satellite position coordinates for earth rotation correction
        svenu               = xyz2enu(svxyzr(svindex,:),estusr);
        
        % tropospheric and ionospheric delay correction	
        counter_corr(svindex) = counter_corr(svindex) + 1;
        if counter_corr(svindex) ==  corrUpt;
            svenu           = xyz2enu(svxyzr(svindex,:),estusr(1:3));
            el_rad(svindex) = atan(svenu(3)/norm(svenu(1:2)));
            az_rad(svindex) = (pi/2)-atan2(svenu(1),svenu(2));
            az(svindex)     = az_rad(svindex) * 180/pi;
            el(svindex)     = el_rad(svindex) * 180/pi;
            temp            = xyz2llh(estusr(1:3));
            user_ll         = [temp(1:2) .* 180 /pi temp(3)];
            ionodel(svindex)        = ionocorr(tot_est(svindex),svxyzr(svindex,:), estusr(1:3));
            tropodel_unb3(svindex)  = abs(trop_UNB3(184,user_ll(1),user_ll(3),el(svindex))); % doy should be corrected accordingly
%             fprintf('iono = %f %f tropo = %f %f  sv_clk = %f  grpdel = %f \n',ionodel(svindex),ionn(svindex),tropodel(svindex),tropodel_unb3(svindex),sv_clk(prn),grpdel *cmn.cSpeed);
            counter_corr(svindex)   = 0;
        end
        
        prvec(svindex) = prvec(svindex) - ionodel(svindex) - tropodel_unb3(svindex); 

    end 

    
    %% Position cal using WLS method
%     [estusr_wls, dop]       = olspos(prvec,svxyzr,estusr_wls); % ordinary least square
%     usrenu_wls(msIndex,:)   = xyz2enu(estusr_wls(1:3),cnslxyz);
%     usr_clk_wls             = estusr_wls(4);
%     TOW_wls(msIndex)        = time ;
%     usrenu_wls(msIndex,:)                 	= xyz2enu(estusr_wls(1:3),cnslxyz);
%     usrllh_wls(msIndex,:)                   = xyz2llh(estusr_wls(1:3));
%     usrllh_wls(msIndex,1:2)                 = usrllh_wls(msIndex,1:2)*180/pi;
%     navSolutionsWLS.usrPos(msIndex,:)       = (1:3);
%     navSolutionsWLS.usrPosENU(msIndex,:)    = usrenu_wls(msIndex,:);
%     navSolutionsWLS.usrPosLLH(msIndex,:)    = usrllh_wls(msIndex,:);
%     navSolutionsWLS.usr_clk(msIndex)        = usr_clk_wls;
%     navSolutionsWLS.TOW(msIndex)       = TOW_wls(msIndex);
%     navSolutionsWLS.DOP(msIndex,:)       = dop;
%     
%     fprintf('WLS: index = %4d TOW = %f E = %f N = %f U = %f B = %f\n\n',msIndex,time,usrenu_wls(msIndex,1),usrenu_wls(msIndex,2),usrenu_wls(msIndex,3),usr_clk_wls);
%       
             
    

    %% Position cal using KF
    for svindex = 1 : svlength
        r                       = sqrt(sum((svxyzr(svindex,:) - estusr_kf(1:3)').^2));
        a(svindex,:)            = (svxyzr(svindex,:)-estusr_kf(1:3)')/r;
        H(svindex,:)            = [-a(svindex,:) 0 0 0 1 0];
        H(svindex+svlength,:)   = [0 0 0 -a(svindex,:) 0 1];
        d_p(svindex,1)          = cmn.cSpeed*((TckResultCT(Acquired.sv(svindex)).carrFreq(msIndex) - signal.IF) / fL1 ) - sum( sv_vel(svindex,:) .* a(svindex,:));
        pr_delta(svindex,1)     = prvec(svindex) - norm(svxyzr(svindex,:)'-estusr_kf(1:3)) - clkBias_kf;
    end

    
    Z = [pr_delta; d_p];
    
    % system propogation
%     state       = zeros(8,1);	
    state = [zeros(1,3) estVel' 0 clkDrift]';
    state       = Transistion_Matrix * state;
    state_cov   = Transistion_Matrix * state_cov * transpose(Transistion_Matrix) + process_noise;
    
    kalman_gain = state_cov * transpose(H) * inv(H * state_cov * transpose( H) +  mesurement_noise);
    
    counterUptR             = counterUptR + 1;
    recordR(counterUptR,:)  = (Z - H * state);
    
    state       = state + kalman_gain * (Z - H * state);
    state_cov   = (eye(num_state) - kalman_gain * H) * state_cov ;
    
    % measurement update
    estusr_kf  = estusr_kf + state(1:3);
    clkBias_kf = clkBias_kf + state(7);
    
    estVel  = state(4:6);
    clkDrift = state(8);
    
    usrenu(msIndex,:)                   = xyz2enu(estusr_kf(1:3),cnslxyz);
    usrllh(msIndex,:)                   = xyz2llh(estusr_kf(1:3));
    usrllh(msIndex,1:2)                 = usrllh(msIndex,1:2)*180/pi;
    navSolutionsKF.usrPos(msIndex,:)      = estusr_kf;
    navSolutionsKF.usrVel(msIndex,:)      = estVel;
    navSolutionsKF.usrPosENU(msIndex,:)   = usrenu(msIndex,:);
    navSolutionsKF.usrPosLLH(msIndex,:)   = usrllh(msIndex,:);
    navSolutionsKF.state(msIndex,:)       = state;
    navSolutionsKF.clkDrift(msIndex)       = clkDrift;
    navSolutionsKF.clkBias(msIndex)       = clkBias_kf;
    navSolutionsKF.TOW(msIndex)           = time;
    for svindex = 1 : svlength
        prn                                 = sv(svindex);
        svInform(prn).transTime(msIndex)    = tot_est(svindex);
    end
    
    TOW_USR_CT(msIndex) = time;
    fprintf('KF: index = %4d TOW = %f E = %f N = %f U = %f  B = %f\n',msIndex, ...
        TOW_USR_CT(msIndex),usrenu(msIndex,1),usrenu(msIndex,2),usrenu(msIndex,3), clkBias_kf); %%%
    
        
end % end for msIndex

% navSolutionsCT = navSolutionsWLS;
navSolutionsCT = navSolutionsKF;

save(['navSolCT_',file.fileName,'_',num2str(20000/1000)], 'navSolutionsCT','eph','TOW_USR_CT');
save(['tckRstCT_',file.fileName,'_',num2str(20000/1000)], 'TckResultCT','CN0_CT');

