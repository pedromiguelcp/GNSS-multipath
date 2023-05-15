function [TckResultVT, navSolutionsVT] = trackingVT_POS(file,signal,track,cmn, Acquired,cnslxyz,eph,sbf,TckResultCT,navSolutionsCT)
%Purpose:
%   Vector tracking and positioning
%Inputs: 
%	file        - parameters related to the data file to be processed,a structure 
%	signal   	- parameters related to signals,a structure
%	track     	- parameters related to signal tracking,a structure
%	cmn         - parameters commmonly used,a structure
%	Acquired 	- acquisition results
%	cnslxyz 	- initial position in ECEF coordinate
%	eph         - ephemeris 
%	sbf         - parameters used for pseudorange estimation
%	TckResultCT             - conventional tracking results
%	navSolutionsCT          - navigation solutions in conventional tracking mode
%Outputs:
%	TckResultVT         - vector tracking results
%	navSolutionsVT   	- navigation solutions in vector tracking mode
%--------------------------------------------------------------------------
%                           SoftXXXGPS v1.0
% 
% Copyright (C) X X
% Written by X X

%% 
% parameter for ionospheirc model (from rinex)
% ftp://cddis.gsfc.nasa.gov/gps/data/daily/2012/brdc/
% GPSTime = 86400*7*1676+376351;
global ALPHA BETA
% ALPHA = [0.1397e-07  0.0000e+00 -0.5960e-07  0.5960e-07];
% BETA  = [0.1106e+06 -0.3277e+05 -0.2621e+06  0.1966e+06];

% 2014/07/03 JAPAN
ALPHA = [0.1490e-07   0.2235e-07  -0.1192e-06  -0.1192e-06];
BETA  = [0.1167e+06   0.1802e+06  -0.1311e+06  -0.4588e+06];

Spacing = [-track.CorrelatorSpacing 0 track.CorrelatorSpacing];
datalength = 3000; %20000 - file.skiptimeVT;    % track.msToProcessVT;
navSolutions = navSolutionsCT;
f0  = signal.codeFreqBasis;
fL1 = signal.Fc;
fs  = signal.Fs;

pdi = 1; % unit:ms integration time
t   = 1e-3;
sv          = Acquired.sv;
svlength    = length(Acquired.sv);
sv_clk      = zeros(1,32);
sv_clk_pos  = zeros(1,32);
eph_idx     = ones(1,svlength);

% Kalman Filter Parameter
num_state   = 8;
state       = zeros(num_state,1);
Dynamic_Model       = diag([0,0,0,0,0,0,0,0]);
Dynamic_Model(1,4)  = 1;
Dynamic_Model(2,5)  = 1;
Dynamic_Model(3,6)  = 1;
Dynamic_Model(7,8)  = 1;
Transistion_Matrix  = eye(length(state)) + Dynamic_Model * pdi * t;

state_cov = diag([1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e0,1e0]);

Sv = 1e2; % m^2/s^3,
process_noise(1:3,1:3) = diag(ones(1,3)*(1/3 * Sv * (pdi*t)^3));
process_noise(1:3,4:6) = diag(ones(1,3)*(1/2 * Sv * (pdi*t)^2));
process_noise(4:6,1:3) = diag(ones(1,3)*(1/2 * Sv * (pdi*t)^2));
process_noise(4:6,4:6) = diag(ones(1,3)*(1/1 * Sv * (pdi*t)^1));
h0      = 2e-20;
h2      = 2e-21;
S_phi   = cmn.cSpeed^2 * h0/2;
S_f     = cmn.cSpeed^2 * 2 * pi^2 * h2;
process_noise(7:8,7:8) = [S_phi*pdi*t+S_f*(pdi*t)^3/3    S_f*(pdi*t)^2/2;...
                                S_f*(pdi*t)^2/2                 S_f*(pdi*t)];
process_noise = 1e1*diag([0.1 0.1 0.1 0.01 0.01 0.01 10 1]);
mesurement_noise(1:svlength,1:svlength) = eye(svlength)*1e-1; 
mesurement_noise(svlength+1:2*svlength,svlength+1:2*svlength) = eye(svlength)*1e-2; 

% for experimental variance estimation
counterUptR         = 0;
thresUptR           = 1000;
counter_r           = 0;
jump                = [0 1 1 0 0 0];


file.skiptimeVT = 1000;

estPos      = navSolutions.usrPos(file.skiptimeVT,:);
estVel      = navSolutions.usrVel(file.skiptimeVT,:);
oldestVel   = estVel;
clkBias     = navSolutions.state(file.skiptimeVT,7);
clkDrift     = navSolutions.state(file.skiptimeVT,8); 
oldclkDrift  = clkDrift;

oldCarrError    = zeros(1,svlength);
carrError       = zeros(1,svlength);
codeError       = zeros(1,svlength);
oldCarrNco      = zeros(1,svlength); 

% parameters for C/N0 estimate
snrIndex	= ones(1,svlength);
K           = 100;
flag_snr    = ones(1,svlength); % flag to calculate C/N0
index_int   = zeros(1,svlength);

% parameter for updating the iono and tropo correction.
corrUpdateSec   = 100;
corrUpt         = corrUpdateSec / (pdi * t);
counter_corr    = corrUpt-1 * ones(svlength,1);

for svindex = 1 : length(Acquired.sv)
    codetemp                = generateCAcode(Acquired.sv(svindex));
    Code(svindex,:)         = [codetemp(end) repmat(codetemp,1,pdi) codetemp(1)];
    codeFreq(svindex)       = TckResultCT(Acquired.sv(svindex)).codeFreq(file.skiptimeVT);
    carrFreq(svindex)       = TckResultCT(Acquired.sv(svindex)).carrFreq(file.skiptimeVT);
    oldCarrFreq(svindex)    = TckResultCT(Acquired.sv(svindex)).carrFreq(file.skiptimeVT);
    remChip(svindex)        = TckResultCT(Acquired.sv(svindex)).remChip(file.skiptimeVT);
    codePhaseStep(svindex)  = codeFreq(svindex)/fs;
    remCarrPhase(svindex)   = TckResultCT(Acquired.sv(svindex)).remCarrPhase(file.skiptimeVT);
    file_ptr(svindex)       = TckResultCT(Acquired.sv(svindex)).absoluteSample(file.skiptimeVT );
    carrFreqBasis(svindex)  = TckResultCT(Acquired.sv(svindex)).carrFreq(file.skiptimeVT);  
    old_abs_sample(svindex) = TckResultCT(Acquired.sv(svindex)).codedelay(file.skiptimeVT);
    codedelay_tck(svindex)  = TckResultCT(Acquired.sv(svindex)).codedelay(file.skiptimeVT);
    oldCarrError(svindex)   = TckResultCT(Acquired.sv(svindex)).PLLdiscri(file.skiptimeVT);
    oldCarrNco(svindex)     = TckResultCT(Acquired.sv(svindex)).carrFreq(file.skiptimeVT) - carrFreqBasis(svindex);
end

% Calculate filter coefficient values
[tau1carr, tau2carr] = calcLoopCoef(10,0.707,0.25);

tic
for msIndex = 1: datalength/pdi    
    for svindex = 1 : length(Acquired.sv)
        prn = sv(svindex);
        
        % read data
        numSample(svindex) = ceil((signal.codelength*pdi-remChip(svindex)) /(codeFreq(svindex)/fs));  %%%
%         delayValue(svindex,msIndex) = numSample(svindex) - signal.Sample;
        fseek(file.fid, file_ptr(svindex)*1,'bof');  % file_ptr has already considered the data type of 'int16'.
        if file.dataPrecision == 2
            [rawsignal, count] = fread(file.fid,numSample(svindex)*file.dataType,'int16') ;
            rawsignal = rawsignal';
            sin_rawsignal = rawsignal(1:2:length(rawsignal));
            cos_rawsignal = rawsignal(2:2:length(rawsignal));
            rawsignal = (sin_rawsignal - mean(sin_rawsignal)) + 1i.*(cos_rawsignal-mean(cos_rawsignal));
        else
            rawsignal = fread(file.fid,numSample*file.dataType,'int8')'; %
        end       
        
        time = eph(prn).TOW1(eph_idx(svindex)) - sbf.sfb1(prn)*1/50 - sbf.nav1(prn)*t ...
                    + pdi*msIndex*t + file.skiptimeVT*pdi*t;
        
        % time to subframe 1 in unit of samples, used for track
        diff_of_dat_tck(svindex) = sbf.sfb1(prn)*20*fs*t...
                                        + sbf.nav1(prn)*fs*t...
                                        + codedelay_tck(svindex)-1; %...
                                        %+ mod(Sample(svindex),(fs*t));
        % time of transmition 
        tot_est_tck(svindex) = eph(prn).TOW1(eph_idx(svindex))...
                                    - diff_of_dat_tck(svindex)/fs...
                                    + pdi*msIndex*t ...
                                    + (1/cmn.cSpeed)*(sv_clk(prn)) ;%...
                                    %+ floor(ftell(file.fid)/(fs*t))*t...
                                    % % + pdi*(msIndex-1)*t - ((fs*pdi*t)-fs_new*pdi*t)/fs + file.skiptimeVT*pdi*t + Sample(svindex)/fs;
        
        % Current SV position and Velocity for tck
        [svxyz_tck(svindex,:) sv_vel(svindex,:) sv_clk(prn) sv_clk_vel(prn)  grpdel(prn)] = ...
                                svPosVel(prn,eph,tot_est_tck(svindex),eph_idx(svindex));
        
        %% Iono, trop correction, follow the book of Paul: pp.268 and eq(7.34) (svxyz - estPos).
        counter_corr(svindex) = counter_corr(svindex) + 1;
        if counter_corr(svindex) ==  corrUpt;
            svenu           = xyz2enu(svxyz_tck(svindex,:),estPos(1:3));
            el_rad(svindex) = atan(svenu(3)/norm(svenu(1:2)));
            az_rad(svindex) = (pi/2)-atan2(svenu(1),svenu(2));
            az(svindex)     = az_rad(svindex) * 180/pi;
            el(svindex)     = el_rad(svindex) * 180/pi;
            
            temp = xyz2llh(estPos);
            user_ll	= [temp(1:2) .* 180 /pi temp(3)];
            ionodel(svindex)        = ionocorr(tot_est_tck(svindex),svxyz_tck(svindex,:), estPos(1:3));
            tropodel_unb3(svindex)  = abs(trop_UNB3(184,user_ll(1),user_ll(3),el(svindex)));
            counter_corr(svindex)   = 0;
        end
        
        %% Predict code freq
        r = sqrt(sum((svxyz_tck(svindex,:) - estPos(1:3)).^2));
%         r = r
        predictedPr_tck(svindex) = r + clkBias + sv_clk(prn) - grpdel(prn)*cmn.cSpeed - tropodel_unb3(svindex) - ionodel(svindex);%
%         predictedPr_tck(svindex) = predictedPr_tck(svindex)
        svxyzr_tck(svindex,:) = erotcorr(svxyz_tck(svindex,:),(predictedPr_tck(svindex)));% - clkBias
        
        r = sqrt(sum((svxyzr_tck(svindex,:) - estPos(1:3)).^2));
        predictedPr_tck(svindex) = r + clkBias + sv_clk(prn) - grpdel(prn)*cmn.cSpeed - tropodel_unb3(svindex) - ionodel(svindex);% + + clkBias + grpdel(prn)*cmn.cSpeed
        
        if msIndex == 1
            codeFreq(svindex) = TckResultCT(Acquired.sv(svindex)).codeFreq(file.skiptimeVT+1);
        else
            codeFreq(svindex) = f0 * (1 - (predictedPr_tck(svindex) - predictedPr_last(svindex))/(cmn.cSpeed*pdi*t)); 
        end
        predictedPr_last(svindex) = predictedPr_tck(svindex);
        
        codePhaseStep(svindex) = codeFreq(svindex)/fs;
        t_CodeEarly      = (0 + Spacing(1) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(1) + remChip(svindex));
        t_CodePrompt     = (0 + Spacing(2) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(2) + remChip(svindex));
        t_CodeLate       = (0 + Spacing(3) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(3) + remChip(svindex));
        
        CodeEarly 	= Code(svindex,(ceil(t_CodeEarly    ) + 1));
        CodePrompt	= Code(svindex,(ceil(t_CodePrompt   ) + 1));
        CodeLate	= Code(svindex,(ceil(t_CodeLate     ) + 1));
        
        remChip(svindex) = (t_CodePrompt(numSample(svindex)) + codePhaseStep(svindex)) - 1023*pdi;
        
        CarrTime = (0: numSample(svindex)) ./ signal.Fs;
        Wave = (2*pi*((carrFreq(svindex)) .* CarrTime)) + remCarrPhase(svindex);
        remCarrPhase(svindex) = rem(Wave(numSample(svindex)+1),2*pi);
       
        carrsig = exp(1i.* Wave(1:numSample(svindex)));
        InphaseSignal    = imag(rawsignal .* carrsig);
        QuadratureSignal = real(rawsignal .* carrsig);
        
        E_i  = sum(CodeEarly    .*InphaseSignal);  E_q = sum(CodeEarly    .*QuadratureSignal);
        P_i  = sum(CodePrompt   .*InphaseSignal);  P_q = sum(CodePrompt   .*QuadratureSignal);
        L_i  = sum(CodeLate     .*InphaseSignal);  L_q = sum(CodeLate     .*QuadratureSignal);
        
        IE(svindex) = E_i; IP(svindex) = P_i; IL(svindex) = L_i;
        QE(svindex) = E_q; QP(svindex) = P_q; QL(svindex) = L_q;
        
        % Calculate CN0
        if (flag_snr(svindex) == 1)
            index_int(svindex) = index_int(svindex) + 1;
            Zk(svindex,index_int(svindex)) = P_i^2 + P_q^2;
            if mod(index_int(svindex),K) == 0
                meanZk  = mean(Zk(svindex,:));
                varZk   = var(Zk(svindex,:));
                NA2     = sqrt(meanZk^2-varZk);
                varIQ   = 0.5 * (meanZk - NA2);
                CN0_VT(snrIndex(svindex),svindex) =  abs(10*log10(1/(1*t) * NA2/(2*varIQ)));
                index_int(svindex)  = 0;
                snrIndex(svindex)   = snrIndex(svindex) + 1;
            end
        end
        
        % PLL discriminator
        carrError(svindex)      = atan(QP(svindex)/IP(svindex) )/(2.0 * pi);
        carrNco(svindex)        = oldCarrNco(svindex) + (tau2carr/tau1carr) * (carrError(svindex) - ...
                                        oldCarrError(svindex)) + carrError(svindex) * (pdi*1e-3/tau1carr);
        oldCarrNco(svindex)     = carrNco(svindex);
        oldCarrError(svindex)   = carrError(svindex);
        carrFreq(svindex)       = carrFreqBasis(svindex) + carrNco(svindex);
        oldCarrFreq(svindex)    = carrFreq(svindex);
        
        % DLL discriminator
        E                   = sqrt(IE(svindex)^2+QE(svindex)^2);
        L                   = sqrt(IL(svindex)^2+QL(svindex)^2);
        DLLdiscri           = 0.5*(E-L)/(E+L);
        codeError(svindex)  = DLLdiscri;
        
        % book of Paul: Eq(7.210)
        Z(msIndex,svindex)          = (codeError(svindex))*cmn.cSpeed / codeFreq(svindex);
        delta_Z_ch(msIndex,svindex) = Z(msIndex,svindex);
        oldCodeDiff(svindex)        = codeError(svindex);
        
        % Data Record
        TckResultVT(prn).P_i(msIndex)                = IP(svindex);
        TckResultVT(prn).P_q(msIndex)                = QP(svindex);
        TckResultVT(prn).E_i(msIndex)                = IE(svindex);
        TckResultVT(prn).E_q(msIndex)                = QE(svindex);
        TckResultVT(prn).L_i(msIndex)                = IL(svindex);
        TckResultVT(prn).L_q(msIndex)                = QL(svindex);
        TckResultVT(prn).carrError(msIndex)          = carrError(svindex);
        TckResultVT(prn).codeError(msIndex)          = codeError(svindex);
        TckResultVT(prn).codePhase(msIndex)          = remChip(svindex);
        TckResultVT(prn).carrPhase(msIndex)          = remCarrPhase(svindex);
        TckResultVT(prn).codeFreq(msIndex)           = codeFreq(svindex);
        TckResultVT(prn).carrFreq(msIndex)           = carrFreq(svindex);
        TckResultVT(prn).absoluteSample(msIndex)     = ftell(file.fid); 
        file_ptr(svindex)                           = TckResultVT(prn).absoluteSample(msIndex);
        TckResultVT(prn).sv_vel(msIndex,:)           = sv_vel(svindex,:); % file_ptr(svindex)
        
        % codedelay calculation
%         if (mod(TckResultVT(prn).absoluteSample(msIndex),fs*t) - old_abs_sample(svindex)) >= 1000
%             idx_sample(svindex) = idx_sample(svindex) - 1;
%         elseif (mod(TckResultVT(prn).absoluteSample(msIndex),fs*t) - old_abs_sample(svindex)) <= -1000
%             idx_sample(svindex) = idx_sample(svindex) + 1;
%         end
%         TckResultVT(prn).codedelay(msIndex)  = idx_sample(svindex)*fs*t + mod(TckResultVT(prn).absoluteSample(msIndex),fs*t);
%         codedelay_tck(svindex)              = TckResultVT(prn).codedelay(msIndex);
%         old_abs_sample(svindex)             = mod(TckResultVT(prn).absoluteSample(msIndex),fs*t);
        
        
%         delayValue(svindex,msIndex) = numSample(svindex) - signal.Sample;        
        TckResultVT(prn).codedelay(msIndex)       = mod(TckResultVT(prn).absoluteSample(msIndex)/4,fs*t);
%         TckResultVT(prn).codedelay(msIndex)       = codedelay_tck(svindex)+ sum(delayValue(svindex,(1:msIndex))) ;   % mod(TckResultCT(prn).absoluteSample(msIndex)/4,fs*t);%
%         TckResultVT(prn).codedelay(msIndex)
        codedelay_tck(svindex)              = TckResultVT(prn).codedelay(msIndex);
        
    end % end for svindex in Tracking
    
    %%
    if msIndex == 1
        for svindex = 1 : svlength
            prn = sv(svindex);
            x1 = [(TckResultCT(prn).absoluteSample(file.skiptimeVT-1))/4 (TckResultVT(prn).absoluteSample(msIndex))/4];
            y1 = [0 delta_Z_ch(msIndex,svindex)];
            y2 = [TckResultCT(Acquired.sv(svindex)).codedelay(file.skiptimeVT-1) TckResultVT(prn).codedelay(msIndex)];
            y3 = [TckResultCT(Acquired.sv(svindex)).carrFreq(file.skiptimeVT-1) carrFreq(svindex)];
            % x_s = [(file.skiptimeVT+msIndex)*pdi*fs*t]
            x_s = (file.skip + file.skiptimeVT-1)*fs*t + msIndex*pdi*fs*t;
            Z(msIndex,svindex)         = interp1(x1,y1,x_s);
            codedelay_pos(svindex)     = interp1(x1,y2,x_s);
            new_carrFreq(svindex)      = interp1(x1,y3,x_s);
            
            sv_clk_pos = sv_clk;            
            
            oldabsoluteSample_pos(svindex) = (TckResultVT(prn).absoluteSample(msIndex));
        end
    else
        for svindex = 1 : svlength
            prn = sv(svindex);
            % x1 = [(file.skiptimeVT+msIndex-1)*pdi*fs*t TckResultVT(prn).absoluteSample(msIndex)]
            % x1 = [(file.skip + file.skiptimeVT)*fs*t+(msIndex-1)*pdi*fs*t  TckResultVT(prn).absoluteSample(msIndex)/4]
%             x1 = [(file.skip + file.skiptimeVT-1)*fs*t+(msIndex-1)*pdi*fs*t  TckResultVT(prn).absoluteSample(msIndex)/4];
            x1 = [oldabsoluteSample_pos(svindex)/4,  TckResultVT(prn).absoluteSample(msIndex)/4]; 
            y1 = [oldcode(svindex), delta_Z_ch(msIndex,svindex)] ;
            y2 = [oldcodedelay_pos(svindex), TckResultVT(prn).codedelay(msIndex)];
            y3 = [oldnew_carrFreq(svindex), carrFreq(svindex)];
            % x_s = [(file.skiptimeVT+msIndex)*pdi*fs*t];
%             x_s = (file.skip + file.skiptimeVT-1)*fs*t + msIndex*pdi*fs*t;
            x_s = (file.skip + file.skiptimeVT)*fs*t + msIndex*pdi*fs*t;
            Z(msIndex,svindex)         = interp1(x1,y1,x_s);
            codedelay_pos(svindex)     = interp1(x1,y2,x_s);
            new_carrFreq(svindex)      = interp1(x1,y3,x_s);
        
            oldabsoluteSample_pos(svindex) = TckResultVT(prn).absoluteSample(msIndex);
        end
        
    end
    
    %%
    oldcode = Z(msIndex,1:svlength);
    oldcodedelay_pos = codedelay_pos;
    oldnew_carrFreq =  new_carrFreq;
    for svindex = 1 : svlength
        prn = sv(svindex);
        diff_of_dat_pos(svindex) = sbf.sfb1(prn)*20*fs*t ...
                                        + sbf.nav1(prn)*fs*t ...
                                        + (codedelay_pos(svindex));
        
%         fs_new = (f0/codeFreq(svindex))*fs;        
   
        tot_est_pos(svindex) = eph(prn).TOW1(eph_idx(svindex))...
                                    - diff_of_dat_pos(svindex)/fs...%
                                    + pdi*msIndex*t ...
                                    + (1/cmn.cSpeed)*sv_clk_pos(prn);
%                                     - ((fs*pdi*t)-fs_new*pdi*t)/fs  + file.skiptimeVT*1*t...
%                                     + (1/cmn.cSpeed)*(clkBias + sv_clk_pos(prn) - ionodel(svindex) - tropodel_unb3(svindex)) - grpdel(prn);
        
     
        [svxyz_pos(svindex,:) sv_vel_pos(svindex,:) sv_clk_pos(prn) sv_clk_vel(prn) grpdel(prn)] = ...
                                    svPosVel(prn,eph,tot_est_pos(svindex),eph_idx(svindex));
        
        r = sqrt(sum((svxyz_pos(svindex,:) - estPos(1:3)).^2));
        predictedPr_pos(svindex) = r + clkBias + sv_clk_pos(prn) - grpdel(prn)*cmn.cSpeed - tropodel_unb3(svindex) - ionodel(svindex);%
        svxyzr_pos(svindex,:) = erotcorr(svxyz_pos(svindex,:),predictedPr_pos(svindex));% - clkBias
        r = sqrt(sum((svxyzr_pos(svindex,:) - estPos(1:3)).^2));
        a_pos(svindex,:) = (svxyzr_pos(svindex,:)-estPos(1:3)) / r;
        H_pos(svindex,:) = [-a_pos(svindex,:) 0 0 0 1 0];
        H_pos(svindex+svlength,:) = [0 0 0 -a_pos(svindex,:) 0 1];
        
        Z(msIndex,svlength+svindex) = cmn.cSpeed * ((new_carrFreq(svindex) - signal.IF) / fL1) ...
                                            + sum(sv_vel_pos(svindex,:) .* (-a_pos(svindex,:)));
    end
    
    newZ = Z(msIndex,:);
    
    state = [zeros(3,1);estVel';0;clkDrift];
    state = Transistion_Matrix * state;
    state_cov = Transistion_Matrix * state_cov * transpose(Transistion_Matrix) + process_noise; % state covariance    
    kalman_gain = state_cov * transpose(H_pos) * inv(H_pos * state_cov * transpose(H_pos) + mesurement_noise); % Kalman gain     
    state = state + kalman_gain * (newZ' - H_pos * state); % state udate   
    state_cov = (eye(num_state) - kalman_gain * H_pos) * state_cov; % state covariance update
    
    estPos     = estPos + state(1:3)';
    clkBias    = clkBias + state(7);
    estVel     = state(4:6)';
    clkDrift   = state(8);

%     state = zeros(8,1);	
%     state = Transistion_Matrix * state;
%     state_cov = Transistion_Matrix * state_cov * transpose(Transistion_Matrix) + process_noise; % state covariance    
%     kalman_gain = state_cov * transpose(H_pos) * inv(H_pos * state_cov * transpose(H_pos) + mesurement_noise); % Kalman gain     
%     state = state + kalman_gain * (newZ' - H_pos * state); % state udate   
%     state_cov = (eye(num_state) - kalman_gain * H_pos) * state_cov; % state covariance update
%     
%     estPos     = estPos + state(1:3)';
%     clkBias    = clkBias + state(7);
%     estVel     = state(4:6)';
%     clkDrift   = state(8);
    
    
    %% record results   
    llh     = xyz2llh(estPos);
    L_b     = llh(1);
    lamda_b = llh(2);
    Cne = [ -sin(lamda_b)           cos(lamda_b)         	0;...
                -sin(L_b)*sin(lamda_b) -sin(L_b)*sin(lamda_b)	cos(L_b);...
                cos(L_b)*cos(lamda_b)  cos(L_b)*sin(lamda_b)	sin(L_b);];    
    usrenuvel(msIndex,:) =  Cne * estVel';
    usrenu(msIndex,:) = xyz2enu(estPos(1:3),cnslxyz);    
    navSolutionsVT.usrPos(msIndex,:)         = estPos;
    navSolutionsVT.usrVel(msIndex,:)         = estVel;
    navSolutionsVT.usrPosENU(msIndex,:)      = usrenu(msIndex,:);
    navSolutionsVT.usrVelENU(msIndex,:)      = usrenuvel(msIndex,:);
    navSolutionsVT.clkDrift(msIndex,:)        = clkDrift;
    navSolutionsVT.clkBias(msIndex,:)        = clkBias;
    navSolutionsVT.state(msIndex,:)          = state;
    navSolutionsVT.svxyz_pos(:,:,msIndex)    = svxyz_pos;
    navSolutionsVT.kalman_gain(:,:,msIndex)	 = kalman_gain;
    navSolutionsVT.state_cov(msIndex,:)      = diag(state_cov);
    navSolutionsVT.meas_inno(msIndex,:)      = ((newZ' - H_pos * state));
    navSolutionsVT.newZ(msIndex,:)           = newZ;
    navSolutionsVT.predicted_z(msIndex,:)    = H_pos * state;
    
    
    %% predict postion and clkBias at next epoch
    estPos   	= estPos - (estVel+oldestVel)/2 * pdi * t; % predicted position
    clkBias     = clkBias - (clkDrift+oldclkDrift)/2 * pdi * t;
    oldestVel   = estVel;
    oldclkDrift	= clkDrift;
    
    
    %% update Q and R by measurement variance
% 	counterUptR	= counterUptR + 1;  % conter for Q and R update
%     if flag_corrCovEst2 == 1 && counterUptR == thresUptR
%         %                 disp('update R');
%         
%         tmpR = diag(var(recordR));
%         mesurement_noise(1:svlength,1:svlength)=tmpR(1:svlength,1:svlength) * 50;       % original: *10
%         mesurement_noise(svlength+1:2*svlength,svlength+1:2*svlength)=tmpR(svlength+1:2*svlength,svlength+1:2*svlength) * (50*1);
%         
%         % mesurement_noise = mesurement_noise;
%         for idx = 1 : svlength
%             if mesurement_noise(idx,idx) >= 5000
%                 mesurement_noise(idx,idx) = 5000;
%             elseif mesurement_noise(idx,idx) <= 0.01;%min_thres(idx_diff_min)
%                 mesurement_noise(idx,idx) = 0.01;%min_thres(idx_diff_min);
%             end
%             if mesurement_noise(idx+svlength,idx+svlength) >= 50
%                 mesurement_noise(idx+svlength,idx+svlength) = 50;
%             elseif mesurement_noise(idx+svlength,idx+svlength) <= 0.01;%min_thres_rate(idx_diff_min_rate)
%                 mesurement_noise(idx+svlength,idx+svlength) = 0.01;%min_thres_rate(idx_diff_min_rate);
%             end
%         end
%         counterUptR = 0;
%         counter_r = counter_r + 1;
%         navSolutionsVT.R(counter_r,:) = diag(mesurement_noise);
%     end
    
    %%
    TOW_USR(msIndex) = time;
    %             fprintf('[%4d] %3d  %e TOW = %6.2f  E = %+3.2f  N = %+3.2f  U = %+3.2f  B = %+3.2f  Ve = %+2.2f Vn = %+2.2f Vu = %+2.2f\n',...
    %                 msIndex,min_thres(idx_diff_min),min_thres_rate(idx_diff_min_rate),time,usrenu(msIndex,1),usrenu(msIndex,2),usrenu(msIndex,3),usr_clk ,usrenuvel(msIndex,1),usrenuvel(msIndex,2),usrenuvel(msIndex,3));
    fprintf('[%4d] TOW = %6.3f  E = %+3.2f  N = %+3.2f  U = %+3.2f  B = %+3.2f  \n',...
        msIndex,time,usrenu(msIndex,1),usrenu(msIndex,2),usrenu(msIndex,3),clkBias);
    
    
    if mod(index_int(svindex),K) == 0
        for svindex = 1 : svlength
            fprintf('%3.2f -- ',CN0_VT(snrIndex(svindex)-1,svindex));
        end, fprintf('\n');
    end
    
end % end for msIndex
toc

TrackingPlot_vt(TckResultVT,Acquired,navSolutionsVT,0,cnslxyz);

save(['navSolVT_',file.fileName,'_',num2str(20000/1000)], 'navSolutionsVT','eph','TOW_USR_CT');
save(['tckRstVT_',file.fileName,'_',num2str(20000/1000)], 'TckResultVT','CN0_VT');


