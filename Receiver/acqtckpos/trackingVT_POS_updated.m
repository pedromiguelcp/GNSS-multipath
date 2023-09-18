function [TckResultVT, navSolutionsVT] = trackingVT_POS_updated(file,signal,track,cmn,solu, Acquired,cnslxyz,eph,sbf,TckResult_Eph, TckResultCT,navSolutionsCT)
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
%                           GPSSDR_vt v1.2
% 
% Written by B. XU and L. T. HSU

%%
Spacing = 0.7:-0.05:-0.7;  % [-0.5,0,0.5];%
datalength = track.msToProcessVT;
f0  = signal.codeFreqBasis;
fL1 = signal.Fc;
fs  = signal.Fs;

pdi = 1; % unit:ms integration time
t   = 1e-3;
sv  = Acquired.sv;
svlength    = length(Acquired.sv);
sv_clk      = zeros(1,32);
sv_clk_pos  = zeros(1,32);
eph_idx     = ones(1,svlength);

% Kalman Filter Parameter
num_state = 8;
error_state = zeros(num_state,1);
Dynamic_Model = diag([0,0,0,0,0,0,0,0]);
Dynamic_Model(1,4)  = 1;
Dynamic_Model(2,5)  = 1;
Dynamic_Model(3,6)  = 1;
Dynamic_Model(7,8)  = 0;
Transistion_Matrix  = eye(length(error_state)) + Dynamic_Model * pdi * t;

state_cov = 1e5*diag([1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e0,1e0]);

process_noise(1:3,1:3) = diag(ones(1,3)*1e0);
process_noise(4:6,4:6) = diag(ones(1,3)*1e-1);
process_noise(7,7) = 1e-1;
process_noise(8,8) = 1e-2;
mesurement_noise(1:svlength,1:svlength) = eye(svlength)*3e-1;
mesurement_noise(svlength+1:2*svlength,svlength+1:2*svlength) = eye(svlength)*1e-1;


% for experimental variance estimation
flag_corrCovEst2 = 1;
counterUptR = 0;
counter_r = 0;
thresUptR = 200/pdi;

% PVT initialization based on STL
estPos      = navSolutionsCT.usrPos(file.skiptimeVT/solu.navSolPeriod,:);  
estVel      = navSolutionsCT.usrVel(file.skiptimeVT/solu.navSolPeriod,:);  
clkBias     = navSolutionsCT.clkBias(file.skiptimeVT/solu.navSolPeriod);   
clkDrift    = navSolutionsCT.clkDrift(file.skiptimeVT/solu.navSolPeriod);
total_state = [estPos,estVel,clkBias,clkDrift]';

oldCarrError    = zeros(1,svlength);
carrError       = zeros(1,svlength);
codeError       = zeros(1,svlength);
oldCarrNco      = zeros(1,svlength);

% parameters for C/N0 estimate
snrIndex	= ones(1,svlength);
K           = 20;
flag_snr    = ones(1,svlength); % flag to calculate C/N0
index_int   = zeros(1,svlength);

% parameter for updating the iono and tropo correction.
corrUpdateSec   = 0.1;
corrUpt         = corrUpdateSec/(pdi*t);
counter_corr    = corrUpt-1 * ones(svlength,1);

% Find the start of VTL, in unit of samples  
sampleStart = zeros(1, svlength); % the first subframe 
for svindex = 1:svlength
    prn = sv(svindex);
    sampleStart(svindex) = ...
        TckResult_Eph(prn).absoluteSample(sbf.nav1(prn)+eph(prn).sfb(1)*20);
    msStartTckVT(prn) = sbf.nav1(prn)+eph(prn).sfb(1)*20 + file.skiptimeVT; %file.skiptimeVT/(solu.navSolPeriod/1);

end

% Tracking initialization based on STL
for svindex = 1:svlength
    prn = sv(svindex);    
    codetemp                = generateCAcode(Acquired.sv(svindex));
    Code(svindex,:)         = [codetemp(end) repmat(codetemp,1,pdi) codetemp(1) ];
    codeFreq(svindex)       = TckResultCT(Acquired.sv(svindex)).codeFreq(msStartTckVT(prn)); %fix(sampleStartTckVT/2/signal.Fs)
    codePhaseStep(svindex)  = codeFreq(svindex)/fs;
    remChip(svindex)        = TckResultCT(Acquired.sv(svindex)).remChip(msStartTckVT(prn));
    carrFreq(svindex)       = TckResultCT(Acquired.sv(svindex)).carrFreq(msStartTckVT(prn));
    oldCarrFreq(svindex)    = TckResultCT(Acquired.sv(svindex)).carrFreq(msStartTckVT(prn));
    remCarrPhase(svindex)   = TckResultCT(Acquired.sv(svindex)).remCarrPhase(msStartTckVT(prn));
    file_ptr(svindex)       = TckResultCT(Acquired.sv(svindex)).absoluteSample(msStartTckVT(prn));
    carrFreqBasis(svindex)  = TckResultCT(Acquired.sv(svindex)).carrFreq(msStartTckVT(prn));
    codedelay_tck(svindex)  = TckResultCT(Acquired.sv(svindex)).codedelay(msStartTckVT(prn));
    oldCarrError(svindex)   = TckResultCT(Acquired.sv(svindex)).carrError(msStartTckVT(prn));
    oldCarrNco(svindex)     = TckResultCT(Acquired.sv(svindex)).carrFreq(msStartTckVT(prn)) - carrFreqBasis(svindex);
        
    svxyzr_tck_last(svindex,:) = zeros(1,3);    
    transmitTimeVT(svindex) = eph(prn).TOW(1) + file.skiptimeVT/1000;
end

% PLL filter coefficientS
[tau1carr, tau2carr] = calcLoopCoef(track.PLLBW,track.PLLDamp,track.PLLGain);

% 
amplitude = 0;
navi_data = 0;
navi_dataL035 = 0; 
deltaPr = zeros(1,length(Acquired.sv));
prRate = zeros(1,length(Acquired.sv));  
clkBias_last = clkBias;
estPos_last = estPos;


% 
localTime = min(transmitTimeVT);

% flag for position. When the file pointers in all tracking channels exceed
% the current measurement point (in samples), we do positioning
flg_pos = zeros(1,svlength);

%
posIndex = 0;

%% Start processing
h = waitbar(0,['Vector Tracking, Length: ',num2str(datalength),' ms,', '  Please wait...']);
tic
for msIndex = 1:datalength/pdi % 
    msIndex
    waitbar(msIndex/(datalength/pdi),h)
    for svindex = 1:length(Acquired.sv)
        prn = sv(svindex);
        numSample(svindex) = ceil((signal.codelength*pdi-remChip(svindex))/(codeFreq(svindex)/fs));  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fseek(file.fid, file_ptr(svindex)*1,'bof');  % file_ptr has already considered the data type of 'int16'.
        if file.dataPrecision == 2
            [rawsignal, ~] = fread(file.fid,numSample(svindex)*file.dataType,'int16') ;
            rawsignal = rawsignal';
            sin_rawsignal = rawsignal(1:2:length(rawsignal));
            cos_rawsignal = rawsignal(2:2:length(rawsignal));
            rawsignal = (sin_rawsignal - mean(sin_rawsignal)) + 1i.*(cos_rawsignal-mean(cos_rawsignal));
        else
            [rawsignal, ~]  = fread(file.fid,numSample(svindex)*file.dataType,'int8');
            rawsignal = rawsignal';
            if file.dataType == 2 
            rawsignal = rawsignal(1:2:length(rawsignal)) + 1i*rawsignal(2:2:length(rawsignal));% For NSL STEREO LBand only
            end
        end
                
        % transmitTime of the last sample of this block
        transmitTimeVT(svindex) =  transmitTimeVT(svindex) + numSample(svindex)/signal.Fs;   
        tot_est_tck(svindex) = transmitTimeVT(svindex);
              
        % SV position and Velocity at the transmitTime
        [svxyz_tck(svindex,:), sv_vel(svindex,:), sv_clk(prn), sv_clk_vel(prn), grpdel(prn)] = ...
            svPosVel(prn,eph, tot_est_tck(svindex),eph_idx(svindex));
        
        %% Iono, trop correction, follow the book of Paul: pp.268 and eq(7.34) (svxyz - estPos).
        counter_corr(svindex) = counter_corr(svindex) + 1;
        if counter_corr(svindex) ==  corrUpt
            svenu           = xyz2enu(svxyz_tck(svindex,:), estPos);
            el_rad(svindex) = atan(svenu(3)/norm(svenu(1:2)));
            %             az_rad(svindex) = (pi/2)-atan2(svenu(1),svenu(2));
            az_rad(svindex) = atan2(svenu(1),svenu(2));
            az(svindex)     = az_rad(svindex) * 180/pi;
            el(svindex)     = el_rad(svindex) * 180/pi;
            
            temp = xyz2llh(estPos);
            user_ll	= [temp(1:2).*180/pi temp(3)];
            ionodel(svindex)        = ionocorr(tot_est_tck(svindex), svxyz_tck(svindex,:), cnslxyz); %estPos(1:3)
            tropodel_unb3(svindex)  = abs(trop_UNB3(cmn.doy,user_ll(1),user_ll(3),el(svindex))); %
            
            counter_corr(svindex)   = 0;
        end
        
        
        %% Predict code freq. based on the navigation solution
        r = sqrt(sum((svxyz_tck(svindex,:) - estPos).^2));
        % Apply corrections
        predictedPr_tck(svindex) = r + clkBias + sv_clk(prn) - grpdel(prn)*cmn.cSpeed - tropodel_unb3(svindex) - ionodel(svindex);%
        % earth rotation correction
        svxyzr_tck(svindex,:) = erotcorr(svxyz_tck(svindex,:),predictedPr_tck(svindex));
        % Corrected range
        r = sqrt(sum((svxyzr_tck(svindex,:) - estPos).^2));
        predictedPr_tck(svindex) = r + clkBias + sv_clk(prn) - grpdel(prn)*cmn.cSpeed - tropodel_unb3(svindex) - ionodel(svindex) ;%
        
        % Predicted code freq.
        if msIndex == 1
            codeFreq(svindex) = TckResultCT(Acquired.sv(svindex)).codeFreq(msStartTckVT(prn));%fix(sampleStartTckVT/2/signal.Fs)
        else
            deltaPr(svindex) = (predictedPr_tck(svindex) - predictedPr_last(svindex))/(pdi*t);
            codeFreq(svindex) = f0*(1 - deltaPr(svindex)/cmn.cSpeed); % 
        end 
        predictedPr_last(svindex) = predictedPr_tck(svindex); 
           
        % new code phase step
        codePhaseStep(svindex) = codeFreq(svindex)/fs;
        
        %%  
        t_CodeEarly       = (0 + Spacing(5) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(5) + remChip(svindex));
        t_CodePrompt      = (0 + Spacing(15) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(15) + remChip(svindex));
        t_CodeLate        = (0 + Spacing(25) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(25) + remChip(svindex));
        
        indx = 1;
        CodeEarly      = Code(svindex,(ceil(t_CodeEarly) + indx));
        CodePrompt     = Code(svindex,(ceil(t_CodePrompt) + indx));
        CodeLate       = Code(svindex,(ceil(t_CodeLate) + indx));
        
        
        %% mixed with local carrier replica
        CarrTime = (0: numSample(svindex)) ./ signal.Fs;
        Wave = (2*pi*((carrFreq(svindex)) .* CarrTime)) + remCarrPhase(svindex);
        %         remCarrPhase(svindex) = rem(Wave(numSample(svindex)+1),2*pi);
        
        carrsig = exp(1i.* Wave(1:numSample(svindex)));
        InphaseSignal    = imag(rawsignal .* carrsig); %
        QuadratureSignal = real(rawsignal .* carrsig); 
        
        %%
        E_i  = sum(CodeEarly    .*InphaseSignal);
        E_q = sum(CodeEarly    .*QuadratureSignal);
        P_i  = sum(CodePrompt   .*InphaseSignal);
        P_q = sum(CodePrompt   .*QuadratureSignal);
        L_i  = sum(CodeLate     .*InphaseSignal);
        L_q = sum(CodeLate     .*QuadratureSignal);
        
        
        % remaining code and carrier phase
        remChip(svindex) = (t_CodePrompt(numSample(svindex)) + codePhaseStep(svindex)) - 1023*pdi;  
        remCarrPhase(svindex) = rem(Wave(numSample(svindex)+1),2*pi);
        
        % some parameters at last epcoh
        clkBias_last = clkBias;
        estPos_last = estPos;
        svxyzr_tck_last(svindex,:) = svxyzr_tck(svindex,:); 
        
        %% CN0
        if (flag_snr(svindex) == 1)
            index_int(svindex) = index_int(svindex) + 1;
            Zk(svindex,index_int(svindex)) = P_i^2 + P_q^2;
            if mod(index_int(svindex),K) == 0
                meanZk  = mean(Zk(svindex,:));
                varZk   = var(Zk(svindex,:));
                NA2     = sqrt(meanZk^2-varZk);
                varIQ   = 0.5 * (meanZk - NA2);
                CN0_VT(snrIndex(svindex),svindex) =  abs(10*log10(1/(1*t*pdi) * NA2/(2*varIQ)));
                index_int(svindex)  = 0;
                snrIndex(svindex)   = snrIndex(svindex) + 1;
            end
        end
        
        %% PLL discriminator
        carrError(svindex)      = atan(P_q/P_i)/(2.0 * pi);
        carrNco(svindex)        = oldCarrNco(svindex) + (tau2carr/tau1carr) * (carrError(svindex) - oldCarrError(svindex)) + carrError(svindex)*(pdi*1e-3/tau1carr) ;
        oldCarrNco(svindex)     = carrNco(svindex);
        oldCarrError(svindex)   = carrError(svindex);
        carrFreq(svindex)       = carrFreqBasis(svindex) + carrNco(svindex);%
        oldCarrFreq(svindex)    = carrFreq(svindex);
        
        %% DLL discriminator
        E	= sqrt(E_i^2+E_q^2);
        L	= sqrt(L_i^2+L_q^2);
        codeError(svindex)	= -0.5*(E-L)/(E+L); 
        
        % pseudorange measurements; pseudorange error correction
        Z(msIndex,svindex) = codeError(svindex)*cmn.cSpeed/codeFreq(svindex);
        
        %% Result Record
        TckResultVT(Acquired.sv(svindex)).E_i(msIndex)            = E_i;
        TckResultVT(Acquired.sv(svindex)).E_q(msIndex)            = E_q;
        TckResultVT(prn).P_i(msIndex)                             = P_i;
        TckResultVT(prn).P_q(msIndex)                             = P_q;
        TckResultVT(Acquired.sv(svindex)).L_i(msIndex)            = L_i;
        TckResultVT(Acquired.sv(svindex)).L_q(msIndex)            = L_q;        
        
        TckResultVT(Acquired.sv(svindex)).amplitude(msIndex)            = amplitude;
        TckResultVT(Acquired.sv(svindex)).navi_data(msIndex)            = sum(navi_data);
        TckResultVT(Acquired.sv(svindex)).navi_dataL035(msIndex)            = sum(navi_dataL035);
        
        
        %%
        TckResultVT(prn).carrError(msIndex)          = carrError(svindex);
        TckResultVT(prn).codeError(msIndex)          = codeError(svindex);
        TckResultVT(prn).remChip(msIndex)            = remChip(svindex);
        TckResultVT(prn).remCarrPhase(msIndex)       = remCarrPhase(svindex);
        TckResultVT(prn).codeFreq(msIndex)           = codeFreq(svindex);
        TckResultVT(prn).carrFreq(msIndex)           = carrFreq(svindex);
        TckResultVT(prn).carrNco(msIndex)           = carrNco(svindex);
        TckResultVT(prn).absoluteSample(msIndex)     = ftell(file.fid);
        file_ptr(svindex)                            = TckResultVT(prn).absoluteSample(msIndex);
        TckResultVT(prn).sv_vel(msIndex,:)           = sv_vel(svindex,:);
        TckResultVT(prn).codedelay(msIndex)          = mod(TckResultVT(prn).absoluteSample(msIndex)/(file.dataPrecision*file.dataType),fs*t);
        
        codedelay_tck(svindex)                       = TckResultVT(prn).codedelay(msIndex);
        
        TckResultVT(prn).deltaPr(msIndex)           = deltaPr(svindex);
        TckResultVT(prn).prRate(msIndex)           = prRate(svindex);        
    end % end for svindex in Tracking
     
    
    %% Navigation solution solving  
    numSample_min = min(numSample)-1;
    
    
    for svindex = 1:svlength
        prn = sv(svindex);
         
        tot_est_pos(svindex) = tot_est_tck(svindex) - (numSample(svindex)-numSample_min)/signal.Fs;
           
        localTime = min(tot_est_pos);
    
        [svxyz_pos(svindex,:), sv_vel_pos(svindex,:), sv_clk_pos(prn), sv_clk_vel(prn), grpdel(prn)] = ...
            svPosVel(prn,eph, tot_est_pos(svindex), eph_idx(svindex));
        
        r = sqrt(sum((svxyz_pos(svindex,:) - estPos).^2));
        predictedPr_pos(svindex) = r  + clkBias + sv_clk_pos(prn) - grpdel(prn)*cmn.cSpeed ...
            - tropodel_unb3(svindex) - ionodel(svindex);    %
        svxyzr_pos(svindex,:) = erotcorr(svxyz_pos(svindex,:),(predictedPr_pos(svindex)));%
        r = sqrt(sum((svxyzr_pos(svindex,:) - estPos).^2));
        a_pos(svindex,:) = (svxyzr_pos(svindex,:)-estPos) / r;
        H_pos(svindex,:) = [-a_pos(svindex,:) 0 0 0 1 0];
        H_pos(svindex+svlength,:) = [0 0 0 -a_pos(svindex,:) 0 1];
        
        % Find pseudorange rate error 
        prr_measured(svindex)	= -(carrFreq(svindex) - signal.IF)*cmn.cSpeed/fL1; %% 22 Jun 2021, for Labsat3w
        prr_predicted(svindex)	= (estVel - sv_vel_pos(svindex,:))*a_pos(svindex,:)'; %
        Z(msIndex,svlength+svindex) = prr_predicted(svindex) - prr_measured(svindex) - clkDrift + sv_clk_vel(prn);
    end
    
    newZ = Z(msIndex,:); % complete measurements
    
    % Kalman filter
    error_state = zeros(num_state,1);
    error_state = Transistion_Matrix * error_state;
    
    state_cov = Transistion_Matrix * state_cov * transpose(Transistion_Matrix) + process_noise;
    kalman_gain = state_cov * transpose(H_pos) * inv(H_pos * state_cov * transpose(H_pos) + mesurement_noise);
    
    counterUptR = counterUptR + 1;  % counter for update measurement noise variance, R
    recordR(counterUptR,:) = ((newZ' - H_pos * error_state));
    
    error_state = error_state + kalman_gain * (newZ' - H_pos * error_state);
    state_cov = (eye(num_state) - kalman_gain * H_pos) * state_cov;
     
    
    total_state =  total_state + error_state;
    estPos = total_state(1:3)';
    estVel = total_state(4:6)';  %
    clkBias = total_state(7);
    clkDrift = total_state(8);
    
    %% record results
    llh     = xyz2llh(cnslxyz);  %
    L_b     = llh(1);
    lamda_b = llh(2);
    C_e_n = [ -sin(lamda_b)           cos(lamda_b)         	 0;...
        -sin(L_b)*cos(lamda_b) -sin(L_b)*sin(lamda_b)	 cos(L_b);...
        -cos(L_b)*cos(lamda_b) -cos(L_b)*sin(lamda_b)	-sin(L_b);];
    usrenuvel(msIndex,:) = C_e_n * estVel';
    
    usrenu(msIndex,:) = xyz2enu(estPos,cnslxyz);
    usrllh(msIndex,:) = xyz2llh(estPos);
    usrllh(msIndex,1:2)	= usrllh(msIndex,1:2)*180/pi;
    navSolutionsVT.localTime(msIndex,:)        = localTime;
    navSolutionsVT.usrPos(msIndex,:)         = estPos;
    navSolutionsVT.usrVel(msIndex,:)         = estVel;
    navSolutionsVT.usrPosENU(msIndex,:)      = usrenu(msIndex,:);
    navSolutionsVT.usrVelENU(msIndex,:)      = usrenuvel(msIndex,:);
    navSolutionsVT.usrPosLLH(msIndex,:)   	 = usrllh(msIndex,:);
    navSolutionsVT.clkDrift(msIndex,:)  	 = clkDrift;
    navSolutionsVT.clkBias(msIndex,:)        = clkBias;
    navSolutionsVT.satePos(msIndex,:)   	 = svxyzr_pos(svindex,:);
    navSolutionsVT.sateVel(msIndex,:)   	 = sv_vel_pos(svindex,:);
    navSolutionsVT.state(msIndex,:)          = error_state;
    navSolutionsVT.svxyz_pos(:,:,msIndex)    = svxyz_pos;
    navSolutionsVT.kalman_gain(:,:,msIndex)	 = kalman_gain;
    navSolutionsVT.state_cov(msIndex,:)      = diag(state_cov); 
    navSolutionsVT.meas_inno(msIndex,:)      = ((newZ' - H_pos * error_state));
    navSolutionsVT.newZ(msIndex,:)           = newZ;
    navSolutionsVT.predicted_z(msIndex,:)    = H_pos * error_state;
    navSolutionsVT.satEA(msIndex,:)      = el;
    navSolutionsVT.satAZ(msIndex,:)      = az;  
    
    %% predict postion and clkBias at next epoch
    total_state  = Transistion_Matrix * total_state;
    estPos = total_state(1:3)';
    clkBias = total_state(7)'; 
    
    %% update Q and R by measurement variance
    if flag_corrCovEst2 == 1 && counterUptR == thresUptR
        tmpR = diag(1/counterUptR*(sum(recordR.^2)));
        mesurement_noise(1:svlength,1:svlength) = tmpR(1:svlength,1:svlength)/10; 
        mesurement_noise(svlength+1:2*svlength,svlength+1:2*svlength) = ...
            tmpR(svlength+1:2*svlength,svlength+1:2*svlength)*1;
        
        for idx = 1 : svlength
            if mesurement_noise(idx,idx) >= 12000
                mesurement_noise(idx,idx) = 12000;
            elseif mesurement_noise(idx,idx) <= 0.01
                mesurement_noise(idx,idx) = 0.01;
            end
            if mesurement_noise(idx+svlength,idx+svlength) >= 400
                mesurement_noise(idx+svlength,idx+svlength) = 400;
            elseif mesurement_noise(idx+svlength,idx+svlength) <= 0.01
                mesurement_noise(idx+svlength,idx+svlength) = 0.01;
            end
            
        end
        counterUptR = 0;
        counter_r = counter_r + 1;
        navSolutionsVT.R(counter_r,:) = diag(mesurement_noise);
    end
         
    %% 
    if mod(msIndex, 1) == 0
        disp([msIndex, localTime, usrenu(msIndex,1),usrenu(msIndex,2),usrenu(msIndex,3), clkBias, clkDrift]); 
    end
     
end % end for msIndex
close(h);
toc

%% Save results
save(['navSolVT_',file.fileName,'_updated'], 'navSolutionsVT','eph');
save(['tckRstVT_',file.fileName,'_updated'], 'TckResultVT','CN0_VT');


