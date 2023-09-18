function [TckResultCT_pos, navSolutionsCT] = trackingCT_POS_updated(file,signal,track,cmn, Acquired,TckResult_Eph, cnslxyz,eph,sbf,solu)
%Purpose:
%   Scalar tracking and positioning. Positioning and tracking are
%   implemented at the same time.
%   Conventional tracking and positioning using EKF and WLS
%Inputs:
%	file        - parameters related to the data file to be processed,a structure
%	signal   	- parameters related to signals,a structure
%	track     	- parameters related to signal tracking,a structure
%	cmn         - parameters commmonly used,a structure
%	Acquired 	- acquisition results
%   TckResult_Eph - Tracking results that used for decoding eph., which
%   also contains information like indexes for first nav. bit transition, subframe,
%   absolute sample index in the IF file for each ms, etc
%	cnslxyz 	- initial position in ECEF coordinate
%	eph         - ephemeris
%	sbf         - parameters used for pseudorange estimation
%
%Outputs:
%	TckResultCT         - conventional tracking results
%	navSolutionsCT   	- navigation solutions in conventional tracking mode
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.2
% 
% Written by B. XU and L. T. HSU


%%
sv_clk              = zeros(1,32);
clkBias_kf         	= 0;
usr_clk_wls         = 0;
clkDrift            = 0;
oldclkDrift         = 0;
estusr              = zeros(1,3);
estusr_wls          = cnslxyz;% 
estusr_kf           = cnslxyz;%
estVel              = zeros(1,3);
oldestVel          	= estVel;
num_state           = 8;

Spacing = 0.6:-0.05:-0.6;
Spacing = [-track.CorrelatorSpacing 0 track.CorrelatorSpacing];

sv      = Acquired.sv ;
f0      = signal.codeFreqBasis;
fs      = signal.Fs ;
pdi     = track.pdi;
t       = signal.ms;
svlength    = length(sv);
datalength  = track.msPosCT;

% Kalman Filter Parameter
num_state   = 8;

% error state vector
error_state = zeros(num_state,1);
total_state = [cnslxyz(1:3),zeros(1,5)]';%zeros(num_state,1);

% system transition matrix
Dynamic_Model = diag(zeros(1,num_state));
Dynamic_Model(1,4)  = 1;
Dynamic_Model(2,5)  = 1;
Dynamic_Model(3,6)  = 1;
Dynamic_Model(7,8)  = 1;
Transistion_Matrix  = eye(length(error_state)) + Dynamic_Model*pdi*t;

% error covariance matrix
state_cov = 1e5*diag([1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e0,1e0]);

% process (System) noise noise covariance matrix
process_noise(1:3,1:3) = diag(ones(1,3)*2e-1);
process_noise(4:6,4:6) = diag(ones(1,3)*1e-1);
process_noise(7,7) = 1e-1;
process_noise(8,8) = 1e-2;

% measurement noise covariance matrix
mesurement_noise(1:svlength,1:svlength) = eye(svlength)*3e-1;
mesurement_noise(svlength+1:2*svlength,svlength+1:2*svlength) = eye(svlength)*1e-1;

% parameters for measurement noise variance update
flag_corrCovEst2 = 1;
counterUptR = 0;
counter_r = 0;
thresUptR = 200/pdi;

% Calculate filter coefficient values
[tau1code, tau2code] = calcLoopCoef(track.DLLBW,track.DLLDamp,track.DLLGain);
[tau1carr, tau2carr] = calcLoopCoef(track.PLLBW,track.PLLDamp,track.PLLGain);
%

% initialize tracking parameters using acquisition results
for svindex = 1:length(sv)
    prn                     = sv(svindex);
    codetemp                = generateCAcode(prn);
    Code(svindex,:)         = [codetemp(end) repmat(codetemp,1,pdi) codetemp(1)];
    Code_plot(svindex,:)    = [codetemp(end-6) codetemp(end-5) codetemp(end-4) codetemp(end-3) codetemp(end-2) codetemp(end-1) codetemp(end) repmat(codetemp,1,pdi) codetemp(1) codetemp(2) codetemp(3) codetemp(4) codetemp(5) codetemp(6)];
    AcqCodeDelay(svindex)   = Acquired.codedelay(svindex);
   
    file_ptr(svindex)       = (signal.Sample - AcqCodeDelay(svindex) +1 ...
        + file.skip *fs*t ... 
        )*file.dataPrecision*file.dataType;
    
    carrFreq(svindex)       = Acquired.fineFreq(svindex);
    AcqFreq(svindex)        = Acquired.fineFreq(svindex);
    
    oldcodedelay_pos(svindex) = 0;
    oldabsoluteSample_pos(svindex) = 0;
end

% parameters for C/N0 estimate
snrIndex	= ones(1,svlength);
K           = 20;
flag_snr    = ones(1,svlength); % flag to calculate C/N0
index_int   = zeros(1,svlength);

eph_idx     = ones(1,svlength);

corrUpdateSec   = 0.1;
corrUpt         = corrUpdateSec/(pdi*t);
counter_corr    = corrUpt-1 * ones(svlength,1);

% Tracking parameters
carrNco      = zeros(1,svlength);
oldCarrNco  = zeros(1,svlength);
oldCarrError       = zeros(1,svlength);
codeNco         = zeros(1,svlength);
code_outputLast     = zeros(1,svlength);
DLLdiscriLast       = zeros(1,svlength);
remChip             = zeros(1,svlength);
codeFreq            = ones(1,svlength)*f0;
remCarrPhase        = zeros(1,svlength);
carrError           = zeros(1,svlength);
codeError           = zeros(1,svlength);
delayValue          = zeros(svlength,datalength/pdi);

%%
localTime = inf;



% Find start and end of measurement point locations in IF signal stream with available
% measurements
sampleStart = zeros(1, svlength);
sampleEnd = inf(1, svlength);

for channelNr = 1:svlength
    prn = sv(channelNr);
    sampleStart(channelNr) = ...
        TckResult_Eph(prn).absoluteSample(sbf.nav1(prn)+eph(prn).sfb(1)*20); % first subframe, in unit of ms
    
    sampleEnd(channelNr) = TckResult_Eph(prn).absoluteSample(end);
end
sampleStartMea = max(sampleStart) + 1; % Now, in  unit of sample
sampleEndMea = min(sampleEnd) - 1;

%--- Measurement step in unit of IF samples -------------------------------
measSampleStep = fix(signal.Fs * solu.navSolPeriod/1000)*file.dataType;%file.dataType;

% flag for position. When the file pointers in all tracking channels exceed
% the current measurement point (in samples), we do positioning
flg_pos = zeros(1,svlength);

% Index for positioning
posIndex = 0;

%%
h = waitbar(0,['Conventional Tracking, Length: ',num2str(datalength),' ms,', '  Please wait...']);
tic

%%%%%%%%%%%%%%%%%%% Multicorrelator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Spacing_plot = [-1.5 -1.25 -1 -0.83 -0.66 -0.5 -0.4 -0.3 -0.2 -0.1 -0.05 0 0.05 0.1 0.2 0.3 0.4 0.5 0.66 0.83 1 1.25 1.5];
Spacing_plot = (-1:0.05:1);
total_corr = size(Spacing_plot,2);
one_side_corr = floor(total_corr/2);
corr_out_plot = zeros(1,total_corr);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RLS VARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_sv=zeros(total_corr,total_corr,svlength);
for subindex=1:svlength
    C_sv(:,:,subindex)=eye(total_corr);
end
enable_RLS = 1;
enable_LMS=1;
enable_MAPA=0;
epsylon=.01; %MAPA
%lr=.0000000000001;
lr=.01;
lambda=0.98;
Mo=0;    %momentum flag
C_ant(:,:,svindex)=C_sv(:,:,svindex);
C_ant(:,:,svindex)=C_sv(:,:,svindex);
beta=.9;
ga=.7;
Mo2=0; %2nd momentum flag
C_antant(:,:,subindex)=C_sv(:,:,subindex);

if enable_RLS || Mo || enable_LMS || enable_MAPA
    new_tracking = 1;
else 
    new_tracking = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = zeros(total_corr, total_corr);
fsttime=1;
G_inv = zeros(total_corr, total_corr);
g_los = zeros(total_corr, 1);
y_unfiltered = zeros(total_corr, 1);
y_filtered = zeros(total_corr, 1);
y_desired = zeros(total_corr, svlength);
a_los_aux = zeros(total_corr, 1);
a_los = 0;
t_los = 0;
delta=0.005;
P_rls=eye(total_corr)/delta;
R=1;
pk=lr;
first_pulse=1;

 count=0;
 D=0;
for msIndex = 1: datalength/pdi % Note that for pdi > 1ms, the index is still denoted as msIndex. 30/04/2020, BING XU
    waitbar(msIndex/(datalength/pdi),h)
    for svindex = 1 :svlength
        prn = sv(svindex);
        
        % read raw data file
        codePhaseStep(svindex) = codeFreq(svindex)/signal.Fs;
        numSample = ceil((signal.codelength*pdi-remChip(svindex))/codePhaseStep(svindex));
        
        delayValue(svindex,msIndex) = numSample - signal.Sample*pdi;
        
        fseek(file.fid, file_ptr(svindex),'bof');
        
        if file.dataPrecision == 2
            rawsignal = fread(file.fid,numSample*file.dataType,'int16')';
            sin_rawsignal = rawsignal(1:2:length(rawsignal));
            cos_rawsignal = rawsignal(2:2:length(rawsignal));
            rawsignal = (sin_rawsignal - mean(sin_rawsignal)) + 1i.*(cos_rawsignal-mean(cos_rawsignal));
        else
            rawsignal = fread(file.fid,numSample*file.dataType,'int8')'; %
            if file.dataType == 2            
                rawsignal = rawsignal(1:2:length(rawsignal)) + 1i*rawsignal(2:2:length(rawsignal));% For NSL STEREO LBand only
            end
        end
        
        file_ptr(svindex)   = file_ptr(svindex) + numSample*file.dataType;  %%%%%%  
        
        if new_tracking %correlators
            time_stamps_plot = zeros(total_corr,numSample); 
            code_replicas_plot = zeros(total_corr,numSample); % tem que ser feito aqui, numSample varia!
            for subindex = 1: total_corr
                time_stamps_plot(subindex,:) = (Spacing_plot(subindex) + remChip(svindex)) : codeFreq(svindex)/signal.Fs : ((numSample -1) * (codeFreq(svindex)/signal.Fs) + Spacing_plot(subindex) + remChip(svindex));
                code_replicas_plot(subindex,:) = Code_plot(svindex,ceil(time_stamps_plot(subindex,:)) + 7); % +7 porquê? No clássico é +1!
            end
        else  % traditinal correlators
            
        
        %% spacing = -0.6:0.05:0.6
        t_CodeEarly       = (0 + Spacing(1) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(1) + remChip(svindex));
        t_CodePrompt      = (0 + Spacing(2) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(2) + remChip(svindex));
        t_CodeLate        = (0 + Spacing(3) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(3) + remChip(svindex));
        
        CodeEarly      = Code(svindex,(ceil(t_CodeEarly) + 1));
        CodePrompt     = Code(svindex,(ceil(t_CodePrompt) + 1));
        CodeLate       = Code(svindex,(ceil(t_CodeLate) + 1));
        
        %%
        end
        
        if new_tracking 
            remChip(svindex) =  time_stamps_plot((total_corr +1)/2,numSample) + codePhaseStep(svindex) - signal.codelength*pdi;
        else
            remChip(svindex) = t_CodePrompt(numSample) + codePhaseStep(svindex) - signal.codelength*pdi;
        end
        
        CarrTime = (0:numSample)./signal.Fs;
        Wave = 2*pi*((carrFreq(svindex)).*CarrTime) + remCarrPhase(svindex);
        remCarrPhase(svindex) = rem(Wave(numSample+1), 2*pi);
        carrsig = exp(1i.* Wave(1:numSample));
        InphaseSignal    = imag(rawsignal .* carrsig);
        QuadratureSignal = real(rawsignal .* carrsig);    
        
        if new_tracking %correlations
            for subindex = 1: total_corr
                corr_inphase(subindex) = sum(code_replicas_plot(subindex,:).*InphaseSignal);
                corr_quadrat(subindex) = sum(code_replicas_plot(subindex,:).*QuadratureSignal);
                corr_out_plot(subindex) = sqrt(corr_inphase(subindex)^2 + corr_quadrat(subindex)^2);
            end
        else
        %%
        E_i      = sum(CodeEarly    .*InphaseSignal);
        E_q      = sum(CodeEarly    .*QuadratureSignal);
        P_i      = sum(CodePrompt   .*InphaseSignal);
        P_q      = sum(CodePrompt   .*QuadratureSignal);
        L_i     = sum(CodeLate     .*InphaseSignal);
        L_q     = sum(CodeLate     .*QuadratureSignal);
        end

            %%%%%%%%%%%%%%%%%% INITIAL STATE %%%%%%%%%%%%%%%%%%
            if (new_tracking && fsttime)
                for subindex = 1 : total_corr
                    for subindex1 = 1 : total_corr 
                        G(subindex1,subindex) = sum(code_replicas_plot(subindex,:).*code_replicas_plot(subindex1,:));
                    end
                end
                G_inv = pinv(G);
                %y_desired(:,svindex) = corr_out_plot'; % no LOS estimation in the first iteration
                %y_filtered = y_desired(:,svindex)';
                fsttime=0;
            end          

            
            %%%%%%%%%%%%%%%%%% ADAPTIVE CHANNEL COMPENSATION %%%%%%%%%%%%%%%%%%
            %Para debug retirar
            y_unfiltered = corr_out_plot;
%             if msIndex > 106
%                 chan=[ 0 0 0 1 .9 .8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
%                 convol=conv(chan, y_unfiltered);
%                 y_unfiltered (1:21) = convol (1:21);
%             end
            %[y_filtered] = filter_rls(y_unfiltered, y_desired(:,svindex)', lambda);

%         if svindex==1
%              y_desired(:,svindex)=y_unfiltered;  % z^-1 do artigo
%              %[y_filtered,chan] = filter_rls(y_unfiltered, y_desired(:,svindex), y_unfiltered', chan, lambda);
%              y_filtered=y_unfiltered;
%         else
%             %y_d=y_desired(:,svindex-1)'-y_filtered;
%             [y_filtered,chan] = filter_rls(y_unfiltered, y_desired(:,svindex-1)', chan, lambda);
%         end
            %%%%%%%%%%%%%%%%%% LOS SIGNAL ESTIMATION BLOCK %%%%%%%%%%%%%%%%%%
         if new_tracking
             %%% LOS estimation
            for subindex = 1 : total_corr
                for subindex1 = 1 : total_corr 
                    g_los(subindex1,1) = sum(code_replicas_plot(subindex,:).*code_replicas_plot(subindex1,:));
                end
                a_los_aux(subindex,1)=(g_los'*G_inv*C_sv(:,:,svindex)*y_unfiltered')/(g_los'*G_inv*g_los);        
            end
            [a_los, t_los] = max(a_los_aux(:,1));
            
            %%%%%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               [amp pico] = max (y_unfiltered);
%               [t_los t_los - pico]
%             if (abs(t_los - pico) > 2)
%                 [msIndex svindex pico t_los count+1]
%                 if (abs(t_los - pico) > 3)
%                     3
%                 end
%                 if (abs(t_los - pico) > 4)
%                     4
%                 end
%                 
%                
%                 count=count+1;
%             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %end
          for subindex = 1 : total_corr 
              y_desired(subindex,svindex) = a_los*sum(code_replicas_plot(subindex,:).*code_replicas_plot(t_los,:));
          end 
         end
            
            
         %%%%%%%%%%%%%%%%%%%%%%%% Channel Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
         if enable_LMS
             for i=1:5
                 if Mo
                     if (msIndex==101 && i==1) %1st iteration (hapens once for each satellite. i01 senõ entra para cada i do for
                         C_ant(:,:,svindex)=C_sv(:,:,svindex);
                         C_aux(:,:,svindex)=C_sv(:,:,svindex);
                     else
                          C_aux(:,:,svindex)=C_sv(:,:,svindex);
                     end
                     C_sv(:,:,svindex) = C_sv(:,:,svindex) + lr*y_unfiltered'*(y_desired(:,svindex)'-y_unfiltered*C_sv(:,:,svindex))+beta*(C_sv(:,:,svindex)-C_ant(:,:,svindex));
                     C_ant(:,:,svindex)=C_aux(:,:,svindex);
                 elseif Mo2
                     if msIndex==101 %1st iteration (hapens once for each satellite
                         C_ant(:,:,svindex)=C_sv(:,:,svindex);
                         C_aux(:,:,svindex)=C_sv(:,:,svindex);
                         C_antant(:,:,svindex)=C_sv(:,:,svindex);
                     else
                         C_aux(:,:,svindex)=C_sv(:,:,svindex);
                     end
                     C_sv(:,:,svindex) = C_sv(:,:,svindex) + lr*y_unfiltered'*(y_desired(:,svindex)'-y_unfiltered*C_sv(:,:,svindex))+(beta+ga)*(C_sv(:,:,svindex)-C_ant(:,:,svindex)) -ga*beta*(C_ant(:,:,svindex)-C_antant(:,:,svindex));
                     C_antant(:,:,svindex)=C_ant(:,:,svindex);
                     C_ant(:,:,svindex)=C_aux(:,:,svindex);
              
                 else
                     miu=pk/(y_unfiltered*pk*y_unfiltered' + 1*R); %R
                     pk=(1-miu*y_unfiltered*y_unfiltered'*pk);
                     if (svindex ==5 && msIndex<3000)
                         lr=0.01;
                     else
                         lr=0.1;
                     end
                     if (msIndex == 3000) && (svindex==1)
                         lr=lr/10;
                     elseif (msIndex == 5000) && (svindex == 1)
                         lr=lr/10;
                     elseif (msIndex == 7000) && (svindex ==1)
                         lr=lr/10;
                     end
                    C_sv(:,:,svindex) = C_sv(:,:,svindex) + 1/(y_unfiltered*y_unfiltered')*lr*y_unfiltered'*(y_desired(:,svindex)'-y_unfiltered*C_sv(:,:,svindex));
                 end % 1/(y_unfiltered*y_unfiltered')*lr*1e13 ->este é o NMLS que até decresce mais depressa no início mas não chega tão baixo como o Kalman
             end     % falta tentar explorar o R. Experimentei 100R e deu mais baixo ainda. 
             
         elseif enable_RLS && ~enable_LMS &&~enable_MAPA
         
            P_rls=eye(total_corr)/delta;
          
            for i=1:20  % RLS iteration number 
                
                kappa = lambda^(-1)*P_rls*y_unfiltered'/(1+lambda^(-1)*y_unfiltered*P_rls*y_unfiltered');
                P_rls = lambda^(-1)*P_rls-lambda^(-1)*kappa*y_unfiltered*P_rls;
                C_sv(:,:,svindex) = C_sv(:,:,svindex) + P_rls*y_unfiltered'*(y_desired(:,svindex)'-y_unfiltered*C_sv(:,:,svindex)); %-1? Assim (sem -1) deu melhor
            end
         elseif enable_MAPA      %MAPA por acabar
             D=y_desired(:,svindex:-1:svindex-K+1);
             Y=y_unfiltered(:,svindex:-1:svindex-K+1);
             C = C + lr*Y'*inv(epsylon*eye(total_corr)+YY')*(D-Y*C);
         else
             ;
             
            
         end
         %end
        %end
        %%%%%%%%%%%%%%%%%% ADAPTIVE CHANNEL COMPENSATION %%%%%%%%%%%%%%%%%%
      
        
        
        

        % Implement code loop filter and generate NCO command
        if new_tracking
            codeError(svindex) = 1*(Spacing_plot(one_side_corr+1)-Spacing_plot(t_los));  % DLL discriminator
        else
            E = sqrt(E_i^2+E_q^2);
            L = sqrt(L_i^2+L_q^2);
            codeError(svindex) = 0.5*(E-L)/(E+L);  % DLL discriminator  
        end       
         if new_tracking % No artigo há aqui o loop filter. Analisar + tarde 
              codeNco(svindex) = (tau2code/tau1code)*2.5*codeError(svindex); % precisa de um ganho
              codeFreq(svindex) = signal.codeFreqBasis - codeNco(svindex);
         
              %For PLL discriminator and localization calculations
              
              P_i = corr_inphase((total_corr +1)/2); 
              P_q = corr_quadrat((total_corr +1)/2); 
              E_i = corr_inphase((total_corr -1)/2); 
              E_q = corr_quadrat((total_corr -1)/2);
              L_i = corr_inphase((total_corr +3)/2); 
              L_q = corr_quadrat((total_corr +3)/2);
         else
        %codeError(svindex) = 0.5*(E-L)/(E+L);  % DLL discriminator
        codeNco(svindex) = code_outputLast(svindex) + (tau2code/tau1code)*(codeError(svindex)...
            - DLLdiscriLast(svindex)) + codeError(svindex)* ((pdi*t)/tau1code);
        DLLdiscriLast(svindex) = codeError(svindex);
        code_outputLast(svindex) = codeNco(svindex);
        codeFreq(svindex) = signal.codeFreqBasis - codeNco(svindex);
        %codeFreq(svindex) = signal.codeFreqBasis + codeNco(svindex);
         end
        



        % PLL discriminator
        carrError(svindex) = atan(P_q/P_i)/(2*pi);  % PLL discriminator
        carrNco(svindex) = oldCarrNco(svindex) + (tau2carr/tau1carr)*(carrError(svindex) ...
            - oldCarrError(svindex)) + carrError(svindex) * ((pdi*t)/tau1carr);
        oldCarrNco(svindex) = carrNco(svindex);
        oldCarrError(svindex) = carrError(svindex);
        carrFreq(svindex)  = AcqFreq(svindex) + carrNco(svindex);  % Modify carrier freq

        % Calculate CN0
        if (flag_snr(svindex) == 1)
            index_int(svindex) = index_int(svindex) + 1;
            Zk(svindex,index_int(svindex)) = P_i^2 + P_q^2;
            if mod(index_int(svindex),K) == 0
                meanZk  = mean(Zk(svindex,:));
                varZk   = var(Zk(svindex,:));
                NA2     = sqrt(meanZk^2-varZk);
                varIQ   = 0.5 * (meanZk - NA2);
                CN0_CT(snrIndex(svindex),svindex) =  abs(10*log10(1/(1*t*pdi) * NA2/(2*varIQ)));
                index_int(svindex)  = 0;
                snrIndex(svindex)   = snrIndex(svindex) + 1;
            end
        end
        
        %% Data Recording
        TckResultCT(prn).E_i(msIndex) = E_i;
        TckResultCT(prn).E_q(msIndex) = E_q;
        TckResultCT(prn).P_i(msIndex) = P_i;
        TckResultCT(prn).P_q(msIndex) = P_q;
        TckResultCT(prn).L_i(msIndex) = L_i;
        TckResultCT(prn).L_q(msIndex) = L_q;
        
        
        TckResultCT(prn).carrError(msIndex)       = carrError(svindex);
        TckResultCT(prn).codeError(msIndex)       = codeError(svindex);
        TckResultCT(prn).codeFreq(msIndex)        = codeFreq(svindex);
        TckResultCT(prn).carrFreq(msIndex)        = carrFreq(svindex);
        TckResultCT(prn).numSample(msIndex)       = numSample;
        TckResultCT(prn).remChip(msIndex)         = remChip(svindex);
        TckResultCT(prn).remCarrPhase(msIndex)    = remCarrPhase(svindex);
        TckResultCT(prn).absoluteSample(msIndex)  = ftell(file.fid);
        TckResultCT(prn).absoluteSampleCodedelay(msIndex)  = mod(TckResultCT(prn).absoluteSample(msIndex)/(file.dataPrecision*file.dataType),fs*t );
        TckResultCT(prn).codedelay(msIndex)       = signal.Sample - AcqCodeDelay(svindex) +1 + sum(delayValue(svindex,(1:msIndex)));
        TckResultCT(prn).codedelay2(msIndex)      = mod( TckResultCT(prn).absoluteSample(msIndex)/(file.dataPrecision*file.dataType),fs*t );
        TckResultCT(prn).delayValue(msIndex)      = delayValue(svindex,msIndex); 
    end % end for svindex in Tracking    
    
    %%
    % Position index of current measurement time in IF signal stream
    % (in unit IF signal sample point)
    currMeasSample = sampleStartMea + measSampleStep*posIndex;
    
    for svIndex=1:svlength
        prn = sv(svIndex);
        if TckResultCT(prn).absoluteSample(msIndex) > currMeasSample
            flg_pos(svIndex) = 1;
        else
            flg_pos(svIndex) = 0;
        end
    end
    
    %%
    if sum(flg_pos) == svlength  
        posIndex = posIndex + 1;
        for svIndex=1:svlength
            prn = sv(svIndex);
            
            % Find index of I_P stream whose integration contains current
            % measurment point location
            for index = 1: length(TckResultCT(prn).absoluteSample)
                if(TckResultCT(prn).absoluteSample(index) > currMeasSample )
                    break
                end
            end
            index = index - 1;
            
            % Update the phasestep based on code freq and sampling frequency
            codePhaseStepX = TckResultCT(prn).codeFreq(index)/signal.Fs;
            
            codePhaseMeas(svIndex) = TckResultCT(prn).remChip(index) + ...
                codePhaseStepX*((currMeasSample - TckResultCT(prn).absoluteSample(index))/file.dataType);%file.dataType
            
            transmitTime(svIndex) = codePhaseMeas(svIndex)/signal.codelength/1000 + ... %TckResultCT(prn).codedelay(msIndex)/(signal.Fs/1000) 
                (index - (sbf.nav1(prn)+eph(prn).sfb(1)*20))/1000 + ...
                eph(prn).TOW(1);
        end
      
        % At first time of fix, local time is initialized by transmitTime and
        % settings.startOffset
        if (localTime == inf)
            maxTime   = max(transmitTime);
            localTime = maxTime + 75/1000; % 68 ms is an assumed travel time
        end
        pseudorange = (ones(1,svlength).*localTime - transmitTime)*cmn.cSpeed;
        
        %
        usr_clk = usr_clk_wls ;%%%
        estusr = estusr_wls;%%%
        
        for svindex = 1 : svlength
            prn = sv(svindex);
            
            tot_est_pos(svindex) = transmitTime(svindex);% ...
%                                                 + (1/cmn.cSpeed)*sv_clk(prn);
            
            % find the sv pos in ECEF at the time of transmision
            [svxyz(svindex,:), sv_vel(svindex,:), sv_clk(prn), sv_clk_vel(prn), grpdel] = ...
                svPosVel(prn,eph,tot_est_pos(svindex),eph_idx(svindex));
            
            % C/A-code pseudorange corrected for satellite clock (in meters) and Tgd(in sec)
            prvec(svindex)      = pseudorange(svindex) + sv_clk(prn) - grpdel*cmn.cSpeed;% -sv_clk(prn)?
            
            % Adjust satellite position coordinates for earth rotation correction
            svxyzr(svindex,:)   = erotcorr(svxyz(svindex,:),prvec(svindex));
            
            % tropospheric and ionospheric delay correction
            counter_corr(svindex) = counter_corr(svindex) + 1;
            if counter_corr(svindex) ==  corrUpt
                svenu           = xyz2enu(svxyzr(svindex,:), estusr(1:3));%
                el_rad(svindex) = atan(svenu(3)/norm(svenu(1:2)));
                %             az_rad(svindex) = (pi/2)-atan2(svenu(1),svenu(2));
                az_rad(svindex) =  atan2(svenu(1),svenu(2));
                az(svindex)     = az_rad(svindex)*180/pi;
                el(svindex)     = el_rad(svindex)*180/pi;
                temp            = xyz2llh(estusr(1:3));
                user_ll         = [temp(1:2).*180/pi temp(3)];
                ionodel(svindex)        = ionocorr(tot_est_pos(svindex),svxyzr(svindex,:), estusr(1:3));
                tropodel_unb3(svindex)  = abs(trop_UNB3(cmn.doy,user_ll(1),user_ll(3),el(svindex)));
                counter_corr(svindex)   = 0;
            end
            
            prvec(svindex) = prvec(svindex) - ionodel(svindex) - tropodel_unb3(svindex); % sign of iono. and trop. error?
        end % for svindex=1:svlength
        
        
        %% Record Ppseudorange measurement 
        navSolutionsWLS.rawPseudorange(posIndex,:) = pseudorange ;        
        
        %% Position cal using LS method
        [estusr_wls, dop]       = olspos(prvec,svxyzr,estusr_wls); % ordinary least square
        [VR, dtRV, ~]     = ...
            LS_SA_code_Vel_xubing(estusr_wls(1:3)', svxyzr, sv_vel, carrFreq', 0.190293672798365, sv_clk_vel(sv));
        
        usrenu_wls(posIndex,:)   = xyz2enu(estusr_wls(1:3),cnslxyz);
        usr_clk_wls             = estusr_wls(4);
        
        llh     = xyz2llh(estusr_wls(1:3));
        L_b     = llh(1);
        lamda_b = llh(2);
        C_e_n = [ -sin(lamda_b)           cos(lamda_b)         	 0;...
            -sin(L_b)*cos(lamda_b) -sin(L_b)*sin(lamda_b)	 cos(L_b);...
            -cos(L_b)*cos(lamda_b) -cos(L_b)*sin(lamda_b)	-sin(L_b);];
        usr_velENU(posIndex,:) = C_e_n * VR ;
        
        usrenu_wls(posIndex,:)                 	= xyz2enu(estusr_wls(1:3),cnslxyz);
        usrllh_wls(posIndex,:)                   = xyz2llh(estusr_wls(1:3));
        usrllh_wls(posIndex,1:2)                 = usrllh_wls(posIndex,1:2)*180/pi;
        
        %     navSolutionsWLS.RxTime(msIndex)       = RxTime(msIndex);
        navSolutionsWLS.usrPos(posIndex,:)       = estusr_wls(1:3);
        navSolutionsWLS.usrVel(posIndex,:)       = VR;
        navSolutionsWLS.usrPosENU(posIndex,:)    = usrenu_wls(posIndex,:);
        navSolutionsWLS.usrPosLLH(posIndex,:)    = usrllh_wls(posIndex,:);
        navSolutionsWLS.clkBias(posIndex)        = usr_clk_wls;
        navSolutionsWLS.usrVelENU(posIndex,:)        = usr_velENU(posIndex,:);
        navSolutionsWLS.clkDrift(posIndex)   = dtRV; % m/s
        navSolutionsWLS.DOP(posIndex,:)       = dop;
        navSolutionsWLS.satEA(posIndex,:)      = el;
        navSolutionsWLS.satAZ(posIndex,:)      = az;
        
        navSolutionsWLS.codePhaseMeas(posIndex,:) = codePhaseMeas;
        %     navSolutionsWLS.test(msIndex,:)      = test;
        
        
        %=== Correct local time by clock error estimation =================
        localTime = localTime - navSolutionsWLS.clkBias(posIndex)/cmn.cSpeed;
        navSolutionsWLS.localTime(posIndex,:) = localTime;
        
        %=== Update local time by measurement  step  ====================================
        localTime = localTime + measSampleStep/signal.Fs;
        
        navSolutionsCT = navSolutionsWLS;
        
        if mod(posIndex, 1) == 0
            fprintf('WLS: msIndex = %4d; PosIndex = %4d; localTime: %f;  2D Err = %f VE = %f B = %f D = %f\n\n', ...
                msIndex, posIndex, localTime, sqrt(usrenu_wls(posIndex,1)^2+usrenu_wls(posIndex,2)^2),usr_velENU(posIndex,1), usr_clk_wls, dtRV);
        end
    end % end for positioning at current measurement epoch
    
end % end for msIndex

close(h);

TckResultCT_pos = TckResultCT;

save(['navSolCT_',num2str(pdi),'ms_',file.fileName], 'navSolutionsCT' );
save(['tckRstCT_',num2str(pdi),'ms_',file.fileName], 'TckResultCT_pos','CN0_CT');

function [h] = filter_ls_to_rls(u, d, lambda)
  
    N=length(u);%number of correlators
    delta=0.005;
    
    P = eye(N)/delta;
    uvec=zeros(N,1);
    h=zeros(N,1);

    for i=1:N
        uvec=u(i,:)';
        kappa = (lambda^(-1)*P*uvec)/(1+lambda^(-1)*uvec'*P*uvec);
        e=d(i)-h'*uvec;
        h=h+kappa*conj(e);
        P = lambda^(-1)*P-lambda^(-1)*kappa*uvec'*P;
    end
end


function [y,w] = filter_rls(u, d,w, lambda)
 
    delta=0.005;
    N=length(d);%number of correlators
    d=d(:);
    uvec=zeros(N,1);
    %w=zeros(N,1);
    P_rls=eye(N)/delta;
    
    for i=1:N
        uvec(1:1:i)=u(i:-1:1);
        
        kappa = lambda^(-1)*P_rls*uvec/(1+lambda^(-1)*uvec'*P_rls*uvec);
        e=d(i)-w'*uvec;
        %e=d(i)-d1(i);
        w=w+kappa*conj(e);
        P_rls = lambda^(-1)*P_rls-lambda^(-1)*kappa*uvec'*P_rls;
    end

    y=filter(w,1,u);      %filtered signal
end

%%%%%   Matrix triangular decimation   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  if i>=6
%                      Kap=12;
%                      for lin=total_corr-Kap+1:total_corr
%                          for col=1:total_corr-Kap-1
%                              if col<lin-Kap
%                                  C(lin,col)=0;
%                              end
%                          end
%                      end
%                      for col=total_corr-Kap+1:total_corr
%                          for lin=1:total_corr-Kap-1
%                              if lin<col-Kap
%                                  C(lin,col)=0;
%                              end
%                          end
%                      end
%                  end


end