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




%%%%%%%% RLS VARS %%%%%%%%
corr_per_chip_plot = 10; % correlators per chip
one_side_chip = 1.5; % total correlator spacing (e.g. -1 to 1, -1.5 to 1.5)
one_side_corr = floor(corr_per_chip_plot*one_side_chip);
total_corr = corr_per_chip_plot*one_side_chip*2+1;
lambda=0.98;
d = zeros(1,total_corr);
u = zeros(1,total_corr);
Spacing_plot = zeros(1,total_corr);
for Index = 1: total_corr
    Spacing_plot(Index) = -one_side_chip + (Index-1)/corr_per_chip_plot;
end
corr_out_plot = zeros(1,total_corr);

enable_RLS = 1;
enable_RLS_compensation_plot = 0;

if enable_RLS & enable_RLS_compensation_plot
    colororder({'[0.4660 0.6740 0.1880]','#D95319'})
    corr_func = figure(1);
    subplot(5,3,1:3);
    yyaxis left
    p0=plot(Spacing_plot, zeros(1,total_corr),'*-','Color',"#77AC30",'DisplayName','Raw Correlation');
    ylabel('Amplitude');
    xlabel('Delay [chips]');
    title('Correlation function');
    %ylim([0 10000])
    xlim([-1.5 1.5])
    hold on;
    yyaxis right
    p0_0=plot(Spacing_plot, zeros(1,total_corr),'*-','Color',"#D95319",'DisplayName','Ideal Correlation');
    legend

    p1sb=subplot(5,3,4);
    p1=plot(zeros(1,1),zeros(1,1),'-*');
    hold on;
    p1_1=plot(zeros(1,1),zeros(1,1),'-*','Color',"#D95319");
    ylabel('Amplitude');
    xlabel('Delay [chips]');
    title('Channel Impulse Response');
    %ylim([-1 1])

    subplot(5,3,5);
    p2=plot(Spacing_plot,zeros(1,total_corr),'-*');
    ylabel('Amplitude');
    xlabel('Delay [chips]');
    title('Correlation Function from Impulse Response');
    %ylim([0 10000])
    xlim([-1.5 1.5])  

    p3sb=subplot(5,3,7);
    p3=plot(zeros(1,1),zeros(1,1),'-*');
    ylabel('Amplitude');
    xlabel('Delay [chips]');
    title('New Channel Impulse Response');
    %ylim([-1 1])

    subplot(5,3,8);
    p4=plot(Spacing_plot,zeros(1,total_corr),'-*');
    ylabel('Amplitude');
    xlabel('Delay [chips]');
    title('Multipath Correlation Function');
    %ylim([0 10000])
    xlim([-1.5 1.5])  
    
    subplot(5,3,10:12);
    p5=plot(Spacing_plot,zeros(1,total_corr),'-*');
    ylabel('Amplitude');
    xlabel('Delay [chips]');
    title('Multipath-free Correlation Function');
    %ylim([0 20000])
    xlim([-1.5 1.5])  

    subplot(5,3,[6 9]);
    p6=plot(zeros(1,1),zeros(1,1),'-*');
    ylabel('Error');
    xlabel('Nº multipaths');
    title('Impulse Response Size - RMSE');

    p7sb=subplot(5,3,13:15);
    p7=plot(Spacing_plot,zeros(1,total_corr),'-*');
    ylabel('Amplitude');
    xlabel('Delay [chips]');
    title('Multipath-free Correlation Function v2');
    %ylim([0 20000])
    xlim([-1.5 1.5])
end

%%
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
        
        %% spacing = -0.6:0.05:0.6
        t_CodeEarly       = (0 + Spacing(1) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(1) + remChip(svindex));
        t_CodePrompt      = (0 + Spacing(2) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(2) + remChip(svindex));
        t_CodeLate        = (0 + Spacing(3) + remChip(svindex)) : codePhaseStep(svindex) : ((numSample -1) * codePhaseStep(svindex) + Spacing(3) + remChip(svindex));
        
        CodeEarly      = Code(svindex,(ceil(t_CodeEarly) + 1));
        CodePrompt     = Code(svindex,(ceil(t_CodePrompt) + 1));
        CodeLate       = Code(svindex,(ceil(t_CodeLate) + 1));
        
        %%
        remChip(svindex) = t_CodePrompt(numSample) + codePhaseStep(svindex) - signal.codelength*pdi;
        
        CarrTime = (0:numSample)./signal.Fs;
        Wave = 2*pi*((carrFreq(svindex)).*CarrTime) + remCarrPhase(svindex);
        remCarrPhase(svindex) = rem(Wave(numSample+1), 2*pi);
        carrsig = exp(1i.* Wave(1:numSample));
        InphaseSignal    = imag(rawsignal .* carrsig);
        QuadratureSignal = real(rawsignal .* carrsig);        
        
        %%
        E_i      = sum(CodeEarly    .*InphaseSignal);
        E_q      = sum(CodeEarly    .*QuadratureSignal);
        P_i      = sum(CodePrompt   .*InphaseSignal);
        P_q      = sum(CodePrompt   .*QuadratureSignal);
        L_i     = sum(CodeLate     .*InphaseSignal);
        L_q     = sum(CodeLate     .*QuadratureSignal);



        %%%%%%%% RLS TRACKING %%%%%%%%
        if enable_RLS & (msIndex > 100) % wait for tracking pull-in

            %%%%%%%% HIGH RESOLUTION CORRELATOR %%%%%%%%
            time_stamps_plot = zeros(total_corr,numSample);
            code_replicas_plot = zeros(total_corr,numSample);
            for subindex = 1: total_corr
                time_stamps_plot(subindex,:) = (Spacing_plot(subindex) + remChip(svindex)) : codeFreq(svindex)/signal.Fs : ((numSample -1) * (codeFreq(svindex)/signal.Fs) + Spacing_plot(subindex) + remChip(svindex));
                code_replicas_plot(subindex,:) = Code_plot(svindex,ceil(time_stamps_plot(subindex,:)) + 7);
                corr_out_plot(subindex) = sqrt(sum(code_replicas_plot(subindex,:).*InphaseSignal)^2 + sum(code_replicas_plot(subindex,:).*QuadratureSignal)^2); %pll should bring signal to inphase
            end
%{


            %%%%%%%% IDEAL CORRELATION FUNCTION %%%%%%%% 
            u=abs(corr_out_plot(:));
            for subindex = 1 : total_corr
                d(1,subindex) = sum(code_replicas_plot(subindex,:).*code_replicas_plot(one_side_corr+1,:));
            end

            %%%%%%%% INVERSE CHANNEL %%%%%%%%
            [ch_imp_resp] = filter_rls2(u, d, lambda);
            %[ch_imp_resp] = filter_rls3(u, d, lambda, 30); % use larger (more precise) impulse response

            %%%%%%%% FIND LOS IMPULSE %%%%%%%%
            [a_los, t_los] = max(ch_imp_resp);
            %[a_los, t_los] = max(ch_imp_resp(1:10)); % but search for LOS only in the first 10 delays (avoid side peaks)
            
            %%%%%%%% ARTIFITIAL CORRELATION FUNCTION %%%%%%%%
            for subindex = 1 : total_corr
                f_corr_out(1,subindex) = sum(code_replicas_plot(subindex,:).*code_replicas_plot(one_side_corr+1+t_los-2,:));
            end

            E_i = f_corr_out(1,one_side_corr); % new correlator values for DLL corrections
            L_i = f_corr_out(1,one_side_corr+2);
            %E_i = corr_out_plot(one_side_corr);
            %L_i = corr_out_plot(one_side_corr+2);

            %%%%%%%% PLOT %%%%%%%%
            if (svindex == 4) & enable_RLS_compensation_plot
                set(p0,'XData',time_stamps_plot(:,1),'YData',abs(corr_out_plot));
                set(p0_0,'XData',time_stamps_plot(:,1),'YData',d);
                set(p1,'XData',(0:1/corr_per_chip_plot:(1/corr_per_chip_plot)*(size(ch_imp_resp,1)-1)),'YData',ch_imp_resp);
                set(p1_1,'XData',(1/corr_per_chip_plot)*(t_los-1),'YData',a_los);
                set(p5,'XData',time_stamps_plot(:,1),'YData',f_corr_out);
                pause(0.01);
            end
%}           

            
            %%%%%%%% IDEAL CORRELATION FUNCTION %%%%%%%% 
            u=abs(corr_out_plot);
            for subindex = 1 : total_corr
                d(1,subindex) = sum(code_replicas_plot(subindex,:).*code_replicas_plot(one_side_corr+1,:));
            end

            %%%%%%%% SEARCH BEST CHANNEL IMPULSE RESPONSE %%%%%%%% 
            for multp_num = 1 : 40
                [ch_imp_resp] = filter_rls(u, d, lambda, multp_num);
                x_corr_out=conv(d,ch_imp_resp);
                x_corr_out=x_corr_out(1:1:total_corr);
                if size(x_corr_out,1) ~= size(u,1)
                    x_corr_out=x_corr_out';
                end
                rmse(1,multp_num) = sqrt(mean((u(1:1:end) - x_corr_out(1:1:end)).^2));
            end
            [a0, t0] = min(rmse);
            
            %%%%%%%% MULTIPATH MITIGATION 1 %%%%%%%%
            [ch_imp_resp] = filter_rls(u, d, lambda, t0); % optimize to avoid running rls again
            t0=ceil(corr_per_chip_plot/2); % force a specific impulse response size to ignore side peaks
            ch_imp_resp = ch_imp_resp(1:1:t0);
            ch_imp_resp1 = ch_imp_resp;
            ch_imp_resp1(1) = 0; % preserve the LOS impulse and remove all others
            x_corr_out=conv(d,ch_imp_resp1(1:1:end)); % get correlation function from interferences
            x_corr_out=x_corr_out(1:1:total_corr);
            for subindex = 1: total_corr
                new_corr(subindex) = u(subindex) - x_corr_out(subindex);
            end

            %%%%%%%% MULTIPATH MITIGATION 2 %%%%%%%%
            multp_max_delay = t0; % remove multipaths up to a certain delay (max would be half correlators)
            new_corr1 = u';
            for subindex = 2 : 2
                if (ch_imp_resp1(subindex) > 0) % only for positive impulses (!?)
                    for subindex1 = 1: total_corr
                        xd(subindex1,1) = sum(code_replicas_plot(subindex1,:).*code_replicas_plot(corr_per_chip_plot+subindex,:));
                    end
                    new_corr1 = new_corr1 - xd.*ch_imp_resp1(subindex);
                end
            end

            E_i = new_corr1(one_side_corr); % new correlator values for DLL corrections
            L_i = new_corr1(one_side_corr+2);

            %%%%%%%% PLOT %%%%%%%%
            if (svindex == 6) & enable_RLS_compensation_plot
                set(p0,'XData',time_stamps_plot(:,1),'YData',abs(corr_out_plot));
                set(p0_0,'XData',time_stamps_plot(:,1),'YData',d);
                figure(corr_func);
                subplot(p1sb);
                cla
                stem((0:1/corr_per_chip_plot:(1/corr_per_chip_plot)*(t0-1)),ch_imp_resp,'filled');
                axis padded
                %set(p1,'XData',(0:1/corr_per_chip_plot:(1/corr_per_chip_plot)*(t0-1)),'YData',ch_imp_resp);
                corr0=conv(d,ch_imp_resp(1:1:end));
                set(p2,'XData',time_stamps_plot(:,1),'YData',corr0(1:1:total_corr));
                subplot(p3sb);
                cla
                stem((0:1/corr_per_chip_plot:(1/corr_per_chip_plot)*(t0-1)),ch_imp_resp1,'filled');
                axis padded
                %set(p3,'XData',(0:1/corr_per_chip_plot:(1/corr_per_chip_plot)*(t0-1)),'YData',ch_imp_resp1);
                set(p4,'XData',time_stamps_plot(:,1),'YData',x_corr_out);
                set(p5,'XData',time_stamps_plot(:,1),'YData',new_corr);
                set(p6,'XData',(1:1:40),'YData',rmse);
                set(p7,'XData',time_stamps_plot(:,1),'YData',new_corr1);
                pause(0.01);
            end
            %{
            %} 
            
        end


        
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
        
        % Implement code loop filter and generate NCO command
        E = sqrt(E_i^2);
        L = sqrt(L_i^2);
        %E = sqrt(E_i^2+E_q^2);
        %L = sqrt(L_i^2+L_q^2);
        codeError(svindex) = 1*(E-L)/(E+L);  % DLL discriminator
        codeNco(svindex) = code_outputLast(svindex) + (tau2code/tau1code)*(codeError(svindex)...
            - DLLdiscriLast(svindex)) + codeError(svindex)* ((pdi*t)/tau1code);
        DLLdiscriLast(svindex) = codeError(svindex);
        code_outputLast(svindex) = codeNco(svindex);
                 codeFreq(svindex) = signal.codeFreqBasis - codeNco(svindex);
        %codeFreq(svindex) = signal.codeFreqBasis + codeNco(svindex);
        
        % PLL discriminator
        carrError(svindex) = atan(P_q/P_i)/(2*pi);  % PLL discriminator
        carrNco(svindex) = oldCarrNco(svindex) + (tau2carr/tau1carr)*(carrError(svindex) ...
            - oldCarrError(svindex)) + carrError(svindex) * ((pdi*t)/tau1carr);
        oldCarrNco(svindex) = carrNco(svindex);
        oldCarrError(svindex) = carrError(svindex);
        carrFreq(svindex)  = AcqFreq(svindex) + carrNco(svindex);  % Modify carrier freq
        
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
            LS_SA_code_Vel_xubing(estusr_wls(1:3)', svxyzr, sv_vel, -carrFreq', 0.190293672798365, sv_clk_vel(sv));
        
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

function [h] = filter_rls(u, d, lambda, M)
  
    N=length(d);%number of correlators
    d=d(:);
    
    %M=ceil(length(d)/2);
    delta=0.005;
    P = eye(M)/delta;
    uvec=zeros(M,1);
    h=zeros(M,1);
    for i=1:N
        if i<M
            uvec(1:1:i)=d(i:-1:1);
        else
            uvec=d(i:-1:i-M+1);
        end
        kappa = lambda^(-1)*P*uvec/(1+lambda^(-1)*uvec'*P*uvec);
        e=u(i)-h'*uvec;
        h=h+kappa*conj(e);
        P = lambda^(-1)*P-lambda^(-1)*kappa*uvec'*P;
    end
end

function [h] = filter_rls2(u, d, lambda)
      
    N=length(d);%number of correlators
    d=d(:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Multipath is just the inverse system
    for i=1:1
        d=[d; 0];
        u=[0; u];
    end
    M=ceil(length(d)/2);
    delta=0.005;
    P = eye(M)/delta;
    uvec=zeros(M,1);
    h=zeros(M,1);
    for i=1:N+1
        if i<M
            uvec(1:1:i)=d(i:-1:1);
        else
            uvec=d(i:-1:i-M+1);
        end
        kappa = lambda^(-1)*P*uvec/(1+lambda^(-1)*uvec'*P*uvec);
        e=u(i)-h'*uvec;
        h=h+kappa*conj(e);
        P = lambda^(-1)*P-lambda^(-1)*kappa*uvec'*P;
    end
end

function [h] = filter_rls3(u, d, lambda, M)
  
    N=length(d);%number of correlators
    d=d(:);
    for i=1:1
        d=[d; 0];
        u=[0; u];
    end
    %M=ceil(length(d)/2);
    delta=0.005;
    P = eye(M)/delta;
    uvec=zeros(M,1);
    h=zeros(M,1);
    for i=1:N
        if i<M
            uvec(1:1:i)=d(i:-1:1);
        else
            uvec=d(i:-1:i-M+1);
        end
        kappa = lambda^(-1)*P*uvec/(1+lambda^(-1)*uvec'*P*uvec);
        e=u(i)-h'*uvec;
        h=h+kappa*conj(e);
        P = lambda^(-1)*P-lambda^(-1)*kappa*uvec'*P;
    end
end
end