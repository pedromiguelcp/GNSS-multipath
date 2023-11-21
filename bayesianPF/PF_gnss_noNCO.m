clear all;
% close all;
clf
clc;

% Define the parameters for the simulation
c = 299792458; % m/s
fs = 20000000; % Sample rate (Hz)
codeFreqBasis = 1.023e6; % C/A code rate (chips/s)
codeFreq = codeFreqBasis;
codelength = 1023;

% correlators will always be aligned with samples
disp(['Maximum number of correlators per chip: ',num2str(ceil(fs/codeFreqBasis)+1)]);

N_corr = 41; % number of correlators
numSample=fs*0.001;
signals = zeros(numSample, N_corr);
signals_deriv = zeros(numSample, N_corr);
shiftidx= zeros(N_corr,1);
% Generate the C/A code
CodeVec = generateCAcode(10); % Generate C/A code for satellite 1
CodeVec = [CodeVec(end) CodeVec];
time_stamps = 0 : codeFreq/fs : ((numSample -1) * (codeFreq/fs));
CodeVecSampled = CodeVec(ceil(time_stamps)+1);
extended_CodeVec=[CodeVecSampled(end) CodeVecSampled CodeVecSampled(1)];
for iCorr = 1:N_corr
    shiftidx(iCorr,1)=(1/((ceil(N_corr/2)-1)/ceil(fs/codeFreqBasis)))*(iCorr-1)-ceil(fs/codeFreqBasis);
    signals(:, iCorr) = circshift(CodeVecSampled', shiftidx(iCorr,1));
end


%Spacing = (-1:0.1:1); % change here the number of correlators (0.1=21correlators, 0.05=41correlators ...)
delayRes = 1/((ceil(N_corr/2)-1)/ceil(fs/codeFreqBasis))*(codeFreqBasis/fs);
Spacing = (codeFreqBasis/fs).*shiftidx';
corr_per_chip = floor(N_corr/2);
code_replicas = zeros(N_corr,numSample);
corr_out = zeros(1,N_corr)';
VE = corr_per_chip-1;
E = corr_per_chip;
P = corr_per_chip+1;
L = corr_per_chip+2;
VL = corr_per_chip+3;

% Specify here the paths delay (chips) and attenuation (linear factor), including LOS (first)
% Channel 3
chan3.delays=[0 2*delayRes 3*delayRes];
chan3.attenuation=[1 3 4];
% Channel 4
chan4.delays=0;
chan4.attenuation=1;


epochs = 5000; %signal tracking epochs

DLLdiscri = zeros(1,epochs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PF %%
T_c=0.001/codelength; % chip duration
L_h=corr_per_chip; % correlators per chip
delta_t=delayRes; % correlator spacing
n_part = 100;%Number of particles
Rww_fil = 0 ; %Process noise
Rvv_fil = 2500*((N_corr-1)/20)*(n_part/250); %Measurement noise (21 corr - 2500, 41 corr - 5000, ...)
M_pf = 3; % number of paths
x_est_bpf = zeros(epochs,2*M_pf);
C_s = zeros(epochs,M_pf);
corr_outPF = zeros(1,N_corr)';




%%%%%%%%%%%%%% SIMULATION CONFIG %%%%%%%%%%%%%%
LOS_delay=0;
dynamic_LOS=0; % time-varying LOS
dynamic_multipath=1; % time-varying multipath
en_plots=1; % enable plots
% initial multipath free channel
path_delays=chan4.delays;
path_attenuation=chan4.attenuation;


rng(10);


for Index = 1: epochs
    
    %%%% varying multipath %%%%
    if dynamic_multipath
        if Index > 30
            path_delays=chan3.delays;
            path_attenuation=chan3.attenuation;
            
            %% testing
            %path_delays(2) = path_delays(2) + (Index-30)/500;
        end
    end
    
    %%%% varying LOS %%%%
    if dynamic_LOS
        if Index>50 && Index<300
            LOS_delay=LOS_delay+delayRes/4;
            path_delays=chan4.delays+LOS_delay; % all signal paths shifted
        end
    end

    %interference injection - input simulated signal
    INCode = 0;
    for idx = 1: size(path_delays, 2)
        INCode = INCode + subsampleSignalShifter(CodeVecSampled, path_delays(idx)*T_c, fs)./path_attenuation(idx);
    end

    %correlator engine
    for iCorr = 1:N_corr   
        code_replicas(iCorr,:) = circshift(CodeVecSampled, shiftidx(iCorr,1));         
        % loop closure
        if Index > 1
            code_replicas(iCorr,:) = subsampleSignalShifter(code_replicas(iCorr,:), -DLLdiscri(Index-1), fs);
        end
    end
    corr_out = code_replicas*INCode;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Initialization
     if (Index == 1)
        particle_pred = zeros(n_part, 2*M_pf);
        particle = zeros(n_part, 2*M_pf);
        particle(:,1) = -0.5 + rand(n_part, 1); % LOS ambiguity
        for idx=1:M_pf-1
            particle(:,2*idx+1) = particle(:,1) + rand(n_part, 1); % Multipath ambiguity
        end
        weight = ones(n_part,1)/n_part;
        phi_ss = signals' * signals;
    end

    for idx=1:n_part
        % Importance Sampling
        % Generate path delays according to the Importance Density Function
        if Index>1
            %particle(idx,1) = particle_pred(idx,1) + randn(1,1)*sqrt(C_s(Index-1, 1));
            pd = makedist('Normal','mu',particle_pred(idx,1),'sigma',sqrt(C_s(Index-1, 1)));
            t = truncate(pd,-0.5,0.5);
            particle(idx,1) = random(t,1,1);
            for subindex=1:M_pf-1
                pd = makedist('Normal','mu',particle_pred(idx,2*subindex+1),'sigma',sqrt(C_s(Index-1,subindex+1)));
                t = truncate(pd,particle(idx,1),particle(idx,1)+1);%random(t,1,1)
                particle(idx,2*subindex+1) = random(t,1,1);
                %particle(idx,2*subindex+1) = particle(idx,1) + abs((particle(idx,2*subindex+1)-particle(idx,1)) + random(t,1,1));%randn(1,1)*sqrt(C_s(Index-1,subindex+1)));
            end
        end

        % Estimate amplitudes according to path delays
        LS_H = phi_ss*corr_out;
        for subindex=0:M_pf-1
            particle(idx,2*subindex+2)=interp1(Spacing,LS_H,particle(idx,2*subindex+1));
            if(isnan(particle(idx,2*subindex+2)))
                particle(idx,2*subindex+2)=0;
            end
        end
        
        % Weight update
        % transform the particle into the measurements domain - correlation output
        INCodePF=0;
        for subindex=0:M_pf-1
            %CodePF = Spacing(P) + particle(idx,2*subindex+1) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(P) + particle(idx,2*subindex+1));
            %INCodePF = INCodePF + Code(ceil(CodePF) + 4).*particle(idx,2*subindex+2);
            INCodePF = INCodePF + subsampleSignalShifter(CodeVecSampled, particle(idx,2*subindex+1)*T_c, fs).*particle(idx,2*subindex+2);
        end


        for subindex = 1: N_corr
            time_stampsPF = (Spacing(subindex)) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(subindex));
            code_replicasPF = Code(ceil(time_stampsPF) + 4);
            corr_outPF(subindex) = sum(code_replicasPF.*INCodePF);
        end
        %corr_outPF=phi_ss*particle(idx,2:2:end);%%%%%%%%%%%%%%%%%%%
        %porquê estimar só numero de path em vez do canal todo?

        innov = norm(corr_out - corr_outPF);
        %weight(idx) = exp( -log(sqrt(2*pi*Rvv_fil)) -(( innov/10000 )^2)/(2*Rvv_fil) );
        weight(idx) = exp( -0.5 * ((innov / Rvv_fil)^2) );
    end
    % Weigth normalization
    weight = weight/sum(weight);
        
    % Estimation
    % Maximum a Posteriori (MAP) estimation
    [~, MAP] = max(weight);
    x_est_bpf(Index,:) = particle(MAP,:);
    
    % Update the covariance of the state error estimation
    for idx=1:n_part
        for subindex=0:M_pf-1
            C_s(Index,subindex+1) = C_s(Index,subindex+1) + weight(idx)*(   (particle(idx,subindex*2+1:subindex*2+2)-x_est_bpf(Index,subindex*2+1:subindex*2+2))*(particle(idx,subindex*2+1:subindex*2+2)-x_est_bpf(Index,subindex*2+1:subindex*2+2))');
        end
    end
    
    % Resampling
    % Initialize next iteration
    % Every particle equal to MAP
    for idx=1:n_part
        particle_pred(idx,:) = x_est_bpf(Index,:);
    end


    % PLOTS
    DLLdiscri(1,Index)=x_est_bpf(Index,1);
    %DLLdiscri(1,Index) = 0; % open loop

    PR_error(Index,1)=(LOS_delay*T_c+x_est_bpf(Index,1))*c;
    if (en_plots)
        if(Index == 1) 
            % plot initial correlation output
            response = figure(1);
            subplot(222);
            p0=plot(Spacing, corr_out./numSample,'-','Color',"#77AC30");
            %ylim([-0.1 1.5])
            xlim([-1.1 1.1])
            ylabel('Normalized Amplitude')
            xlabel('Delay [chips]')
            title('ACF')
            grid on;
            hold on;
            % plot initial position of central correlators and filter peak
            p1=plot(Spacing((VE : 1 : VL)),corr_out((VE : 1 : VL))/numSample,'*','Color',"#77AC30");
    
            % plot multipath information
            drawnow
            subplot(221);
            xlabel('Delay in chips')
            ylabel('Amplitude')
            title('Channel Impulse Response (CIR)')
            drawnow
    
            subplot(223);
            discri=plot((1 : 1 : Index), DLLdiscri(1,(1 : 1 : Index)),'-','Color',"#7E2F8E");
            ylabel('Seconds')
            xlabel('Epochs (ms)')
            title('LOS delay state estimate')
            subplot(224);
            pr=plot((1 : 1 : Index), PR_error((1 : 1 : Index),1),'-','Color',"#7E2F8E");
            ylabel('Meters')
            xlabel('Epochs (ms)')
            title('Pseudorange error')
        else
            %%%% REFRESH PLOTS %%%%
            set(p0,'XData',Spacing,'YData',corr_out./numSample);
            set(p1,'XData',Spacing((VE : 1 : VL)),'YData',corr_out(VE : 1 : VL)/numSample);
            drawnow
            set(discri,'XData',(1 : 1 : Index),'YData',DLLdiscri(1,(1 : 1 : Index)));
            drawnow
            set(pr,'XData',(1 : 1 : Index),'YData',PR_error((1 : 1 : Index),1));
            drawnow
            figure(response);
            subplot(221)
            hold on
            cla
            title('Channel Impulse Response (CIR)');
            stem(path_delays(1)-path_delays(1),1/path_attenuation(1),'b');
            if size(path_delays,2)>1
                stem(path_delays(2:end)-path_delays(1),1./path_attenuation(2:end),'r');
            end
            hold on
            stem((-L_h*delta_t:delta_t:L_h*delta_t),X_k(2:end),'--og');
            drawnow
            pause(0.005);
        end
    end
end

%{

%}
