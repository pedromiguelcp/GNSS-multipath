clear all;
close all;
clc;
% Define the parameters for the simulation
c = 299792458; % m/s
fs = 26000000; % Sample rate (Hz)
codeFreqBasis = 1.023e6; % C/A code rate (chips/s)
codeFreq = codeFreqBasis;
codelength = 1023;
numSample = round(codelength/(codeFreq/fs));
remChip = 0;
code_outputLast = 0;
DLLdiscriLast = 0;
cross_corr = 0;

Spacing = (-1:0.1:1); % change here the number of correlators (0.1=21correlators, 0.05=41correlators ...)
total_corr = size(Spacing,2);
corr_per_chip = floor(total_corr/2);
time_stamps = zeros(corr_per_chip*2+1,numSample);
code_replicas = zeros(corr_per_chip*2+1,numSample);
corr_out = zeros(1,corr_per_chip*2+1)';
VE = corr_per_chip-1;
E = corr_per_chip;
P = corr_per_chip+1;
L = corr_per_chip+2;
VL = corr_per_chip+3;

% Specify here the multipath delay (chips) and attenuation (linear factor), compared to the LOS signal
chan1.delays=[0.1 0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.9];
chan1.attenuation=[2 2.5 3 2 1.5 1.6 1.7 4 3];
%Channel 2
chan2.delays=[0.1 0.25 0.35 0.45 0.55];
chan2.attenuation=[1.8 2.2 3 3.5 4];
% Channel 3
chan3.delays=[0.1 0.2 0.3];
chan3.attenuation=[2 4 8];
% Channel 4
chan4.delays=[];
chan4.attenuation=[];


epochs = 5000; %signal tracking epochs

% Generate the C/A code
Code = generateCAcode(1); % Generate C/A code for satellite 1
Code = [Code(end-1) Code(end) Code Code(1) Code(2) Code(3)];%enable replica slidding

DLLdiscri = zeros(1,epochs);

%%%%% SELECT THE ALGORITHM %%%%%
enable_LS = 0;
enable_EKF = 1;
if enable_LS || enable_EKF
    new_tracking = 1;
else
    new_tracking = 0; % the traditional DLL will be used
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMMON %%
LS_G = zeros(total_corr, total_corr);
a_los = 0;
t_los = 0;
%% EKF %%
L_ekf=corr_per_chip;
P_ekf=2*L_ekf+1; % total correlators
X_k=zeros(L_ekf+2,1); % state vector
Z_k=zeros(P_ekf,1); % measurements
Z_K_predict=zeros(P_ekf,1); % predicted measurements
tau_k=0; % tau LOS estimate
h_k=zeros(L_ekf+1,1); % channel estimate
A=eye(L_ekf+2); % state transition matrix
P_k=eye(L_ekf+2); % state estimate covariance
K_k=zeros(L_ekf+2,P_ekf); % kalman gain
F_k=zeros(P_ekf,L_ekf+2); % linearization matrix
S_k=zeros(P_ekf,P_ekf); % residual covariance matrix
R_k=zeros(P_ekf,P_ekf); % measurement noise covariance matrix
E_k=zeros(P_ekf,1); % residual between measurement and predicted measurements
IDENT=eye(L_ekf+2); % identity matrix
B=eye(L_ekf+2)*0.0001; % process noise covariance matrix
B(1,1)=0.000001;
drv_fi_ss=zeros(P_ekf,L_ekf+1); % derivative

LOS_delay=0;
dynamic_LOS=0; % time-varying LOS
dynamic_multipath=1; % time-varying multipath
en_plots=1; % enable plots

for Index = 1: epochs
    
    %%%% varying multipath %%%%
    if Index < 30
        mltpth_delays=chan4.delays;
        mltpth_attenuation=chan4.attenuation;
    elseif dynamic_multipath
        if Index<60
            %LOS_delay=0.1;
            mltpth_delays=chan3.delays;
            mltpth_attenuation=chan3.attenuation;
        elseif Index<80
            mltpth_delays=chan2.delays;
            mltpth_attenuation=chan2.attenuation;
        elseif Index<100
            mltpth_delays=chan1.delays;
            mltpth_attenuation=chan1.attenuation;
        end
    end
    
    %%%% varying LOS %%%%
    if dynamic_LOS
        if Index<200
            LOS_delay=LOS_delay+0.001; % simulate negative doppler
        else
            LOS_delay=LOS_delay-0.000;
        end
    end
    mltpth_delays=mltpth_delays+LOS_delay; % multipaths are relative to the LOS

    %numSample = round((codelength-remChip)/(codeFreq/fs)); 
    
    LOS_Code = Spacing(P)+LOS_delay : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(P)+LOS_delay);
    INCode   = Code(ceil(LOS_Code) + 2);

    %interference injection
    for subindex = 1: size(mltpth_delays, 2)
        multipath = Spacing(P) + mltpth_delays(subindex) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(P) + mltpth_delays(subindex));
        INCode = INCode + Code(ceil(multipath) + 2)./mltpth_attenuation(subindex);
    end

    %correlator engine
    for subindex = 1: total_corr
        time_stamps(subindex,:) = (Spacing(subindex) + remChip) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(subindex) + remChip);
        code_replicas(subindex,:) = Code(ceil(time_stamps(subindex,:)) + 2);
        corr_out(subindex) = sum(code_replicas(subindex,:)     .*INCode);
    end

    remChip   = (time_stamps(P,numSample) + codeFreq/fs) - codeFreqBasis*0.001;



    if new_tracking
        if (Index == 1) % compute correlation matrix
            for subindex = 1 : total_corr
                for subindex1 = 1 : total_corr 
                    LS_G(subindex1,subindex) = sum(code_replicas(subindex,:).*code_replicas(subindex1,:));
                end
            end
        end 

        %%%%%%%%% ADAPTIVE CHANNEL COMPENSATION %%%%%%%%%
        % LS
        if enable_LS
            if (Index==1)
                LS_S=pinv(LS_G'*LS_G)*LS_G';
            end
            LS_H = LS_S*corr_out;
            [a_los, t_los] = max(LS_H);

        % EKF    
        elseif enable_EKF 
            % Initialization
            if (Index == 1)
                X_k = zeros(L_ekf+2,1);
                X_k(2,1) = 1; % only LOS
                fi_ss = LS_G(:,L_ekf+1:P_ekf);
                F_k(:,2:L_ekf+2) = fi_ss; % equation 46

                drv_fi_ss(:,1)=(LS_G(:,L_ekf+1)-LS_G(:,L_ekf))/(1/L_ekf);
                for drv_idx=2:L_ekf+1 % this should be wrong!!!!!!!
                    drv_fi_ss(:,drv_idx)=(fi_ss(:,drv_idx)-LS_G(:,L_ekf))./((drv_idx-1)*(1/L_ekf));
                end
            end
            % prediction
            X_k = A*X_k; % equation 33
            P_k = A*P_k*A'+B; % equation 34

            % update
            Z_k = corr_out';
            Z_K_predict = X_k(2:L_ekf+2,1)'*fi_ss'; % equation 39
            E_k = Z_k - Z_K_predict; % equation 36

            F_k(:,1) = (X_k(2:L_ekf+2,1)'*drv_fi_ss')'; % equation 45
            
            h_r = 0;
            for idx=1:L_ekf
                h_r = h_r + X_k(idx+2)*X_k(idx+2);
            end
            S_k = F_k*P_k*F_k' + LS_G.*(1+h_r); % equation 43
            K_k = P_k*F_k'*pinv(S_k); % equation 42
            X_k = X_k+K_k*E_k'; % equation 40
            P_k = (IDENT-K_k*F_k)*P_k; % equation 41
            t_los = X_k(1,1);
        end
    end

    % DLL discriminator and loop filter
    if new_tracking
        if enable_EKF
            DLLdiscri(1,Index) = -t_los;
        else
            DLLdiscri(1,Index) = Spacing(P) - Spacing(t_los);
        end
        code_output = (0.3749245/0.007030542258775)*2.5*DLLdiscri(1,Index); % this can be wrong!!!!!!!
        code_output = 0; %uncomment to disable DLL corrections
    else % this is the traditional DLL loop filter and nco
        DLL_E = sqrt(corr_out(E)^2);
        DLL_L = sqrt(corr_out(L)^2);
        DLLdiscri(1,Index) = 0.5 * (DLL_E-DLL_L)/(DLL_E+DLL_L);
        code_output = code_outputLast + (0.3749245/0.007030542258775)*(DLLdiscri(1,Index) - DLLdiscriLast) + DLLdiscri(1,Index)* (0.001/0.007030542258775);
        DLLdiscriLast = DLLdiscri(1,Index);
        code_outputLast = code_output;
    end

    codeFreq = codeFreqBasis - code_output;
    PR_error(Index,1)=(time_stamps(P)-LOS_delay)*(0.001/1023)*c;%compute pseudorange error - 1 PRN period is 1ms (1023 chips)

    % PLOTS
    if (en_plots)
        if(Index == 1) 
            % plot initial correlation output
            response = figure(1);
            subplot(222);
            p0=plot(Spacing(2 : 1 : corr_per_chip*2), corr_out(2 : 1 : corr_per_chip*2)./numSample,'-','Color',"#77AC30");
            ylim([-0.1 5])
            xlim([-1.5 1.5])
            ylabel('Normalized Amplitude')
            xlabel('Delay [chips]')
            title('ACF')
            hold on;
            % plot initial position of central correlators and filter peak
            p1=plot(Spacing((VE : 1 : VL)),corr_out((VE : 1 : VL))/numSample,'*','Color',"#77AC30");
            p2=plot([0 0],[0 1],'-*','Color',"#D95319");
            % plot initial desired response
    
            % plot multipath information
            drawnow
            subplot(221);
            xlabel('Delay in chips')
            ylabel('Amplitude')
            title(['Channel Impulse Response (CIR)'])
            drawnow
    
            subplot(223);
            discri=plot((1 : 1 : Index), DLLdiscri(1,(1 : 1 : Index)),'-','Color',"#7E2F8E");
            ylabel('DLL discriminator')
            xlabel('Epochs (ms)')
            title('DLL corrections')
            subplot(224);
            pr=plot((1 : 1 : Index), PR_error((1 : 1 : Index),1),'-','Color',"#7E2F8E");
            ylabel('Meters')
            xlabel('Epochs (ms)')
            title('Pseudorange error')
        else
            %%%% REFRESH PLOTS %%%%
            set(p0,'XData',time_stamps(:, 1),'YData',corr_out./numSample);
            set(p1,'XData',time_stamps((VE : 1 : VL), 1),'YData',corr_out(VE : 1 : VL)/numSample);
            set(p2,'XData',[t_los t_los],'YData',[0 1]);
            drawnow
            set(discri,'XData',(1 : 1 : Index),'YData',DLLdiscri(1,(1 : 1 : Index)));
            drawnow
            set(pr,'XData',(1 : 1 : Index),'YData',PR_error((1 : 1 : Index),1));
            drawnow
            figure(response);
            subplot(221)
            hold on
            cla
            title(['Channel Impulse Response (CIR)']);
            stem(LOS_delay,1,'b');
            stem(mltpth_delays+LOS_delay,1./mltpth_attenuation,'r');
            hold on
            if enable_EKF
                stem((0:1/L_ekf:1),X_k(2:end),'--og');
            end
            drawnow
            pause(0.005);
        end
    end
end

%{

%}