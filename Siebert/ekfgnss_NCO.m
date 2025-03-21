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
d_CodeVecSampled=(extended_CodeVec(3:1:end)-extended_CodeVec(1:1:end-2))/(2*(1/fs));
for iCorr = 1:N_corr
    shiftidx(iCorr,1)=(1/((ceil(N_corr/2)-1)/ceil(fs/codeFreqBasis)))*(iCorr-1)-ceil(fs/codeFreqBasis);
    signals(:, iCorr) = circshift(CodeVecSampled', shiftidx(iCorr,1));
    signals_deriv(:, iCorr) = circshift(d_CodeVecSampled', shiftidx(iCorr,1));
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
 chan3.attenuation=[1 4 8];
%chan3.delays=[0 5*delayRes];
%chan3.attenuation=[1 2];
% Channel 4
chan4.delays=0;
chan4.attenuation=1;


epochs = 5000; %signal tracking epochs

DLLdiscri = zeros(1,epochs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EKF %%
T_c=0.001/codelength; % chip duration
L_h=corr_per_chip; % correlators per chip
N_meas=N_corr+1; % number of measurement variables
N_tap=N_corr; % total correlators
N_st=N_corr+1; % total states
delta_t=delayRes; % correlator spacing
A=eye(N_st); % process matrix

R=zeros(N_meas,N_meas); % measurement noise covariance matrix
Q=eye(N_st); % process covariance matrix
std_dev_Q_t=0.04*T_c;
std_dev_Q_h=0.0025;
std_dev_P_t=0.1*T_c;
std_dev_P_h=0.025;
std_dev_constr=0.001; % constraint variance

X_k=zeros(N_st,1); % state vector
P_k=eye(N_st); % state estimate covariance
Z_k=zeros(N_meas,1); % measurements
Z_k_predict=zeros(N_meas,1); % predicted measurements
J_k=zeros(N_meas,N_st); % jacobian of measurement model
h_k=zeros(N_corr,1); % channel estimate
K_k=zeros(N_st,N_meas); % kalman gain
S_k=zeros(N_meas,N_meas); % residual covariance matrix
E_k=zeros(N_meas,1); % residual between measurement and predicted measurements

z_constr=0; % constraining measurement
phi_ss = zeros(N_corr, N_corr); % correlation matrix
d_phi_ss = zeros(N_corr, N_corr); % derivative of correlation matrix
IDENT=eye(N_st); % identity matrix
t_los = 0; % LOS delay estimate

%%%%%%%%%%%%%% SIMULATION CONFIG %%%%%%%%%%%%%%
LOS_delay=0;
dynamic_LOS=1; % time-varying LOS
dynamic_multipath=0; % time-varying multipath
en_plots=1; % enable plots
% initial multipath free channel
path_delays=chan4.delays;
path_attenuation=chan4.attenuation;


rng(10);

enable_shifter = true

for Index = 1: epochs
    
    %%%% varying multipath %%%%
    if dynamic_multipath
        if Index > 30
            path_delays=chan3.delays;
            path_attenuation=chan3.attenuation;
            
            %% testing
            path_delays(2) = path_delays(2) + (Index-30)/500*enable_shifter;
        end
    end
    
    %%%% varying LOS %%%%
    if dynamic_LOS
        if Index>50 && Index<300
            LOS_delay=delayRes;
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EKF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% INITIALIZATION %%%%%%
    if (Index == 1)
        phi_ss = signals' * signals;
        d_phi_ss = signals' * signals_deriv;

        Q(1,1)=std_dev_Q_t^2;
        Q(2:end,2:end)=Q(2:end,2:end).*(std_dev_Q_h^2);

        R(1:N_meas-1,1:N_meas-1)=phi_ss.*0.001;
        R(N_meas,N_meas)=std_dev_constr^2;

        X_k(1,1)=0; % initial delay 0
        X_k(1+L_h+1,1)=1; % initial multipath free scenario (only LOS pulse)
        
        P_k(1,1)=std_dev_P_t^2;
        P_k(2:end,2:end)=P_k(2:end,2:end).*(std_dev_P_h^2);

        J_k(1:end-1,2:end)=phi_ss;
    end


    %%%%%% PREDICTION %%%%%%
    X_k = A*X_k;
    P_k = A*P_k*A' + Q;

    %%%%%% UPDATE %%%%%%
    % measurements
    z_constr=(1/(2*L_h-1))*(abs(X_k(1+L_h+1,1))^(-2));
    z_constr_aux=0; 
    for idx=1:N_tap
        if(idx~=(L_h+1))
            z_constr_aux=z_constr_aux+abs(X_k(idx+1,1))^2;
        end
    end
    z_constr=z_constr*z_constr_aux;

    Z_k(1:end-1,1) = corr_out'; % actual correlation output
    Z_k_predict(1:end-1,1) = phi_ss*X_k(2:end,1); % predicted correlation output 
    Z_k(N_meas,1) = 0; % constraining measurement
    Z_k_predict(N_meas,1) = z_constr; % predicted constraining measurement
    E_k = Z_k - Z_k_predict;

    % jacobian 
    J_k(1:end-1,1)=d_phi_ss*X_k(2:end,1);    
    J_k(N_meas,1+L_h+1)=-(2/(2*L_h-1))*(X_k(1+L_h+1,1)/(abs(X_k(1+L_h+1,1))^4))*z_constr_aux;
    for idx=1:N_tap
        if(idx~=(L_h+1))
            J_k(N_meas,1+idx)=(2/(2*L_h-1))*(X_k(1+idx,1)/(abs(X_k(1+L_h+1,1))^2));
        end
    end
    J_k(N_meas,1)=0;
    %J_k(N_meas,:)=zeros(1,N_st); % uncomment to remove constraining measurement

    % kalman equations
    S_k = J_k*P_k*J_k' + R;
    K_k = P_k*J_k'*pinv(S_k);
    X_k = X_k+K_k*E_k;
    P_k = (IDENT-K_k*J_k)*P_k;

    t_los = X_k(1,1); % feedback for local replicas
    %t_los = 0; % open loop


    % PLOTS
    DLLdiscri(1,Index)=X_k(1,1);
    %PR_error(Index,1)=(ceil(t_los/delayRes)-LOS_delay)*(0.001/1023)*c;%compute pseudorange error - 1 PRN period is 1ms (1023 chips)
    PR_error(Index,1)=(LOS_delay*T_c+t_los)*c;
    if (en_plots)
        if(Index == 1) 
            % plot initial correlation output
            response = figure(1);
            subplot(222);
            p0=plot(Spacing, corr_out./numSample,'-','Color',"#77AC30");
            ylim([-0.1 1.5])
            xlim([-1.1 1.1])
            ylabel('Normalized Amplitude')
            xlabel('Delay [chips]')
            title('ACF')
            grid on;
            hold on;
            % plot initial position of central correlators and filter peak
            p1=plot(Spacing((VE : 1 : VL)),corr_out((VE : 1 : VL))/numSample,'*','Color',"#77AC30");
            % plot initial desired response
            % plot predicted measurements
            p2=plot(Spacing,Z_k_predict(1:end-1,1)/numSample);
            p3=plot(Spacing,(Z_k_predict(1:end-1,1)-corr_out)/numSample);
    
            % plot multipath information
            drawnow
            subplot(221);
            xlabel('Delay in chips')
            ylabel('Amplitude')
            title('Channel Impulse Response (CIR)')
            drawnow
    
            subplot(223);
            discri=plot((1 : 1 : Index), DLLdiscri(1,(1 : 1 : Index)),'-','Color',"#7E2F8E");
            ylabel('Chips')
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
            set(p2,'YData',Z_k_predict(1:end-1,1)/numSample);
            set(p3,'YData',(Z_k_predict(1:end-1,1)-corr_out)/numSample);
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
