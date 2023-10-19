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
Code = [Code(end-1) Code(end) Code Code(1) Code(2) Code(3) Code(4)];%enable replica slidding

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
s=1;
T_int=0.001; % integration time
T_c=0.001/codelength; % chip duration
L_h=corr_per_chip; % correlators per chip
L_c=L_h;
N_corr=total_corr;
N_meas=N_corr+1; % number of measurement variables
N_tap=N_corr; % total correlators
N_st=N_corr+1; % total states
delta_t=1/L_h; % correlator spacing
T_h=delta_t;
W_CIR=L_h*T_h; % CIR ring
W_bank=W_CIR;
X_k=zeros(N_st,1); % state vector
h_k=zeros(N_tap,1); % channel taps
A=eye(N_st); % process matrix
Q=eye(N_st); % process covariance matrix
cov_Q_t=0.02*T_c/(s^2);
cov_Q_h=5*(10^(-4));
v_k=zeros(N_st,1);% process noise
%v_K=mvnrnd(0,Q);
z_constr=0; % constraining measurement
Z_k=zeros(N_meas,1); % measurements
Z_k_predict=zeros(N_meas,1); % predicted measurements
R=zeros(N_meas,N_meas); % measurement noise covariance matrix
W=W_bank+delta_t;
alpha=0.9; % correlator output covariance scaling
R_wn_tilde_c=zeros(N_corr,N_corr);
w=zeros(N_corr,1); % weighting function
R_wn_c=eye(N_corr,N_corr);
constr_dev=0.001; % standard deviation of constraint
phi_ss = zeros(N_corr, N_corr); % correlation matrix
d_phi_ss = zeros(N_corr, N_corr); % derivative of correlation matrix
P_k=eye(N_st); % state estimate covariance
state_cov_t=0.01*T_c;
state_cov_d=0.01*T_c/s;
state_cov_h=0.5;
J_k=zeros(N_meas,N_st); % jacobian of measurement model
h_k=zeros(N_corr,1); % channel estimate
K_k=zeros(N_st,N_meas); % kalman gain
S_k=zeros(N_meas,N_meas); % residual covariance matrix
E_k=zeros(N_meas,1); % residual between measurement and predicted measurements
IDENT=eye(N_st); % identity matrix
tau_k=0; % LOS delay estimate



LOS_delay=0;
dynamic_LOS=0; % time-varying LOS
dynamic_multipath=1; % time-varying multipath
en_plots=1; % enable plots


rng(10);

for Index = 1: epochs
    
    %%%% varying multipath %%%%
    if Index < 30
        mltpth_delays=chan4.delays;
        mltpth_attenuation=chan4.attenuation;
    elseif dynamic_multipath
        if Index<60
            %LOS_delay=0.2;
            mltpth_delays=chan3.delays;
            mltpth_attenuation=chan3.attenuation;
        elseif Index<80
            %LOS_delay=0;
            mltpth_delays=chan2.delays;
            mltpth_attenuation=chan2.attenuation;
        elseif Index<100
            %LOS_delay=-0.2;
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
            %%%%%% INITIALIZATION %%%%%%
            if (Index == 1)
                
                for idx = 1 : N_corr
                    for idx1 = 1 : N_corr 
                        phi_ss(idx1,idx) = sum(code_replicas(idx,:).*code_replicas(idx1,:));
                    end
                end

                for idx = 1 : N_corr
                    for idx1 = 1 : N_corr
                        extended_code_replica=[code_replicas(idx,end) code_replicas(idx,:) code_replicas(idx,1)];
                        d_s1=(extended_code_replica(3:1:end)-extended_code_replica(1:1:end-2))/(2*(1/fs));
                        s2=code_replicas(idx1,:);
                        d_phi_ss(idx,idx1) = sum(d_s1 .* s2);
                    end
                end


                Q(1,1)=(cov_Q_t^2)*(T_int^4)/4; % equation 20
                Q(2:end,2:end)=(cov_Q_h^2).*Q(2:end,2:end);
                idx_aux=1;
                for idx=(-L_c*delta_t):(L_c*delta_t)
                    % f_Tukey window (equation 26)
                    if(abs(idx)<((1-alpha)*W))
                        w(idx_aux,1)=1;
                    elseif ((abs(idx)>=((1-alpha)*W)) && (abs(idx)<=W))
                        w(idx_aux,1)=0.5+0.5*cos((pi/alpha)*((abs(idx)/W)+alpha-1));
                    else
                        w(idx_aux,1)=0;
                    end
                    w(idx_aux,1)=w(idx_aux,1)^(-1); % equation 25
                    idx_aux=idx_aux+1;
                end
                R_wn_tilde_c=w*w'.*phi_ss; % equation 24
                R_wn_c=R_wn_tilde_c./2; % equation 28
                R(1:N_meas-1,1:N_meas-1)=R_wn_c; % equation 27
                R(N_meas,N_meas)=constr_dev^2;

                X_k(1,1)=0; % initial delay 0
                X_k(1+L_h+1,1)=1; % initial multipath free scenario
                
                P_k(1,1)=state_cov_t^2;
                P_k(2:end,2:end)=P_k(2:end,2:end).*state_cov_h; % equation 32
                

                J_k(1:end-1,2:end)=phi_ss; % appendix A2
            end


            %%%%%% PREDICTION %%%%%%
            X_k = A*X_k; % equation 29
            P_k = A*P_k*A' + Q; % equation 30


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
            J_k(1:end-1,1)=d_phi_ss*X_k(2:end,1); % appendix A1       
            J_k(N_meas,1+L_h+1)=-(2/(2*L_h-1))*(X_k(1+L_h+1,1)/(abs(X_k(1+L_h+1,1))^4))*z_constr_aux; % equation A3-a
            for idx=1:N_tap
                if(idx~=(L_h+1))
                    J_k(N_meas,1+idx)=(2/(2*L_h-1))*(X_k(1+idx,1)/(abs(X_k(1+L_h+1,1))^2)); % equation A3-b
                end
            end
            J_k(N_meas,1)=0;
            %J_k(N_meas,:)=zeros(1,N_st);

            % kalman equations
            S_k = J_k*P_k*J_k' + R; % equation 34-b
            K_k = P_k*J_k'*pinv(S_k); % equation 33-b
            X_k = X_k+K_k*E_k; % equation 33-a
            P_k = (IDENT-K_k*J_k)*P_k; % equation 34-a

            t_los = X_k(1,1);
        end
    end

    % DLL discriminator and loop filter
    if new_tracking
        if enable_EKF
            DLLdiscri(1,Index) = t_los;
        else
            DLLdiscri(1,Index) = Spacing(P) - Spacing(t_los);
        end
        code_output = (0.3749245/0.007030542258775)*2.5*DLLdiscri(1,Index)*1000000; % this can be wrong!!!!!!!
        %code_output = code_outputLast + (0.3749245/0.000007030542258775)*(DLLdiscri(1,Index) - DLLdiscriLast) + DLLdiscri(1,Index)* (0.001/0.000007030542258775);
        DLLdiscriLast = DLLdiscri(1,Index);
        code_outputLast = code_output;
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
                stem((-L_h*delta_t:delta_t:L_h*delta_t),X_k(2:end),'--og');
            end
            drawnow
            pause(0.005);
        end
    end
end

%{

%}