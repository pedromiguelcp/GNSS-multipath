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

Spacing = (-1:0.05:1);
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

svlength=1;
svindex=1;

%%%%% SELECT THE ALGORITHM %%%%%
enable_LS = 0;
enable_LMS = 0;
enable_NLMS = 0;
enable_EKF = 1;
if enable_LS || enable_LMS || enable_NLMS || enable_EKF
    new_tracking = 1;
else
    new_tracking = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMMON %%
G = zeros(total_corr, total_corr);
G_inv = zeros(total_corr, total_corr);
g_los = zeros(total_corr, 1);
a_los_aux = zeros(total_corr, 1);
a_los = 0;
t_los = 0;
y_unfiltered = zeros(total_corr, 1);
y_filtered = zeros(total_corr, svlength);
y_desired = zeros(total_corr, svlength);
desired = zeros(total_corr, 1);
C_sv=zeros(total_corr,total_corr,svlength);
gradient=zeros(total_corr,total_corr,svlength);
for subindex=1:svlength
    C_sv(:,:,subindex)=eye(total_corr);
    gradient(:,:,subindex)=eye(total_corr).*0;
end
error=0;
weights=eye(total_corr);
weights_aux=eye(total_corr);
%% LMS %%
lms_lr=.0000000000001;
LMS_iterations=5;
autocorr_matrix = zeros(total_corr,total_corr);
corr_matrix = zeros(total_corr,total_corr);
%% NLMS %%
nlms_lr=.001;
gamma = 1e-12; % small positive constant to avoid singularity
NLMS_iterations=1;
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
B=eye(L_ekf+2)*0.0001; % process noise covariance
B(1,1)=0.1;
drv_fi_ss=zeros(P_ekf,L_ekf+1);
drv_fi_ss1=zeros(P_ekf,L_ekf+1);
drv_fi_ss2=zeros(P_ekf,L_ekf+1);


LOS_delay=0;
dynamic_LOS=0; % time-varying LOS
dynamic_multipath=1; % time-varying multipath
en_plots=1;

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
        elseif Index<150
            mltpth_delays=chan2.delays;
            mltpth_attenuation=chan2.attenuation;
        end
    end
    
    %%%% varying LOS %%%%
    if dynamic_LOS
        if Index<200
            LOS_delay=LOS_delay-0.001; % simulate positive doppler
        else
            LOS_delay=LOS_delay-0.000; % simulate negative doppler
        end
    end
    mltpth_delays=mltpth_delays+LOS_delay; % multipaths are relative to the LOS

    %numSample = round((codelength-remChip)/(codeFreq/fs)); 
    
    LOS_Code = Spacing(P)+LOS_delay : codeFreqBasis/fs : ((numSample -1) * (codeFreqBasis/fs) + Spacing(P)+LOS_delay);
    INCode   = Code(ceil(LOS_Code) + 2);

    %interference injection
    for subindex = 1: size(mltpth_delays, 2)
        multipath = Spacing(P) + mltpth_delays(subindex) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(P) + mltpth_delays(subindex));
        INCode = INCode + Code(ceil(multipath) + 2)./mltpth_attenuation(subindex);
    end

    %correlator engine
    for subindex = 1: total_corr
        time_stamps(subindex,:) = (Spacing(subindex) + remChip) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(subindex) + remChip);
    end

    for subindex = 1: total_corr
        code_replicas(subindex,:) = Code(ceil(time_stamps(subindex,:)) + 2);
    end
    
    for subindex = 1: total_corr
        corr_out(subindex) = sum(code_replicas(subindex,:)     .*INCode);
    end
    %corr_out=corr_out+std(corr_out)*randn(1,21)/10;

    remChip   = (time_stamps(P,numSample) + codeFreq/fs) - codeFreqBasis*0.001;



    if new_tracking
        %%% INITIAL STATE %%%
        if (Index == 1)% same for all sv
            LS_G = zeros(total_corr, total_corr);
            for subindex = 1 : total_corr
                for subindex1 = 1 : total_corr 
                    G(subindex1,subindex) = sum(code_replicas(subindex,:).*code_replicas(subindex1,:));
                    LS_G(subindex1,subindex) = sum(code_replicas(subindex,:).*code_replicas(subindex1,:));
                end
            end
            G_inv = pinv(G);
            LS_S=pinv(LS_G'*LS_G)*LS_G';
        end 
        y_unfiltered = corr_out';
        weights=C_sv(:,:,svindex);

        %%%%%%%%% LOS SIGNAL ESTIMATION %%%%%%%%%
        for subindex = 1 : total_corr
            for subindex1 = 1 : total_corr 
                g_los(subindex1,1) = sum(code_replicas(subindex,:).*code_replicas(subindex1,:));
            end
            a_los_aux(subindex,1)=(g_los'*G_inv*weights*y_unfiltered')/(g_los'*G_inv*g_los);        
        end
        [a_los, t_los] = max(a_los_aux(:,1));
        for subindex = 1 : total_corr 
            y_desired(subindex,svindex) = a_los*sum(code_replicas(subindex,:).*code_replicas(t_los,:));
        end
        desired=y_desired(:,svindex);


        
        %%%%%%%%% ADAPTIVE CHANNEL COMPENSATION %%%%%%%%%
        % LS
        if enable_LS
            LS_H = LS_S*y_unfiltered';
            [a_los, t_los] = max(LS_H);

        % LMS
        elseif enable_LMS 
            for i=1:LMS_iterations
                error=(desired'-y_unfiltered*weights);
                weights = weights + lms_lr*y_unfiltered'*error;
            end

        % NLMS
        elseif enable_NLMS 
            for i=1:NLMS_iterations
                error=(desired'-y_unfiltered*weights);
                weights = weights + (nlms_lr/(gamma+(abs(y_unfiltered)*abs(y_unfiltered)')))*y_unfiltered'*error;
            end  

        % EKF    
        elseif enable_EKF 
            if (Index == 1)
                X_k=zeros(L_ekf+2,1);
                X_k(2,1)=1;
                fi_ss=LS_G(:,L_ekf+1:P_ekf);
                F_k(:,2:L_ekf+2)=fi_ss;

                drv_x=0:(1/L_ekf):(2);
                for drv_idx=1:L_ekf+1
                    drv_fi_ss1(1:P_ekf-1,drv_idx)=diff(fi_ss(:,drv_idx)')./diff(drv_x);
                end
               
                drv_fi_ss1(P_ekf,:)=drv_fi_ss1(P_ekf-1,:);
                drv_fi_ss1(P_ekf,L_ekf+1)=drv_fi_ss1(P_ekf,L_ekf);

                drv_fi_ss1(L_ekf+1,1)=0;
                drv_fi_ss1(L_ekf+2,2)=0;
                drv_fi_ss1(L_ekf+3,3)=0;
                drv_fi_ss1(L_ekf+4,4)=0;
                drv_fi_ss1(L_ekf+5,5)=0;
                drv_fi_ss1(L_ekf+6,6)=0;
                drv_fi_ss1(L_ekf+7,7)=0;
                drv_fi_ss1(L_ekf+8,8)=0;
                drv_fi_ss1(L_ekf+9,9)=0;
                drv_fi_ss1(L_ekf+10,10)=0;
                drv_fi_ss1(L_ekf+11,11)=0;
                for drv_idx=1:L_ekf+1
                    aux=drv_fi_ss1(drv_idx,:);
                    drv_fi_ss1(drv_idx,:)=drv_fi_ss1(P_ekf-drv_idx+1,:);
                    drv_fi_ss1(P_ekf-drv_idx+1,:)=aux;
                end

                drv_x=0:(1/L_ekf):(1);
                for drv_idx=1:P_ekf
                    drv_fi_ss(drv_idx,1:L_ekf)=diff(fi_ss(drv_idx,:))./diff(drv_x);
                end
                drv_fi_ss(:,L_ekf+1)=drv_fi_ss(:,L_ekf);
                drv_fi_ss(L_ekf+1,L_ekf+1)=drv_fi_ss(L_ekf,L_ekf+1);
                %drv_fi_ss(L_ekf+1,1)=0;
                %drv_fi_ss(L_ekf+2,2)=0;
                %drv_fi_ss(L_ekf+3,3)=0;
                %drv_fi_ss(L_ekf+4,4)=0;
                %drv_fi_ss(L_ekf+5,5)=0;
                %drv_fi_ss(L_ekf+6,6)=0;
                %drv_fi_ss(L_ekf+7,7)=0;
                %drv_fi_ss(L_ekf+8,8)=0;
                %drv_fi_ss(L_ekf+9,9)=0;
                %drv_fi_ss(L_ekf+10,10)=0;
                %drv_fi_ss(L_ekf+11,11)=0;
                %for drv_idx=1:P_ekf
                %    for drv_idx1=1:L_ekf+1
                %        if abs(drv_fi_ss(drv_idx,drv_idx1))<500
                %            drv_fi_ss(drv_idx,drv_idx1)=-26000;
                %        end
                %        if drv_fi_ss(drv_idx,drv_idx1)<0
                %            drv_fi_ss(drv_idx,drv_idx1)=-26000;
                %        end
                %        if drv_fi_ss(drv_idx,drv_idx1)>0
                %            drv_fi_ss(drv_idx,drv_idx1)=26000;
                %        end
                %    end
                %end
                drv_fi_ss2(:,1)=(LS_G(:,L_ekf+1)-LS_G(:,L_ekf))./(1/L_ekf);
                for drv_idx=2:L_ekf+1
                    drv_fi_ss2(:,drv_idx)=(fi_ss(:,drv_idx)-fi_ss(:,1))./((drv_idx-1)*(1/L_ekf));
                end
            end
            % prediction
            X_k=A*X_k;
            P_k=A*P_k*A'+B;

            % update
            Z_k = corr_out';
            Z_K_predict = X_k(2:L_ekf+2,1)'*fi_ss';
            E_k = Z_k - Z_K_predict;

            F_k(:,1) = (X_k(2:L_ekf+2,1)'*drv_fi_ss2')';

            S_k = F_k*P_k*F_k'+ LS_G;%
            K_k = P_k*F_k'*pinv(S_k);
            X_k = X_k+K_k*E_k';
            P_k = (IDENT-K_k*F_k)*P_k;
            t_los=X_k(1,1);
        end



        C_sv(:,:,svindex)=weights; % update weights
    end





    % DLL discriminator and loop filter
    if new_tracking
        if enable_EKF
            DLLdiscri(1,Index) = -t_los;
            code_output= (0.3749245/0.007030542258775)*2.5*DLLdiscri(1,Index);
        end
    else
        DLL_E           = sqrt(corr_out(E)^2);
        DLL_L           = sqrt(corr_out(L)^2);
        DLLdiscri(1,Index)       = 0.5 * (DLL_E-DLL_L)/(DLL_E+DLL_L);
        code_output     = code_outputLast + (0.3749245/0.007030542258775)*(DLLdiscri(1,Index) - DLLdiscriLast) + DLLdiscri(1,Index)* (0.001/0.007030542258775);
        DLLdiscriLast   = DLLdiscri(1,Index);
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
            %p3=plot(time_stamps(:, 1), y_desired(:,svindex)/numSample,'-','Color',"#D95319"); % tirei d coloquei y
    
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
            %set(p3,'XData',time_stamps(:, 1),'YData',y_desired(:,svindex)/numSample); 
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