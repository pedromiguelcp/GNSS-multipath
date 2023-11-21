clear all;
close all;
clc;
SimpleKalman(GetVolt());
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

Spacing = (-1:0.1:1);
total_corr = size(Spacing,2);
corr_per_chip = floor(total_corr/2);
time_stamps = zeros(total_corr,numSample);
code_replicas = zeros(total_corr,numSample);
corr_out = zeros(1,total_corr)';
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
chan4_1.delays=[0.2 0.4];
chan4_1.attenuation=[2 3];


epochs = 5000; %signal tracking epochs

% Generate the C/A code
Code = generateCAcode(1); % Generate C/A code for satellite 1
Code = [Code(end-3) Code(end-2) Code(end-1) Code(end) Code Code(1) Code(2) Code(3) Code(4) Code(5) Code(6)];%enable replica slidding

DLLdiscri = zeros(1,epochs);

svlength=1;
svindex=1;

%%%%% SELECT THE ALGORITHM %%%%%
enable_LS = 0;
enable_LMS = 0;
enable_NLMS = 0;
enable_PF = 0;
enable_RBPF = 1;
if enable_LS || enable_LMS || enable_NLMS || enable_EKF || enable_PF || enable_RBPF
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
%% PF %%
n_part = 100;%Number of particles
n_iter = 10;
Rww_fil = 0 ; %Process noise
Rvv_fil = 2500*((total_corr-1)/20)*(n_part/250); %Measurement noise (21 corr - 2500, 41 corr - 5000, ...)
M_pf = 1; % number of paths
x_est_bpf = zeros(epochs,2*M_pf);
C_s = zeros(epochs,M_pf);
corr_outPF = zeros(1,total_corr)';


LOS_delay=0;
dynamic_LOS=0; % time-varying LOS
dynamic_multipath=0; % time-varying multipath
en_plots=1;

for Index = 1: epochs
    
    %%%% varying multipath %%%%
    if Index < 30
        mltpth_delays=chan4.delays;
        mltpth_attenuation=chan4.attenuation;
    elseif dynamic_multipath
        if Index<40
            LOS_delay=0.1;
            %mltpth_delays=chan4_1.delays;
            %mltpth_attenuation=chan4_1.attenuation;
        elseif Index<50
            LOS_delay=0.2;
            %mltpth_delays=chan4.delays;
            %mltpth_attenuation=chan4.attenuation;
        elseif Index<60
            LOS_delay=0.3;
            %mltpth_delays=chan4.delays;
            %mltpth_attenuation=chan4.attenuation;
        elseif Index<70
            LOS_delay=0.4;
            %mltpth_delays=chan4.delays;
            %mltpth_attenuation=chan4.attenuation;
        end
    end
    
    %%%% varying LOS %%%%
    if dynamic_LOS
        if Index<200
            LOS_delay=LOS_delay+0.001; % simulate positive doppler
        elseif Index<400
            LOS_delay=LOS_delay-0.001; % simulate negative doppler
        else
            LOS_delay=LOS_delay;
        end
    end
    mltpth_delays=mltpth_delays+LOS_delay; % multipaths are relative to the LOS

    %numSample = round((codelength-remChip)/(codeFreq/fs)); 
    
    LOS_Code = Spacing(P)+LOS_delay : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(P)+LOS_delay);
    INCode   = Code(ceil(LOS_Code) + 4);

    %interference injection
    for subindex = 1: size(mltpth_delays, 2)
        multipath = Spacing(P) + mltpth_delays(subindex) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(P) + mltpth_delays(subindex));
        INCode = INCode + Code(ceil(multipath) + 4)./mltpth_attenuation(subindex);
    end

    %correlator engine
    for subindex = 1: total_corr
        time_stamps(subindex,:) = (Spacing(subindex) + remChip) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(subindex) + remChip);
        code_replicas(subindex,:) = Code(ceil(time_stamps(subindex,:)) + 4);
        corr_out(subindex) = sum(code_replicas(subindex,:)     .*INCode);
    end

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

        if enable_LMS | enable_NLMS
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
        end


        
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
            
        % RBPF    
        elseif enable_RBPF
            % Initialization
            if (Index == 1)
                particle_pred = zeros(n_part, 2*M_pf);
                particle = zeros(n_part, 2*M_pf);
                particle(:,1) = -0.5 + rand(n_part, 1); % LOS ambiguity
                for idx=1:M_pf-1
                    particle(:,2*idx+1) = particle(:,1) + rand(n_part, 1); % Multipath ambiguity
                end
                weight = ones(n_part,1)/n_part;
            end
            % run newton-raphson algorithm

            % form the gaussian importance density given the previous likelihood approx.

            for idx=1:n_part
                % kalman prediction
                xp = F_k*x;
                Pp = F_k*P*F_k' + Q;





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
                LS_H = LS_S*corr_out;
                for subindex=0:M_pf-1
                    particle(idx,2*subindex+2)=interp1(Spacing,LS_H,particle(idx,2*subindex+1));
                    if(isnan(particle(idx,2*subindex+2)))
                        particle(idx,2*subindex+2)=0;
                    end
                end


                dt = 0.2;
                t  = 0:dt:100;
                Nsamples = length(t);
                Xsaved = zeros(Nsamples, 1);
                Zsaved = zeros(Nsamples, 1);
                for k=1:Nsamples
                    z = GetVolt();  
                    volt = SimpleKalman(z);
                    
                    Xsaved(k) = volt;
                    Zsaved(k) = z;
                end
                figure
                plot(t, Xsaved, 'o-')
                hold on
                plot(t, Zsaved, 'r:*') 


                
                % Weight update
                % transform the particle into the measurements domain - correlation output
                INCodePF=0;
                for subindex=0:M_pf-1
                    if particle(idx,2*subindex+2) ~= 0
                        CodePF = Spacing(P) + particle(idx,2*subindex+1) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(P) + particle(idx,2*subindex+1));
                        INCodePF = INCodePF + Code(ceil(CodePF) + 4).*particle(idx,2*subindex+2);
                    end
                end

                for subindex = 1: total_corr
                    time_stampsPF = (Spacing(subindex)) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(subindex));
                    code_replicasPF = Code(ceil(time_stampsPF) + 4);
                    corr_outPF(subindex) = sum(code_replicasPF.*INCodePF);
                end
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
                    C_s(Index,subindex+1) = C_s(Index,subindex+1) + weight(idx)*((particle(idx,subindex*2+1:subindex*2+2)-x_est_bpf(Index,subindex*2+1:subindex*2+2))*(particle(idx,subindex*2+1:subindex*2+2)-x_est_bpf(Index,subindex*2+1:subindex*2+2))');
                end
            end
            
            % Resampling
            % effective sample size to indicate degeneracy
            N_eff = 1/( sum(weight.^2) );
            if N_eff < (2/3)*n_part
                %Resampling
                cdf = cumsum(weight);
                %Systematic resampling
                sam = rand/n_part;
                for i=1:n_part
                    samInd = sam + (i-1)/n_part;
                    ind = find( samInd<=cdf ,1);
                    particle_pred(i) = particle(ind);
                end
            else
                for idx=1:n_part
                    particle_pred(idx,:) = particle(Index,:);
                end
            end
            


        % PF    
        elseif enable_PF
            % Initialization
            if (Index == 1)
                particle_pred = zeros(n_part, 2*M_pf);
                particle = zeros(n_part, 2*M_pf);
                particle(:,1) = -0.5 + rand(n_part, 1); % LOS ambiguity
                for idx=1:M_pf-1
                    particle(:,2*idx+1) = particle(:,1) + rand(n_part, 1); % Multipath ambiguity
                end
                weight = ones(n_part,1)/n_part;
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
                LS_H = LS_S*corr_out;
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
                    if particle(idx,2*subindex+2) ~= 0
                        CodePF = Spacing(P) + particle(idx,2*subindex+1) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(P) + particle(idx,2*subindex+1));
                        INCodePF = INCodePF + Code(ceil(CodePF) + 4).*particle(idx,2*subindex+2);
                    end
                end

                for subindex = 1: total_corr
                    time_stampsPF = (Spacing(subindex)) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(subindex));
                    code_replicasPF = Code(ceil(time_stampsPF) + 4);
                    corr_outPF(subindex) = sum(code_replicasPF.*INCodePF);
                end
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
                    C_s(Index,subindex+1) = C_s(Index,subindex+1) + weight(idx)*((particle(idx,subindex*2+1:subindex*2+2)-x_est_bpf(Index,subindex*2+1:subindex*2+2))*(particle(idx,subindex*2+1:subindex*2+2)-x_est_bpf(Index,subindex*2+1:subindex*2+2))');
                end
            end
            
            % Resampling
            % Initialize next iteration
            % Every particle equal to MAP
            for idx=1:n_part
                particle_pred(idx,:) = x_est_bpf(Index,:);
            end
        end



        C_sv(:,:,svindex)=weights; % update weights
    end



    % DLL discriminator and loop filter
    if new_tracking
        if enable_PF
            DLLdiscri(1,Index) = -x_est_bpf(Index,1);
            code_output= (0.3749245/0.007030542258775)*2.5*DLLdiscri(1,Index);
            %code_output= 0;
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
            ylim([-0.1 1.2])
            xlim([-0.1 1])
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
            if enable_PF
                for subindex=1:M_pf
                    if subindex==1
                        stem(x_est_bpf(Index,1),x_est_bpf(Index,2),'--om');
                    else
                        stem(x_est_bpf(Index,(3:2:end)),x_est_bpf(Index,(4:2:end)),'--og');
                    end
                end
            end
            drawnow
            pause(0.005);
        end
    end
end

function volt = SimpleKalman(z)
    persistent A H Q R 
    persistent x P
    persistent firstRun
        response = figure(1);
        subplot(222);
        p0=plot((1), (1),'-','Color',"#77AC30");
        title('estimate x')
        hold on;
        % plot initial position of central correlators and filter peak
        
        % plot initial desired response

        % plot multipath information
        drawnow
        subplot(221);
        p1=plot((1), (1),'-','Color',"#77AC30");
        title('error cov P')
        drawnow

        subplot(223);
        p2=plot((1), (1),'-','Color',"#77AC30");
        title('measurement')
        subplot(224);
        p3=plot((1), (1),'-','Color',"#77AC30");
        title('gain')

    if isempty(firstRun)


        A = 1;
        H = 1;
        
        Q = 1;
        R = 4;

        x=zeros(6001,1);
        P=zeros(6001,1);
        z=zeros(6001,1);
        k=zeros(6001,1);
        x(1,1) = 14;
        P(1,1) =  6;
        
        firstRun = 1;  
    end

    for i=1:6000  
        if i>=250
            z(i,1) = 14.4 + 10 +  40*randn(1,1);
            R=40;
        else
            z(i,1) = 14.4 + 4*randn(1,1);
        end
        xp = A*x(i,1);
        Pp = A*P(i,1)*A' + Q;

        K(i,1) = Pp*H'*inv(H*Pp*H' + R);

        x(i+1,1) = xp + K(i,1)*(z(i,1) - H*xp);
        P(i+1,1) = Pp - K(i,1)*H*Pp;


        volt = x;


                    %%%% REFRESH PLOTS %%%%
                    set(p0,'XData',(1:1:i),'YData',x(1:i,1));
                    drawnow
                    set(p1,'XData',(1:1:i),'YData',P(1:i,1));
                    drawnow
                    set(p2,'XData',(1:1:i),'YData',z(1:i,1));
                    drawnow
                    set(p3,'XData',(1:1:i),'YData',K(1:i,1));
                    drawnow
                    %pause(0.001);

    end
end

function z = GetVolt()
    persistent count
    persistent firstRun1
    if isempty(firstRun1)
        
        count = 0;
        firstRun1=1;
    end
    count = count+1;
    w = 0 + 4*randn(1,1);
    if count>=50
        z = 14.4 + 50 + w;
    else
        z = 14.4 + w;
    end
end
%{

%}