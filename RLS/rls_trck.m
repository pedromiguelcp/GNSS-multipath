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

corr_per_chip = 10;
Spacing = zeros(corr_per_chip*2+1);
time_stamps = zeros(corr_per_chip*2+1,numSample);
code_replicas = zeros(corr_per_chip*2+1,numSample);
corr_out = zeros(1,corr_per_chip*2+1)';
f_corr_out = zeros(1,corr_per_chip*2+1)';


for Index = 1: (corr_per_chip*2+1)
    Spacing(Index) = -1 + (Index-1)/corr_per_chip;
end

VE = corr_per_chip-1;
E = corr_per_chip;
P = corr_per_chip+1;
L = corr_per_chip+2;
VL = corr_per_chip+3;

% Specify here the multipath delay (chips) and attenuation (linear factor), compared to the LOS signal
chan1.delays=[];
chan1.attenuation=[];
chan2.delays=[0.05 0.1 0.15];
chan2.attenuation=[2 4 2];
chan3.delays=[0.1 0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.9];
chan3.attenuation=[2 2.5 3 2 1.5 1.6 1.7 4 3];

epochs = 5000; %signal tracking epochs

DLLdiscri = zeros(1,epochs);

% Generate the C/A code
Code = generateCAcode(1); % Generate C/A code for satellite 1
Code = [Code(end-1) Code(end) Code Code(1) Code(2) Code(3)];%enable replica slidding

%rls params
d = zeros(1,(corr_per_chip*2+1));
channel3 = zeros(ceil((length(d)+1)/2),1);
delta=0.005;
P_rls = eye(ceil((length(d)+1)/2))/delta;
lambda=0.95;
a_los=zeros(1,epochs);
t_los=zeros(1,epochs);
trigger_new_corr_func = 0;
trigger_new_los = 0;
move_los=0;

noise_power=0;
rmse = zeros(1,corr_per_chip*2+1);
load(['ibaseband.mat']);

for Index = 1: epochs
    % section used to change the multipath and LOS position during simulation
    if Index <= 500
        mltpth_delays=chan3.delays;
        mltpth_attenuation=chan3.attenuation;
    elseif Index <= 100
        mltpth_delays=chan2.delays;
        mltpth_attenuation=chan2.attenuation;
        trigger_new_corr_func = 1;
    elseif Index <= 500 % stop moving the LOS at 500 ms
        move_los = move_los+0.001;% change here to where the LOS will be moving (+/-) and the step size (e.g. +0.001 = LOS will be delayed 0.001 chips per ms)
        trigger_new_los = 1;
    end

    %%%%%%% CODE TRACKING %%%%%%%%
    %numSample = round((codelength-remChip)/(codeFreq/fs)); 
    LOS_Code = Spacing(P)+move_los: codeFreqBasis/fs : ((numSample -1) * (codeFreqBasis/fs) + Spacing(P))+move_los;
    INCode   = Code(ceil(LOS_Code) + 2) + noise_power.*randn(1,numSample);

    %interference injection
    for subindex = 1: size(mltpth_delays, 2)
        multipath = Spacing(P) + mltpth_delays(subindex) +move_los: codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(P) + mltpth_delays(subindex))+move_los;
        INCode = INCode + Code(ceil(multipath) + 2)./mltpth_attenuation(subindex) + noise_power.*randn(1,numSample);
    end
    %INCode = iBasebandSignal; % from FGI receiver real data where the RLS performs poorly

    %correlator engine
    for subindex = 1: (corr_per_chip*2)+1
        time_stamps(subindex,:) = (Spacing(subindex) + remChip) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(subindex) + remChip);
    end
    for subindex = 1: (corr_per_chip*2)+1
        code_replicas(subindex,:) = Code(ceil(time_stamps(subindex,:)) + 2);
    end 
    for subindex = 1: (corr_per_chip*2)+1
        corr_out(subindex) = sum(code_replicas(subindex,:)     .*INCode);
    end
    %corr_out=corr_out+std(corr_out)*randn(1,21)/10;

    remChip   = (time_stamps(P,numSample) + codeFreq/fs) - codeFreqBasis*0.001;
    u=corr_out';
    for ind = 1 : corr_per_chip*2+1
        for subindex = 1 : corr_per_chip*2+1
            d(1,subindex) = sum(code_replicas(subindex,:).*code_replicas(ind,:));
        end
        rmse(1,ind) = sqrt(mean((u - d).^2));
    end
    [a_los, t_los] = min(rmse);
    for subindex = 1 : corr_per_chip*2+1
        f_corr_out(subindex) = sum(code_replicas(subindex,:).*code_replicas(t_los,:));
    end



    %%%%%%% MULTIPATH COMPENSATION %%%%%%%%
    % due to doppler variations we do not know where the LOS peak is and to where it is moving
    % we assume one position for the LOS peak
    % and through the inverse channel determine the true LOS position relative to the one we assumed

    % first we produce the pseudo LOS, centered at the current position of the prompt correlator (corr_per_chip+1)
        %for subindex = 1: (corr_per_chip*2)+1
        %    d(subindex) = sum(code_replicas(subindex,:).*code_replicas(corr_per_chip+1,:));
        %end

    % compute the inverse channel
    % input signal - pseudo LOS, desired signal - current correlation function (distorted)
        %[channel3, P_rls] = filter_rls(corr_out, d, lambda, channel3, P_rls);

    % find the true LOS delay and amplitude
        %[a_los(Index), t_los(Index)] = max(channel3);

    % at this point we know where the LOS is (from the inverse channel)
    % now we produce an artifitial LOS (perfect triangle) located in that LOS position
    % if the inverse channel indicates that the true LOS peak is in the prompt correlator position, the DLL does not need to do corrections
    % but if the true LOS is to the left (early) or right (late) of the prompt correlator
    % the artifitial correlation function should have the peak shifted accordingly
    % the DLL will compare the Early Late correlator values and they will not be equal, so it will produce a correction
        %for subindex = 1: (corr_per_chip*2)+1
        %    f_corr_out(subindex) = sum(code_replicas(subindex,:).*code_replicas(corr_per_chip+1+t_los(Index)-2,:));
        %end



    % the DLL is the same, but we substitute the real correlator values by the artifitial correlator values
    %DLL_E  = sqrt(corr_out(E)^2);
    %DLL_L  = sqrt(corr_out(L)^2);
    DLL_E   = sqrt(f_corr_out(E)^2);
    DLL_L   = sqrt(f_corr_out(L)^2);
    DLLdiscri(1,Index)  = 0.5 * (DLL_E-DLL_L)/(DLL_E+DLL_L);
    %DLLdiscri(1,Index) = 0.5 * (time_stamps(corr_per_chip+1,1)-time_stamps(corr_per_chip+1+t_los(Index)-2,1)); % this is another discriminator implementation
    code_output = code_outputLast + (0.3749245/0.007030542258775)*(DLLdiscri(1,Index) - DLLdiscriLast) + DLLdiscri(1,Index)* (0.001/0.007030542258775);
    
    %code_output = 0; % uncomment to disable the DLL correction
    DLLdiscriLast   = DLLdiscri(1,Index);
    code_outputLast = code_output;
    codeFreq        = codeFreqBasis - code_output;
    

    

    %%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%
    if(Index == 1) 
        % plot initial correlation output
        colororder({'[0.4660 0.6740 0.1880]','#D95319'})
        response = figure(1);
        subplot(2,2,1:2);
        yyaxis left
        p0=plot(time_stamps(:,1), corr_out,'-','DisplayName','Raw corr. func.' ,'Color',"#77AC30");
        %ylim([-0.1 8])
        xlim([-1.5 1.5])
        ylabel('Amplitude')
        xlabel('Delay [chips]')
        title('Raw Correlation function')
        hold on;
        % plot raw corr func - central correlators
        p1=plot(time_stamps(VE : 1 : VL,1),corr_out(VE : 1 : VL),'*','HandleVisibility','off' ,'Color',"#77AC30");
        % plot filtered corr func
        yyaxis right
        p2=plot(time_stamps(:, 1), f_corr_out,'-','DisplayName','Artificial corr. func.' ,'Color',"#D95319");
        % plot filtered corr func - central correlators
        p3=plot(time_stamps(VE : 1 : VL,1), f_corr_out(VE : 1 : VL),'*','HandleVisibility','off','Color',"#D95319"); 
        legend
        % plot multipath information
        p4=plot([0 0],[0 1],'-bo','HandleVisibility','off');
        drawnow
        subplot(2,2,3);
        p5=plot((1 : 1 : Index), DLLdiscri(1,(1 : 1 : Index)),'-','HandleVisibility','off','Color',"#7E2F8E");
        ylabel('DLL discriminator')
        xlabel('Epochs (ms)')
        title('DLL corrections')
        subplot(2,2,4);
        p7=plot(Spacing(E : 1 : end-1,1), channel3,'-*','HandleVisibility','off','Color',"#7E2F8E");
        ylim([-0.5 2])
        ylabel('Magnitude')
        xlabel('H')
        title('Inverse Channel')
        drawnow
    else
        if trigger_new_corr_func
            set(p0,'XData',time_stamps(:,1),'YData',corr_out);
        end
        set(p1,'XData',time_stamps(VE : 1 : VL,1),'YData',corr_out(VE : 1 : VL));
        set(p2,'XData',time_stamps(:,1),'YData',f_corr_out);
        set(p3,'XData',time_stamps(VE : 1 : VL,1),'YData',f_corr_out(VE : 1 : VL));
        if trigger_new_los
            set(p4,'XData',[move_los move_los],'YData',[1 0]);
        end
        set(p5,'XData',(1 : 1 : Index),'YData',DLLdiscri(1,(1 : 1 : Index)));
        set(p7,'XData',Spacing(E : 1 : end-1,1),'YData',channel3);
    end 
    pause(0.01);  
end



function [h, P] = filter_rls(u, d, lambda, h, P)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Multipath is just the inverse system
    N=length(d);%number of correlators
    d=d(:);
    for i=1:1
        d=[d; 0];
        u=[0; u];
    end
    M=ceil(length(d)/2);
    %delta=0.005;
    %P = eye(M)/delta;
    uvec=zeros(M,1);
    %h=zeros(M,1);
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

