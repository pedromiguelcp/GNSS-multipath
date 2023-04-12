clear all;
close all;
clc;
% Define the parameters for the simulation
c = 299792458; % m/s
fs = 12500000; % Sample rate (Hz)
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
weights3 = zeros(corr_per_chip*2+1,1);
channel3 = zeros(1,corr_per_chip*2+1);
y3 = zeros(1,(corr_per_chip*2+1));
d = zeros(1,(corr_per_chip*2+1));
delta=0.1;
P_rls3 = eye(corr_per_chip*2+1)/delta;
lambda=0.95;
a_los=zeros(1,epochs);
t_los=zeros(1,epochs);
trigger_new_corr_func = 0;
trigger_new_los = 0;
move_los=0;

for Index = 1: epochs
    if Index <= 50
        mltpth_delays=chan3.delays;
        mltpth_attenuation=chan3.attenuation;
    elseif Index <= 100
        mltpth_delays=chan2.delays;
        mltpth_attenuation=chan2.attenuation;
        trigger_new_corr_func = 1;
    elseif Index <= 500
        move_los = move_los-0.001;
        trigger_new_los = 1;
    end

    %numSample = round((codelength-remChip)/(codeFreq/fs)); 
    if Index == 1
        LOS_Code = Spacing(P) : codeFreqBasis/fs : ((numSample -1) * (codeFreqBasis/fs) + Spacing(P));
    else
        LOS_Code = Spacing(P)+move_los: codeFreqBasis/fs : ((numSample -1) * (codeFreqBasis/fs) + Spacing(P))+move_los;
    end
    INCode   = Code(ceil(LOS_Code) + 2);

    %interference injection
    for subindex = 1: size(mltpth_delays, 2)
        multipath = Spacing(P) + mltpth_delays(subindex) +move_los: codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(P) + mltpth_delays(subindex))+move_los;
        INCode = INCode + Code(ceil(multipath) + 2)./mltpth_attenuation(subindex);
    end

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



    %%%%%%% MULTIPATH COMPENSATION %%%%%%%%


    for subindex = 1: (corr_per_chip*2)+1
        d(subindex) = sum(code_replicas(subindex,:).*code_replicas(corr_per_chip+1,:));
    end

    % filter size = total correlators, full slidding
    % imp response = half correlators, delayed
    [channel3, weights3, y3, P_rls3] = filter_rls2(corr_out, d, weights3, P_rls3, lambda);

    RMSE3 = sqrt(mean((y3 - d).^2));

    [a_los(Index), t_los(Index)] = max(channel3);

    % filtered correlation function
    % the traditional DLL should correct the disalignment
    % the RLS is useless here because we already know where the LOS is
    for subindex = 1: (corr_per_chip*2)+1
        f_corr_out(subindex) = sum(code_replicas(subindex,:).*code_replicas(corr_per_chip+1+t_los(Index)-2,:));
    end


    



    %DLL
    %DLL_E           = sqrt(corr_out(E)^2);
    DLL_E           = sqrt(f_corr_out(E)^2);
    %DLL_L           = sqrt(corr_out(L)^2);
    DLL_L           = sqrt(f_corr_out(L)^2);
    DLLdiscri(1,Index)       = 0.5 * (DLL_E-DLL_L)/(DLL_E+DLL_L);
    %DLLdiscri(1,Index)       = 0.5 * (time_stamps(corr_per_chip+1,1)-time_stamps(corr_per_chip+1+t_los(Index)-2,1));
    code_output     = code_outputLast + (0.3749245/0.007030542258775)*(DLLdiscri(1,Index) - DLLdiscriLast) + DLLdiscri(1,Index)* (0.001/0.007030542258775);
    %if t_los(Index) == 2
    %    code_output = 0;
    %end
    DLLdiscriLast   = DLLdiscri(1,Index);
    code_outputLast = code_output;
    codeFreq        = codeFreqBasis - code_output;
    

    

    %%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%
    if(Index == 1) 
        % plot initial correlation output
        response = figure(1);
        subplot(2,2,1:2);
        p0=plot(time_stamps(:,1), corr_out./numSample,'-','DisplayName','Raw corr. func.' ,'Color',"#77AC30");
        ylim([-0.1 8])
        xlim([-1.5 1.5])
        ylabel('Normalized Amplitude')
        xlabel('Delay [chips]')
        title('Raw Correlation function')
        hold on;
        % plot raw corr func - central correlators
        p1=plot(time_stamps(VE : 1 : VL,1),corr_out(VE : 1 : VL)./numSample,'*','HandleVisibility','off' ,'Color',"#77AC30");
        % plot filtered corr func
        p2=plot(time_stamps(:, 1), f_corr_out./numSample,'-','DisplayName','Filtered corr. func.' ,'Color',"#D95319");
        % plot filtered corr func - central correlators
        p3=plot(time_stamps(VE : 1 : VL,1), f_corr_out(VE : 1 : VL)./numSample,'*','HandleVisibility','off','Color',"#D95319"); 
        legend
        % plot multipath information
        p4=plot([0 0],[1 0],'-bo','HandleVisibility','off');
        %for subindex = 1: size(mltpth_delays, 2)
        %    plot([mltpth_delays(subindex) mltpth_delays(subindex)],[1/mltpth_attenuation(subindex) 0],'-bo','HandleVisibility','off');
        %end
        drawnow
        subplot(2,2,3);
        p5=plot((1 : 1 : Index), DLLdiscri(1,(1 : 1 : Index)),'-','HandleVisibility','off','Color',"#7E2F8E");
        ylabel('DLL discriminator')
        xlabel('Epochs (ms)')
        title('DLL corrections')
        %subplot(2,2,3);
        %p6=plot(weights3,'-*','Color',"#7E2F8E");
        %ylim([-2.5 2.5])
        %ylabel('Magnitude')
        %xlabel('W')
        %title('Filter weights')
        %drawnow
        subplot(2,2,4);
        p7=plot(Spacing(E : 1 : end-1,1), channel3,'-*','HandleVisibility','off','Color',"#7E2F8E");
        ylim([-0.5 2])
        ylabel('Magnitude')
        xlabel('H')
        title('Inverse Channel')
        drawnow
    else
        if trigger_new_corr_func
            set(p0,'XData',time_stamps(:,1),'YData',corr_out./numSample);
        end
        set(p1,'XData',time_stamps(VE : 1 : VL,1),'YData',corr_out(VE : 1 : VL)./numSample);
        set(p2,'XData',time_stamps(:,1),'YData',f_corr_out./numSample);
        set(p3,'XData',time_stamps(VE : 1 : VL,1),'YData',f_corr_out(VE : 1 : VL)./numSample);
        if trigger_new_los
            set(p4,'XData',[move_los move_los],'YData',[1 0]);
        end
        set(p5,'XData',(1 : 1 : Index),'YData',DLLdiscri(1,(1 : 1 : Index)));
        %set(p6,'XData',(1 : 1 : length(weights3)),'YData',weights3);
        set(p7,'XData',Spacing(E : 1 : end-1,1),'YData',channel3);
    end 
    pause(0.01);  
end

%{

%}


function [h, w, xp, P_rls] = filter_rls2(u, d, w, P_rls, lambda)
  
    N=length(d);%number of correlators
    d=d(:);
    uvec=zeros(N,1);
    
    for i=1:N
        uvec(1:1:i)=u(i:-1:1);
        
        kappa = lambda^(-1)*P_rls*uvec/(1+lambda^(-1)*uvec'*P_rls*uvec);
        e=d(i)-w'*uvec;
        w=w+kappa*conj(e);
        P_rls = lambda^(-1)*P_rls-lambda^(-1)*kappa*uvec'*P_rls;
    end

    xp=filter(w,1,u);      %filtered signal
    
    
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

    %imp=conv(w,h);%must be (approaches) de impulse function
end

