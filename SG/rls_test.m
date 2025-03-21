%clear all;
%close all;
%clc;
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
chan2.delays=[0.05 0.1];
chan2.attenuation=[2 4];
chan3.delays=[0.1 0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.9];
chan3.attenuation=[2 2.5 3 2 1.5 1.6 1.7 4 3];

epochs = 2; %signal tracking epochs

% Generate the C/A code
Code = generateCAcode(1); % Generate C/A code for satellite 1
Code = [Code(end-1) Code(end) Code Code(1) Code(2) Code(3)];%enable replica slidding

%rls params
Path_len=size(chan2.delays,2)+1;
weights = zeros(Path_len,1);
weights1 = zeros(corr_per_chip+1,1);
weights2 = zeros(corr_per_chip+1,1);
weights3 = zeros(corr_per_chip*2+1,1);
channel = zeros(1,Path_len);
channel1 = zeros(1,corr_per_chip+1);
channel2 = zeros(1,corr_per_chip+1);
channel3 = zeros(1,corr_per_chip*2+1);
channel4 = zeros(1,corr_per_chip*2+1);
channel5 = zeros(1,3);
y = zeros(1,(corr_per_chip*2+1));
y1 = zeros(1,(corr_per_chip*2+1));
y2 = zeros(1,(corr_per_chip*2+1));
y3 = zeros(1,(corr_per_chip*2+1));
d = zeros(1,(corr_per_chip*2+1));
delta=0.1;
P_rls = eye(Path_len)/delta;
P_rls1 = eye(corr_per_chip+1)/delta;
P_rls2 = eye(corr_per_chip+1)/delta;
P_rls3 = eye(corr_per_chip*2+1)/delta;
lambda=0.98;
los_delay=0;

for Index = 1: epochs
    if Index == 1
        mltpth_delays=chan1.delays;
        mltpth_attenuation=chan1.attenuation;
    else
        mltpth_delays=chan3.delays;
        mltpth_attenuation=chan3.attenuation;
    end

    %numSample = round((codelength-remChip)/(codeFreq/fs)); 
    if Index == 1
        LOS_Code = Spacing(P) : codeFreqBasis/fs : ((numSample -1) * (codeFreqBasis/fs) + Spacing(P));
    else
        LOS_Code = Spacing(P)+los_delay: codeFreqBasis/fs : ((numSample -1) * (codeFreqBasis/fs) + Spacing(P))+los_delay;
    end
    INCode   = Code(ceil(LOS_Code) + 2);

    %interference injection
    for subindex = 1: size(mltpth_delays, 2)
        multipath = Spacing(P) + mltpth_delays(subindex) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(P) + mltpth_delays(subindex));
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


    if Index ==1
        d=corr_out;
    else

        % original (path length)
        % imp response original (path length)
        [channel, weights, y, P_rls] = filter_rls(corr_out, d, weights, P_rls, lambda);
        % filter size = half correlators
        % imp response = half correlators
        [channel1, weights1, y1, P_rls1] = filter_rls(corr_out, d, weights1, P_rls1, lambda);
        % filter size = half correlators, full slidding
        % imp response = half correlators, full slidding
        [channel2, weights2, y2, P_rls2] = filter_rls1(corr_out, d, weights2, P_rls2, lambda);
        % filter size = total correlators, full slidding
        % imp response = half correlators, delayed
        [channel3, weights3, y3, P_rls3] = filter_rls2(corr_out, d, weights3, P_rls3, lambda);

        [channel_new] = filter_rls_new1(corr_out, d, lambda);

        
        
        [channel4] = filter_rls3(corr_out, d, lambda);
        
        [channel5] = filter_rls4(corr_out, d, lambda);

        figure
        plot(d,'DisplayName','LOS');
        hold on;
        plot(corr_out,'DisplayName','Interference');
        plot(y,'DisplayName','Original');
        RMSE = sqrt(mean((y - d).^2));
        plot(y1,'DisplayName','Filtered v2');
        RMSE1 = sqrt(mean((y1 - d).^2));
        plot(y2,'DisplayName','Filtered v3');
        RMSE2 = sqrt(mean((y2 - d).^2));
        plot(y3,'DisplayName','Filtered v4');
        RMSE3 = sqrt(mean((y3 - d).^2));
        legend

        figure(2)
        subplot(4,1,1);
        xlim([-2/corr_per_chip 1.1])
        plot([los_delay los_delay],[0 1],'-bo','Color',"#D95319",'DisplayName','LOS');
        xlim([-2/corr_per_chip 1.1])
        hold on;
        for mltpath_idx=1:size(mltpth_delays,2)
            plot([mltpth_delays(mltpath_idx) mltpth_delays(mltpath_idx)],[0 1/mltpth_attenuation(mltpath_idx)],'-bo','Color',"#D95319");
        end
        subplot(4,1,2);
        for mltpath_idx=0:size(channel,1)-1
            plot([1/corr_per_chip*mltpath_idx 1/corr_per_chip*mltpath_idx],[0 abs(channel(mltpath_idx+1))],'-bo','Color',"#D95319");
            hold on;
        end
        xlim([-2/corr_per_chip 1.1])
        subplot(4,1,3);
        for mltpath_idx=0:size(channel3,1)-1
            plot([(-1/corr_per_chip + 1/corr_per_chip*mltpath_idx) (-1/corr_per_chip + 1/corr_per_chip*mltpath_idx)],[0 abs(channel3(mltpath_idx+1))],'-bo','Color',"#D95319");
            hold on;
        end
        xlim([-2/corr_per_chip 1.1])
        subplot(4,1,4);
        for mltpath_idx=0:size(channel4,1)-1
            plot([(0 + 1/corr_per_chip*mltpath_idx) (0 + 1/corr_per_chip*mltpath_idx)],[0 abs(channel4(mltpath_idx+1))],'-bo','Color',"#D95319");
            hold on;
        end
        xlim([-2/corr_per_chip 1.1])
        done=1;
    end
    
  
end

%{

%}
function [h, w, xp, P_rls] = filter_rls(u, d, w, P_rls, lambda)
 
   
    xi=d;  %error  
    M=length(w);%number of multipaths
    N=length(d);%number of correlators
    d=d(:);
    
    for i=M:N
        uvec=u(i:-1:i-M+1);
        kappa = lambda^(-1)*P_rls*uvec/(1+lambda^(-1)*uvec'*P_rls*uvec);
        xi(i)=d(i)-w'*uvec;
        w=w+kappa*conj(xi(i));
        P_rls = lambda^(-1)*P_rls-lambda^(-1)*kappa*uvec'*P_rls;
    end
    xp=filter(w,1,u);      %filtered signal
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Multipath is just the inverse system
    delta=0.005;
    P = eye(M)/delta;

    h=zeros(M,1);
    for i=M:N
        uvec=d(i:-1:i-M+1);
        kappa = lambda^(-1)*P*uvec/(1+lambda^(-1)*uvec'*P*uvec);
        xi(i)=u(i)-h'*uvec;
        h=h+kappa*conj(xi(i));
        P = lambda^(-1)*P-lambda^(-1)*kappa*uvec'*P;
    end

    %imp=conv(w,h);%must be (approaches) de impulse function

end

function [h, w, xp, P_rls] = filter_rls1(u, d, w, P_rls, lambda)
 
   
    e=0;  %error  
    M=length(w);%number of multipaths
    N=length(d);%number of correlators
    d=d(:);
    uvec=zeros(M,1);
    
    for i=1:N
        if i<M
            uvec(1:1:i)=u(i:-1:1);
        else
            uvec=u(i:-1:i-M+1);
        end
        
        kappa = lambda^(-1)*P_rls*uvec/(1+lambda^(-1)*uvec'*P_rls*uvec);
        e=d(i)-w'*uvec;
        w=w+kappa*conj(e);
        P_rls = lambda^(-1)*P_rls-lambda^(-1)*kappa*uvec'*P_rls;
    end

    xp=filter(w,1,u);      %filtered signal
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Multipath is just the inverse system
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
        xi(i)=u(i)-h'*uvec;
        h=h+kappa*conj(xi(i));
        P = lambda^(-1)*P-lambda^(-1)*kappa*uvec'*P;
    end

    %imp=conv(w,h);%must be (approaches) de impulse function

end

function [h, w, xp, P_rls] = filter_rls2(u, d, w, P_rls, lambda)
 
   
    e=0;  %error  
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
        xi(i)=u(i)-h'*uvec;
        h=h+kappa*conj(xi(i));
        P = lambda^(-1)*P-lambda^(-1)*kappa*uvec'*P;
    end

    %imp=conv(w,h);%must be (approaches) de impulse function
end



function [h] = filter_rls3(u, d, lambda)

    d=d(:);
    N=length(d);%number of correlators

    e=0;  %error  
    delta=0.005;
    P = eye(N)/delta;
    uvec=zeros(N,1);
    h=zeros(N,1);
    for i=1:N
        uvec(1:1:i)=d(i:-1:1);
        kappa = lambda^(-1)*P*uvec/(1+lambda^(-1)*uvec'*P*uvec);
        xi(i)=u(i)-h'*uvec;
        h=h+kappa*conj(xi(i));
        P = lambda^(-1)*P-lambda^(-1)*kappa*uvec'*P;
    end
end


function [h] = filter_rls4(u, d, lambda)
    
    d=d(:);
    N=length(d);%number of correlators

    e=0;  %error  
    delta=0.005;
    P = eye(N)/delta;
    uvec=zeros(N,1);
    h=zeros(N,1);
    for i=1:N
        uvec(1:1:i)=d(i:-1:1);
        kappa = lambda^(-1)*P*uvec/(1+lambda^(-1)*uvec'*P*uvec);
        e=u(i)-h'*uvec;
        h=h+kappa*conj(e);
        P = lambda^(-1)*P-lambda^(-1)*kappa*uvec'*P;
    end
end


function [h] = filter_rls_new(u, d, lambda)
  
    N=length(d);%number of correlators
    d=d(:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Multipath is just the inverse system
    for i=1:floor(N/2)*2
        d=[d; 0];
        u=[0; u];
    end
    M=ceil(length(d)/2);
    delta=0.005;
    P = eye(M)/delta;
    uvec=zeros(M,1);
    h=zeros(M,1);
    for i=1:N+floor(N/2)*2
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
function [h] = filter_rls_new1(u, d, lambda)
  
    N=length(d);%number of correlators
    d=d(:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Multipath is just the inverse system
    for i=1:4
        d=[d; 0];
        u=[0; u];
    end
    M=ceil(length(d)/2);
    delta=0.005;
    P = eye(M)/delta;
    uvec=zeros(M,1);
    h=zeros(M,1);
    for i=1:N+4
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
