%clear all;
%close all;
%clc;
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

corr_per_chip = 20;
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
load(['ibaseband.mat']);

for Index = 1: 1
    % section used to change the multipath and LOS position during simulation
    if Index <= 50
        mltpth_delays=chan3.delays;
        mltpth_attenuation=chan3.attenuation;
    elseif Index <= 100
        mltpth_delays=chan3.delays;
        mltpth_attenuation=chan3.attenuation;
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
    
    
    %%%%%%% RLS MULTIPATH COMPENSATION %%%%%%%%
    figure;
    subplot(5,2,1:2);
    plot(corr_out,'-*');
    title('Original correlation function');
    % PSEUDO LOS
    for subindex = 1: (corr_per_chip*2)+1
        d(subindex) = sum(code_replicas(subindex,:).*code_replicas(corr_per_chip+1,:));
    end
    % BEST CHANNEL IMPULSE RESPONSE
    for multp_num = 1 : 50
        [ch_imp_resp] = filter_rls(corr_out, d, lambda, multp_num);
        corr=conv(d,ch_imp_resp);
        corr=corr(:);
        x_corr_out=corr(1:1:(corr_per_chip*2)+1);
        rmse(1,multp_num) = sqrt(mean((corr_out - x_corr_out).^2));
    end
    [a0, t0] = min(rmse);
    [ch_imp_resp] = filter_rls(corr_out, d, lambda, t0);
    subplot(5,2,3);
    plot(ch_imp_resp,'-*');
    title('Channel Impulse Response');
             
    [a1, t1] = max(ch_imp_resp);
    corr=conv(d,ch_imp_resp(1:1:end));
    subplot(5,2,4);
    plot(corr(1:1:(corr_per_chip*2)+1),'-*');
    title('Correlation Function from Impulse Response'); 

    % MULTIPATH INTERFERENCE COMPONENT
    ch_imp_resp(1) = 0;
    subplot(5,2,5);
    plot(ch_imp_resp,'-*');
    title('New Channel Impulse Response');
    corr=conv(d,ch_imp_resp(1:1:end));
    corr=corr(:);
    x_corr_out=corr(1:1:(corr_per_chip*2)+1);
    subplot(5,2,6);
    plot(x_corr_out,'-*');
    title('Multipath Correlation Function');
        
    % MULTIPATH COMPENSATION v1
    for subindex = 1: (corr_per_chip*2)+1
        new_corr(subindex) = corr_out(subindex) - x_corr_out(subindex);
    end
    subplot(5,2,7:8);
    plot(new_corr,'-*');
    title('Multipath-free Correlation Function');
        
    % MULTIPATH COMPENSATION v2
    new_corr1 = corr_out(:);
    for subindex = 2 : length(ch_imp_resp)
        if (ch_imp_resp(subindex) > 0) & (subindex <= corr_per_chip+1)
            for subindex1 = 1: (corr_per_chip*2)+1
                xd(subindex1,1) = sum(code_replicas(subindex1,:).*code_replicas(corr_per_chip+subindex,:));
            end
            new_corr1 = new_corr1 - xd.*ch_imp_resp(subindex);
        end
    end

    subplot(5,2,9:10);
    plot(new_corr1,'-*');
    title('Multipath-free Correlation Function v2');

    figure;
    plot(rmse,'-*');
    title('Multipath number RMSE');
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

    pause(0.01);  
end



function [h] = filter_rls(u, d, lambda, M)
  
    N=length(d);%number of correlators
    d=d(:);
    
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

