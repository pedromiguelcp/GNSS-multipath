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

corr_per_chip = 30;
Spacing = zeros(corr_per_chip*2+1);
time_stamps = zeros(corr_per_chip*2+1,numSample);
code_replicas = zeros(corr_per_chip*2+1,numSample);
initial_results = zeros(1,corr_per_chip*2+1);
corr_out = zeros(1,corr_per_chip*2+1);
corr_out_MEDLL = zeros(1,corr_per_chip*2+1);
corr_out_MEDLL1 = zeros(1,corr_per_chip*2+1);


for Index = 1: (corr_per_chip*2+1)
    Spacing(Index) = -1 + (Index-1)/corr_per_chip;
end

VE = corr_per_chip-1;
E = corr_per_chip;
P = corr_per_chip+1;
L = corr_per_chip+2;
VL = corr_per_chip+3;

% Specify here the multipath delay (chips) and attenuation (linear factor), compared to the LOS signal
mltpth_delays = [0.2 0.3 0.35];
mltpth_attenuation = [2 2.5 3];

epochs = 500; %signal tracking epochs
convergence_it = 0;

% Generate the C/A code
Code = generateCAcode(1); % Generate C/A code for satellite 1
Code = [Code(end) Code Code(1) Code(2) Code(3)];%enable replica slidding


% RLS data
weights = zeros(5,5);
filter_output = zeros(epochs,5);
desired_response = zeros(epochs,5);
filter_error = zeros(epochs+1,5);
val = zeros(1,epochs);
idx = zeros(1,epochs);
P_rls=eye(5);

%tracking
for Index = 1: epochs

    %numSample = round((codelength-remChip)/(codeFreq/fs)); 
    
    LOS_Code    = Spacing(P) : codeFreqBasis/fs : ((numSample -1) * (codeFreqBasis/fs) + Spacing(P));
    INCode      = Code(ceil(LOS_Code) + 1);

    %interference injection
    for subindex = 1: size(mltpth_delays, 2)
        multipath = Spacing(P) + mltpth_delays(subindex) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(P) + mltpth_delays(subindex));
        INCode = INCode + Code(ceil(multipath) + 1)./mltpth_attenuation(subindex);
    end

    %correlator engine
    for subindex = 2: (corr_per_chip*2)
        time_stamps(subindex,:) = (Spacing(subindex) + remChip) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(subindex) + remChip);
    end

    for subindex = 2: (corr_per_chip*2)
        code_replicas(subindex,:) = Code(ceil(time_stamps(subindex,:)) + 1);
    end

    for subindex = 2: (corr_per_chip*2)
        corr_out(subindex) = sum(code_replicas(subindex,:)     .*INCode);
    end

    remChip   = (time_stamps(P,numSample) + codeFreq/fs) - codeFreqBasis*0.001;
    


    % ADAPTIVE CHANNEL COMPENSATION
    [weights, filter_output(Index,:), P_rls] = filter_rls(corr_out(VE : 1 : VL), filter_error(Index,:), weights, P_rls, 0.95);
    
    % LOS SIGNAL ESTIMATION BLOCK
    [val(Index), idx(Index)] = max(filter_output(Index,:));
    cnt=1;
    for subindex = VE : VL
        desired_response(Index,cnt) = 1*sum(code_replicas(subindex,:).*code_replicas(VE + idx(Index) -1,:));
        cnt=cnt+1;
    end
    filter_error(Index+1,:) = desired_response(Index,:) - filter_output(Index,:);
    
    
    
    %y=filter_output(Index,:);%1x5
    %g=desired_response(Index,:);%1x5
    %G=;%
    %a=(g*y)/(g*g);%1x1
    %desired_response(Index,:)=desired_response(Index,:).*a;
    

    


    %DLL
    DLL_E           = sqrt(corr_out(E)^2);
    DLL_L           = sqrt(corr_out(L)^2);

    DLLdiscri       = 0.5 * (DLL_E-DLL_L)/(DLL_E+DLL_L);
    newDLLdiscri    = Spacing(P) - Spacing(idx(Index));
    %DLLdiscri = newDLLdiscri;
    code_output     = code_outputLast + (0.3749245/0.007030542258775)*(DLLdiscri - DLLdiscriLast) + DLLdiscri* (0.001/0.007030542258775);
    DLLdiscriLast   = DLLdiscri;
    code_outputLast = code_output;
    codeFreq        = codeFreqBasis - code_output;

    if(Index == 1) 
        initial_results = corr_out./numSample;
    end
end

% Plot the results
plot(Spacing(2 : 1 : corr_per_chip*2), initial_results(2 : 1 : corr_per_chip*2),'-');
hold on;
for subindex = VE : VL
    plot(Spacing(subindex),initial_results(subindex),'-g*');
    plot(time_stamps(subindex, 1),corr_out(subindex)/numSample,'r*');
end
plot([0 0],[1 0],'-bo');
for subindex = 1: size(mltpth_delays, 2)
    plot([mltpth_delays(subindex) mltpth_delays(subindex)],[1/mltpth_attenuation(subindex) 0],'-bo');
end

ylabel('Normalized Amplitude')
xlabel('Delay [chips]')
title('Initial/Final correlation result')


function [w, xp, P_rls] = filter_rls(u, e, w, P_rls, lambda)

    xp(1,:) = u(1,:)*w';
    % update kappa as perRLS
    kappa = lambda^(-1)*P_rls*u(1,:)'/(1+lambda^(-1)*u(1,:)*P_rls*u(1,:)');
    % update weights
    w = w+kappa*e(1,:); 
    % update as per R
    P_rls = lambda^(-1)*P_rls-lambda^(-1)*u(1,:)*kappa*P_rls; 
end


%{

%}