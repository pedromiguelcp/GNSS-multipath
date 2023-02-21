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


mltpth_delays = [0.2];
mltpth_attenuation = [2];

epochs = 10; %signal tracking epochs
convergence_it = 0;

% Generate the C/A code
Code = generateCAcode(1); % Generate C/A code for satellite 1
Code = [Code(end) Code Code(1) Code(2) Code(3)];%enable replica slidding



%tracking
for Index = 1: epochs

    %numSample = round((codelength-remChip)/(codeFreq/fs)); 
    
    LOS_Code    = Spacing(P) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(P));
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
    
    %DLL
    DLL_E           = sqrt(corr_out(E)^2);
    DLL_L           = sqrt(corr_out(L)^2);

    DLLdiscri       = 0.5 * (DLL_E-DLL_L)/(DLL_E+DLL_L);
    code_output     = code_outputLast + (0.3749245/0.007030542258775)*(DLLdiscri - DLLdiscriLast) + DLLdiscri* (0.001/0.007030542258775);
    DLLdiscriLast   = DLLdiscri;
    code_outputLast = code_output;
    codeFreq        = codeFreqBasis - code_output;

    if(Index == 1) 
        initial_results = corr_out./numSample;
    end
end

%RLS
it=6;
val = zeros(1,size(mltpth_delays, 2)+1);
idx = zeros(1,size(mltpth_delays, 2)+1);
u = zeros(it,5);
d = zeros(it,5);

for i = 0 : 20
      
    [val(i+1), idx(i+1)] = max(corr_out);

    num=1;
    
    for index = 1 : it
        num=1;
        for  subindex = VE: VL
            u(index,num) = corr_out(subindex)/numSample;
            d(index,num) = sum(code_replicas(subindex,:).*code_replicas(idx(i+1),:))/numSample;
            num=num+1;
        end
    end
    
    lambda = 0.9;
    mu = 0.1;

    %filter data
    [W,y] = rls(u, d, lambda);

    [W,y] = lms(u, d, mu);

end


function [W, xp] = rls(u, d, lambda)
    % Maximum number of time step that can be predicted
    N = min(size(u, 1),size(d, 1));
    Nin = size(u,2);
    Nout = size(d,2);
    % Intializatize weight matrix and associated parameters for RLS predictor
    w = zeros(Nout,Nin);
    W = [];
    delta=0.5;
    % reset filter variable between monte carlo runs
    P=eye(Nin)*delta;
    for n = 1:N
        W = [W;w];
        xp(n,:) = u(n,:)*w';
        e(n,:) = d(n,:)-xp(n,:);
        % update kappa as perRLS
        kappa = lambda^(-1)*P*u(n,:)'/(1+lambda^(-1)*u(n,:)*P*u(n,:)');
        % update weights
        w = w+kappa'*e(n); 
        % update as per R
        P = lambda^(-1)*P-lambda^(-1)*u(n,:)*kappa*P; 
    end
end

function [W, e] = lms(u, d, mu)
    % Maximum number of time step that can be predicted
    N = min(size(u, 1),size(d, 1));
    Nin = size(u,2);
    Nout = size(d,2);
    % Intializatize weight matrix and associated parameters for LMS predictor
    w = zeros(Nout, Nin);
    W = [];
    for n = 1:N
        W = [W;w];
        % Predict next sample and error
        xp(n, :) = u(n,:)*w';
        e(n,:) = d(n,:)-xp(n,:);
        % Adapt weight matrix ans step size
        w = w + mu * e(n,:)' * u(n,:);
    end
end
