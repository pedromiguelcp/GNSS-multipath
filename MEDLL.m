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


mltpth_delays = [0.2 0.35 0.45];
mltpth_attenuation = [2 2.5 3];

epochs = 100; %signal tracking epochs
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

%MEDLL
val = zeros(1,size(mltpth_delays, 2)+1);
idx = zeros(1,size(mltpth_delays, 2)+1);

for i = 0 : size(mltpth_delays, 2)

    while convergence_it < 20

        %step 1
        corr_out_MEDLL = corr_out;
        for k = 0 : i-1
            for corridx = 2: (corr_per_chip*2)
                cross_corr = sum(code_replicas(corridx,:).*code_replicas(idx(k+1),:));
                corr_out_MEDLL(corridx) = corr_out_MEDLL(corridx) - cross_corr*val(k+1)/numSample;
            end
        end
        
        %step 2
        [val(i+1), idx(i+1)] = max(corr_out_MEDLL);

        %step 3
        for k = 0 : i-1
            corr_out_MEDLL = corr_out;
            for l = 0 : i
                if l ~= k
                    for corridx = 2: (corr_per_chip*2)
                        cross_corr = sum(code_replicas(corridx,:).*code_replicas(idx(l+1),:));
                        corr_out_MEDLL(corridx) = corr_out_MEDLL(corridx) - cross_corr*val(l+1)/numSample;
                    end
                end
            end
            [val(k+1), idx(k+1)] = max(corr_out_MEDLL);
        end
       
        %step 4, change to find convergence
        convergence_it = convergence_it + 1;
    end

    %convergence epoch results
    for j = 0 : i
        X = sprintf('MEDLL -> Peak %d: Amplitude = %.2f  Delay = %.2f',j, val(j+1)/numSample,(-1+idx(j+1)*1/30));
        disp(X) 
    end
    convergence_it = 0;
end


% Plot the results
plot(Spacing(2 : 1 : corr_per_chip*2), initial_results(2 : 1 : corr_per_chip*2),'-');
hold on;
for subindex = VE : VL
    plot(Spacing(subindex),initial_results(subindex),'-g*');
    plot(time_stamps(subindex, 1),corr_out(subindex)/numSample,'r*');
end

ylabel('Normalized Amplitude')
xlabel('Delay [chips]')
title('Initial/Final correlation result')

%{
if the first and second peaks are switched for some reason, 
    the LOS peak is the one with the lowest delay, even if it has a lower amplitude compared to other peaks
number of multipath must be chosen somehow.. or rather estimated on-the-fly
    the stop criterion seems to bring unnecessary overhead
The phase component is missing from both estimation and modulation
%}
