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
corr_out1 = zeros(numSample,corr_per_chip*2+1);
corr_out_SML = zeros(1,corr_per_chip*2+1);
corr_out_SML1 = zeros(numSample,corr_per_chip*2+1);

full_correlation = zeros(1,numSample);


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

epochs = 10; %signal tracking epochs
convergence_it = 0;

% Generate the C/A code
Code = generateCAcode(1); % Generate C/A code for satellite 1
Code = [Code(end) Code Code(1) Code(2) Code(3)];%enable replica slidding

Y = zeros(numSample,corr_per_chip*2+1);

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
        corr_out(subindex) = sum(code_replicas(subindex,:).*INCode);
        corr_out1(:,subindex) = code_replicas(subindex,:).*INCode;
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

%SML
val = zeros(1,size(mltpth_delays, 2)+1);
idx = zeros(1,size(mltpth_delays, 2)+1);
G = zeros(numSample,size(mltpth_delays, 2)+1);
A = zeros(size(mltpth_delays, 2)+1,size(mltpth_delays, 2)+1);
D = zeros(size(mltpth_delays, 2)+1,size(mltpth_delays, 2)+1);


for i = 0 : size(mltpth_delays, 2)

    while convergence_it < 20

        %step 1
        corr_out_SML = corr_out;
        for k = 0 : i-1
            for corridx = 2: (corr_per_chip*2)
                cross_corr = sum(code_replicas(corridx,:).*code_replicas(idx(k+1),:));
                corr_out_SML(corridx) = corr_out_SML(corridx) - cross_corr*A(k+1,1);

                corr_out_SML1(:,corridx) = corr_out1(:,corridx) - ((code_replicas(corridx,:).*code_replicas(idx(k+1),:))*A(k+1,1))';
            end
        end
        
        %step 2
        [val(i+1), idx(i+1)] = max(corr_out_SML);
        [val1, ] = min(idx(idx>0));
        %Y(:,i+1) = code_replicas(idx(i+1),:).*INCode; % keep feeding the MLE algorithm
        G(:,i+1) = code_replicas(idx(i+1),:).*code_replicas(val1,:);
        %A = pinv(G'*G)*(G'*Y);
        A(i+1,1) = val(i+1)/numSample;
        

        %step 3
        for k = 0 : i-1
            corr_out_SML = corr_out;
            for l = 0 : i
                if l ~= k
                    for corridx = 2: (corr_per_chip*2)
                        cross_corr = sum(code_replicas(corridx,:).*code_replicas(idx(l+1),:));
                        corr_out_SML(corridx) = corr_out_SML(corridx) - cross_corr*A(l+1,1);
                        
                        Y(:,corridx) = corr_out1(:,corridx) - ((code_replicas(corridx,:).*code_replicas(idx(l+1),:))*A(l+1,1))';
                    end
                end
            end

            %{
            perceber como calcular o arg max
            tem tudo de ficar em matrizes como o Y para no final o vetor de amplitudes ter os resultados
            pq dá com LOS + 2 multipaths mas não mais?!
            como injetar no step 3 o y na equação? como está na linha 135?
            %}



            [val(k+1), idx(k+1)] = max(corr_out_SML);
            [val1, ] = min(idx(idx>0));

            %Y(:,k+1) = code_replicas(idx(k+1),:).*INCode;
            G(:,k+1) = code_replicas(idx(k+1),:).*code_replicas(val1,:);
            A = pinv(G'*G)*(G'*Y);
        end
       
        %step 4, change to find convergence
        convergence_it = convergence_it + 1;
    end

    %convergence epoch results
    for j = 0 : i
        X = sprintf('SML -> Peak %d: Amplitude = %.2f  Delay = %.2f',j, A(j+1,1),(-1+idx(j+1)*1/30));
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
