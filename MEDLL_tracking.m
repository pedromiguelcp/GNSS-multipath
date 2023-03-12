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


mltpth_delays = [0.066 0.133 0.233 0.4];
mltpth_attenuation = [1.1 1.2 1.3 1.35];

epochs = 400; %signal tracking epochs
convergence_it = 0;

% Generate the C/A code
Code = generateCAcode(1); % Generate C/A code for satellite 1
Code = [Code(end) Code Code(1) Code(2) Code(3)];%enable replica slidding

%MEDLL
val = zeros(1,size(mltpth_delays, 2)+1);
idx = zeros(1,size(mltpth_delays, 2)+1);

DLLdiscri = zeros(1,epochs);
MEDLL_track = 0;

%tracking
for Index = 1: epochs

    %numSample = round((codelength-remChip)/(codeFreq/fs)); 
    
    LOS_Code    = Spacing(P): codeFreqBasis/fs : ((numSample -1) * (codeFreqBasis/fs) + Spacing(P));
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


    if(Index == 1) 
        % plot initial correlation output
        response = figure(1);
        [val_max, idx_max] = max(corr_out);
        plot(Spacing(2 : 1 : corr_per_chip*2), corr_out(2 : 1 : corr_per_chip*2)./val_max,'-','Color',"#77AC30");
        ylabel('Normalized Amplitude')
        xlabel('Delay [chips]')
        title('Correlation function')
        hold on;
        % plot initial position of central correlators and filter peak
        p1=plot(Spacing((E : 1 : L)),corr_out((E : 1 : L))/numSample,'*','Color',"#77AC30");
        % plot multipath information
        plot([0 0],[1 0],'-bo');
        for subindex = 1: size(mltpth_delays, 2)
            p2=plot([mltpth_delays(subindex) mltpth_delays(subindex)],[1/mltpth_attenuation(subindex) 0],'-bo');
        end
        % one bar for each signal component
        m0=plot([0 0],[0 0],'-','Color',"#D95319");

        if size(mltpth_delays, 2) > 0
            m1=plot([0 0],[0 0],'-','Color',"#D95319");
        end
        if size(mltpth_delays, 2) > 1
            m2=plot([0 0],[0 0],'-','Color',"#D95319");
        end
        if size(mltpth_delays, 2) > 2
            m3=plot([0 0],[0 0],'-','Color',"#D95319");
        end
        if size(mltpth_delays, 2) > 3
            m4=plot([0 0],[0 0],'-','Color',"#D95319");
        end
        if size(mltpth_delays, 2) > 4
            m5=plot([0 0],[0 0],'-','Color',"#D95319");
        end
        drawnow
        dll_fig = figure(2);
        discri=plot((1 : 1 : Index), DLLdiscri(1,(1 : 1 : Index)),'-','Color',"#7E2F8E");
        ylabel('DLL discriminator')
        xlabel('Epochs (ms)')
        title('DLL corrections')
        drawnow
    end





    if MEDLL_track == 1
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
            x_data=[];
            y_data=[];
            %convergence epoch results
            for j = 0 : i
                X = sprintf('MEDLL -> Peak %d: Amplitude = %.2f  Delay = %.2f',j, val(j+1)/numSample,(-1+idx(j+1)*1/corr_per_chip));
                disp(X) 
                
                %plot([(-1+idx(subindex)*1/corr_per_chip) (-1+idx(j+1)*1/corr_per_chip)],[val(subindex)/numSample 0],'-bo');
                x_data=[(time_stamps(P, 1)-1+idx(j+1)*1/corr_per_chip) (time_stamps(P, 1)-1+idx(j+1)*1/corr_per_chip)];
                y_data=[val(j+1)/numSample 0];

                if j == 0
                    set(m0,'XData',x_data,'YData',y_data);
                end
                if j == 1
                    set(m1,'XData',x_data,'YData',y_data);
                end
                if j == 2
                    set(m2,'XData',x_data,'YData',y_data);
                end
                if j == 3
                    set(m3,'XData',x_data,'YData',y_data);
                end
                if j == 4
                    set(m4,'XData',x_data,'YData',y_data);
                end
                if j == 5
                    set(m5,'XData',x_data,'YData',y_data);
                end
                drawnow
            end
            
            convergence_it = 0;
        end
    end







    
    %DLL
    DLL_E           = sqrt(corr_out(E)^2);
    DLL_L           = sqrt(corr_out(L)^2);

    %DLLdiscri(1,Index)       = 0.5 * (DLL_E-DLL_L)/(DLL_E+DLL_L);
    DLLdiscri(1,Index)       = (corr_out(E) - corr_out(L))/numSample;
    if MEDLL_track==1
        [amplitude, idx_delay] = max(val);
        newDLLdiscri    = (-1+idx(idx_delay)*1/corr_per_chip)*(-1);
        DLLdiscri(1,Index)       = newDLLdiscri;
    end
    code_output     = code_outputLast + (0.3749245/0.007030542258775)*(DLLdiscri(1,Index) - DLLdiscriLast) + DLLdiscri(1,Index)* (0.001/0.007030542258775);
    DLLdiscriLast   = DLLdiscri(1,Index);
    code_outputLast = code_output;
    codeFreq        = codeFreqBasis - code_output;

    [val_max, idx_max] = max(corr_out);
    set(p1,'XData',time_stamps((VE : 1 : VL), 1),'YData',corr_out(VE : 1 : VL)/val_max);
    drawnow
    set(discri,'XData',(1 : 1 : Index),'YData',DLLdiscri(1,(1 : 1 : Index)));
    drawnow
end








%{
if the first and second peaks are switched for some reason, 
    the LOS peak is the one with the lowest delay, even if it has a lower amplitude compared to other peaks
number of multipath must be chosen somehow.. or rather estimated on-the-fly
    the stop criterion seems to bring unnecessary overhead
The phase component is missing from both estimation and modulation
%}
