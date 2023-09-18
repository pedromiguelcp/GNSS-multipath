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
remChip_plot = 0;
code_outputLast = 0;
DLLdiscriLast = 0;
cross_corr = 0;

corr_per_chip = 5;
corr_per_chip_plot = 20;

Spacing = zeros(corr_per_chip*2+1);
Spacing_plot = zeros(corr_per_chip_plot*2+1);
time_stamps = zeros(corr_per_chip*2+1,numSample);
time_stamps_plot = zeros(corr_per_chip_plot*2+1,numSample);
code_replicas = zeros(corr_per_chip*2+1,numSample);
code_replicas_plot = zeros(corr_per_chip_plot*2+1,numSample);
corr_out = zeros(1,corr_per_chip*2+1);
corr_out_plot = zeros(1,corr_per_chip_plot*2+1);


for Index = 1: (corr_per_chip*2+1)
    Spacing(Index) = -1 + (Index-1)/corr_per_chip;
end
for Index = 1: (corr_per_chip_plot*2+1)
    Spacing_plot(Index) = -1 + (Index-1)/corr_per_chip_plot;
end

VE = corr_per_chip-1;
E = corr_per_chip;
P = corr_per_chip+1;
L = corr_per_chip+2;
VL = corr_per_chip+3;

% Specify here the multipath delay (chips) and attenuation (linear factor), compared to the LOS signal
mltpth_delays = [];
mltpth_attenuation = [];

epochs = 2; %signal tracking epochs

% Generate the C/A code
Code = generateCAcode(1); % Generate C/A code for satellite 1
Code = [Code(end-1) Code(end) Code Code(1) Code(2) Code(3)];%enable replica slidding

DLLdiscri = zeros(1,epochs);

tap_weights = zeros(1, (corr_per_chip_plot*2)+1);

Filter_ON = 0;
Dynamic_simulation = 1;



%tracking
for Index = 1: epochs

    %numSample = round((codelength-remChip)/(codeFreq/fs)); 
    if (Dynamic_simulation == 1) && (Index == 2)
        mltpth_delays = [0.1 0.25 0.35];
        mltpth_attenuation = [2 2.5 3];
    end

    if (Index == 1) || (Index == 2)
        LOS_Code = Spacing(P) : codeFreqBasis/fs : ((numSample -1) * (codeFreqBasis/fs) + Spacing(P));
        INCode   = Code(ceil(LOS_Code) + 2);
        LOS_Code_plot = Spacing_plot(corr_per_chip_plot+1) : codeFreqBasis/fs : ((numSample -1) * (codeFreqBasis/fs) + Spacing_plot(corr_per_chip_plot+1));
        INCode_plot   = Code(ceil(LOS_Code_plot) + 2) + randn(numSample);
    
        %interference injection
        for subindex = 1: size(mltpth_delays, 2)
            multipath = Spacing(P) + mltpth_delays(subindex) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing(P) + mltpth_delays(subindex));
            INCode = INCode + Code(ceil(multipath) + 2)./mltpth_attenuation(subindex);
            multipath_plot = Spacing_plot(corr_per_chip_plot+1) + mltpth_delays(subindex) : codeFreq/fs : ((numSample -1) * (codeFreq/fs) + Spacing_plot(corr_per_chip_plot+1) + mltpth_delays(subindex));
            INCode_plot = INCode_plot + Code(ceil(multipath_plot) + 2)./mltpth_attenuation(subindex);
        end
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
    remChip   = (time_stamps(P,numSample) + codeFreq/fs) - codeFreqBasis*0.001;

    % plot purposes - more resolution in the correlation function plot
    for subindex = 1: (corr_per_chip_plot*2)+1%remove
        time_stamps_plot(subindex,:) = Spacing_plot(subindex) : codeFreqBasis/fs : (numSample -1) * (codeFreqBasis/fs) + Spacing_plot(subindex);
    end
    for subindex = 1: (corr_per_chip_plot*2)+1%remove
        code_replicas_plot(subindex,:) = Code(ceil(time_stamps_plot(subindex,:)) + 2);
    end
    for subindex = 1: (corr_per_chip_plot*2)+1%remove
        corr_out_plot(subindex) = sum(code_replicas_plot(subindex,:)     .*INCode_plot);
    end

    % INPUT SIGNAL 

    for subindex = 1 : (corr_per_chip_plot*2)+1
        u(1,subindex) = sum(code_replicas_plot(subindex,:).*code_replicas_plot(corr_per_chip_plot+1,:));
    end
    d(1,:) = corr_out_plot;

    mu = 0.1;
    lambda = 0.98;
    it = 50;
    [y] = modified_lms(u, d, it, mu);
    %[y] = modified_rls(u, d, lambda);



    %% Channel compensation %%
    % remove contribution of multipath components on the correlation function
    tap_weights = abs(d./u);






    %DLL
    DLL_E           = sqrt(corr_out(E)^2);
    DLL_L           = sqrt(corr_out(L)^2);

    DLLdiscri(1,Index)       = 0.5 * (DLL_E-DLL_L)/(DLL_E+DLL_L);
    if Filter_ON==1
        newDLLdiscri    = Spacing(P) - Spacing(tau_los(Index)); %difference between delay of prompt and delay of max peak
        DLLdiscri(1,Index)       = newDLLdiscri;
    end
    if (Filter_ON==1) && (Index > 1) && (DLLdiscri(1,Index) == 0) && (DLLdiscri(1,Index-1) == 0)
        code_output = 0; % normal DLL always keeps sliding the correlators (drifting) to search near peaks
    else
        code_output     = code_outputLast + (0.3749245/0.007030542258775)*(DLLdiscri(1,Index) - DLLdiscriLast) + DLLdiscri(1,Index)* (0.001/0.007030542258775);
    end

    
    DLLdiscriLast   = DLLdiscri(1,Index);
    code_outputLast = code_output;
    codeFreq        = codeFreqBasis - code_output;


    if(Index == 1) 
        % plot initial correlation output
        [val_max, idx_max] = max(corr_out_plot);
        colororder({'[0.4660 0.6740 0.1880]','[0 0.4470 0.7410]'})
        response = figure(1);
        subplot(1,2,1);
        yyaxis left
        corr_function=plot(Spacing_plot(2 : 1 : corr_per_chip_plot*2), corr_out_plot(2 : 1 : corr_per_chip_plot*2)./12500,'-','Color',"#77AC30",'DisplayName','Correlation output');
        ylabel('Normalized Amplitude')
        xlabel('Delay [chips]')
        title('Correlation function')
        xlim([-1.2 1.2])
        ylim([0 2])
        hold on;
        % plot initial desired response
        if Filter_ON==1
            p3=plot(time_stamps(:, 1), d(1,:)/val_max,'-','Color',"#D95319",'DisplayName','Desired response');
            p4=plot(time_stamps(:, 1), y(1,:)/val_max,'-','Color',"#7E2F8E",'DisplayName','Filter output');
        end
        % plot initial position of central correlators and filter peak
        p1=plot(Spacing((VE : 1 : VL)),corr_out((VE : 1 : VL))/12500,'*','Color',"#77AC30",'DisplayName','VE-VL');
        if Filter_ON==1
            p2=plot(time_stamps((VE : 1 : VL), 1),d(1,(VE : 1 : VL))/val_max,'*','Color',"#D95319",'DisplayName','VE-VL');
            p5=plot(time_stamps((VE : 1 : VL), 1),y(1,(VE : 1 : VL))/val_max,'*','Color',"#7E2F8E",'DisplayName','VE-VL');
        end
        % plot multipath information
        legend('Location','northwest')
        drawnow
        plot([-1 0],[0 1],'-bo','Color',"#D95319",'DisplayName','LOS');
        legend('Location','northwest')
        plot([0 1],[1 0],'-bo','HandleVisibility','off','Color',"#D95319");
        %for subindex = 1: size(mltpth_delays, 2)
        %    plot([mltpth_delays(subindex) mltpth_delays(subindex)],[1/mltpth_attenuation(subindex) 0],'-bo','HandleVisibility','off');
        %end
        yyaxis right
        ylim([0 10])
        for subindex = 2: corr_per_chip_plot*2
            if(subindex == 2)
                plot([Spacing_plot(subindex) Spacing_plot(subindex)],[tap_weights(subindex) 0],'-o','Color',"#0072BD" ,'DisplayName','Tap weights');
                legend('Location','northwest')
            else
                plot([Spacing_plot(subindex) Spacing_plot(subindex)],[tap_weights(subindex) 0],'-o','Color',"#0072BD" ,'HandleVisibility','off');
            end
        end
        
        drawnow
        
    else
        subplot(1,2,2);
        yyaxis left
        corr_function=plot(Spacing_plot(2 : 1 : corr_per_chip_plot*2), corr_out_plot(2 : 1 : corr_per_chip_plot*2)./12500,'-','Color',"#77AC30",'DisplayName','Correlation output');
        ylabel('Normalized Amplitude')
        xlabel('Delay [chips]')
        title('Correlation function')
        xlim([-1.2 1.2])
        ylim([0 2])
        hold on;
        % plot initial desired response
        if Filter_ON==1
            p3=plot(time_stamps(:, 1), d(1,:)/val_max,'-','Color',"#D95319",'DisplayName','Desired response');
            p4=plot(time_stamps(:, 1), y(1,:)/val_max,'-','Color',"#7E2F8E",'DisplayName','Filter output');
        end
        % plot initial position of central correlators and filter peak
        p1=plot(Spacing((VE : 1 : VL)),corr_out((VE : 1 : VL))/12500,'*','Color',"#77AC30",'DisplayName','VE-VL');
        if Filter_ON==1
            p2=plot(time_stamps((VE : 1 : VL), 1),d(1,(VE : 1 : VL))/val_max,'*','Color',"#D95319",'DisplayName','VE-VL');
            p5=plot(time_stamps((VE : 1 : VL), 1),y(1,(VE : 1 : VL))/val_max,'*','Color',"#7E2F8E",'DisplayName','VE-VL');
        end
        % plot multipath information
        legend('Location','northwest')
        drawnow
        plot([-1 0],[0 1],'-bo','Color',"#D95319",'DisplayName','LOS');
        legend('Location','northwest')
        plot([0 1],[1 0],'-bo','HandleVisibility','off','Color',"#D95319");
        %for subindex = 1: size(mltpth_delays, 2)
        %    plot([mltpth_delays(subindex) mltpth_delays(subindex)],[1/mltpth_attenuation(subindex) 0],'-bo','HandleVisibility','off');
        %end
        yyaxis right
        ylim([0 10])
        for subindex = 2: corr_per_chip_plot*2
            if(subindex == 2)
                plot([Spacing_plot(subindex) Spacing_plot(subindex)],[tap_weights(subindex) 0],'-o','Color',"#0072BD" ,'DisplayName','Tap weights');
                legend('Location','northwest')
            else
                plot([Spacing_plot(subindex) Spacing_plot(subindex)],[tap_weights(subindex) 0],'-o','Color',"#0072BD" ,'HandleVisibility','off');
            end
        end
        
        drawnow
    end

    %{
    %%%% REFRESH PLOTS %%%%
    [val_max, idx_max] = max(corr_out_plot);
    if Filter_ON==1
        % full curve - desired response
        set(p3,'XData',time_stamps(:, 1),'YData',d(Index,:)/val_max); 
        % full curve - filter output
        set(p4,'XData',time_stamps(:, 1),'YData',y(Index,:)/val_max); 
        % central peaks
        set(p2,'XData',time_stamps((VE : 1 : VL), 1),'YData',d(Index,(VE : 1 : VL))/val_max);
        set(p5,'XData',time_stamps((VE : 1 : VL), 1),'YData',y(Index,(VE : 1 : VL))/val_max);
        
    end
    set(p1,'XData',time_stamps((VE : 1 : VL), 1),'YData',corr_out(VE : 1 : VL)/12500);
    set(corr_function,'XData',Spacing_plot(2 : 1 : corr_per_chip_plot*2),'YData',corr_out_plot(2 : 1 : corr_per_chip_plot*2)/12500);
    drawnow
    set(discri,'XData',(1 : 1 : Index),'YData',DLLdiscri(1,(1 : 1 : Index)));
    drawnow
%}
    pause(0.2);
end

%{

%}

function [e] = modified_lms(u, d, it, mu)
    % Maximum number of time step that can be predicted
    N = min(size(u, 1),size(d, 1));
    Nin = size(u,2);
    Nout = size(d,2);
    % Intializatize weight matrix and associated parameters for LMS predictor
    w = ones(1, Nout);
    for n = 1:it
        % Predict next sample and error
        xp(n, :) = u.*w;
        e(n,:) = d-xp(n,:);
        % Adapt weight matrix ans step size
        w = w + (0.0001/norm(u)).*e(n,:).*u;
    end
end


function [xp] = modified_rls(u, d, lambda)
    % Maximum number of time step that can be predicted
    N = min(size(u, 1),size(d, 1));
    Nin = size(u,2);
    Nout = size(d,2);
    % Intializatize weight matrix and associated parameters for RLS predictor
    w = ones(1,Nin);
    % regularization parameter
    delta=1;
    % reset filter variable between monte carlo runs
    P=eye(Nin)*delta;
    for n = 1:N
        xp(n,:) = u(n,:).*w;
        e(n,:) = d(n,:)-xp(n,:);
        % update kappa as perRLS
        kappa = lambda^(-1)*P*u(n,:)'/(1+lambda^(-1)*u(n,:)*P*u(n,:)');
        % update weights
        w = w+kappa'.*e(n,:); 
        % update as per R
        P = lambda^(-1)*P-lambda^(-1)*u(n,:)*kappa*P; 
    end
end