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


for Index = 1: (corr_per_chip*2+1)
    Spacing(Index) = -1 + (Index-1)/corr_per_chip;
end

VE = corr_per_chip-1;
E = corr_per_chip;
P = corr_per_chip+1;
L = corr_per_chip+2;
VL = corr_per_chip+3;

% Specify here the multipath delay (chips) and attenuation (linear factor), compared to the LOS signal
mltpth_delays = [0.1 0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.9];
%mltpth_delays = [0.1 0.25 0.35 0.4 0.6 0.7 0.8];
mltpth_attenuation = [2 2.5 3 2 1 0.5 0.3 4 3];
%mltpth_attenuation = [2 2.5 3 3.5 4 4.5 5];
chan1.delays=mltpth_delays;
chan1.attenuation=mltpth_attenuation;
%Channel 2
mltpth_delays = [0.1 0.25 0.35 0.45 0.55];
mltpth_attenuation = [1.8 2.2 3 3.5 4];
chan2.delays=mltpth_delays;
chan2.attenuation=mltpth_attenuation;
% Channel 3
mltpth_delays = [0.1 0.25 0.35];
mltpth_attenuation = [2 4 8];
chan3.delays=mltpth_delays;
chan3.attenuation=mltpth_attenuation;


epochs = 300; %signal tracking epochs

% Generate the C/A code
Code = generateCAcode(1); % Generate C/A code for satellite 1
Code = [Code(end-1) Code(end) Code Code(1) Code(2) Code(3)];%enable replica slidding

DLLdiscri = zeros(1,epochs);

% RLS data
%weights = eye((corr_per_chip*2+1),(corr_per_chip*2+1));
Path_len=4; % needs to be continuously estimated therefore must be computed in the for cycle
weights = zeros(epochs,Path_len);
channel = zeros(epochs,Path_len); % only for debug purposes
y = zeros(epochs,(corr_per_chip*2+1));
d = zeros(epochs,(corr_per_chip*2+1));
e = zeros(epochs+1,(corr_per_chip*2+1));
%P_rls = eye((corr_per_chip*2+1));
delta=0.005;
P_rls = eye(Path_len)/delta;

% convergence related variables
wlast=weights;

% LOS estimation data
G = zeros((corr_per_chip*2+1),(corr_per_chip*2+1));
g = zeros((corr_per_chip*2+1),1);
amplitude = zeros((corr_per_chip*2+1),1);
a_los = zeros(1,epochs);
tau_los = zeros(1,epochs);
val = zeros(1,epochs);
idx = zeros(1,epochs);

%dynamic channel tracking

dedt_max=0;
lambda_max=0.97;
lambda_min=0.1;
lambda=0.95; %lambda_max; % first time (a value is needed)
lambdaIndex=0;

Filter_ON = 1;

%tracking
channelTracking=0;

sem_rls=0; %rls test flag
for Index = 1: epochs
    mltpth_delays=chan3.delays;
    mltpth_attenuation=chan3.attenuation;
    if channelTracking  %Change channel
        if Index < 51
            mltpth_delays=chan2.delays;
            mltpth_attenuation=chan2.attenuation;
        elseif Index<101
                mltpth_delays=chan2.delays;
                mltpth_attenuation=chan2.attenuation;
        elseif Index<151
                mltpth_delays=chan3.delays;
                mltpth_attenuation=chan3.attenuation;
        elseif Index<201
                mltpth_delays=chan1.delays;
                mltpth_attenuation=chan1.attenuation;
        end
    end



    %numSample = round((codelength-remChip)/(codeFreq/fs)); 
    
    LOS_Code = Spacing(P) : codeFreqBasis/fs : ((numSample -1) * (codeFreqBasis/fs) + Spacing(P));
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
    


    %%%%%%%% ADAPTIVE CHANNEL COMPENSATION %%%%%%%%
      if Index==1
        y(Index,:)=corr_out(:,length(corr_out(1,:)));
    else

       [channel(Index+1,:),weights(Index +1,:), wlast, y(Index,:), P_rls] = filter_rls(corr_out, d(Index-1,:), weights(Index,:)', wlast, P_rls, lambda);
       
    end
  
    %[val(Index), idx(Index)] = max(y(Index,:));
  if sem_rls
      y(Index,:)=corr_out(:,length(corr_out(1,:)));
  end

    %%%%%%%% LOS SIGNAL ESTIMATION BLOCK %%%%%%%%  
    if(Index == 1) %this only needs to be computed once
        % compute known matrix G (LxL)
        for subindex = 1 : (corr_per_chip*2+1)
            for subsubindex = 1 : (corr_per_chip*2+1)
                G(subindex,subsubindex) = sum(code_replicas(subindex,:).*code_replicas(subsubindex,:));
            end
        end
        G = inv(G);
    end
    
    %ML equation - search LOS delay to maximize LOS amplitude - eq42
    for search_index = 1 : (corr_per_chip*2+1)
        for subindex = 1 : (corr_per_chip*2+1)
            g(subindex,1) = sum(code_replicas(search_index,:).*code_replicas(subindex,:));
        end    
        amplitude(search_index,1) = (g'*G*y(Index,:)')/(g'*G*g);
    end  
    % get the max amplitude instead of doing equation again (line 117)
    % equivalent to eq43
    [a_los(Index), tau_los(Index)] = max(amplitude); % the return is the same as the commented line 98
    
    % compute desired response
    for subindex = 1 : (corr_per_chip*2)
        %auto-correlation between LOS replica and the other replicas (perfect triangle)
        d(Index,subindex) = sum(code_replicas(subindex,:).*code_replicas(tau_los(Index),:))*a_los(Index);
    end

    % compute error
    e(Index+1,:) = d(Index,:) - y(Index,:);

     %lambda computation for non-stationary channel tracking
    if Index > 1
     if channelTracking
        if abs(norm(e(Index+1,:))-norm(e(Index,:))) > dedt_max
            dedt_max=abs(norm(e(Index+1,:))-norm(e(Index,:))); %update dedt_max if necessary
        end
        deltae = (norm(e(Index+1,:))-norm(e(Index,:)))/norm(e(Index+1,:)); %delta_dedt, avoids false positives from non-ideal convergence 
         %deltae = norm(e(Index+1,:)-e(Index,:))/norm(e(Index+1,:));

        if abs(deltae) > 0.40  %error increased; channel changed or convergence issues
            lambda1=-((lambda_max-lambda_min)/dedt_max)*abs((norm(e(Index+1,:))-norm(e(Index,:))))+1;
            if lambda1<lambda %avoid successive entries that increase lambda just after a channel change was detected
                lambda=lambda1;
                lambdaIndex=Index;
            end
        else
            if Index > (lambdaIndex +8)
                lambda=lambda_max;
            else
                lambda=min(lambda + 0.1, lambda_max); %gradual inclusion of correlation estimates in the computation of the correlation Matrix estimate
            end
        end
     end
    end

    %DLL
    DLL_E           = sqrt(corr_out(E)^2);
    DLL_L           = sqrt(corr_out(L)^2);

    DLLdiscri(1,Index)       = 0.5 * (DLL_E-DLL_L)/(DLL_E+DLL_L);
    if Filter_ON==1
        newDLLdiscri    = Spacing(P) - Spacing(tau_los(Index)); %difference between delay of prompt and delay of max peak
        DLLdiscri(1,Index)       = newDLLdiscri;
    end
    code_output     = code_outputLast + (0.3749245/0.007030542258775)*(DLLdiscri(1,Index) - DLLdiscriLast) + DLLdiscri(1,Index)* (0.001/0.007030542258775);
    DLLdiscriLast   = DLLdiscri(1,Index);
    code_outputLast = code_output;
    codeFreq        = codeFreqBasis - code_output;


    if(Index == 1) 
        % plot initial correlation output
        response = figure(1);
        subplot(2,2,1);
        plot(Spacing(2 : 1 : corr_per_chip*2), corr_out(2 : 1 : corr_per_chip*2)./numSample,'-','Color',"#77AC30");
        ylim([-0.1 8])
        xlim([-1.5 1.5])
        ylabel('Normalized Amplitude')
        xlabel('Delay [chips]')
        title('Initial/Final correlation result')
        hold on;
        % plot initial position of central correlators and filter peak
        p1=plot(Spacing((VE : 1 : VL)),corr_out((VE : 1 : VL))/numSample,'*','Color',"#77AC30");
        p2=plot(time_stamps((VE : 1 : VL), 1),y(epochs,(VE : 1 : VL))/numSample,'*','Color',"#D95319"); % tirei d coloquei y
        % plot initial desired response
        if Filter_ON==1
            p3=plot(time_stamps(:, 1), y(epochs,:)/numSample,'-','Color',"#D95319"); % tirei d coloquei y
        end
        % plot multipath information
        plot([0 0],[1 0],'-bo');
        for subindex = 1: size(mltpth_delays, 2)
            plot([mltpth_delays(subindex) mltpth_delays(subindex)],[1/mltpth_attenuation(subindex) 0],'-bo');
        end
        drawnow
        subplot(2,2,2);
        discri=plot((1 : 1 : Index), DLLdiscri(1,(1 : 1 : Index)),'-','Color',"#7E2F8E");
        ylabel('DLL discriminator')
        xlabel('Epochs (ms)')
        title('DLL corrections')
        subplot(2,2,3);
        filtw=plot((1 : 1 : Path_len), weights(Index+1,:),'-*','Color',"#7E2F8E");
        ylim([-2.5 2.5])
        ylabel('Magnitude')
        xlabel('W')
        title('Filter weights')
        drawnow
        subplot(2,2,4);
        impresp=plot((1 : 1 : Path_len), channel(Index+1,:),'-*','Color',"#7E2F8E");
        ylim([-0.5 2])
        ylabel('Magnitude')
        xlabel('H')
        title('Inverse Channel')
        drawnow
    end

    %%%% REFRESH PLOTS %%%%
    if Filter_ON==1
        % full curve
        set(p3,'XData',time_stamps(:, 1),'YData',y(Index,:)/numSample); 
        % central peaks
        set(p2,'XData',time_stamps((VE : 1 : VL), 1),'YData',y(Index,(VE : 1 : VL))/numSample);
        
    end
    set(p1,'XData',time_stamps((VE : 1 : VL), 1),'YData',corr_out(VE : 1 : VL)/numSample);
    drawnow
    set(discri,'XData',(1 : 1 : Index),'YData',DLLdiscri(1,(1 : 1 : Index)));
    drawnow
    set(filtw,'XData',(1 : 1 : Path_len),'YData',weights(Index+1,:));
    drawnow
    set(impresp,'XData',(1 : 1 : Path_len),'YData',channel(Index+1,:));
    drawnow
    %pause(0.1);
end

%{

%}
function [h, w, wlast, xp, P_rls] = filter_rls(u, d, w, wlast, P_rls, lambda)
 
   
    xi=d;  %error  
    M=length(w);%number of multipaths
    N=length(d);%number of correlators
    d=d(:),
    
    for i=M:N
        uvec=u(i:-1:i-M+1);
        kappa = lambda^(-1)*P_rls*uvec/(1+lambda^(-1)*uvec'*P_rls*uvec);
        xi(i)=d(i)-w'*uvec;
        w=w+kappa*conj(xi(i));
        P_rls = lambda^(-1)*P_rls-lambda^(-1)*kappa*uvec'*P_rls;
    end
    xp=filter(w,1,u);      %filtered signal
    %xp1=w'*u;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Multipath is just the inverse system
    % computation only for debug purposes
    % the role of u and d must be changed
    delta=0.005;
    P = eye(M)/delta;

    h=zeros(M,1);
    for i=M:N
        uvec=d(i:-1:i-M+1);
        kappa = lambda^(-1)*P*uvec/(1+lambda^(-1)*uvec'*P*uvec);
        xi(i)=u(i)-h'*uvec;
        h=h+kappa*conj(xi(i));
        P_rls = lambda^(-1)*P-lambda^(-1)*kappa*uvec'*P;
    end

    %imp=conv(w,h);%must be (approaches) de impulse function

end





