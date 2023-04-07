function  TckResultCT = trackingCT(file,signal,track,Acquired)
%Purpose:
%   Perform signal tracking using conventional DLL and PLL
%Inputs:
%	file        - parameters related to the data file to be processed
%	signal      - parameters related to signals,a structure
%	track       - parameters related to signal tracking 
%	Acquired    - acquisition results
%Outputs:
%	TckResultCT	- conventional tracking results, e.g. correlation values in 
%                   inphase prompt (P_i) channel and in qudrature prompt 
%                   channel (P_q), etc.
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
% 
% Written by B. XU and L. T. HSU


%%
Spacing = [-track.CorrelatorSpacing 0 track.CorrelatorSpacing -2];

%% *************** %%
corr_per_chip_plot = 10;
Spacing_plot = zeros(1,corr_per_chip_plot*2+1);
for Index = 1: (corr_per_chip_plot*2+1)
    Spacing_plot(Index) = -1 + (Index-1)/corr_per_chip_plot;
end
%% *************** %%


[tau1code, tau2code] = calcLoopCoef(track.DLLBW,track.DLLDamp,track.DLLGain);
[tau1carr, tau2carr] = calcLoopCoef(track.PLLBW,track.PLLDamp,track.PLLGain);

datalength = track.msEph;
delayValue = zeros(length(Acquired.sv),datalength);
noiseCNOfromSNR = zeros(1,datalength);
CN0fromSNR = zeros(1,datalength);
G = zeros((corr_per_chip_plot*2+1),(corr_per_chip_plot*2+1));
g = zeros((corr_per_chip_plot*2+1),1);

detail_corr_fun = 1;
just_multipath = 0;
no_plots = 0;

%MEDLL
MEDLL_track = 0;
multipaths=0;
val = zeros(1,multipaths+1);
idx = zeros(1,multipaths+1);
convergence_it = 0;

for svindex = 1:length(Acquired.sv)
    remChip = 0;
    remPhase=0;
    remSample = 0;
    carrier_output=0;
    carrier_outputLast=0;
    PLLdiscriLast=0;
    code_output=0;
    code_outputLast=0;
    DLLdiscriLast=0;
    Index = 0;
    AcqDoppler = Acquired.fineFreq(svindex)-signal.IF;
    AcqCodeDelay = Acquired.codedelay(svindex);
    
    Codedelay = AcqCodeDelay;
    codeFreq = signal.codeFreqBasis;
    carrierFreqBasis = Acquired.fineFreq(svindex);
    carrierFreq = Acquired.fineFreq(svindex);
    
    % set the file position indicator according to the acquired code delay
    fseek(file.fid,(signal.Sample-AcqCodeDelay-1+file.skip*signal.Sample)*file.dataPrecision*file.dataType,'bof');  % 
    
    Code = generateCAcode(Acquired.sv(svindex));
    Code = [Code(end-2) Code(end-1) Code(end) Code Code(1)];
    
    h = waitbar(0,['Channel:',num2str(svindex),'  Tracking, please wait...']);
    index_int = 0;
    snrIndex = 1;
    Zk = zeros(1,20);
    CN0_CT = zeros(1,datalength/20);
    CN0_CT_plot = zeros(1,datalength);

    
    for IndexSmall = 1: datalength        
        waitbar(IndexSmall/datalength)
        Index = Index + 1;
        
        
        remSample = ((signal.codelength-remChip) / (codeFreq/signal.Fs));
        numSample = round((signal.codelength-remChip)/(codeFreq/signal.Fs)); 
        delayValue(svindex,IndexSmall) = numSample - signal.Sample;
        
        if file.dataPrecision == 2
            rawsignal = fread(file.fid,numSample*file.dataType,'int16')'; 
            sin_rawsignal = rawsignal(1:2:length(rawsignal));
            cos_rawsignal = rawsignal(2:2:length(rawsignal));
            rawsignal0DC = sin_rawsignal - mean(sin_rawsignal) + 1i*(cos_rawsignal-mean(cos_rawsignal));
        else
            rawsignal0DC = fread(file.fid,numSample*file.dataType,'int8')';  
        end
        
        t_CodeEarly    = (0 + Spacing(1) + remChip) : codeFreq/signal.Fs : ((numSample -1) * (codeFreq/signal.Fs) + Spacing(1) + remChip);
        t_CodePrompt   = (0 + Spacing(2) + remChip) : codeFreq/signal.Fs : ((numSample -1) * (codeFreq/signal.Fs) + Spacing(2) + remChip);
        t_CodeLate     = (0 + Spacing(3) + remChip) : codeFreq/signal.Fs : ((numSample -1) * (codeFreq/signal.Fs) + Spacing(3) + remChip);
        t_CodeVeryEarly     = (0 + Spacing(4) + remChip) : codeFreq/signal.Fs : ((numSample -1) * (codeFreq/signal.Fs) + Spacing(4) + remChip);
        CodeEarly      = Code(ceil(t_CodeEarly) + 3);
        CodePrompt     = Code(ceil(t_CodePrompt) + 3);
        CodeLate       = Code(ceil(t_CodeLate) + 3);
        CodeVeryEarly       = Code(ceil(t_CodeVeryEarly) + 3);

        
        CarrTime = (0 : numSample)./signal.Fs;
        Wave     = (2*pi*(carrierFreq .* CarrTime)) + remPhase ;  
        remPhase =  rem( Wave(numSample+1), 2*pi); 
        carrsig = exp(1i.* Wave(1:numSample));
        InphaseSignal    = imag(rawsignal0DC .* carrsig);
        QuadratureSignal = real(rawsignal0DC .* carrsig);
        
        E_i  = sum(CodeEarly    .*InphaseSignal);  E_q = sum(CodeEarly    .*QuadratureSignal);
        P_i  = sum(CodePrompt   .*InphaseSignal);  P_q = sum(CodePrompt   .*QuadratureSignal);
        L_i  = sum(CodeLate     .*InphaseSignal);  L_q = sum(CodeLate     .*QuadratureSignal);
        VE_i  = sum(CodeLate     .*InphaseSignal);  VE_q = sum(CodeLate     .*QuadratureSignal);




        %% *********** %%
        if (detail_corr_fun == 1)
            if ((abs(P_i) < abs(E_i)) || (abs(P_i) < abs(L_i)) || (IndexSmall == 1) || (just_multipath == 0) )
                time_stamps_plot = zeros(corr_per_chip_plot*2+1,numSample);
                code_replicas_plot = zeros(corr_per_chip_plot*2+1,numSample);
                corr_out_plot = zeros(1,corr_per_chip_plot*2+1);
                for subindex = 1: (corr_per_chip_plot*2)+1
                    time_stamps_plot(subindex,:) = (Spacing_plot(subindex) + remChip) : codeFreq/signal.Fs : ((numSample -1) * (codeFreq/signal.Fs) + Spacing_plot(subindex) + remChip);
                    code_replicas_plot(subindex,:) = Code(ceil(time_stamps_plot(subindex,:)) + 3);
                    corr_out_plot(subindex) = sum(code_replicas_plot(subindex,:)     .*InphaseSignal);
                end
            end
            
            %SNR
            noiseCNOfromSNR(1, Index) = VE_i + VE_q;
            noiseLevel = noiseCNOfromSNR(1, 1:1:Index);
            noiseVariance = sum((noiseLevel-mean(noiseLevel)).^2)/length(noiseLevel); % Variance of noise level
            signalPower = P_i.^2 + P_q.^2; % Signal power
            CN0fromSNR(1, Index)=10*log10(((signalPower)/noiseVariance)/0.001);
            if Index == 2
                CN0fromSNR(1, 1)=CN0fromSNR(1, 2);
            end

            index_int = index_int + 1;
            Zk(1,index_int) = P_i^2 + P_q^2;
            if mod(index_int,20) == 0
                meanZk  = mean(Zk(1,:));
                varZk   = var(Zk(1,:));
                NA2     = sqrt(meanZk^2-varZk);
                varIQ   = 0.5 * (meanZk - NA2);
                CN0_CT(1,snrIndex) =  abs(10*log10(1/(1*1e-3*1) * NA2/(2*varIQ)));
                index_int  = 0;
                snrIndex   = snrIndex + 1;
            end
            if snrIndex > 1
                CN0_CT_plot(1,Index) = CN0_CT(1,snrIndex-1);
            else
                CN0_CT_plot(1,Index) = 0;
            end
            

            LS_G = zeros(numSample, (corr_per_chip_plot*2)+1);
            for subindex = 1: (corr_per_chip_plot*2)+1%remove
                LS_G(:, subindex) = code_replicas_plot(subindex,:);
            end
            LS_D = InphaseSignal';
            LS_H = inv(LS_G'*LS_G)*(LS_G'*LS_D);
            [a_los, tau_los] = max(abs(LS_H));
        end
        %% *********** %%
        %derivatives = zeros(1,corr_per_chip_plot*2);
        %for subindex = 2 : corr_per_chip_plot*2
        %    derivatives(1,subindex) = corr_out_plot(subindex)-corr_out_plot(subindex-1);
        %end
        %figure;
        %plot(derivatives);








        %{
            if MEDLL_track == 1
            for i = 0 : multipaths
    
                while convergence_it < 20
            
                    %step 1
                    corr_out_MEDLL = abs(corr_out_plot);
                    for k = 0 : i-1
                        for corridx = 2: (corr_per_chip_plot*2)
                            cross_corr = sum(code_replicas_plot(corridx,:).*code_replicas_plot(idx(k+1),:));
                            corr_out_MEDLL(corridx) = corr_out_MEDLL(corridx) - cross_corr*val(k+1)/numSample;
                        end
                    end
                    
                    %step 2
                    [val(i+1), idx(i+1)] = max(corr_out_MEDLL);
            
                    %step 3
                    for k = 0 : i-1
                        corr_out_MEDLL = abs(corr_out_plot);
                        for l = 0 : i
                            if l ~= k
                                for corridx = 2: (corr_per_chip_plot*2)
                                    cross_corr = sum(code_replicas_plot(corridx,:).*code_replicas_plot(idx(l+1),:));
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
                    X = sprintf('MEDLL -> Peak %d: Amplitude = %.2f  Delay = %.2f',j, val(j+1)/numSample,(-1+idx(j+1)*1/corr_per_chip_plot));
                    disp(X) 
                end
                
                convergence_it = 0;
            end
        end
%}





        
        remChip   = (t_CodePrompt(numSample) + codeFreq/signal.Fs) - signal.codeFreqBasis*signal.ms;
        
        
        % DLL
        E               = sqrt(E_i^2+E_q^2);
        L               = sqrt(L_i^2+L_q^2);
        DLLdiscri       = 0.5 * (E-L)/(E+L);
        %[a_los, tau_los] = max(abs(corr_out_plot));
        %DLLdiscri       = t_CodePrompt(1) - time_stamps_plot(tau_los);
        
        code_output     = code_outputLast + (tau2code/tau1code)*(DLLdiscri - DLLdiscriLast) + DLLdiscri* (0.001/tau1code);
        DLLdiscriLast   = DLLdiscri;
        code_outputLast = code_output;
        codeFreq        = signal.codeFreqBasis - code_output;
        
        % PLL
        PLLdiscri           = atan(P_q/P_i) / (2*pi);
        carrier_output      = carrier_outputLast + (tau2carr/tau1carr)*(PLLdiscri - PLLdiscriLast) + PLLdiscri * (0.001/tau1carr);
        carrier_outputLast  = carrier_output;  
        PLLdiscriLast       = PLLdiscri;
        carrierFreq         = carrierFreqBasis + carrier_output;  % Modify carrier freq based on NCO command
        
        % Data Record
        TckResultCT(Acquired.sv(svindex)).P_i(Index)            = P_i;
        TckResultCT(Acquired.sv(svindex)).P_q(Index)            = P_q;
        TckResultCT(Acquired.sv(svindex)).E_i(Index)            = E_i;
        TckResultCT(Acquired.sv(svindex)).E_q(Index)            = E_q;
        TckResultCT(Acquired.sv(svindex)).L_i(Index)            = L_i;
        TckResultCT(Acquired.sv(svindex)).L_q(Index)            = L_q;
        TckResultCT(Acquired.sv(svindex)).PLLdiscri(Index)      = PLLdiscri;
        TckResultCT(Acquired.sv(svindex)).DLLdiscri(Index)      = DLLdiscri;
        TckResultCT(Acquired.sv(svindex)).codedelay(Index)      = Codedelay + sum(delayValue(1:Index));
        TckResultCT(Acquired.sv(svindex)).remChip(Index)        = remChip;
        TckResultCT(Acquired.sv(svindex)).codeFreq(Index)       = codeFreq;  
        TckResultCT(Acquired.sv(svindex)).carrierFreq(Index)    = carrierFreq;  
        TckResultCT(Acquired.sv(svindex)).remPhase(Index)       = remPhase;
        TckResultCT(Acquired.sv(svindex)).remSample(Index)      = remSample;
        TckResultCT(Acquired.sv(svindex)).numSample(Index)      = numSample;
        TckResultCT(Acquired.sv(svindex)).delayValue(Index)     = delayValue(svindex,IndexSmall);


        % PLOTS
        %{
        %}
        if (detail_corr_fun==1) && (no_plots==0)
            if(IndexSmall == 1)
         
                % plot initial correlation output
                response = figure(1);
                subplot(3,2,1:2);
                %corr_function=plot(Spacing_plot(2 : 1 : corr_per_chip_plot*2), abs(corr_out_plot(2 : 1 : corr_per_chip_plot*2)),'-','Color',"#77AC30",'DisplayName','Correlation output');
                corr_function=plot(Spacing_plot(2 : 1 : corr_per_chip_plot*2), corr_out_plot(2 : 1 : corr_per_chip_plot*2)/max(abs(corr_out_plot)),'-','Color',"#77AC30",'DisplayName','Correlation output');
                ylabel('Normalized Amplitude')
                xlabel('Delay [chips]')
                title('Correlation function')
                xlim([-1 1])
                %ylim([0 10000])
                hold on;
                % plot initial position of central correlators and filter peak
                %p1=plot(Spacing(1:1:3),[abs(E_i) abs(P_i) abs(L_i)],'*','Color',"#77AC30",'DisplayName','Central correlators (E, P, L)');
                p1=plot(Spacing(1:1:3),[E_i P_i L_i]./max(abs(corr_out_plot)),'*','Color',"#77AC30",'DisplayName','Central correlators (E, P, L)');
                yyaxis right
                ls=plot(Spacing_plot(1 : 1 : corr_per_chip_plot*2+1), LS_H,'-','DisplayName','Least Squares');
                peak_ls=plot([Spacing_plot(tau_los) Spacing_plot(tau_los)], [LS_H(tau_los) 0],'r-','DisplayName','Estimated Peak');
                % plot multipath information
                legend
                drawnow
                subplot(3,2,3);
                dll_discri=plot((1 : 1 : Index), TckResultCT(Acquired.sv(svindex)).DLLdiscri(1 : 1 : Index),'-','Color',"#7E2F8E");
                ylabel('DLL discriminator')
                xlabel('Epochs (ms)')
                title('DLL corrections')
                subplot(3,2,4);
                pll_discri=plot((1 : 1 : Index), TckResultCT(Acquired.sv(svindex)).PLLdiscri(1 : 1 : Index),'-','Color',"#7E2F8E");
                ylabel('PLL discriminator')
                xlabel('Epochs (ms)')
                title('PLL corrections')

                subplot(3,2,5);
                ylabel('Promp Correlator magnitude','Color','b')
                yyaxis left
                i_corr=plot((1 : 1 : Index), TckResultCT(Acquired.sv(svindex)).P_i(1 : 1 : Index),'g-','DisplayName','Ip');
                hold on;
                q_corr=plot((1 : 1 : Index), TckResultCT(Acquired.sv(svindex)).P_q(1 : 1 : Index),'b-','DisplayName','Qp');
                yyaxis right
                ylabel('CN0','Color','y')
                cn0=plot((1 : 1 : Index), CN0fromSNR(1, 1 : 1 : Index),'b-','Color',"#EDB120",'DisplayName','CN0');
                cn01=plot((1 : 1 : Index), CN0_CT_plot(1, 1 : 1 : Index),'b-','Color',"#EDB120");
                xlabel('Epochs (ms)')
                title('In-phase Quatradure magnitude')

                subplot(3,2,6);
                code_delay_plot=plot((1 : 1 : Index), TckResultCT(Acquired.sv(svindex)).codedelay(1 : 1 : Index),'-','Color',"#7E2F8E");
                ylabel('Sample Delay')
                xlabel('Epochs (ms)')
                title('Code Delay')
                
                drawnow
            elseif (just_multipath == 0) || (abs(P_i) < abs(E_i)) || (abs(P_i) < abs(L_i))
                
                %set(p1,'XData',[t_CodeEarly(1) t_CodePrompt(1) t_CodeLate(1)],'YData',[abs(E_i) abs(P_i) abs(L_i)]);
                set(p1,'XData',[t_CodeEarly(1) t_CodePrompt(1) t_CodeLate(1)],'YData',[E_i P_i L_i]./max(abs(corr_out_plot)));
                %set(corr_function,'XData',time_stamps_plot(2 : 1 : corr_per_chip_plot*2,1),'YData',abs(corr_out_plot(2 : 1 : corr_per_chip_plot*2))/max(corr_out_plot));
                set(corr_function,'XData',time_stamps_plot(2 : 1 : corr_per_chip_plot*2,1),'YData',corr_out_plot(2 : 1 : corr_per_chip_plot*2)/max(abs(corr_out_plot)));
                set(ls,'XData',time_stamps_plot(1 : 1 : corr_per_chip_plot*2+1,1),'YData',LS_H);
                set(peak_ls,'XData',[time_stamps_plot(tau_los) time_stamps_plot(tau_los)],'YData',[LS_H(tau_los) 0]);
                
                set(dll_discri,'XData',(1 : 1 : Index),'YData',TckResultCT(Acquired.sv(svindex)).DLLdiscri(1 : 1 : Index));
                
                set(pll_discri,'XData',(1 : 1 : Index),'YData',TckResultCT(Acquired.sv(svindex)).PLLdiscri(1 : 1 : Index));
                drawnow
                set(i_corr,'XData',(1 : 1 : Index),'YData',TckResultCT(Acquired.sv(svindex)).P_i(1 : 1 : Index));
                set(q_corr,'XData',(1 : 1 : Index),'YData',TckResultCT(Acquired.sv(svindex)).P_q(1 : 1 : Index));
                set(cn0,'XData',(1 : 1 : Index),'YData',CN0fromSNR(1, 1 : 1 : Index));
                set(cn01,'XData',(1 : 1 : Index),'YData',CN0_CT_plot(1, 1 : 1 : Index));

                set(code_delay_plot,'XData',(1 : 1 : Index),'YData',TckResultCT(Acquired.sv(svindex)).codedelay(1 : 1 : Index));
                drawnow
            end
        end
        
        disp("Iteration: " + IndexSmall + "  Signal: " + svindex + "/" + length(Acquired.sv));
    end
    close(h);
end % end for