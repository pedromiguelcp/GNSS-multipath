function TrackingPlot(TckResult,prnList, acquired)
    for i=1:1
        prnNumber = prnList(i);
        figure(prnNumber);

        subplot(4,3,1);
        plot(TckResult(prnNumber).P_i(:),TckResult(prnNumber).P_q(:),'b*');
        title('Discrete-Time Scatter Plot');
        xlabel('I prompt');
        ylabel('Q prompt');
        grid on;

        subplot(4,3,2:3);
        plot(TckResult(prnNumber).P_i(:),'r-');
        hold on;
        plot(TckResult(prnNumber).P_q(:),'b-');
        title('Ip and Qp of Tracking Result');
        legend('I_P','Q_P');
        grid on;

        subplot(4,3,4);
        plot(TckResult(prnNumber).carrError(:),'b-');
        title('PLL Discriminator');
        xlabel('Time(ms)');
        ylabel('Amplitude');
        grid on;

        subplot(4,3,5:6)
        plot(sqrt(TckResult(prnNumber).E_i(:).^2+TckResult(prnNumber).E_q(:).^2),'b*-');hold on
        plot(sqrt(TckResult(prnNumber).L_i(:).^2+TckResult(prnNumber).L_q(:).^2),'g*-');hold on
        plot(sqrt(TckResult(prnNumber).P_i(:).^2+TckResult(prnNumber).P_q(:).^2),'r*-');hold on
        legend('Early','Late','Prompt');
        title('Correlation');
        xlabel('Time(ms)');
        grid on;

        subplot(4,3,7);
        plot(TckResult(prnNumber).carrFreq,'b');
        title('Carrier frequency');
        xlabel('Time(ms)');
        ylabel('Hz');
        grid on;

        subplot(4,3,8);
        plot(TckResult(prnNumber).codeError,'r-');
        title('DLL Discriminator');
        xlabel('Time(ms)');
        ylabel('Amplitude');
        grid on;

        subplot(4,3,9);
        plot(TckResult(prnNumber).codedelay,'r*-');
        title('Code Delay Sample');
        xlabel('Time(ms)');
        ylabel('Sample');
        grid on;

        RawNavigationData = TckResult(prnNumber).P_i(:);
        NaviData(find(RawNavigationData>=0)) = 1;
        NaviData(find(RawNavigationData<0)) = -1; 
        datalength = length(TckResult(prnNumber).P_i);
        subplot(4,3,10:12);
        stairs(NaviData);
        axis([1 datalength -1.2 1.2])
        title('Navigation Data');
        xlabel('Time(ms)');
        grid on; 
         











        plot_size=5000;
        %+/- 5kHz possible variation
        pll_correction = zeros(1,plot_size);
        for it=1:plot_size
            pll_correction(it) = TckResult(prnNumber).carrFreq(it) - 6.5e6;
        end

        figure(prnNumber+1);
        subplot(5,1,1);
        plot(pll_correction(1:plot_size),'b');
        title("Carrier frequency variation from IF (6.5 MHz)");
        xlabel('Time(ms)');
        ylabel('Hz');
        grid on;



        N=TckResult(prnNumber).pseudorange(1) / (299792458/1575.42e6);
        ambN = zeros(1,plot_size);
        for it=1:plot_size
            ambN(it) = N;
        end

        subplot(5,1,2);
        plot(ambN(1:plot_size),'b');
        title("N");
        xlabel('Time(ms)');
        ylabel('Phase cycles ambiguity');
        grid on;




        beta = zeros(1,plot_size);
        for it=1:plot_size
            beta(it) = TckResult(prnNumber).carrNco(it)*0.001;%freq to phase cycles
        end        

        subplot(5,1,3);
        plot(beta(1:plot_size),'b');
        title('Beta evolution');
        xlabel('Time(ms)');
        ylabel('Phase cycles');
        grid on;




        subplot(5,1,4);
        plot(TckResult(prnNumber).pseudorange(1:plot_size),'b');
        title('code measurement');
        xlabel('Time(ms)');
        ylabel('meters');
        grid on;



        phase_mesurements = zeros(1,plot_size);
        for it=1:plot_size
            if it==1
                %phase_mesurements(it) = TckResult(prnNumber).pseudorange(1) - beta(it)*(299792458/TckResult(prnNumber).carrFreq(it));
                %phase_mesurements(it) = TckResult(prnNumber).pseudorange(1) - (TckResult(prnNumber).carrFreq(it)-6.5e6)*0.001*(299792458/6.5e6);
                phase_mesurements(it) = TckResult(prnNumber).pseudorange(1)/(299792458/1575.42e6) + (TckResult(prnNumber).carrFreq(it)-6.5e6)*0.001*(1575.42e6/TckResult(prnNumber).carrFreq(it));
            else
                %phase_mesurements(it) = phase_mesurements(it-1) - (beta(it)-beta(it-1))*(299792458/TckResult(prnNumber).carrFreq(it)); 
                %phase_mesurements(it) = phase_mesurements(it-1) - (TckResult(prnNumber).carrFreq(it)-6.5e6)*0.001*(299792458/6.5e6);
                phase_mesurements(it) = phase_mesurements(it-1) + (TckResult(prnNumber).carrFreq(it)-6.5e6)*0.001*(1575.42e6/TckResult(prnNumber).carrFreq(it));
            end        
        end
        subplot(5,1,5);
        plot(phase_mesurements(1:plot_size),'b');
        title('phase cycles measurement');
        xlabel('Time(ms)');
        ylabel('meters');
        grid on;




        figure(prnNumber+2);
        code_minus_phase =  zeros(1,plot_size);
        for it=1:plot_size
            code_minus_phase(it) = TckResult(prnNumber).pseudorange(it) - phase_mesurements(it)*(299792458/1575.42e6);     
        end

        subplot(1,1,1);
        plot(code_minus_phase(1:plot_size),'b');
        title('code minus phase');
        xlabel('Time(ms)');
        ylabel('meters');
        grid on;

    end