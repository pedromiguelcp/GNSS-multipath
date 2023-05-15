function TrackingPlot3(TckResult,Acquired,navSolutions,flag_pt)

svnum = length(Acquired.sv);

datalength= length(TckResult(Acquired.sv(1)).P_i(:))-1;
% period = 1 : datalength;
% period = 70 : datalength;
period = 1 : datalength;
% % period = 2001 : 2001 + 10000;

if (flag_pt == 1)
for svindex = 1 : svnum
    figure(svindex);

    subplot(5,3,1);
    plot(TckResult(Acquired.sv(svindex)).P_i(period),TckResult(Acquired.sv(svindex)).P_q(period),'b*');
    title('Discrete-Time Scatter Plot');
    xlabel('I prompt');
    ylabel('Q prompt');
    grid on;
    subplot(5,3,2:3);
    plot(TckResult(Acquired.sv(svindex)).P_i(period),'r-');
    hold on;
    plot(TckResult(Acquired.sv(svindex)).P_q(period),'b-');
    title('Ip and Qp of Tracking Result');
    legend('I_P','Q_P');
    grid on;

    subplot(5,3,11:12);
    plot(TckResult(Acquired.sv(svindex)).carrError(period).* 2*pi*180/pi,'b*-');
    title('Raw PLL Discriminator');
%     xlabel('Time(ms)');
    ylabel('Degree');
    grid on;

    subplot(5,3,5:6)
    plot(sqrt(TckResult(Acquired.sv(svindex)).E_i(period).^2+TckResult(Acquired.sv(svindex)).E_q(period).^2),'b*-');hold on
    plot(sqrt(TckResult(Acquired.sv(svindex)).L_i(period).^2+TckResult(Acquired.sv(svindex)).L_q(period).^2),'g*-');hold on
    plot(sqrt(TckResult(Acquired.sv(svindex)).P_i(period).^2+TckResult(Acquired.sv(svindex)).P_q(period).^2),'r*-');hold on
    legend('Early','Late','Prompt');
    title('Correlation');
%     xlabel('Time(ms)');
    grid on;


    subplot(5,3,8:9);
    plot(TckResult(Acquired.sv(svindex)).codeError(period),'b*-');
    title('Raw DLL Discriminator');
%     xlabel('Time(ms)');
    ylabel('Amplitude');
    grid on;

%     subplot(5,3,8);
%     plot(TckResult(Acquired.sv(svindex)).codedelay(period),'r*-');
%     title('Code Delay Sample');
%     ylabel('Sample');
%     grid on;

    subplot(5,3,10);
    plot(TckResult(Acquired.sv(svindex)).remCodePhase(period));
    title('Code Phase')
    xlabel('Time(ms)');
    ylabel('Chip');
    grid on;    
    
    subplot(5,3,4);
    plot(TckResult(Acquired.sv(svindex)).carrFreq(period) - 4123968);    
    title('Carr Freq');
    grid on;
          
    subplot(5,3,7);
    plot(TckResult(Acquired.sv(svindex)).codeFreq(period) - 1.023e6);    
    title('Code Freq');
    grid on;
    
    subplot(5,3,13);
    plot(TckResult(Acquired.sv(svindex)).remCarrPhase(period) .* 180/pi);
    title('Carr Phase')
    xlabel('Time(ms)');
    ylabel('Chip');
    grid on;       
end % end for
end
figure(100);
subplot(4,3,1:6);
plot(navSolutions.usrPosENU(period,1),navSolutions.usrPosENU(period,2),'b*');
xlabel('East(m)');
ylabel('North(m)');
grid on;
subplot(4,3,7:9);
plot(navSolutions.usrPosENU(period,3));
xlabel('epoch(1ms)');
ylabel('Up(m)');
subplot(4,3,10:12);
plot(navSolutions.state(period,7));
xlabel('epoch(1ms)');
ylabel('Time(m)');

grid on;

figure(101);
subplot(411);
plot(navSolutions.usrVel(period,1));
title('East direction');
ylabel('m/s');
grid on;
subplot(412);
plot(navSolutions.usrVel(period,2));
title('North direction');
ylabel('m/s');
grid on;
subplot(413);
plot(navSolutions.usrVel(period,3));
title('Up direction');
ylabel('m/s');
grid on;
subplot(414);
plot(navSolutions.state(period,8));
title('Time direction');
ylabel('m/s');
grid on;


figure(103);
subplot(411);
plot(navSolutions.usrPosENU(period,1));
title('E direction');
ylabel('m');
grid on;
subplot(412);
plot(navSolutions.usrPosENU(period,2));
title('N direction');
ylabel('m');
grid on;
subplot(413);
plot(navSolutions.usrPosENU(period,3));
title('U direction');
ylabel('m');
grid on;
subplot(414);
plot(navSolutions.usr_clk(period));
title('Time direction');
ylabel('m');
grid on;


figure(104);
subplot(411);
plot(-1*navSolutions.usrVelENU(period,1));
% plot(usr_vel_enu(period,1));
title('E direction');
ylabel('m/s');
grid on;
subplot(412);
plot(-1*navSolutions.usrVelENU(period,2));
% plot(usr_vel_enu(period,2));
title('N direction');
ylabel('m/s');
grid on;
subplot(413);
plot(-1*navSolutions.usrVelENU(period,3));
% plot(usr_vel_enu(period,3));
title('U direction');
ylabel('m/s');
grid on;
subplot(414);
plot(navSolutions.clkBias(period));
title('Time direction');
ylabel('m/s');
grid on;   

end