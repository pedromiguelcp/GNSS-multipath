function TrackingPlot_vt(TckResult,Acquired,navSolutions,flag_pt,cnslxyz)

svnum = length(Acquired.sv);

datalength= length(TckResult(Acquired.sv(1)).P_i(:))-1;


if flag_pt == 1
for svindex = 1 : svnum
    figure(svindex);
%     subplot(311);
%     plot(TckResult(Acquired.sv(svindex)).absoluteSample(:));
%     subplot(312);
%     plot(TckResult(Acquired.sv(svindex)).file_ptr(:));
%     subplot(313);
%     plot(TckResult(Acquired.sv(svindex)).absoluteSample(:)-TckResult(Acquired.sv(svindex)).file_ptr(:));
    

    subplot(5,3,1);
    plot(TckResult(Acquired.sv(svindex)).P_i(1:datalength),TckResult(Acquired.sv(svindex)).P_q(1:datalength),'b*');
    title('Discrete-Time Scatter Plot');
    xlabel('I prompt');
    ylabel('Q prompt');
    grid on;
    hold on;
    plot(0,0,'r*')
    
    subplot(5,3,2:3);
    plot(TckResult(Acquired.sv(svindex)).P_i(1:datalength),'r-');
    hold on;
    plot(TckResult(Acquired.sv(svindex)).P_q(1:datalength),'b-');
    title('Ip and Qp of Tracking Result');
    legend('I_P','Q_P');
    grid on;

    subplot(5,3,11:12);
%     plot(TckResult(Acquired.sv(svindex)).carrError(1:datalength),'b*-');
%     ylabel('Amplitude');    
    plot(TckResult(Acquired.sv(svindex)).carrError(1:datalength) .* 2*pi*180/pi,'b*-');
    ylabel('degree');
    title('Raw PLL Discriminator');
    xlabel('Epoch');    
    grid on;

    subplot(5,3,5:6)
    plot(sqrt(TckResult(Acquired.sv(svindex)).E_i(1:datalength).^2+TckResult(Acquired.sv(svindex)).E_q(1:datalength).^2),'b*-');hold on
    plot(sqrt(TckResult(Acquired.sv(svindex)).L_i(1:datalength).^2+TckResult(Acquired.sv(svindex)).L_q(1:datalength).^2),'g*-');hold on
    plot(sqrt(TckResult(Acquired.sv(svindex)).P_i(1:datalength).^2+TckResult(Acquired.sv(svindex)).P_q(1:datalength).^2),'r*-');hold on
    legend('Early','Late','Prompt');
    title('Correlation');
%     xlabel('Time');
    grid on;


    subplot(5,3,8:9);
    plot(TckResult(Acquired.sv(svindex)).codeError,'b*-');
    title('Raw DLL Discriminator');
%     xlabel('Time');
    ylabel('Chip');
    grid on;

    subplot(5,3,10);
%     plot(TckResult(Acquired.sv(svindex)).carrPhase(1:datalength));
    plot(TckResult(Acquired.sv(svindex)).codePhase(1:datalength));
    ylabel('chip');
    title('Code Phase')
    xlabel('Time');    
    grid on;     
    
    subplot(5,3,13);
%     plot(TckResult(Acquired.sv(svindex)).carrPhase(1:datalength));
    plot(TckResult(Acquired.sv(svindex)).carrPhase(1:datalength).*180/pi);
    ylabel('degree');
    title('Carrier Phase')
    xlabel('Time');    
    grid on;        
    
    subplot(5,3,14:15);
%     plot(TckResult(Acquired.sv(svindex)).codePhase(1:datalength));
%     title('Code Phase')
%     xlabel('Time');
%     ylabel('Chip');
    plot(sqrt(TckResult(Acquired.sv(svindex)).sv_vel(1:datalength,1).^2+TckResult(Acquired.sv(svindex)).sv_vel(1:datalength,2).^2+TckResult(Acquired.sv(svindex)).sv_vel(1:datalength,3).^2));
    title('SV Velocity')
    xlabel('Epoch');
    ylabel('m/s');    
    grid on;    
    
    subplot(5,3,4);
    plot(TckResult(Acquired.sv(svindex)).carrFreq(1:datalength) - 4123968);    
    title('Carr Freq');
    xlabel('Epoch');
    ylabel('Hz');
    grid on;
          
    subplot(5,3,7);
    plot(TckResult(Acquired.sv(svindex)).codeFreq(1:datalength) - 1.023e6);    
    title('Code Freq');
    xlabel('Epoch');
    ylabel('Hz');
    grid on;    
    
end % end for
end

figure(100);
subplot(4,3,1:6);
plot(navSolutions.usrPosENU(1:datalength,1),navSolutions.usrPosENU(1:datalength,2),'b*');
title('User EN Position');
xlabel('East(m)');
ylabel('North(m)');
grid on;
subplot(4,3,7:9);
plot(navSolutions.usrPosENU(1:datalength,3));
xlabel('epoch(1ms)');
ylabel('Up(m)');
subplot(4,3,10:12);
plot(navSolutions.clkBias(1:datalength));
xlabel('epoch(1ms)');
ylabel('Epoch(m)');


figure(101);
subplot(411);
plot(navSolutions.usrVel(1:datalength,1));
title('X direction');
ylabel('m/s');
grid on;
subplot(412);
plot(navSolutions.usrVel(1:datalength,2));
title('Y direction');
ylabel('m/s');
grid on;
subplot(413);
plot(navSolutions.usrVel(1:datalength,3));
title('Z direction');
ylabel('m/s');
grid on;
subplot(414);
plot(navSolutions.clkBias(1:datalength));
title('Time direction');
ylabel('m/s');
grid on;

figure(102);
subplot(411);
plot(navSolutions.usrPos(1:datalength,1)-cnslxyz(1));
title('X direction');
ylabel('m');
grid on;
subplot(412);
plot(navSolutions.usrPos(1:datalength,2)-cnslxyz(2));
title('Y direction');
ylabel('m');
grid on;
subplot(413);
plot(navSolutions.usrPos(1:datalength,3)-cnslxyz(3));
title('Z direction');
ylabel('m');
grid on;
subplot(414);
plot(navSolutions.clkBias(1:datalength));
title('Time direction');
ylabel('m');
grid on;

figure(103);
subplot(411);
plot(navSolutions.usrPosENU(1:datalength,1));
title('E direction');
ylabel('m');
grid on;
subplot(412);
plot(navSolutions.usrPosENU(1:datalength,2));
title('N direction');
ylabel('m');
grid on;
subplot(413);
plot(navSolutions.usrPosENU(1:datalength,3));
title('U direction');
ylabel('m');
grid on;
subplot(414);
plot(navSolutions.clkBias(1:datalength));
title('Time direction');
ylabel('m');
grid on;


[cnslllh] = xyz2llh(cnslxyz);
L_b = cnslllh(1);
lamda_b = cnslllh(2);
% Cen = [-sin(L_b)*cos(lamda_b)          -sin(lamda_b) -cos(L_b)*sin(lamda_b);...
%        -sin(L_b)*sin(lamda_b)           cos(lamda_b) -cos(L_b)*sin(lamda_b);...
%                  cos(lamda_b)                      0              -sin(L_b);];

% Cne = [-sin(L_b)*cos(lamda_b) -sin(L_b)*sin(lamda_b)               cos(L_b);...
%                 -sin(lamda_b)           cos(lamda_b)                      0;...
%        -cos(L_b)*cos(lamda_b) -cos(L_b)*sin(lamda_b)              -sin(L_b);];

      
% Cne = [         -sin(lamda_b)           cos(lamda_b)                      0;...
%        -sin(L_b)*sin(lamda_b) -sin(L_b)*sin(lamda_b)               cos(L_b);...
%         cos(L_b)*cos(lamda_b)  cos(L_b)*sin(lamda_b)               sin(L_b);];
%    
% for idx = 1 : datalength
%     usr_vel_enu(idx,:) =  -1*Cne * navSolutions.usrVel(idx,:)';
% end
figure(104);
subplot(411);
plot(-1*navSolutions.usrVelENU(1:datalength,1));
% plot(usr_vel_enu(1:datalength,1));
title('E direction');
ylabel('m/s');
grid on;
subplot(412);
plot(-1*navSolutions.usrVelENU(1:datalength,2));
% plot(usr_vel_enu(1:datalength,2));
title('N direction');
ylabel('m/s');
grid on;
subplot(413);
plot(-1*navSolutions.usrVelENU(1:datalength,3));
% plot(usr_vel_enu(1:datalength,3));
title('U direction');
ylabel('m/s');
grid on;
subplot(414);
plot(navSolutions.clkBias(1:datalength));
title('Time direction');
ylabel('m/s');
grid on;   


end