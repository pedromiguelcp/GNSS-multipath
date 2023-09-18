% the ground truth has an update rate of 5 HZ (time increments of 0.2 s)
% but the receiver0 output has an update rate of 1kHZ so we interpolate
% time is repeated 13 times (maybe because of higher update rate of other
% sensors)
load("groundTruth\llh_NovAtel.mat");
load("CT_pos_RLS.mat");
[x,y]=unique(llh(:,1));
lat_interp = interp1(x,llh(y,14),navSolutionsCT.usrTime);
long_interp = interp1(x,llh(y,15),navSolutionsCT.usrTime);
heigh_interp = interp1(x,llh(y,16),navSolutionsCT.usrTime);

% see horizontal velocity like the article figure
% matches the zone B & C of the kinematic test
figure;
plot(sqrt(navSolutionsCT.usrVel(:,1).^2+navSolutionsCT.usrVel(:,2).^2));
ylabel('meters/s');
xlabel('Epoch (ms)');
title('Receiver Horizontal velocity');
% see the horizontal positioning error
receiverXYZ = llh2xyz(navSolutionsCT.usrPosLLH);
groundTruthXYZ = llh2xyz([lat_interp,long_interp,heigh_interp]);

usrPosENU = zeros(size(groundTruthXYZ,1),3);
for epoch = 1 : size(groundTruthXYZ,1)
    usrPosENU(epoch,:)  = xyz2enu(groundTruthXYZ(epoch,:),receiverXYZ(epoch,:))';
end
usrPosError3D = sqrt(usrPosENU(:,1).^2+usrPosENU(:,2).^2+usrPosENU(:,3).^2);
usrPosError2D = sqrt(usrPosENU(:,1).^2+usrPosENU(:,2).^2);

figure;
plot(usrPosError3D,'*');
ylabel('meters');
xlabel('Epoch (ms)');
title('Receiver 3D Position error');
figure;
plot(usrPosError2D,'*');
ylabel('meters');
xlabel('Epoch (ms)');
title('Receiver 2D Position error');
figure;
plot(usrPosError3D(1:75000),'*');
ylabel('meters');
xlabel('Epoch (ms)');
title('Receiver 3D Position error');
figure;
plot(usrPosError2D(1:75000),'*');
ylabel('meters');
xlabel('Epoch (ms)');
title('Receiver 2D Position error');

% difference in errors can be due to the initial interpolation
% or the Ground Truth user position is not only computed using GNSS