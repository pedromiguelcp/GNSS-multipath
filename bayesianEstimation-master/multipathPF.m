close all; 
clear all; 
clc;

%Number of particles
n_part = 100;

%Parameters for the simulation
duration = 1; %seconds
dt = 0.01; %seconds
n_iter = floor( duration/dt );
t = dt*(0:(n_iter-1));

%Noise
Rww_fil = 0 ; %Process noise
Rvv_fil = 0.01; %Measurement noise

% Initial covariance of state error estimation
ambiguity = 0.5; % chips
covariance = zeros(n_iter,1);
covariance(1) = ambiguity^2;

%Resampling
N_t = n_part;

% (1) Initialize particles from -0.5 to 0.5 (ambiguity)
particle = -0.5 + rand(n_part, 1);
particle_pred = zeros(n_part,1);
weight = ones(n_part,1)/n_part;

%Estimation
x_est_bpf = zeros(n_iter,1);
x_est_bpf(1) = mean(particle);

%True state (debug)
y_true = zeros(n_iter,1);
y = zeros(n_iter,1);
for k=1:n_iter
    y_true(k) = k*0.01;
    y(k) = y_true(k) + sqrt(Rvv_fil)*randn;
end



for k=2:n_iter
    for i=1:n_part
        % (2) Importance sampling
        particle(i) = particle_pred(i) + sqrt(covariance(k-1)) * randn;
        
        % (3) Weight
        % (3.1) Weight update
        y(k) = corr_out';
        Z_K_predict = X_k(2:L_ekf+2,1)'*fi_ss';
        innov = y(k) - particle(i);
        weight(i) = exp( -log(sqrt(2*pi*Rvv_fil)) -(( innov )^2)/(2*Rvv_fil) );
    end

    
    % (3.2) Weigth normalization
    weight = weight/sum(weight);

    
    % (4) Estimation
    % (4.1) Maximum a Posteriori delay estimation
    [val,idx] = max(weight);
    x_est_bpf(k) = particle(idx);
    % (4.2) Update covariances of state error estimation
    for i=1:n_part
        covariance(k) = covariance(k) + weight(i)*(particle(i)-x_est_bpf(k))*(particle(i)-x_est_bpf(k));
    end


    % (5) Resampling
    for i=1:n_part
        particle_pred(i) = x_est_bpf(k);
    end
end

disp(x_est_bpf(n_iter));

figure
subplot(3,1,1)
plot(t,y_true);
hold on
plot(t,x_est_bpf);
grid, xlabel('Time [s]'), ylabel('State');
%title(sprintf('State Estimation with $$x_0=$$ %1.1f and $$P_0=$$ %1.0d',x0, P0));
title('State Estimation');
legend('True State','BPF Estimate');


subplot(3,1,2)
error = y_true-x_est_bpf;
plot(t,error);
grid, xlabel('Time [s]'), ylabel('Error');
title('State Estimation Error');

subplot(3,1,3)
plot(t,covariance);
grid, xlabel('Time [s]'), ylabel('Covariance');
title('State Estimation Error Covariance');