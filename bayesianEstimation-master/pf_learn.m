% Particle Filter Demo with Plots

% Parameters
numParticles = 100;    % Number of particles
numIterations = 10;    % Number of iterations
motionModelNoise = 0.1;    % Noise in motion model
measurementNoise = 0.5;    % Noise in measurement model

% True state trajectory (for demonstration purposes)
trueStates = [0:0.1:1, 1:-0.1:0, zeros(1,9)];

% Initialize particles
particles = struct();
particles.state = rand(1, numParticles);  % Random initial states
particles.weight = ones(1, numParticles) / numParticles;  % Uniform weights

% Initialize plots
figure;
statePlot = subplot(2, 1, 1);
weightsPlot = subplot(2, 1, 2);

% Main loop
for t = 1:numIterations
    % Motion update (random walk)
    particles.state = particles.state + motionModelNoise * randn(1, numParticles);
    
    % Measurement update (simulated measurement)
    measurements = trueStates(t) + measurementNoise * randn(1, numParticles);
    likelihoods = normpdf(measurements, particles.state, measurementNoise);
    
    % Update weights
    particles.weight = particles.weight .* likelihoods;
    particles.weight = particles.weight / sum(particles.weight);
    
    % Resampling
    resampledIndices = randsample(1:numParticles, numParticles, true, particles.weight);
    particles.state = particles.state(resampledIndices);
    particles.weight = ones(1, numParticles) / numParticles;
    
    % Estimate state (mean of particles)
    estimatedState = mean(particles.state);
    
    % Display current iteration information
    fprintf('Iteration: %d, True State: %.2f, Estimated State: %.2f\n', t, trueStates(t), estimatedState);
    
    % Update plots
    plot(statePlot, 1:t, trueStates(1:t), 'b-', 'LineWidth', 2); hold on;
    plot(statePlot, 1:t, ones(1, t) * estimatedState, 'r--', 'LineWidth', 2);
    xlim(statePlot, [1, numIterations]);
    xlabel(statePlot, 'Iteration');
    ylabel(statePlot, 'State');
    title(statePlot, 'True State vs Estimated State');
    legend(statePlot, 'True State', 'Estimated State');
    
    bar(weightsPlot, particles.state, particles.weight);
    xlim(weightsPlot, [0, 1]);
    xlabel(weightsPlot, 'State');
    ylabel(weightsPlot, 'Weight');
    title(weightsPlot, 'Particle Weights');
    
    drawnow;
    hold off;
    
    % Pause to observe the plots (optional, adjust as needed)
    pause(0.5);
end