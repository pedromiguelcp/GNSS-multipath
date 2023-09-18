% Particle Filter Implementation

% Define the number of particles
numParticles = 100;

% Initialize particles with random positions
particles = rand(numParticles, 2); % Assuming a 2D state space

% Define the motion model parameters
motionNoise = 0.1; % Motion noise standard deviation

% Define the measurement model parameters
measurementNoise = 0.2; % Measurement noise standard deviation

% Generate a true state trajectory (for demonstration purposes)
trueStates = generateTrueStates();

% Perform the particle filter updates
for t = 1:numel(trueStates)
    % Motion update
    particles = motionUpdate(particles, motionNoise);
    
    % Measurement update
    weights = measurementUpdate(particles, trueStates(t), measurementNoise);
    
    % Resampling
    particles = resampling(particles, weights);
    
    % Estimate the current state
    currentState = estimateState(particles, weights);
    
    % Display the results
    disp(['True state: ', num2str(trueStates(t)), ' Estimated state: ', num2str(currentState)]);
end

% Helper functions

function states = generateTrueStates()
    % Generate a true state trajectory (for demonstration purposes)
    states = [1 2 3 4 5]; % Example true state trajectory
end

function particles = motionUpdate(particles, motionNoise)
    % Perform the motion update step by adding random noise to particles
    numParticles = size(particles, 1);
    particles = particles + randn(numParticles, 2) * motionNoise; % Assuming a 2D state space
end

function weights = measurementUpdate(particles, trueState, measurementNoise)
    % Compute the likelihood of each particle given the measurement
    % Here, we assume a simple measurement model where the measurement is the same as the true state
    % with additive Gaussian noise
    distances = sqrt(sum((particles - trueState).^2, 2));
    likelihoods = exp(-0.5 * (distances.^2) / (measurementNoise^2));
    
    % Normalize the likelihoods to obtain the weights
    weights = likelihoods / sum(likelihoods);
end

function particles = resampling(particles, weights)
    % Perform resampling to select a new set of particles based on the weights
    
    numParticles = size(particles, 1);
    
    % Compute cumulative sum of weights
    cumulativeWeights = cumsum(weights);
    
    % Generate random indices for resampling
    resamplingIndices = rand(numParticles, 1);
    
    % Initialize variables for resampled particles and index
    resampledParticles = zeros(size(particles));
    index = 1;
    
    % Perform resampling
    for i = 1:numParticles
        while cumulativeWeights(index) < resamplingIndices(i)
            index = index + 1;
        end
        resampledParticles(i, :) = particles(index, :);
    end
    
    % Assign the resampled particles
    particles = resampledParticles;
end

function currentState = estimateState(particles, weights)
    % Estimate the current state by computing the weighted mean of the particles
    currentState = sum(particles .* weights, 1);
end