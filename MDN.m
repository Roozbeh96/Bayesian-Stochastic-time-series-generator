clear; clc; close all;

%% Step 1: Generate Synthetic Training Data
n_samples = 1000;
X = linspace(-10, 10, n_samples)'; % Input X values
Y = sin(X) + 0.3 * X + 0.5 * randn(n_samples, 1); % Y = sin(X) + noise

% Normalize X and Y for better training
X = (X - mean(X)) / std(X);
Y = (Y - mean(Y)) / std(Y);

% Plot Initial Data (X, Y)
figure;
scatter(X, Y, 10, 'b', 'filled');
xlabel('X'); ylabel('Y');
title('Scatter Plot of X vs Y (Before Training)');
grid on;

% Define the number of Gaussian components
numGaussians = 3;

%% Step 2: Define the MDN Network
layers = [
    featureInputLayer(1, 'Name', 'input')  % X input layer
    fullyConnectedLayer(64, 'Name', 'fc1')
    reluLayer('Name', 'relu1')
    fullyConnectedLayer(64, 'Name', 'fc2')
    reluLayer('Name', 'relu2')
    
    fullyConnectedLayer(numGaussians * 3, 'Name', 'output') % Output for μ, σ, π
];

mdnNet = dlnetwork(layers);


%% Step 4: Train the MDN Model

numEpochs = 1000;
learningRate = 0.001;
miniBatchSize = 32;

% Convert data to dlarray format for deep learning
X_dl = dlarray(X', 'CB'); % Convert X to dlarray (Channel, Batch)
Y_dl = dlarray(Y', 'CB'); % Convert Y to dlarray (Channel, Batch)

% Training loop
for epoch = 1:numEpochs
    % Compute loss using dlfeval() to enable differentiation
    [loss, gradients] = dlfeval(@computeLoss, mdnNet, X_dl, Y_dl, numGaussians);

    % Update network parameters
    mdnNet = dlupdate(@(W, dW) W - learningRate * dW, mdnNet, gradients);

    % Display loss every 100 epochs
    if mod(epoch, 100) == 0
        disp(['Epoch ', num2str(epoch), ', Loss: ', num2str(extractdata(loss))]);
    end
end





% Generate predictions for test X values
X_test = linspace(-10, 10, 10000)';
X_test = (X_test - mean(X_test)) / std(X_test); % Normalize
X_test_dl = dlarray(X_test', 'CB');
Y_pred = forward(mdnNet, X_test_dl);

% Sample Y from P(Y | X)
Y_samples = sampleFromMDN(Y_pred, numGaussians);

%% Step 6: Plot Results
figure;
scatter(X, Y, 10, 'b', 'filled'); hold on;
scatter(X_test, Y_samples, 10, 'r', 'filled');
xlabel('X');
ylabel('Y');
title('MDN: True Data (Blue) vs Sampled Data (Red)');
grid on;


%% Step 3: Define the Custom Loss Function


function loss = mdnLoss(Y_true, Y_pred, numGaussians)
    % Ensure Y_pred is a dlarray
    Y_pred = dlarray(Y_pred);  

    % Extract components from predicted values
    params = Y_pred;  

    % Get batch size
    batchSize = size(Y_true, 2);

    % Reshape into mixture components (Ensure dlarray)
    pi_raw = exp(params(1:numGaussians, :)); 
    pi = pi_raw ./ sum(pi_raw, 1);  % Normalize mixing coefficients
    
    mu = params(numGaussians + 1 : 2 * numGaussians, :); % Means
    sigma = exp(params(2 * numGaussians + 1 : 3 * numGaussians, :)); % Std-devs

    % Ensure Y_true is properly shaped as a dlarray
    Y_true = reshape(Y_true, 1, batchSize);

    % Compute Gaussian probabilities (Ensure elementwise computations)
    gaussians = (1 ./ (sigma .* sqrt(2 * pi))) .* exp(-0.5 * ((Y_true - mu) ./ sigma) .^ 2);

    % Compute weighted sum of Gaussians
    weighted_gaussians = sum(pi .* gaussians, 1);

    % Compute Negative Log-Likelihood Loss (Ensure scalar dlarray)
    loss = -mean(log(weighted_gaussians + 1e-8), 'all');

    % Ensure loss is a dlarray scalar (needed for dlgradient)
    loss = dlarray(loss);
end




%% Step 5: Sampling from P(Y | X)
function Y_sampled = sampleFromMDN(Y_pred, numGaussians)
    % Extract parameters from predicted values
    params = extractdata(Y_pred);  % Convert dlarray to standard array

    % Extract mixture components
    pi = exp(params(1:numGaussians, :)); % Convert logits to probabilities
    pi = pi ./ sum(pi, 1);  % Normalize so that sum(pi) = 1 per sample
    
    mu = params(numGaussians + 1 : 2 * numGaussians, :); % Means
    sigma = exp(params(2 * numGaussians + 1 : 3 * numGaussians, :)); % Ensure sigma > 0

    % Number of samples (batch size)
    numSamples = size(pi, 2);

    % Sample a Gaussian component per sample
    rand_vals = rand(1, numSamples);  % Generate random values for each sample
    cumsum_pi = cumsum(pi, 1);  % Compute cumulative sum over Gaussians
    component = sum(rand_vals > cumsum_pi, 1) + 1;  % Select Gaussian per sample

    % Extract the corresponding mean and std-dev
    selected_mu = mu(sub2ind(size(mu), component, 1:numSamples));
    selected_sigma = sigma(sub2ind(size(sigma), component, 1:numSamples));

    % Sample from the chosen Gaussian component
    Y_sampled = selected_mu + selected_sigma .* randn(size(selected_mu));
end

%%
function [loss, gradients] = computeLoss(mdnNet, X, Y, numGaussians)
    % Ensure X and Y are dlarray objects
    X = dlarray(X, 'CB');  
    Y = dlarray(Y, 'CB');  

    % Forward pass (Y_pred should be a dlarray)
    Y_pred = forward(mdnNet, X);  

    % Compute MDN loss (Ensure output is dlarray)
    loss = mdnLoss(Y, Y_pred, numGaussians);

    % Compute gradients (Loss must be a traced dlarray)
    gradients = dlgradient(loss, mdnNet.Learnables);
end

