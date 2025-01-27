% Clear workspace
close all;
clc;
clear;

% Parameters for Bouc-Wen SDOF system
m = 1;          % Mass
c = 0.02;       % Damping
k = 1;          % Stiffness
a = 0.5;        % Post-yield ratio
A = 1;          % Bouc-Wen parameter
beta = 0.5;     % Bouc-Wen parameter
gamma = 0.5;    % Bouc-Wen parameter
n = 5;          % Bouc-Wen parameter

% Time parameters
dt = 0.02;      % Time step
T = 10.0;       % Total duration

% External force parameters
F_std = 3.0;    % Standard deviation of white noise
numSamples = 50; % Number of samples

% Create an instance of the LSTMModel class
lstmModel = LSTMModel(m, c, k, a, A, beta, gamma, n, dt, T, F_std, numSamples);

% Generate data
lstmModel.generateData();

% Train the LSTM model
lstmModel.trainModel();

% Evaluate the LSTM model
lstmModel.evaluateModel();

% Save the figure in the data folder
lstmModel.saveFigure();