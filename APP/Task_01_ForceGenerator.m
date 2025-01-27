close all 
clc
clear

% Parameters
dt = 0.02; % Time step
T = 10.0; % Total duration
F_std = 3.0; % Standard deviation of the white noise

% Create an instance of the WhiteNoiseLoadGeneratorModel class
whiteNoiseGenerator = WhiteNoiseLoadGeneratorModel(dt, T, F_std);

% Plot the generated white noise excitation
whiteNoiseGenerator.plotForce();

% Save the figure in the data folder
whiteNoiseGenerator.saveFigure();