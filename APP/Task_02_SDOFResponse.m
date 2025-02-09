close all 
clc
clear

% Parameters for white noise generation
dt = 0.02;      % Time step
T = 10.0;       % Total duration
F_std = 3.0;    % Standard deviation of white noise

% Parameters for Bouc-Wen SDOF system
m = 1;          % Mass
c = 0.02;       % Damping
k = 1;          % Stiffness
a = 0.5;        % Post-yield ratio
A = 1;          % Bouc-Wen parameter
beta = 0.5;     % Bouc-Wen parameter
gamma = 0.5;    % Bouc-Wen parameter
n = 5;          % Bouc-Wen parameter

% Number of samples
num_samples = 50;

% Initialize arrays to store responses and forces
u_samples = cell(num_samples, 1); % Store displacement responses
t_samples = cell(num_samples, 1); % Store time vectors
F_samples = cell(num_samples, 1); % Store force vectors

% Loop to generate samples and solve the system
for i = 1:num_samples
    % Generate white noise excitation
    whiteNoiseGenerator = WhiteNoiseLoadGeneratorModel(dt, T, F_std);
    F = whiteNoiseGenerator.F; % White noise excitation
    
    % Store the force
    F_samples{i} = F;
    
    % Create an instance of the BoucWenSDOFModel class
    boucWenSystem = BoucWenSDOFModel(m, c, k, a, A, beta, gamma, n, dt, T, F);
    
    % Simulate the system
    boucWenSystem.simulate();
    
    % Store the response and time vector
    u_samples{i} = boucWenSystem.u;
    t_samples{i} = boucWenSystem.t;
end

% Plot 6 subplots in 2 columns: Force F(t) on the left, Response u(t) on the right
figure;
for i = 1:3
    % Plot Force F(t) on the left column
    subplot(3, 2, 2*i-1); % Left column (1st, 3rd, 5th subplot)
    plot(t_samples{i}, F_samples{i}, 'r', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel(['Force F_' num2str(i) '(t)']);
    title(['Force Sample ' num2str(i)]);
    grid on;
    
    % Plot Response u(t) on the right column
    subplot(3, 2, 2*i); % Right column (2nd, 4th, 6th subplot)
    plot(t_samples{i}, u_samples{i}, 'b', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel(['Displacement u_' num2str(i) '(t)']);
    title(['Response Sample ' num2str(i)]);
    grid on;
end

% Get current date in 'yyyymmdd_HHMM' format
dateString = datestr(now, 'yyyymmdd_HHMM');

% Construct the base file name using the parameters
Ofile = sprintf('%s - Task02_SDOF_RK4_Response', dateString);

% Get the project's root directory and move one step back
currentFile = mfilename('fullpath'); % Full path of the current script/function
projectRoot = fileparts(fileparts(currentFile)); % Move one step back

% Define the 'Data' folder path within the project
dataFolder = fullfile(projectRoot, 'Data');

% Ensure the 'Data' folder exists
if ~isfolder(dataFolder)
    mkdir(dataFolder);
end

% Define the subfolder path inside the 'Data' folder using Ofile
subFolder = fullfile(dataFolder, Ofile);

% Ensure the subfolder exists
if ~isfolder(subFolder)
    mkdir(subFolder);
end

% Create full file paths for .fig and .png
figName = fullfile(subFolder, sprintf('%s.fig', Ofile));
pngName = fullfile(subFolder, sprintf('%s.png', Ofile));

% Save the current figure as .fig file
savefig(gcf, figName);

% Export the current figure as .png with high resolution
exportgraphics(gcf, pngName, 'Resolution', 600);