classdef WhiteNoiseLoadGeneratorModel < handle
    properties
        dt          % Time step
        T           % Total duration
        t           % Time vector
        F_std       % Standard deviation of the white noise
        F           % White noise excitation
    end
    
    methods
        % Constructor
        function obj = WhiteNoiseLoadGeneratorModel(dt, T, F_std)
            % Initialize properties
            obj.dt = dt;
            obj.T = T;
            obj.t = 0:dt:T; % Time vector
            obj.F_std = F_std;
            
            % Generate white noise excitation
            obj.generateWhiteNoise();
        end
        
        % Method to generate white noise
        function generateWhiteNoise(obj)
            obj.F = obj.F_std * randn(size(obj.t)); % White noise with mean 0 and std F_std
        end
        
        % Method to plot the white noise excitation
        function plotForce(obj)
            figure;
            plot(obj.t, obj.F, 'b', 'LineWidth', 1.5);
            xlabel('Time (s)');
            ylabel('Force F(t)');
            title('White Noise Excitation F(t)');
            grid on;
        end

        function saveFigure(obj)
            % Get current date in 'yyyymmdd_HHMM' format
            dateString = datestr(now, 'yyyymmdd_HHMM');
            
            % Construct the base file name using the parameters
            Ofile = sprintf('%s - Task01_White_Noise_Force', dateString);
            
            % Get the project's root directory and move one step back
            currentFile = mfilename('fullpath'); % Full path of the current script/function
            projectRoot = fileparts(fileparts(fileparts(currentFile))); % Move one step back

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
        end
    end
end