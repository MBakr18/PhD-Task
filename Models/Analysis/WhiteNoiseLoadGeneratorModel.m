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
    end
end