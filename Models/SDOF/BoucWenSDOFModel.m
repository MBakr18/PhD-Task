classdef BoucWenSDOFModel < handle
    properties
        % System parameters
        m       % Mass
        c       % Damping
        k       % Stiffness
        a       % Post-yield ratio
        A       % Bouc-Wen parameter
        beta    % Bouc-Wen parameter
        gamma   % Bouc-Wen parameter
        n       % Bouc-Wen parameter
        
        % Time parameters
        dt      % Time step
        T       % Total duration
        t       % Time vector
        
        % External force
        F       % White noise excitation
        
        % State variables
        u       % Displacement
        v       % Velocity
        z       % Hidden state
    end
    
    methods
        % Constructor
        function obj = BoucWenSDOFModel(m, c, k, a, A, beta, gamma, n, dt, T, F)
            % Initialize system parameters
            obj.m = m;
            obj.c = c;
            obj.k = k;
            obj.a = a;
            obj.A = A;
            obj.beta = beta;
            obj.gamma = gamma;
            obj.n = n;
            
            % Initialize time parameters
            obj.dt = dt;
            obj.T = T;
            obj.t = 0:dt:T; % Time vector
            
            % Initialize external force
            obj.F = F; % White noise excitation
            
            % Initialize state variables
            obj.u = zeros(size(obj.t)); % Displacement
            obj.v = zeros(size(obj.t)); % Velocity
            obj.z = zeros(size(obj.t)); % Hidden state
        end
        
        % Method to define the system of ODEs for ode45
        function dydt = boucWenODEs(obj, t, y)
            % Unpack state variables
            u = y(1); % Displacement
            v = y(2); % Velocity
            z = y(3); % Hidden state
            
            % Interpolate the external force F(t) at the current time t
            F_t = interp1(obj.t, obj.F, t);
            
            % Compute derivatives
            du_dt = v;
            dz_dt = obj.A * v - obj.beta * abs(v) * abs(z)^(obj.n-1) * z - obj.gamma * v * abs(z)^obj.n;
            dv_dt = (F_t - obj.c * v - obj.a * obj.k * u - (1 - obj.a) * obj.k * z) / obj.m;
            
            % Pack derivatives
            dydt = [du_dt; dv_dt; dz_dt];
        end
        
        % Method to simulate the Bouc-Wen system using ode45
        function simulate(obj)
            % Initial conditions
            y0 = [obj.u(1); obj.v(1); obj.z(1)]; % [u(0); v(0); z(0)]
            
            % Solve the system using ode45
            [t_sol, y_sol] = ode45(@(t, y) obj.boucWenODEs(t, y), obj.t, y0);
            
            % Store results
            obj.u = y_sol(:, 1); % Displacement
            obj.v = y_sol(:, 2); % Velocity
            obj.z = y_sol(:, 3); % Hidden state
            obj.t = t_sol;       % Time vector
        end
        
        % Method to plot the results
        function plotResults(obj)
            figure;
            
            % Plot displacement
            subplot(3, 1, 1);
            plot(obj.t, obj.u, 'b', 'LineWidth', 1.5);
            xlabel('Time (s)');
            ylabel('Displacement u(t)');
            title('Displacement Response');
            grid on;
            
            % Plot velocity
            subplot(3, 1, 2);
            plot(obj.t, obj.v, 'r', 'LineWidth', 1.5);
            xlabel('Time (s)');
            ylabel('Velocity v(t)');
            title('Velocity Response');
            grid on;
            
            % Plot hidden state
            subplot(3, 1, 3);
            plot(obj.t, obj.z, 'g', 'LineWidth', 1.5);
            xlabel('Time (s)');
            ylabel('Hidden State z(t)');
            title('Hidden State Evolution');
            grid on;
        end

        function saveFigure(obj)
            % Get current date in 'yyyymmdd' format
            dateString = datestr(now, 'yyyymmdd_HHMM');
            
            % Retrieve parameter names from object properties
            firstParameter = obj.a; % Check property name for typos if errors occur
            
            % Construct the base file name using the parameters
            Ofile = sprintf('SDOF_%s_RK4_Response', firstParameter);
            
            % Get the project's root directory
            projectRoot = fileparts(mfilename('fullpath')); % Assumes this function is in the project root
            
            % Define the 'Data' folder path within the project
            dataFolder = fullfile(projectRoot, 'Data');
            
            % Ensure the 'Data' folder exists
            if ~isfolder(dataFolder)
                mkdir(dataFolder);
            end
            
            % Create full file paths for .fig and .png
            figName = fullfile(dataFolder, sprintf('%s-%s.fig', Ofile, dateString));
            pngName = fullfile(dataFolder, sprintf('%s-%s.png', Ofile, dateString));
            
            % Save the current figure as .fig file
            savefig(gcf, figName);
            
            % Export the current figure as .png with high resolution
            exportgraphics(gcf, pngName, 'Resolution', 600);
        end
    end
end