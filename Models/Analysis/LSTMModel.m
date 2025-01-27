classdef LSTMModel < handle
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
        
        % External force parameters
        F_std   % Standard deviation of white noise
        numSamples % Number of samples
        
        % Data
        inputData % Input data (force F(t))
        outputData % Output data (displacement u(t))
        
        % LSTM model
        net     % Trained LSTM model
    end
    
    methods
        % Constructor
        function obj = LSTMModel(m, c, k, a, A, beta, gamma, n, dt, T, F_std, numSamples)
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
            
            % Initialize external force parameters
            obj.F_std = F_std;
            obj.numSamples = numSamples;
        end
        
       % Method to generate data
        function generateData(obj)
            % Initialize cell arrays to store responses and forces
            obj.inputData = cell(obj.numSamples, 1); % Force F(t) as cell array
            obj.outputData = cell(obj.numSamples, 1); % Displacement u(t) as cell array
            
            % Loop to generate samples and solve the system
            for i = 1:obj.numSamples
                % Generate white noise excitation
                whiteNoiseGenerator = WhiteNoiseLoadGeneratorModel(obj.dt, obj.T, obj.F_std);
                F = whiteNoiseGenerator.F; % White noise excitation (501×1)
                
                % Store the force as a 1×501 numeric vector in a cell
                obj.inputData{i} = F; % Transpose to row vector [1×501]
                
                % Solve the system using Runge-Kutta
                boucWenSystem = BoucWenSDOFModel(obj.m, obj.c, obj.k, obj.a, obj.A, obj.beta, obj.gamma, obj.n, obj.dt, obj.T, F);
                boucWenSystem.simulate();
                
                % Store the response as a 1×501 numeric vector in a cell
                obj.outputData{i} = boucWenSystem.u'; % Transpose to row vector [1×501]
            end
        end

        % Method to train the LSTM model
        function trainModel(obj)
            % Normalize the data across all samples
            allInputs = cat(2, obj.inputData{:}); % Concatenate all input samples
            allOutputs = cat(2, obj.outputData{:}); % Concatenate all output samples
            
            inputMean = mean(allInputs(:));
            inputStd = std(allInputs(:));
            outputMean = mean(allOutputs(:));
            outputStd = std(allOutputs(:));
            
            % Normalize each sample
            inputData = cellfun(@(x) (x - inputMean) / inputStd, obj.inputData, 'UniformOutput', false);
            outputData = cellfun(@(x) (x - outputMean) / outputStd, obj.outputData, 'UniformOutput', false);
            
            % Split into training and validation sets (80% training, 20% validation)
            trainRatio = 0.8;
            numTrain = floor(trainRatio * obj.numSamples);
            
            trainInput = inputData(1:numTrain); % Training input cell array
            trainOutput = outputData(1:numTrain); % Training output cell array
            valInput = inputData(numTrain+1:end); % Validation input cell array
            valOutput = outputData(numTrain+1:end); % Validation output cell array
            
            % Define the LSTM model
            layers = [
                sequenceInputLayer(1) % Input layer (1 feature: F(t))
                lstmLayer(10, 'OutputMode', 'sequence') % LSTM layer with 10 hidden units
                dropoutLayer(0.2) % Dropout layer to prevent overfitting
                fullyConnectedLayer(1) % Fully connected output layer (1 output: u(t))
                regressionLayer % Regression layer for training
            ];
            
            % Training options
            options = trainingOptions('adam', ...
                'MaxEpochs', 100, ...
                'MiniBatchSize', 16, ...
                'InitialLearnRate', 0.01, ...
                'ValidationData', {valInput, valOutput}, ...
                'ValidationFrequency', 10, ...
                'Verbose', false, ...
                'Plots', 'training-progress');
            
            % Train the LSTM model
            obj.net = trainNetwork(trainInput, trainOutput, layers, options);
        end
        
       % Method to evaluate the LSTM model
        function evaluateModel(obj)
            % Generate a new load sample
            newWhiteNoiseGenerator = WhiteNoiseLoadGeneratorModel(obj.dt, obj.T, obj.F_std);
            newF = newWhiteNoiseGenerator.F; % New force F(t) (501×1)
            
            % Normalize the new force using the same parameters as training data
            allInputs = cell2mat(obj.inputData);
            inputMean = mean(allInputs(:));
            inputStd = std(allInputs(:));
            newF = (newF' - inputMean) / inputStd; % Normalize and transpose to [1×501]
            
            % Solve the system using Runge-Kutta
            newBoucWenSystem = BoucWenSDOFModel(obj.m, obj.c, obj.k, obj.a, obj.A, obj.beta, obj.gamma, obj.n, obj.dt, obj.T, newF');
            newBoucWenSystem.simulate();
            trueResponse = newBoucWenSystem.u; % True response from Runge-Kutta
            
            % Predict the response using the trained LSTM model
            newFCell = {newF'}; % Wrap in a cell array {1×501}
            predictedResponse = predict(obj.net, newFCell); % Output is a cell array
            
            % Denormalize the responses
            allOutputs = cell2mat(obj.outputData);
            outputMean = mean(allOutputs(:));
            outputStd = std(allOutputs(:));
            trueResponse = trueResponse * outputStd + outputMean;
            predictedResponse = predictedResponse{1} * outputStd + outputMean; % Extract from cell
            
            % Plot the comparison
            figure;
            subplot(2, 1, 1);
            plot(obj.t, trueResponse, 'b', 'LineWidth', 1.5);
            hold on;
            plot(obj.t, predictedResponse, 'r--', 'LineWidth', 1.5);
            xlabel('Time (s)');
            ylabel('Displacement u(t)');
            title('Response Time Histories');
            legend('Runge-Kutta', 'LSTM Prediction');
            grid on;
            
            % Plot peak responses
            subplot(2, 1, 2);
            plot(obj.t, max(trueResponse) * ones(size(obj.t)), 'b', 'LineWidth', 1.5);
            hold on;
            plot(obj.t, max(predictedResponse) * ones(size(obj.t)), 'r--', 'LineWidth', 1.5);
            xlabel('Time (s)');
            ylabel('Peak Displacement');
            title('Peak Responses');
            legend('Runge-Kutta', 'LSTM Prediction');
            grid on;
        end

        function saveFigure(obj)
            % Get current date in 'yyyymmdd_HHMM' format
            dateString = datestr(now, 'yyyymmdd_HHMM');
            
            % Construct the base file name using the parameters
            Ofile = sprintf('%s - Task03_RK4_vs_LSTM', dateString);
            
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