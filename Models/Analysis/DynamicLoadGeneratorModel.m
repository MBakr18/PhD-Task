classdef DynamicLoadGeneratorModel < handle
    properties
        tt
        dt
        NTotal
        Ae
        we
        acc % Acceleration matrixNOTMD
    end
    
    methods
        % Constructor
        function obj = DynamicLoadGeneratorModel(tt, dt, NTotal, Ae, we)
            obj.tt = tt;
            obj.dt = dt;
            obj.NTotal = NTotal;
            obj.Ae = Ae;
            obj.we = we;
            obj.acc = zeros(NTotal, (tt / dt + 1));

            generateSinusoidalLoads(obj);
        end
        
        function generateZeroLoads(obj)
            obj.acc = zeros(obj.NTotal, (obj.tt / obj.dt + 1));
        end
        
        function generateSinusoidalLoads(obj)
            for ip = 1:obj.NTotal
                for iip = 1:(obj.tt / obj.dt + 1)
                    t = (iip - 1) * obj.dt;
                    obj.acc(ip, iip) = obj.Ae * sin(obj.we * t);
                end
            end
        end
        
        function plotLoad(obj, story)
            figure;
            plot(0:obj.dt:obj.tt, obj.acc(story, :), 'k');
            grid on;
            axis tight;
            xlabel('time');
            ylabel('Load');
            title(['Story ' num2str(story) ' Load']);
        end
    end
end
