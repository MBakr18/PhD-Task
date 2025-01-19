classdef DynamicAnalysisModel < handle
    properties
        acc % Acceleration data
        dt % Time step
        Ns % Number of stories
        NTMD % Number of TMDs
        NTotal
        Ks
        K % Stiffness coefficients
        C
        M
        Hs % Heights of stories
        HsCummulative

        % ... other properties will be used in calculation
        t % Time vector
        u % Displacement vector
        ud % Velocity vector
        udd % Acceleration vector
        Qb % Base shear vector
        Mb % Base moment vector
        beta = 1/4;                     
        gama = 1/2;

        ntt
        MaxNORM % Placeholder for max norm
        MaxDisp % Placeholder for max displacement
        MaxAcc % Placeholder for max acceleration
        MaxBaseQ % Placeholder for max base shear
        MaxBaseM % Placeholder for max base moment
        topStoryDrift % Placeholder for max top story drift
        
        TMDStroke
        MaxStroke

        Drift
        MaxDrift
        KhStatus
        PhStatus
    end
    
    methods
        function obj = DynamicAnalysisModel(acc, dt, Ns, NTotal, NTMD, Ks, K, M, C, Hs)
            obj.acc = acc;
            obj.dt = dt;
            obj.Ns = Ns;
            obj.NTMD = NTMD;
            obj.NTotal = NTotal;
            obj.Ks = Ks;
            obj.K = K;
            obj.M = M;
            obj.C = C;
            obj.Hs = ones(1, Ns(1));
            obj.Hs = Hs * obj.Hs;
            performAnalysis(obj);
            % plotDisplacement(obj)
        end
        
        function performAnalysis(obj)
            [~, nt] = size(obj.acc);
            obj.ntt = nt;
            te = (nt - 1) * obj.dt;
            obj.t = 0:obj.dt:te;
            n = obj.Ns + obj.NTMD;
           
            % Temp = 0;
            % obj.HsCummulative = zeros(size(obj.Hs));
            % 
            % for i = 1 : obj.Ns
            %     obj.HsCummulative(i) = Temp + obj.Hs(i);
            %     Temp = obj.HsCummulative(i);
            % end

            % Perform numerical integration and force calculation

            % Prepare factors 
            a1 = (obj.gama/(obj.beta*obj.dt)) * obj.C + (1 /(obj.beta*obj.dt*obj.dt)) * obj.M;
            a2 = (1/(obj.beta*obj.dt)) * obj.M +(( obj.gama/obj.beta)-1)* obj.C;
            a3 = ((1/(2*obj.beta))-1) * obj.M + obj.dt*((obj.gama/(2*obj.beta)-1)) * obj.C;
            I = ones(n,1);
            Kh = obj.K + a1;
            obj.KhStatus = Kh;

            obj.u   =  zeros  (n, nt);
            obj.ud  =  zeros  (n, nt);
            obj.udd =  zeros  (n, nt);
            
            obj.t(1)=0;
            % Initialize matrices for iteration
            ph = zeros(n, nt);
            obj.PhStatus = zeros(n, nt);
            R1 = zeros(1, nt); 
            
            % ensure that the forces at TMDs are zeros
            if obj.NTMD  > 0
                for i = obj.Ns+1 :obj.Ns+obj.NTMD     
                    for ii = 1 :nt
                        obj.acc(i,ii)  = 0;
                    end
                end
            end

            kl = zeros(1, obj.Ns);
            for  ix = 1:obj.Ns
                kl(ix) = obj.Ks;
            end
            % start iteration at each time step
            for ts = 2 :nt                            %ts is the time step 
                obj.t(ts)= (ts-1)*obj.dt;
                ut (1:n,1) = obj.u(1:n,ts-1);
                udt (1:n,1) = obj.ud(1:n,ts-1);
                uddt (1:n,1) = obj.udd(1:n,ts-1);    
                  
                % the problem is here in ph
                ph(1:n, ts) = obj.acc(1:n,ts) + a1 *obj.u(1:n ,ts-1) + a2 *obj.ud(1:n ,ts-1) + a3 * obj.udd(1:n ,ts-1) ;    
                obj.PhStatus(1:n, ts) = ph(1:n, ts);
                
                obj.u(1:n,ts) = Kh\ph(1:n, ts);
                obj.ud(1:n,ts) = obj.gama/(obj.beta*obj.dt)*(obj.u(1:n ,ts)-obj.u(1:n ,ts-1))   +(1-(obj.gama/obj.beta))*obj.ud(1:n,ts-1) +  obj.dt * (1-(obj.gama/(2*obj.beta))) * obj.udd(1:n,ts);
                obj.udd(1:n,ts) = (1/(obj.beta*obj.dt*obj.dt))*(obj.u(1:n ,ts)-obj.u(1:n ,ts-1)) - (1/(obj.beta*obj.dt))*obj.ud(1:n,ts-1) - (1/(2*obj.beta)-1)*obj.udd(1:n,ts-1);                      
                   
                % storke
                for iTMD = 1 : obj.NTMD
                    obj.TMDStroke(iTMD, ts) = obj.u(obj.Ns+iTMD, ts) - obj.u(obj.Ns,ts);
                end

                % Drift
                if obj.Ns > 1
                    for i = 1 : obj.Ns
                        if (i==1)
                            obj.Drift(i,ts) = obj.u(i, ts) - 0;
                        else
                            obj.Drift(i,ts) = obj.u(i, ts) - obj.u(i-1, ts);
                        end
                    end
                end
                % Element Forces
            
                obj.Qb(ts)=0;
                obj.Mb(ts)=0;
                HH = 0;
                Q = zeros(obj.Ns-1, nt);
                for i =1 : obj.Ns-1
                    if i==1 && i<obj.Ns 
                        Q(i,ts) = kl(i) * obj.u(i,ts) ; 
                        Q(i,ts) = kl(i)  * (obj.u(i,ts)) - kl(i+1) * (obj.u(i+1,ts) - obj.u(i,ts));
                    end
                    if i > 1 && i < obj.Ns        
                        Q(i,ts) = kl(i) * (obj.u(i,ts) - obj.u(i-1,ts) ) - kl(i+1) * (obj.u(i+1,ts) - obj.u(i,ts));   
                    end
                    if i  == obj.Ns && i < obj.Ns  
                        Q(i,ts) = kl (i) * (obj.u(i,ts) - obj.u(i-1,ts));   
                    end
                    HH = HH + obj.Hs(i);
                    obj.Qb(ts) = obj.Qb(ts) + Q(i, ts);
                    obj.Mb(ts)  = obj.Mb(ts) + Q(i, ts) * HH;
                end
               
                %%Response Norm
            
                R1 (ts) = 0;
                for ii = 1:obj.Ns
                    R1 (ts) = R1 (ts) + (obj.u(ii,ts)^2);        %Displacement Norm at ts       
                end    
                R1 (ts) = sqrt(R1 (ts));    
            end
           
            
            % Update properties accordingly:
            obj.MaxNORM = max(R1);
            obj.MaxDisp = max(abs(obj.u(obj.Ns,1500:nt))); 
            obj.MaxAcc = max(abs(obj.udd(obj.Ns,1500:nt)));
            obj.MaxBaseQ = max(abs(obj.Qb))  ; 
            obj.MaxBaseM = max(abs(obj.Mb));

            % Max Stroke 
            for iTMD = 1 : obj.NTMD
                obj.MaxStroke(iTMD,1) = max(abs(obj.TMDStroke(iTMD, 1500:nt)));
            end
            % Max Drift at each story
            for iStory = 1 : obj.Ns
                obj.MaxDrift(iStory, 1) = max(abs(obj.Drift(iStory, 1500:nt)));
            end
            % Calculate the maximum value and its index
            [maxValue, maxIndex] = max(abs(obj.u(obj.Ns, 1500:nt)));
            index_of_maxValue = maxIndex + 1499; % Adjust index to the original array

            % [maxValue, maxIndex] = max(abs(obj.u(obj.Ns, 1:nt)));
            % index_of_maxValue = maxIndex; % Adjust index to the original array
            if obj.Ns > 1
                obj.topStoryDrift = abs((obj.u(obj.Ns,index_of_maxValue)) - abs(obj.u(obj.Ns-1,index_of_maxValue)));
            end        
        end
        
        function plotDisplacement(obj)
            % Plot displacement
            figure; 
            plot(obj.t,zeros(1,obj.ntt), obj.t, obj.u(obj.Ns,1:obj.ntt)*1000  ,'k'); 
            grid on; axis tight; 
            xlabel('Time [sec]'); ylabel('Disp [mm]'); 
        end    
                
        function plotDrift(obj)
            % Plot displacement
            figure; 

            x = obj.MaxDrift * 1000;
            y = obj.HsCummulative;
            plot(x, y); 
            grid on; axis tight; 
            xlabel('Story Drift [mm]'); ylabel('Height[m]'); 
        end
    end
end
