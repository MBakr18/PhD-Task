classdef SDOFModel < handle
    properties
        Ns % Number of stories
        Ks % Story Stiffness
        Ms % Story mass
        Cs % Structural damping
        Tzeta_s
        
        NTMD % Number of Tuned Mass Dampers
        Tzeta % TMD damping ratio
        TMR % Mass Ratio for Tuned Mass Dampers
        Ntotal % Totla degree of freedom

        K % Structure stiffness
        M % Structure mass
        C % Structure damping
        Omegas % Strcture Omegas
        ModeShapes  % Structure mode shapes
        TimePeriods % structure Time Periods
        wstr0 % min omega of structure

        % Additional properties for auxillary damper
        damperAuxillaryMassRatio % auxillary mass ratio
        kd % damper stiffness
        cd % damper damping
        
        ka % auxillary stiffness
        ca % auxillary damping
    end
    
    methods
        % Constructor
        function obj = SDOFModel(Ns, Ks, Ms, Cs, NTMD, TMR, Tzeta, damperAuxillaryMassRatio)
            obj.Ns = Ns;
            obj.Ks = Ks;
            obj.Ms = Ms;
            obj.Cs = Cs;
            obj.NTMD = NTMD;
            obj.TMR = TMR;
            obj.Tzeta = Tzeta;
            obj.Ntotal = obj.Ns + obj.NTMD;

            % Additional properties for auxillary damper
            obj.damperAuxillaryMassRatio = damperAuxillaryMassRatio;

            computeStructureProperties(obj);
            computeTMDProperties(obj);
        end

        function computeStructureProperties(obj)      
            % Initialize matrices and parameters
            obj.K = zeros(obj.Ns);
            obj.M = zeros(obj.Ns);
            obj.C = zeros(obj.Ns);

            % Compute Structure Stiffness Matrix (K)
            for ix = 1:obj.Ns
                if ix > 1
                    obj.K(ix-1, ix-1) = obj.K(ix-1, ix-1) + obj.Ks;
                    obj.K(ix-1, ix) = -obj.Ks;
                    obj.K(ix, ix-1) = -obj.Ks;
                end    
                obj.K(ix, ix) = obj.K(ix, ix) + obj.Ks;  
            end
            % Compute Structure Mass Matrix (M)
            obj.M = obj.Ms * eye(obj.Ns);
            
            % Eigenvalue problem for Structure
            [phi, lam] = eig(obj.K, obj.M);
            w = sqrt(diag(lam));
            obj.Omegas = w;
            obj.TimePeriods = 1 * 2 * pi ./ w;
            obj.wstr0 = min(obj.Omegas);

            % Compute Structure Damping (C)
            if (obj.Ns > 1)
                ccc = [obj.Tzeta_s * 2*(w(1) * w(2)) / (w(1) + w(2)), obj.Tzeta_s * 2/(w(1) + w(2))];
                obj.C = ccc(1) * obj.M  + ccc(2) * obj.K;
                obj.C(obj.Ns, obj.Ns) = -sum(obj.C(obj.Ns, :));
            else
                % obj.C = obj.Tzeta_s * 2 * w * obj.M;
                obj.C = obj.Cs;
            end
        end

        function computeTMDProperties(obj)
            if obj.NTMD > 0
                n = obj.Ns + obj.NTMD;
        
                % Expand matrices K and M to accommodate TMDs
                obj.K = obj.K * eye(n);
                obj.M = obj.M * eye(n);
                obj.C = obj.C * eye(n);
                Wd = zeros(obj.NTMD, 1);
   
                % Loop through each Tuned Mass Damper
                for ix = 1:obj.NTMD
                    if ix == 1
                        obj.M(obj.Ns + ix, obj.Ns + ix) = (obj.TMR * obj.damperAuxillaryMassRatio) * (obj.Ms * obj.Ns);
                        obj.kd = obj.M(obj.Ns + ix, obj.Ns + ix) * (obj.wstr0)^2;
                    else
                        obj.M(obj.Ns + ix, obj.Ns+ix) = (obj.TMR * (1 - obj.damperAuxillaryMassRatio)) * (obj.Ms * obj.Ns);
                        obj.ka = obj.M(obj.Ns+ix, obj.Ns+ix) * (obj.wstr0)^2;
                    end
                end
                
                for ix = 1:obj.NTMD
     
                    if ix == 1
                        obj.K(obj.Ns+ix, obj.Ns) = -obj.kd;  % Uncoupled TMDs
                        obj.K(obj.Ns, obj.Ns+ix) = -obj.kd;  % Uncoupled TMDs
                        K_k = obj.kd+obj.ka;
                        obj.K(obj.Ns, obj.Ns) = obj.K(obj.Ns, obj.Ns) + obj.kd;
                    else 
                        K_k = obj.ka;

                        obj.K(obj.Ns+ix, ix) = -obj.ka;  % Uncoupled TMDs
                        obj.K(ix, obj.Ns+ix) = -obj.ka;  % Uncoupled TMDs
                    end
                             
                    % Update the stiffness matrix K
                    obj.K(obj.Ns + ix, obj.Ns + ix) = K_k;
                end
                % Compute Eigenvalues and Eigenvectors
                [phi, lam] = eig(obj.K, obj.M);
                w = sqrt(diag(lam));
                obj.Omegas = w;
                obj.TimePeriods = 1 * 2 * pi ./ w;
                obj.ModeShapes = phi(1:obj.Ntotal, :) ./ phi(1, :);
        
                % Compute damping matrix C using Non-Classical Damping
                for ix = 1:obj.NTMD
                    ccc = 2 * obj.Tzeta * obj.wstr0 * obj.M(obj.Ns + ix, obj.Ns + ix);
                    if ix == 1
                        obj.cd = ccc;
                    else
                        obj.ca = ccc;
                    end
                end
                

                for ix = 1:obj.NTMD                    
                    if ix == 1
                        obj.C(obj.Ns, obj.Ns+ix) = -obj.cd;
                        obj.C(obj.Ns+ix, obj.Ns) = -obj.cd;
                        obj.C(obj.Ns + ix, obj.Ns + ix) = obj.cd+obj.ca;
                        obj.C(obj.Ns, obj.Ns) = obj.Cs + obj.cd;
                    else
                        obj.C(obj.Ns+ix, ix) = -obj.ca;  % Uncoupled TMDs
                        obj.C(ix, obj.Ns+ix) = -obj.ca;  % Uncoupled TMDs
                        obj.C(obj.Ns + ix, obj.Ns + ix) = obj.ca;
                    end
                end 
            end
        end
    end
end
