classdef classQuadrocopter < classModell
    %CLASSQUADROCOPTER Konfiguration eines Quadrocopters, Bereitstellung
    %der Massenmatrix und Coriolis Kräfte
    %   Der Quadrokopter in + Konfiguration hat folgende Spezifikationen:   
    %   Gross Weight: m = 1.022 kg
    %   Trägheitsmoment
    %   I = [Ix  0   0 ] Ix = 0.0093886 kg m^2
    %       [0   Iy  0 ] Iy = 0.0093886 kg m^2
    %       [0   0   Iz] Iz = 0.018406  kg m^2
    %   Berechnet mit GUI_Modeling.m
    
    properties(Constant)
        I = [0.0093886, 0.0093886, 0.018406];
        % Trägheitsmoment gesamt in kg m^2
        m   = 1.022;                                % Gesamtgewicht des Quadrokopters in kg
        I_M = 10;                                   % Trägheitsmoment der Motoren/Rotoren in kg m^2
        kT  = 1.5e-07;                              % N/(RPM^2)
        kQ  = 3e-09;                                % N*m/(RPM)^2
        d   = 0.22;                                 % Abstand in m
    end
    
    properties
        I_1;
    end
    
    properties(...
        SetAccess = protected, ... 
        GetAccess = public ...
        )
        state;      % the state matrix [r, q, v, w, dot v, dot w] in R^19 for all time instances
        control;    % the control matrix [w_M1, w_M2, w_M3, w_M4]
        M;          % matrix of mass matrix entries for all time instances
        MD;         % matrix of Jacobian entries of mass matrix for all time instances
        MDD;        % matrix of Hessian entries of mass matrix for all time instances
        theta;      % matrix of theta entries for all time instances
        thetaD;     % matrix of Jacobian entries of theta for all time instances
        thetaDD;    % matrix of Hessian entries of theta for all time instances
        env;
    end
    
    methods(Access=public) 
        function cQ = classQuadrocopter(vargin)
            
            if nargin > 0
                obj1 = vargin(1); 
                
                if ~isobject(obj1)
                    error('Parameter 1 ist kein Object');
                end
                mc = metaclass(obj1);
                
                if strcmp(mc.SuperclassList(1).Name, 'classEnv')
                    cQ.env = obj1;
                else
                    error('Objekt erbt nicht von classEnv')
                end
            end
        end
    end
    
    methods(Access=private, Static)
        %R Berechnung der Rotationsmatrix aus den Quaternionen
        function ret = R(q)
            
            if length(q) ~=4 
                error('Länge des Quaternionvektors ist falsch: Geforderte Länge 4');
            end
            
            q = quatnormalize(q);
                        
            % Rotationsmatrix: Geometriekalküle Seite 201/Wikipedia:
            % Quaternionen
            ret = [
                    [1-2*(q(3)^2+q(4)^2), -2 *q(1) * q(4) + 2 * q(2) * q(3), 2 * q(1) * q(3) + 2 * q(2) * q(4)];
                    [2 *q(1) * q(4) + 2 * q(2) * q(3), 1-2*(q(2)^2+q(4)^2), -2 * q(1) * q(2) + 2 * q(3) * q(4)];
                    [-2 *q(1) * q(3) + 2 * q(2) * q(4), 2 * q(1) * q(2) + 2 * q(3) * q(4), 1-2*(q(2)^2+q(3)^2)]
                ];
        end
        function ret = R3(q)
            
            q = quatnormalize(q);
            % -2 *q(1) * q(3) + 2 * q(2) * q(4), 2 * q(1) * q(2) + 2 * q(3) * q(4), 1-2*(q(2)^2+q(3)^2)
            n   = size(q, 2);
            ret = zeros(3, n);
            ret(1, :) = -2 *q(1, :) * q(3, :) + 2 * q(2, :) * q(4, :);
            ret(2, :) = 2 * q(1, :) * q(2, :) + 2 * q(3, :) * q(4, :);
            ret(3, :) = 1- 2* (q(2, :).^2+q(3, :).^2);
        end
    end
    
    methods
        function ret = get.I_1(cq)
            ret = [
                [cq.I(1),         0,        0];
                [0,         cq.I(2),        0];
                [0,               0,   cq.I(3)];
            ];
        end
        %M Bestimmung der Massenmatrix des Quadrocopters
        function ret = get.M(cq)
            if isempty(cq.M)
                cq.M = zeros(6, 6);
                cq.M(1:3, 1:3) = eye(3) .* cq.m;
                cq.M(4:6, 4:6) = cq.I_1;
            end
            ret = cq.M;
        end
        %
        function ret = get.theta(cq)
            
            if isempty(cq.theta)
                %state = [r, q, v, w, dot v, dot w]^T in R^19
                
                q   = cq.state(4:7    , :);
                r_g = R3(q) * cq.env.g; %
                v   = cq.state(8:10   , :);
                w   = cq.state(11:12  , :);
                
                n = size(q, 2); % n Nur einmal setzen, oben?? 
                
                cq.theta = zeros(6, n);
                cq.theta(1, :) = cq.m * (w(2, :) * v(3, :) - w(3, :) * v(2, :) + r_g(1, :));
                cq.theta(2, :) = cq.m * (w(1, :) * v(3, :) - w(3, :) * v(1, :) + r_g(2, :));
                cq.theta(3, :) = cq.m * (w(1, :) * v(2, :) - w(2, :) * v(1, :) + r_g(3, :));
                cq.theta(4, :) = cq.I(3) *  w(2, :) * v(3, :) - cq.I(2) * w(3, :) * v(2, :);
                cq.theta(5, :) = -cq.I(3) *  w(1, :) * v(3, :) + cq.I(1) * w(3, :) * v(1, :);
                cq.theta(6, :) = cq.I(2) *  w(1, :) * v(2, :) - cq.I(1, :) * w(2, :) * v(1, :);
                
            end
            ret = cq.theta;
        end
        
        function ret = get.thetaD(cq)
            if isempty(cq.theta)
                q   = cq.state(4:7    , :);
                r_g = R3(q) * cq.env.g; %
                v   = cq.state(8:10   , :);
                w   = cq.state(11:12  , :);
                %   state = [r, q, v, w, dot v, dot w]^T in R^19
                
                %%   Ableitung (groß) Theta 1 nach ...
                cq.thetaD(1, :) = cq.m * cq.env.g *  -2 * q(3); % q0
                cq.thetaD(2, :) = cq.m * cq.env.g *   2 * q(4); % q1
                cq.thetaD(3, :) = cq.m * cq.env.g *  -2 * q(1); % q2
                cq.thetaD(4, :) = cq.m * cq.env.g *   2 * q(2); % q3
                cq.thetaD(5, :) = cq.m * cq.env.g *   2 * q(2); % v_1
                    
                
                
            end
            
        end
        
    end
    
    methods(Test)
        % Test der Rotationsmatrix
        function testR(obj)
            alpha = pi/2;
            alpha = alpha/2;
            r = [0 0 1]; %Drehachse
            r = r/norm(r);  
            r = sin(alpha) .* r;
            q = [cos(alpha) r(1) r(2) r(3)];
            R = classQuadrocopter.R(q);
            obj.verifyEqual(R * [1; 0; 0], [0; 1; 0], 'AbsTol', 1e-15);
        end
        
        function testTheta(obj)
            
            
        end
    end
    
end

