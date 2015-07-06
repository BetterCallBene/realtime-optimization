classdef Quadrocopter < Model
    %% QUADROCOPTER repraesentiert einen bestimmten Quadrocopter und seine physikalischen Eigenschaften.
    % * $$I = [0.0093886, 0.0093886, 0.018406]$$
    % 
    
    properties(Constant, GetAccess=public)
		n_state = 13;
		n_contr = 4;

        I = [0.0093886, 0.0093886, 0.018406];
        % Tr�gheitsmoment gesamt in kg m^2
        m   = 1.022;                                % Gesamtgewicht des Quadrokopters in kg
        %I_M = 0.0001;     %ToDo                              % Tr�gheitsmoment der Motoren/Rotoren in kg m^2
        kT  = 1.5e1;        %1.5e-7                      % N/(RPM^2)
        kQ  = 3e-01;        %3e-9                        % N*m/(RPM)^2
        d   = 0.22;                                 % Abstand in m
        c_w = 0.5;                                  % Luftwiderstand
        A   = 0.04;                                 % Stroemungsrelevanteflaeche
        motor_m = 0.075;                            % Masse Motor
        motor_r = 0.015;                            % Motor R 1.5 cm
        u_min = 0.01;                                  % Minimale Umdrehungen pro Minute
        u_max = 1;                                % Maximale Umdrehungen pro Minute
        rho = 1.2;                               % Luftdichte
    end
    
    %
    % Jm (Motor Rotation Inertia for Rotating Component only)
    % Mass of rotating component is 52.7% of the total mass of the motor + prop
    %mRC = (motor_m)*(0.527);
    %Jm = ((mRC)*(motor_r)^2)/2; % Jm = mr^2/2
    %
    
    
    properties
        I_M;
        F_w;                                      % Luftwiderstand
        flagI_M;
        flagF_w;
    end

	properties(Dependent)
		n_var;
        I_1;
    end
    
    methods
        function QC = Quadrocopter(varargin)
            QC.flagI_M = true;
            QC.flagF_w = true;
        end
        % Jm (Motor Rotation Inertia for Rotating Component only)
        % Mass of rotating component is 52.7% of the total mass of the motor + prop
        %mRC = (motor_m)*(0.527);
        %Jm = ((mRC)*(motor_r)^2)/2; % Jm = mr^2/2
        %
        function IM = get.I_M(cq)
            %if cq.flagI_M
            %    mRC = (cq.motor_m)*(0.527);
            %    cq.I_M = ((mRC)*(cq.motor_r)^2)/2;
            %    cq.flagI_M = false;
            %end
            IM = 4.4466e-06; %BB: Performance wgn 
        end
        function ret = get.I_1(cq)
            ret = [
                [cq.I(1),         0,        0];
                [0,         cq.I(2),        0];
                [0,               0,   cq.I(3)];
            ];
        end

		function ret = get.n_var(obj)
			%ret = obj.n_state + obj.n_contr;
            ret = 17; % Performance wgn
        end
        
        function res = get.F_w(obj)
            if obj.flagF_w 
                obj.F_w = 1/2 * obj.rho * obj.A *  obj.c_w;
                obj.flagF_w = false;
            end
            res = obj.F_w;
        end
        
        function res = getF_w(obj, v)
            res = obj.F_w * v.^2;
        end
        
        
    end
    
end
