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
        kT  = 1.5e-07;                              % N/(RPM^2)
        kQ  = 3e-09;                                % N*m/(RPM)^2
        d   = 0.22;                                 % Abstand in m
        motor_m = 0.075;                            % Masse Motor
        motor_r = 0.015;                            % Motor R 1.5 cm
    end
    
    %
    % Jm (Motor Rotation Inertia for Rotating Component only)
    % Mass of rotating component is 52.7% of the total mass of the motor + prop
    %mRC = (motor_m)*(0.527);
    %Jm = ((mRC)*(motor_r)^2)/2; % Jm = mr^2/2
    %
    
    

	properties(Dependent)
		n_var;
        I_1;
        I_M;
    end
    
    methods
        function res = Quadrocopter(varargin)
            
        end
        % Jm (Motor Rotation Inertia for Rotating Component only)
        % Mass of rotating component is 52.7% of the total mass of the motor + prop
        %mRC = (motor_m)*(0.527);
        %Jm = ((mRC)*(motor_r)^2)/2; % Jm = mr^2/2
        %
        function IM = get.I_M(cq)
            mRC = (cq.motor_m)*(0.527);
            IM = ((mRC)*(cq.motor_r)^2)/2;
        end
        function ret = get.I_1(cq)
            ret = [
                [cq.I(1),         0,        0];
                [0,         cq.I(2),        0];
                [0,               0,   cq.I(3)];
            ];
        end

		function ret = get.n_var(obj)
			ret = obj.n_state + obj.n_contr;
		end
    end
    
end
