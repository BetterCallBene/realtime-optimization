classdef Quadrocopter < Model
    % QUADROCOPTER repräsentiert einen bestimmten Quadrocopter und seine
    % physikalischen Eigenschaften.
    
    properties(Constant, GetAccess=public)
		n_state = 13;
		n_contr = 4;

        I = [0.0093886, 0.0093886, 0.018406];
        % Tr�gheitsmoment gesamt in kg m^2
        m   = 1.022;                                % Gesamtgewicht des Quadrokopters in kg
        I_M = 0.0001;     %ToDo                              % Tr�gheitsmoment der Motoren/Rotoren in kg m^2
        kT  = 1.5e-07;                              % N/(RPM^2)
        kQ  = 3e-09;                                % N*m/(RPM)^2
        d   = 0.22;                                 % Abstand in m
        
    end

	properties(Dependent)
		n_var;
        I_1;
    end
    
    methods
        function ret = get.I_1(cq)
            ret = [
                [cq.I(1),         0,        0];
                [0,         cq.I(2),        0];
                [0,               0,   cq.I(3)];
            ];
        end

		function ret = n_var.get(obj)
			ret = obj.n_state + obj.n_contr;
		end
    end
    
end
