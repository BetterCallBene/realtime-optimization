classdef(Abstract) classDyn < handle
% classDyn providing the right hand side of the ODE
%   Reformulation of the second-order ODE as a first-order ODE

    properties(Abstract)
        state; % the state vector combining the joint angles theta and their first time derivates dot_theta for all time instances
        contr; % the control vector combining the two joint torque values for all time instances
        robot; % handle for a classRobot element providing the mass matrix and the coriolis terms
    end
    
    methods(Abstract)
        %constructor
        % other functions
        dot(obj,ind);
        dotD(obj,ind);
        dotDD(obj,ind);
    end
end