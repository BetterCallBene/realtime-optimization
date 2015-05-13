classdef(Abstract) classDyn < handle  & classTestEnv
% classDyn providing the right hand side of the ODE
%   Reformulation of the second-order ODE as a first-order ODE

    properties(Abstract)
        environment; % in this object all outside parameters are stored, like gravity, wind, time mesh,...
        state; % the state vector combining the joint angles theta and their first time derivates dot_theta for all time instances
        contr; % the control vector combining the two joint torque values for all time instances
        robot; % handle for a classRobot element providing the mass matrix and the coriolis terms
    end
    
    properties(SetAccess = private, GetAccess = private)
        tol = 1e-4;
    end
    
    methods(Abstract)
        %constructor
        % other functions
        dot(obj,ind);
        dotD(obj,ind);
        dotDD(obj,ind);
    end
     
    methods(Test)
        
        function testdotD(obj)
           %TODO: Debug do we have to init timemesh?
           %Initialize state with random values
           obj.state=rand(13,1);
           obj.contr=rand(4,1);
           
           stco = [obj.state; obj.contr];
           
           eps = 1e-6;
           
           num_dotD = zeros(13,17);
           
           %Test at time step 1
           %TODO: Decide, if we also want to loop over time
           for j=1:17
               
               %plus epsilon shift
               stco_p = stco;
               stco_p(j) = stco_p(j) + eps; 
               obj.state=stco_p(1:13);
               obj.contr=stco_p(14:17);
               dot_p = obj.dot(1);
               
               %minus epsilon shift
               stco_n = stco;
               stco_n(j) = stco_n(j) - eps; 
               obj.state=stco_n(1:13);
               obj.contr=stco_n(14:17);
               dot_n = obj.dot(1);
               
               %central difference
               dotD_xj = (dot_p - dot_n)/(2*eps);
               num_dotD(:,j) = dotD_xj;
               
               obj.state=state_old;
               
           end
           
           %Calculate analytic solution
           ana_dotD = obj.dotD(obj,1);
           
           %Assert that the result has the correct form
           obj.assertSize(ana_dotD, [13,17]);
           
           %Compare the results
           obj.assertLessThan(max(abs(ana_dotD - num_dotD)), obj.tol);
           
        end
          
        function testdotDD(obj)
            %TODO: Debug do we have to init timemesh? Is the data structure
            %correct?
            %Initialize state with random values
            obj.state=rand(13,1);
            obj.contr=rand(4,1);
            
            stco = [obj.state; obj.contr];
            
            eps = 1e-6;
            
            num_dotDD = cell(17,1);
            
            %TODO: Decide if we also want to loop over time
            for j = 1:17
                %plus epsilon shift
                stco_p = stco;
                stco_p(j) = stco_p(j) + eps;
                obj.state=stco_p(1:13);
                obj.contr=stco_p(14:17);
                dotD_p = obj.dotD(1);
                
                %minus epsilon shift
                stco_n = stco;
                stco_n(j) = stco_n(j) - eps;
                obj.state=stco_n(1:13);
                obj.contr=stco_n(14:17);
                dotD_n = obj.dotD(1);
                
                %central difference
                dotDD_xj = (dotD_p - dotD_n)/(2*eps);
                num_dotDD{j} = dotDD_xj;
                
                obj.state=state_old;
            end
                ana_dotDD = obj.dotDD(obj,1);
            for j = 1:length(num_dotDD)
                obj.assertLessThan(max(abs(ana_dotDD{j} - num_dotDD{j})),obj.tol);
            end    
        end
    end
end