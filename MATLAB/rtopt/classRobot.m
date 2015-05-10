classdef classRobot < handle
% classRobot providing mass matrix and coriolis terms of 2D robot arm
%   Use the analytic versions for mass matrix and coriolis terms for the
%   2D robot arm to compute these values and their derivatives

    properties
        state;      % the state matrix combining the joint angles theta and their first time derivates dot_theta for all time instances
        param;      % handle for a classDynParam element providing the parameters for the dynamics
    end
        
    properties (SetAccess = private)
        sine;       % vector of sine values for all joint angles at all time instances
        cosine;     % vector of cosine values for all joint angles at all time instances

        M_add;      % summands in mass matrix computation not depending on the state vector
        M_fac;      % factors in mass matrix computation not depending on the state vector
        theta_fac;  % factors in theta computation not depending on the state vector
        
        M;          % matrix of mass matrix entries for all time instances
        MD;         % matrix of Jacobian entries of mass matrix for all time instances
        MDD;        % matrix of Hessian entries of mass matrix for all time instances
        theta;      % matrix of theta entries for all time instances
        thetaD;     % matrix of Jacobian entries of theta for all time instances
        thetaDD;    % matrix of Hessian entries of theta for all time instances
    end
    
    methods
        %constructor
        function cR = classRobot(varargin)
            % void constructor
        end
        
        %set methods
        function set.param(obj,param)
            % set param element
            if (isa(param,'classDynParam'))
                obj.param = param;
                obj.emptyResults();
            else
                error('Wrong type for param');
            end
        end
        
        function set.state(obj,state)
            % set current state matrix
            obj.state = state;
            obj.emptyResults();
        end

        %get methods
        function res = get.sine(obj)
            % get sine values corresponding to current state
            if (isempty(obj.sine))
                obj.sine = sin(obj.state(1:2,:));
            end
            res = obj.sine;
        end   
        
        function res = get.cosine(obj)
            % get cosine values corresponding to current state
            if (isempty(obj.cosine))
                obj.cosine = cos(obj.state(1:2,:));
            end
            res = obj.cosine;
        end 
        
        function res = get.M_add(obj) 
            % get the 'constant' summands for mass matrix computations
            % (if not yet stored, computed them)
            if (isempty(obj.M_add))
                length_val  = obj.param.Length;
                mass_val    = obj.param.Mass;
                
                obj.M_add = zeros(2,1);
                obj.M_add(1) = sum(obj.param.Inertia) ...
                    + mass_val(1)*(.5*length_val(1))^2 ...
                    + mass_val(2)*(length_val(1)^2 + ...
                                               (.5*length_val(2))^2);
                obj.M_add(2) = .25*mass_val(2)*length_val(2)^2 ...
                    + obj.param.Inertia(2);
            end
            res = obj.M_add;
        end
        
        function res = get.M_fac(obj) 
            % get the 'constant' factos for mass matrix computations
            % (if not yet stored, computed them)
            if (isempty(obj.M_fac))
                length_val  = obj.param.Length;
                mass_val    = obj.param.Mass;
                
                obj.M_fac = zeros(2,1);
                obj.M_fac(1) = mass_val(2)*length_val(1) ...
                                    *length_val(2);
                obj.M_fac(2) = .5*mass_val(2)*length_val(1) ...
                                    *length_val(2);
            end
            res = obj.M_fac;
        end   
        
        function res = get.theta_fac(obj) 
            % get the 'constant' factors for theta computations
            % (if not yet stored, computed them)
            if (isempty(obj.theta_fac))
                obj.theta_fac  = .5*obj.param.Mass(2)*obj.param.Length(1) ...
                                    *obj.param.Length(2);
            end
            res = obj.theta_fac;
        end         
        
        function res = get.M(obj)
            % get the matrix of mass matrix entries
            % (if not yet stored, computed them)
            if (isempty(obj.M))
                cosine_val          = obj.cosine;
                M_fac_val           = obj.M_fac;
                M_add_val           = obj.M_add;
                n                   = size(cosine_val,2);
          
                obj.M       = zeros(3,n);
                obj.M(1,:)  = M_add_val(1) + M_fac_val(1)*cosine_val(2,:);
                obj.M(2,:)  = M_add_val(2) + M_fac_val(2)*cosine_val(2,:);
                obj.M(3,:)  = M_add_val(2);                
            end
            res     = obj.M;      
        end
        
        function res = get.MD(obj)
            % get the matrix of Jacobian entries of mass matrix
            % (if not yet stored, computed them)            
            if (isempty(obj.MD))
                sine_val          = obj.sine;
                M_fac_val           = obj.M_fac;
                n                   = size(sine_val,2);
          
                obj.MD       = zeros(2,n);
                obj.MD(1,:)  = -M_fac_val(1)*sine_val(2,:);
                obj.MD(2,:)  = -M_fac_val(2)*sine_val(2,:);               
            end
            res     = obj.MD;      
        end
        
        function res = get.MDD(obj)
             % get the matrix of Hessian entries of mass matrix
            % (if not yet stored, computed them)   
            if (isempty(obj.MDD))
                cosine_val          = obj.cosine;
                M_fac_val           = obj.M_fac;
                n                   = size(cosine_val,2);
          
                obj.MDD       = zeros(2,n);
                obj.MDD(1,:)  = -M_fac_val(1)*cosine_val(2,:);
                obj.MDD(2,:)  = -M_fac_val(2)*cosine_val(2,:);             
            end
            res     = obj.MDD;      
        end     
        
        function res = get.theta(obj)
            % get the matrix of theta entries
            % (if not yet stored, computed them)               
            if (isempty(obj.theta))
                sine_val            = obj.sine;
                theta_fac_val       = obj.theta_fac;
                n                   = size(sine_val,2);
                
                obj.theta       = zeros(2,n);
                obj.theta(1,:)  = -theta_fac_val*obj.state(4,:).*...
                    (2*obj.state(3,:) + obj.state(4,:)).*sine_val(2,:);
                obj.theta(2,:)  = theta_fac_val* ...
                   obj.state(3,:).^2.*sine_val(2,:);
            end
            res = obj.theta;            
        end
        
        function res = get.thetaD(obj)
            % get the matrix of Jacobian entries of theta
            % (if not yet stored, computed them)               
            if (isempty(obj.thetaD))

                sine_val            = obj.sine;
                cosine_val          = obj.cosine;
                theta_fac_val       = obj.theta_fac;
                n                   = size(sine_val,2);
                
                obj.thetaD          = zeros(5,n);
                obj.thetaD(1,:)     = -theta_fac_val*obj.state(4,:).*...
                    (2*obj.state(3,:) + obj.state(4,:)).*cosine_val(2,:); %Ableitung Theta1 nach theta2
                obj.thetaD(2,:)     = -2*theta_fac_val* ...
                    obj.state(4,:).*sine_val(2,:); % Ableitung Theta1 nach theta1 punkt
                obj.thetaD(3,:)     = -2*theta_fac_val* ...
                    (obj.state(3,:) + obj.state(4,:)).*sine_val(2,:);% Ableitung Theta1 nach theta2 punkt
                
                obj.thetaD(4,:)     = theta_fac_val* ...
                   obj.state(3,:).^2.*cosine_val(2,:);
                obj.thetaD(5,:)     = 2*theta_fac_val* ...
                   obj.state(3,:).*sine_val(2,:);               
            end
            res = obj.thetaD;            
        end   
        
         function res = get.thetaDD(obj)
            % get the matrix of Hessian entries of theta
            % (if not yet stored, computed them)                
            if (isempty(obj.thetaDD))
                sine_val            = obj.sine;
                cosine_val          = obj.cosine;
                theta_fac_val       = obj.theta_fac;
                
                n                    = size(sine_val,2);
                obj.thetaDD          = zeros(6,n);
                obj.thetaDD(1,:)     = theta_fac_val*obj.state(4,:).*...
                    (2*obj.state(3,:) + obj.state(4,:)).*sine_val(2,:);
                obj.thetaDD(2,:)     = -theta_fac_val* ...
                    obj.state(3,:).^2 .*sine_val(2,:);                
                obj.thetaDD(3,:)     = -2*theta_fac_val* ...
                    obj.state(4,:).*cosine_val(2,:);
                obj.thetaDD(4,:)     = 2*theta_fac_val* ...
                    obj.state(3,:).*cosine_val(2,:);
                obj.thetaDD(5,:)     = -2*theta_fac_val* ...
                    (obj.state(3,:) + obj.state(4,:)).*cosine_val(2,:);
                obj.thetaDD(6,:)     = 2*theta_fac_val* sine_val(2,:);            
            end
            res = obj.thetaDD;            
        end          
        
        % other methods
        function emptyResults(obj)
            % if state or param are change all these values
            % have to be computed again.
            obj.M       = [];
            obj.MD      = [];
            obj.MDD     = [];
            
            obj.theta   = [];
            obj.thetaD  = [];
            obj.thetaDD = [];
            
            obj.sine    = [];
            obj.cosine  = [];
        end
        
        function res = getM(obj,ind)
            % the actual mass matrix at time instance ind
            M_val       = obj.M;
            
            res         = zeros(2,2);
            res(1:2,1)  = M_val(1:2,ind);
            res(1,2)    = res(2,1);
            res(2,2)    = M_val(3,ind);
        end
        
        function res = getMD(obj,ind)
            % the actual 1st derivative of mass matrix at time instance ind
            res             = cell(1,4);           
            res{2}          = zeros(2,2);
            res{2}(1:2,1)   = obj.MD(1:2,ind);
            res{2}(1,2)     = res{2}(2,1);
        end
        
        function res = getMDD(obj,ind)
            % the actual 2nd derivative of mass matrix at time instance ind
            res             = cell(4,4);
            res{2,2}        = zeros(2,2);
            res{2,2}(1:2,1) = obj.MDD(1:2,ind);
            res{2,2}(1,2)   = res{2,2}(2,1);
        end       
        
        function res = getTheta(obj,ind)
            % the actual coriolis terms
            res = obj.theta(:,ind);
        end
        
        function res = getThetaD(obj,ind)
            % the actual 1st derivative of coriolis terms
            res = zeros(2,4);
            
            thetaD_val = obj.thetaD;
            
            res(1,2:4) = thetaD_val(1:3,ind)';
            res(2,2:3) = thetaD_val(4:5,ind)';
        end  

        function res = getThetaDD(obj,ind)
            % the actual 2nd derivative of coriolis terms
            res         = cell(2,4);

            thetaDD_val = obj.thetaDD;

            res{2,2}    = thetaDD_val(1:2,ind);
            res{2,3}    = thetaDD_val(3:4,ind);
         
            res{2,4}    = zeros(2,1);
            res{2,4}(1) = thetaDD_val(5,ind);
            
            res{3,3}    = zeros(2,1);
            res{3,3}(2) = thetaDD_val(6,ind);
            
            res{3,4}    = zeros(2,1);
            res{3,4}(1) = -thetaDD_val(6,ind);
            
            res{4,4}    = zeros(2,1);
            res{4,4}(1) = -thetaDD_val(6,ind);  
        end  
        
    end
end