classdef classQuadrocopterDyn < classDyn & classTestEnv
% classDyn providing the right hand side of the ODE
%   Reformulation of the second-order ODE as a first-order ODE

    properties
        state; % the state vector combining the joint angles theta and their first time derivates dot_theta for all time instances
        contr; % the control vector combining the two joint torque values for all time instances
        robot; % handle for a classRobot element providing the mass matrix and the coriolis terms
    end
    
    properties(SetAccess = private, GetAccess = private)
    end
    
    properties(Dependent, GetAccess= public)
        n_var;
        n_contr;
    end
    
    
    
    methods
        %constructor
        function cR = classQuadrocopterDyn(varargin)
            % a classRobot element needed as input for constructor
            if (nargin == 0) %BasisModell
                cR.robot = classQuadrocopter();
            elseif (nargin == 1)
                if (isa(varargin{1},'classQuadrocopter'))
                    cR.robot = varargin{1};
                else
                    error('wrong class type for robot');
                end
            else
                error('wrong number of inputs');     
            end
        end
        
        function ret =  get.n_var(obj)
            if (~isempty(obj.robot))
                ret = obj.robot.n_var;
            end
        end
        
        function ret = get.n_contr(obj)
            if (~isempty(obj.robot))
                ret = obj.robot.n_contr;
            end
        end
        
        %set methods
        function set.state(obj,state)
            % set the current values for the states at all time instances
            % and pass it on the robot
            obj.state = state;
            if (~isempty(obj.robot))
                obj.robot.state = state;
            end
        end
        
        %set methods
        function set.contr(obj,contr)
            % set the current values for the states at all time instances
            % and pass it on the robot
            obj.contr = contr;
            if (~isempty(obj.robot))
                obj.robot.contr = contr;
            end
        end
        
        % other functions
        % TODO: Mapleskript einpflegen, Dynamik 
        function val = dot(obj,ind) % Rechte Seite der ODE aus (Quadrocoptermodell.pdf: (1.46))
            % compute the right hand side of the ode for a given time
            % instance ind
            if ((size(obj.contr,2)==size(obj.state,2))&&(ind <= size(obj.contr,2)))
                
                M           = obj.robot.getM(ind);
                Theta       = obj.robot.getTheta(ind);
                T           = obj.robot.getT(ind);
                Rv          = obj.robot.getRv(ind);
                Q           = obj.robot.getQ(ind);
                
                
                val          = zeros(13,1);
                val(1:3,1)   = Rv;
                val(4:7,1)   = Q;
                val(8:13,1)  = M\(T-Theta);
            else
                error('wrong state and control lengths wrt index.');
            end
        end
        % Testen
        
        function val = dotD(obj,ind)
            % compute the Jacobian of the right hand side of the ode for 
            % a given time instance ind            
            if ((size(obj.contr,2)==size(obj.state,2))&&(ind <= size(obj.contr,2)))
                
                M           = obj.robot.getM(ind);
                RvD         = obj.robot.getRvD(ind);
                QD          = obj.robot.getQD(ind);
                TD           = obj.robot.getTD(ind);
                %MD          = obj.robot.getMD(ind);
                %theta       = obj.robot.getTheta(ind);
                thetaD      = obj.robot.getThetaD(ind);
                
                val         = zeros(13,17);
                
                val(1:3, 4:10) = RvD;
                val(4:7, 4:7) = QD(:, 1:4);
                val(4:7, 11:13) = QD(:, 5:7);
                val(8:13, 4:10) = -M\thetaD(:, 1:7);
                val(8:13, 11:13) = M\(TD(:, 1:3) - thetaD(:, 8:end));
                val(8:13, 14:17) = M\TD(:, 4:end);
            else
                error('wrong state and control lengths wrt index.');
            end
        end
        
        function val = dotDD(obj,ind)
            % compute the Hessian of the right hand side of the ode for 
            % a given time instance ind                
            val = cell(6, 6);
            error('dotDD: muss noch implementiert werden')
            if ((size(obj.contr,2)==size(obj.state,2))&&(ind <= size(obj.contr,2)))
            end    
        end  
        
        function res = getJ(obj, ind)
            J_ = obj.J;
            res = zeros(6, 7);
            for i = 1:length(J_)
                res(J_(i, 1), J_(i, 2)) = J_(i, ind + 2);
            end
        end

        function res = getH(obj, ind)
            res = cell(1, 6);

            H_ = obj.H;
            for i = 1:length(H_)
                if isempty(res{H_(i, 1)})
                    res{H_(i, 1)} = zeros(7, 7);
                end
                res{H_(i, 1)}(H_(i, 2), H_(i, 3)) = H_(i, 4);
            end
        end 
  
    end
    
    methods(Test)
        function testdot(obj)
            n_int_ = 50;
            ind = 50;
            
            n_var_ = obj.n_var;
            n_contr_ = obj.n_contr;
            
            obj.state = rand(n_var_, n_int_+1);
            obj.contr = rand(n_contr_, n_int_+1);
            
            val = obj.dot(ind);
            
            obj.verifySize(val, [13, 1]);
        end
        
        function testdotD(obj)
            n_int_ = 50;
            ind = 2;
            
            n_var_ = obj.n_var;
            n_contr_ = obj.n_contr;
            
            obj.state = rand(n_var_, n_int_+1);
            obj.contr = rand(n_contr_, n_int_+1);
            
            val = obj.dotD(ind);
            
            obj.verifySize(val, [13, 17]);
        end
    end
    
    
end