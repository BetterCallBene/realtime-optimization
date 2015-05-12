classdef classQuadrocopterDyn < classQuadrocopterBasisDyn & classTestEnv
% classDyn providing the right hand side of the ODE
%   Reformulation of the second-order ODE as a first-order ODE

    properties
        %contr; % the control vector combining the two joint torque values for all time instances
        robot; % handle for a classRobot element providing the mass matrix and the coriolis terms
    end
    
    properties(SetAccess = private, GetAccess = private)
    end
    
    properties(Dependent, GetAccess= public)
        %n_var;
        %n_contr;
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
        
%         function ret =  get.n_var(obj)
%             if (~isempty(obj.robot))
%                 ret = obj.robot.n_var;
%             end
%         end
%         
%         function ret = get.n_contr(obj)
%             if (~isempty(obj.robot))
%                 ret = obj.robot.n_contr;
%             end
%         end
        
        %set methods
        %function set.state(obj,state)
            % set the current values for the states at all time instances
            % and pass it on the robot
        %    obj.state = state;
        %    if (~isempty(obj.robot))
        %        obj.robot.state = state;
        %    end
        %end
        
        %set methods
        %function set.contr(obj,contr)
            % set the current values for the states at all time instances
            % and pass it on the robot
        %    obj.contr = contr;
        %    if (~isempty(obj.robot))
        %        obj.robot.contr = contr;
        %    end
        %end
        
        % other functions
        % TODO: Mapleskript einpflegen, Dynamik 
        function res = dot(obj,ind) % Rechte Seite der ODE aus (Quadrocoptermodell.pdf: (1.46))
            % compute the right hand side of the ode for a given time
            % instance ind
            if ((size(obj.contr,2)==size(obj.state,2))&&(ind <= size(obj.contr,2)))
                F_ = obj.F;
                % ToDo Funktionsausgabe
            else
                error('wrong state and control lengths wrt index.');
            end
        end
        % Testen
        
        function res = dotD(obj,ind)
            % compute the Jacobian of the right hand side of the ode for 
            % a given time instance ind            
            if ((size(obj.contr,2)==size(obj.state,2))&&(ind <= size(obj.contr,2)))
                
                J_ = obj.J;
                res = zeros(13, 17);
                for i = 1:length(J_)
                    res(J_(i, 1), J_(i, 2)) = J_(i, ind + 2);
                end
            else
                error('wrong state and control lengths wrt index.');
            end
        end
        
        function res = dotDD(obj,ind)
            % compute the Hessian of the right hand side of the ode for 
            % a given time instance ind                
            if ((size(obj.contr,2)==size(obj.state,2))&&(ind <= size(obj.contr,2)))
                res = cell(1, 13);
                H_ = obj.H;
                for i = 1:length(H_)
                    if isempty(res{H_(i, 1)})
                        res{H_(i, 1)} = zeros(17, 17);
                    end
                    res{H_(i, 1)}(H_(i, 2), H_(i, 3)) = H_(i, 4);
                end
            end    
        end  
    end
    
    methods(Test)
%         function testdot(obj)
%             n_int_ = 50;
%             ind = 50;
%             
%             n_var_ = obj.n_var;
%             n_contr_ = obj.n_contr;
%             
%             obj.state = rand(n_var_, n_int_+1);
%             obj.contr = rand(n_contr_, n_int_+1);
%             
%             val = obj.dot(ind);
%             
%             obj.verifySize(val, [13, 1]);
%        end
        
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
        
        function testdotDDs(obj)
            n_int_ = 50;
            ind = 2;
            
            n_var_ = obj.n_var;
            n_contr_ = obj.n_contr;
            
            obj.state = rand(n_var_, n_int_+1);
            obj.contr = rand(n_contr_, n_int_+1);
            
            val = obj.dotDD(ind);
            
            obj.verifySize(val, [1, 13]);
        end
    end
    
    
end