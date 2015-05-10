classdef classRobotDyn < handle
% classDyn providing the right hand side of the ODE
%   Reformulation of the second-order ODE as a first-order ODE

    properties
        state; % the state vector combining the joint angles theta and their first time derivates dot_theta for all time instances
        contr; % the control vector combining the two joint torque values for all time instances
        robot; % handle for a classRobot element providing the mass matrix and the coriolis terms
    end
    
    methods
        %constructor
        function cR = classDyn(varargin)
            % a classRobot element needed as input for constructor
            if (nargin == 1)
                if (isa(varargin{1},'classRobot'))
                    cR.robot = varargin{1};
                else
                    error('wrong class type for robot');
                end
            else
                error('wrong number of inputs');     
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
        
        % other functions
        function val = dot(obj,ind)
            % compute the right hand side of the ode for a given time
            % instance ind
            if ((size(obj.contr,2)==size(obj.state,2))&&(ind <= size(obj.contr,2)))
                
                M           = obj.robot.getM(ind);
                theta       = obj.robot.getTheta(ind);
                
                val         = zeros(4,1);
                val(1:2,1)  = obj.state(3:4,ind);
                val(3:4,1)  = M\(obj.contr(:,ind)-theta);
            else
                error('wrong state and control lengths wrt index.');
            end
        end
        
        function val = dotD(obj,ind)
            % compute the Jacobian of the right hand side of the ode for 
            % a given time instance ind            
            if ((size(obj.contr,2)==size(obj.state,2))&&(ind <= size(obj.contr,2)))
                
                M           = obj.robot.getM(ind);
                MD          = obj.robot.getMD(ind);
                theta       = obj.robot.getTheta(ind);
                thetaD      = obj.robot.getThetaD(ind);
                
                val         = zeros(4,6);
                val(1,3)    = 1;
                val(2,4)    = 1;
                
                for i=1:4
                    if (~isempty(MD{i}))
                        val(3:4,i) = -M\(MD{i}*...
                            (M\(obj.contr(:,ind)-theta))+thetaD(:,i));
                    else
                        val(3:4,i) = -M\thetaD(:,i);
                    end
                end
                
                val(3:4,5:6) = inv(M);                
            else
                error('wrong state and control lengths wrt index.');
            end
        end
        
        function val = dotDD(obj,ind)
            % compute the Hessian of the right hand side of the ode for 
            % a given time instance ind                
            if ((size(obj.contr,2)==size(obj.state,2))&&(ind <= size(obj.contr,2)))
                
                M           = obj.robot.getM(ind);
                MD          = obj.robot.getMD(ind);
                MDD         = obj.robot.getMDD(ind);
                thetaD      = obj.robot.getThetaD(ind);
                thetaDD     = obj.robot.getThetaDD(ind);
                
                dot         = obj.dot(ind);
                
                val     = cell(4,1);
                val{1}  = zeros(6,6);
                val{2}  = zeros(6,6);
                val{3}  = zeros(6,6);
                val{4}  = zeros(6,6);

                for i=1:4
                    for j=1:4
                        sum         = 0;
                        if (i <= j)
                            if (~isempty(thetaDD{i,j}))
                                sum     = sum -thetaDD{i,j};
                            end
                            if (~isempty(MDD{i,j}))
                                sum     = sum - MDD{i,j}*dot(3:4,1);
                            end
                        else
                            if (~isempty(thetaDD{j,i}))
                                sum     = sum -thetaDD{j,i};
                            end
                            if (~isempty(MDD{j,i}))
                                sum     = sum - MDD{j,i}*dot(3:4,1);
                            end
                        end
                        if (~isempty(MD{i}))
                            sum     = sum + MD{i}*(M\thetaD(:,j));
                        end
                        if (~isempty(MD{j}))
                            sum     = sum + MD{j}*(M\thetaD(:,i));
                        end      
                        if ((~isempty(MD{i}))&&(~isempty(MD{j})))
                            sum     = sum + MD{j}*(M\(MD{i}*dot(3:4,1))) ...
                                        + MD{i}*(M\(MD{j}*dot(3:4,1))); 
                        end
                        
                        if (length(sum)>1) 
                            dotDD_val = M\sum;          

                            val{3}(i,j) = dotDD_val(1);
                            val{3}(j,i) = dotDD_val(1);
                            val{4}(i,j) = dotDD_val(2);
                            val{4}(j,i) = dotDD_val(2);   
                        end
                    end
                end

                for j=1:2
                    if (~isempty(MD{j}))
                        dotDD_val = -inv(M)*MD{j}*inv(M);
                        val{3}(5:6,j) = dotDD_val(1:2,1);
                        val{3}(j,5:6) = dotDD_val(1:2,1)';
                        val{4}(j,5:6) = dotDD_val(1:2,2)';
                        val{4}(5:6,j) = dotDD_val(1:2,2);
                    end
                end               
            else
                error('wrong state and control lengths wrt index.');
            end
        end        
        
    end
end