classdef BasisQDyn < BasisGenQDyn
    
    properties
        state; % the state matrix [r, q, v, w] in R^(13xn) for all time instances
        contr; % the control matrix [u1, u2, u3, u4] in R^(4xn) for all time instances
    end
    
    methods
        function bQDyn = BasisGenQDyn(varargin)
            
            if (nargin >= 1)
                
                mc = metaclass(varargin{1});
                if strcmp(mc.getSuperclassList(1).Name, 'Model')
                    bQDyn.robot = varargin{1};
                else
                    error('First argument must be instance of Model');
                end
            end
            if (nargin >= 2)
                if isa(varargin{2}, 'Environment')
                    bQDyn.environment = varargin{2};
                else
                    error('Second argument must be instance of Environment');
                end
            end
        end
        
        function res = dot(obj,ind) % Rechte Seite der ODE aus (Quadrocoptermodell.pdf: (1.46))
            % compute the right hand side of the ode for a given time
            % instance ind
            if ((size(obj.contr,2)==size(obj.state,2))&&(ind <= size(obj.contr,2)))
                res = obj.F;
            else
                error('wrong state and control lengths wrt index.');
            end
        end
        
        
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
                    res{H_(i, 1)}(H_(i, 2), H_(i, 3)) = H_(i, ind + 3);
                end
            end
        end
        
    end
    
end