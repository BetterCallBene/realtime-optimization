classdef RiccatiManager <  TestEnv
    % RICCATIMANAGER Central class to perform a riccati solution
    
    properties(Access=private)
        robot;
        A;
        B;
        M;
        P;
        Q;
        R;
        
        C;
        D;
        
        nabla_s_star;
        nabla_lambda;
        nabla_q;
        
        nabla_mu;
                        
        delta_lambda;
        delta_s;
        delta_q;
        
        delta_mu;
        
        horizon;
        n_lambda;
        n_state;
        n_contr;
        n_mu;
        n_var;
    end
    
    properties
        delta;
    end
    
       
    methods
        function o = RiccatiManager(varargin)
            
             if(nargin == 0)
                global TEST;
                if ~(~isempty(TEST) && TEST == true)
                    error('wrong number of inputs');
                else
                    o.horizon = 20;
                    o.robot = Quadrocopter();
                    o.n_mu = 0;
                end
            elseif(nargin == 1)
                o.horizon = varargin{1};
                o.robot = Quadrocopter();
                o.n_mu = 0;
            elseif(nargin == 2)
                o.horizon = varargin{1};
                if (isa(varargin{2},'Model'))
                    o.robot = varargin{2};
                end
                o.n_mu = 0;
             elseif(nargin == 3)
                 o.horizon = varargin{1};
                 if (isa(varargin{2}, 'Model'))
                     o.robot = varargin{2};
                 end
                 o.n_mu = varargin{3};
            else
                error('wrong number of inputs');
            end
            
            o.initialize();
        end
        
        function doStep(o,i, LDD_i, LD_i )
            %TODO: implement
        end
        
        function solveStep(o,i)
           %TODO: implement 
        end
        
        
       
    end
    
    methods(Test)
        
        
        
    end
    
    methods
        
        function initialize(o)
            %Initialize storage
            o.n_lambda = o.robot.n_state;
            o.n_state = o.robot.n_state;
            o.n_contr = o.robot.n_contr;
            o.n_var = o.n_lambda + o.n_state + o.n_contr + o.n_mu;
            
            o.A = cell(o.horizon,1);
            o.B = cell(o.horizon,1);
            o.M = cell(o.horizon,1);
            o.P = cell(o.horizon+1,1);
            o.Q = cell(o.horizon+1,1);
            o.R = cell(o.horizon,1);
            o.C = cell(o.horizon,1);
            o.D = cell(o.horizon,1);
            
            o.nabla_s_star = cell(o.horizon+1, 1);
            o.nabla_lambda = cell(o.horizon+1,1);
            o.nabla_q = cell(o.horizon ,1);
            o.nabla_mu = cell(o.horizon, 1);
            
            o.delta = zeros( o.n_var * (o.horizon +1 ) - (o.n_mu +  o.n_contr), 1);
            
            
            o.delta_lambda = cell(o.horizon +1, 1);
            o.delta_s = cell(o.horizon +1, 1);
            o.delta_q = cell(o.horizon, 1);
            o.delta_mu = cell(o.horizon,1);
        end
    end
end