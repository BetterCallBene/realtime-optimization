classdef(Abstract) classDyn < handle  & classTestEnv
    % classDyn providing the right hand side of the ODE
    %   Reformulation of the second-order ODE as a first-order ODE
    
    properties(SetAccess = protected, GetAccess= public)
        environment; % in this object all outside parameters are stored, like gravity, wind, time mesh,...
        robot; % handle for a classRobot element providing the mass matrix and the coriolis terms
        state; % the state matrix [r, q, v, w] in R^(13xn) for all time instances
        contr; % the control matrix [u1, u2, u3, u4] in R^(4xn) for all time instances
        
    end
    
    properties(SetAccess = private, GetAccess = private)
        tol = 1e-4; %Tolerance for the num differentiation
    end
    
    properties(Dependent, GetAccess = public)
        n_state_contr;
    end
    
    properties(Constant, GetAccess = public)
        n_state = 13;   % Count of Variables
        n_contr = 4;    % Count of Controls
        
    end
    
    methods(Abstract)
        dot(obj,ind);
        dotD(obj,ind);
        dotDD(obj,ind);
    end
    
    methods
        
        function ret = get.n_state_contr(obj)
            ret = obj.n_contr + obj.n_state;
        end
        
        function set.state(cq, state)
            
            if size(state, 1) ~= cq.n_state
                error(strcat('Gr��e der State Matrix ist falsch. Erwartete Gr��e  ',int2str(cq.n_state),'xn'));
            end
            n = size(state, 2);
            cq.set_interval(n); %TODO: ist im moment eine Dummyfunktion
            
            for i=1:n
                q = state(4:7, i);
                state(4:7, i) = q./norm(q);
            end
            
            cq.state = state;
            cq.emptyResults();
        end
        
        function set.contr(cq, cntrl)
            
            if size(cntrl, 1) ~= cq.n_contr
                error(strcat('Gr��e der State Matrix ist falsch. Erwartete Gr��e  ',int2str(cq.n_contr),'xn'));
            end
            cq.set_interval(size(cntrl, 2));
            cq.contr = cntrl;
            cq.emptyResults();
        end
        
        
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
                
                obj.state=stco(1:13);
                obj.contr=stco(14:17);
                
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
            
            num_dotDD = zeros(17,13,17);
            
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
                num_dotDD(j,:,:) = dotDD_xj;
                
                obj.state=stco(1:13);
                obj.contr=stco(14:17);
            end
            ana_dotDD = obj.dotDD(obj,1);
            for j = 1:length(num_dotDD)
                %TODO: Adapt data formats
                obj.assertLessThan(max(abs(ana_dotDD{j} - num_dotDD(j,:,:))),obj.tol);
            end
        end
    end
end