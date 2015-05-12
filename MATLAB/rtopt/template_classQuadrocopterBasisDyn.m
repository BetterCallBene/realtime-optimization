classdef classQuadrocopterBasisDyn < handle & classTestEnv 
    %CLASSQUADROCOPTERBASISDYN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        F;
        J;
        H;
        
        state; % the state matrix [r, q, v, w] in R^(13xn) for all time instances
        contr; % the control matrix [u1, u2, u3, u4] in R^(4xn) for all time instances
    end
    
    properties(Constant, GetAccess=private)
        I = [0.0093886, 0.0093886, 0.018406];
        % Trägheitsmoment gesamt in kg m^2
        m   = 1.022;                                % Gesamtgewicht des Quadrokopters in kg
        I_M = 0.0001;     %ToDo                              % Trägheitsmoment der Motoren/Rotoren in kg m^2
        kT  = 1.5e-07;                              % N/(RPM^2)
        kQ  = 3e-09;                                % N*m/(RPM)^2
        d   = 0.22;                                 % Abstand in m
        g   = 9.81;
    end
    
    properties(...
        Dependent, ...
        GetAccess = public ...
    )
        n_state_contr;
    end

    properties(...
            Constant, ...
            GetAccess = public ...
            )
        n_var = 13;        % Count of Variables
        n_contr = 4;       % Count of Controls
    end
    
    properties(GetAccess=private, SetAccess = protected)
         n_int;
         isEmptyF;
         isEmptyJ;
         isEmptyH;
    end
    
    methods(Access=private)
        function set_interval(obj, n_int)
            obj.n_int = n_int;
        end
    end
    
    methods
        function cQ = classQuadrocopterBasisDyn(vargin)
            %cQ@classModell();
        end

        
        function ret = get.n_state_contr(obj)
            ret = obj.n_contr + obj.n_var;
        end
        
        function emptyResults(obj)
            
            %emptyResults@classModell(obj);
            
            obj.F   = [];
            obj.J  = [];
            obj.H = [];
            
            obj.isEmptyF = true;    
            obj.isEmptyJ = true;
            obj.isEmptyH = true;
        end

        function set.state(cq, state)
            
            if size(state, 1) ~= 13
                error('Größe der State Matrix ist falsch. Erwartete Größe 13xn');
            end
            n = size(state, 2);
            cq.set_interval(n);
            
            for i=1:n
                q = state(4:7, i);
                state(4:7, i) = q./norm(q);
            end
            
            cq.state = state;
            cq.emptyResults();
        end
        
        function set.contr(cq, cntrl)
            
            if size(cntrl, 1) ~= 4
                error('Größe der State Matrix ist falsch. Erwartete Größe 4xn');
            end
            cq.set_interval(size(cntrl, 2));
            cq.contr = cntrl;
            cq.emptyResults();
        end
        
        function res = get.F(obj)
            if obj.isEmptyF 
                
                q   = obj.state(4:7    , :);
                v   = obj.state(8:10   , :);
                omega   = obj.state(11:12  , :);
                u   = obj.contr;
            
                I = obj.I;
                IM = obj.I_M;
                m = obj.m;
            
                kT  = obj.kT;
                kQ  = obj.kQ;
                d   = obj.d;
                g   = obj.g;
                
                {0}

                obj.F = cg;
                obj.isEmptyF = false;
            end
            res = obj.F;
        end
        
        function res = get.J(obj)
            if obj.isEmptyJ
                q   = obj.state(4:7    , :);
                v   = obj.state(8:10   , :);
                omega   = obj.state(11:12  , :);
                u   = obj.contr;
            
                I = obj.I;
                IM = obj.I_M;
                m = obj.m;
            
                kT  = obj.kT;
                kQ  = obj.kQ;
                d   = obj.d;
                g   = obj.g;
                
                {1}
                
                obj.J = cg0;
                obj.isEmptyJ = false;
            end
            res = obj.J;
        end
        
        function res = get.H(obj)
            if obj.isEmptyH
                q   = obj.state(4:7    , :);
                v   = obj.state(8:10   , :);
                omega   = obj.state(11:12  , :);
                u   = obj.contr;
            
                I = obj.I;
                IM = obj.I_M;
                m = obj.m;
            
                kT  = obj.kT;
                kQ  = obj.kQ;
                d   = obj.d;
                g   = obj.g;

                {2}
                
                obj.H = cg1;
                obj.isEmptyH = false;
            end
            res = obj.H;
        end
    end
    
    methods(Test)
        function testFJH(obj)

            n_int_ = 3;
            n_state_ = obj.n_var;
            n_contr_ = obj.n_contr;
            
            obj.state = rand(n_state_, n_int_+1);
            obj.contr = rand(n_contr_, n_int_+1);

            F_ = obj.F;
            H_ = obj.H;
            J_ = obj.J;
        end
    end
    
end

