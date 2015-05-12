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
        d   = 0.22;                       % Abstand in m
        g   =9.81;
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
                
                t1 = u(3, :) .^ 2;
t2 = u(2, :) .^ 2;
t3 = u(1, :) .^ 2;
t4 = u(4, :) .^ 2;
t5 = -u(1, :) + u(2, :) - u(3, :) + u(4, :);
t6 = kT .* d;
cg = [0 0 kT .* (t3 + t2 + t1 + t4) t6 .* (t2 - t4) + IM .* omega(2, :) .* t5 t6 .* (t1 - t3) - IM .* omega(1, :) .* t5 kQ .* (-t3 + t2 - t1 + t4)];


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
                
                t1 = 2 .* kT;
t2 = t1 .* u(1, :);
t3 = t1 .* u(2, :);
t4 = t1 .* u(3, :);
t1 = t1 .* u(4, :);
t5 = -u(1, :) + u(2, :) - u(3, :) + u(4, :);
t6 = IM .* omega(2, :);
t7 = IM .* omega(1, :);
t8 = 2 .* kQ;
cg0 = [3 4 t2; 3 5 t3; 3 6 t4; 3 7 t1; 4 2 IM .* t5; 4 4 -t6; 4 5 t3 .* d + t6; 4 6 -t6; 4 7 -t1 .* d + t6; 5 1 -IM .* t5; 5 4 -t2 .* d + t7; 5 5 -t7; 5 6 t4 .* d + t7; 5 7 -t7; 6 4 -t8 .* u(1, :); 6 5 t8 .* u(2, :); 6 6 -t8 .* u(3, :); 6 7 t8 .* u(4, :);];

                
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

                t1 = 2 .* kT;
t2 = t1 .* d;
t3 = 2 .* kQ;
cg1 = [3 4 4 t1; 3 5 5 t1; 3 6 6 t1; 3 7 7 t1; 4 2 4 -IM; 4 2 5 IM; 4 2 6 -IM; 4 2 7 IM; 4 4 2 -IM; 4 5 2 IM; 4 5 5 t2; 4 6 2 -IM; 4 7 2 IM; 4 7 7 -t2; 5 1 4 IM; 5 1 5 -IM; 5 1 6 IM; 5 1 7 -IM; 5 4 1 IM; 5 4 4 -t2; 5 5 1 -IM; 5 6 1 IM; 5 6 6 t2; 5 7 1 -IM; 6 4 4 -t3; 6 5 5 t3; 6 6 6 -t3; 6 7 7 t3;];

                
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

            obj.getH(1)
        end
    end
    
end

