classdef classQuadrocopter < classModell
    %CLASSQUADROCOPTER Konfiguration eines Quadrocopters, Bereitstellung
    %der Massenmatrix und Coriolis Kräfte
    %   Der Quadrokopter in + Konfiguration hat folgende Spezifikationen:   
    %   Gross Weight: m = 1.022 kg
    %   Trägheitsmoment
    %   I = [Ix  0   0 ] Ix = 0.0093886 kg m^2
    %       [0   Iy  0 ] Iy = 0.0093886 kg m^2
    %       [0   0   Iz] Iz = 0.018406  kg m^2
    %   Berechnet mit GUI_Modeling.m
    
    properties(Constant, GetAccess=private)
        I = [0.0093886, 0.0093886, 0.018406];
        % Trägheitsmoment gesamt in kg m^2
        m   = 1.022;                                % Gesamtgewicht des Quadrokopters in kg
        I_M = 0.0001;     %ToDo                              % Trägheitsmoment der Motoren/Rotoren in kg m^2
        kT  = 1.5e-07;                              % N/(RPM^2)
        kQ  = 3e-09;                                % N*m/(RPM)^2
        d   = 0.22;                                 % Abstand in m
        
        startIndexHesseThetaDD1 = 0;
        startIndexHesseThetaDD2 = 8;
        startIndexHesseThetaDD3 = 16;
        startIndexHesseThetaDD4 = 22;
        startIndexHesseThetaDD5 = 24;
        startIndexHesseThetaDD6 = 26;
        
        startIndexHesseTDD1 = 0;
        startIndexHesseTDD2 = 0;
        startIndexHesseTDD3 = 0;
        startIndexHesseTDD4 = 4;
        startIndexHesseTDD5 = 14;
        startIndexHesseTDD6 = 24;
        
        startIndexHesseRDD1 = 0;
        startIndexHesseRDD2 = 30;
        startIndexHesseRDD3 = 60;
        
        startIndexHesseQDD1 = 0;
        startIndexHesseQDD2 = 6;
        startIndexHesseQDD3 = 12;
        startIndexHesseQDD4 = 18;
    end
    
    %EmptyFlags für Matrizen: bessere Performance
    properties(SetAccess=private, ...
              GetAccess=private  ...
    )
       isEmptyM;    
       isEmptyMD;
       isEmptyMDD;
       
       isEmptyTheta;    
       isEmptyThetaD;
       isEmptyThetaDD;
       
       isEmptyT;
       isEmptyTD;
       isEmptyTDD;
       
       isEmptyRges;
       isEmptyRvD;
       isEmptyRvDD;
       
       isEmptyQ;
       isEmptyQD;
       isEmptyQDD;
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
    
    properties(...
            SetAccess = public, ...
            GetAccess = public  ...
    )
        state; % the state matrix [r, q, v, w] in R^(13xn) for all time instances
        contr; % the control matrix [u1, u2, u3, u4] in R^(4xn) for all time instances
    end
    
    properties(...
        SetAccess = protected, ... 
        GetAccess = public ...
        )
             
        M;          % matrix of mass matrix entries for all time instances
        MD;         % matrix of Jacobian entries of mass matrix for all time instances
        MDD;        % matrix of Hessian entries of mass matrix for all time instances
        
        theta;      % matrix of theta entries for all time instances
        thetaD;     % matrix of Jacobian entries of theta for all time instances
        thetaDD;    % matrix of Hessian entries of theta for all time instances
        
        T;          % matrix of T entries for all time instances
        TD;         % matrix of Jacobian entries of T for all time instances
        TDD;        % matrix of Hessian entries of T for all time instances
        
        Rges;
        RvD;
        RvDD;
        
        Q;
        QD;
        QDD;
        
        I_1;        %Trägheitsmatrix des Systems
        
        env;        % pointer to a Instanz of the classEnv
        n_int;
    end
    
    
        
    methods(Access=public) 
        function cQ = classQuadrocopter(vargin)
            cQ@classModell();
            
            if (nargin == 0)
                cQ.env = classEnvConstant();
            elseif nargin == 1
                obj1 = vargin(1); 
                
                if ~isobject(obj1)
                    error('Parameter 1 ist kein Object');
                end
                mc = metaclass(obj1);
                
                if strcmp(mc.SuperclassList(1).Name, 'classEnv')
                    cQ.env = obj1;
                else
                    error('Objekt erbt nicht von classEnv')
                end
            else
                error('Aktuell nur kein oder ein Übergabeparameter erlaubt');
            end
        end
    end
    
    methods(Access=private, Static)
        %R Berechnung der Rotationsmatrix aus den Quaternionen
        function ret = R(q)
            
            if length(q) ~=4 
                error('Länge des Quaternionvektors ist falsch: Geforderte Länge 4');
            end
                        
            % Rotationsmatrix: Geometriekalküle Seite 201/Wikipedia:
            % Quaternionen
            ret = [
                    [1-2*(q(3)^2+q(4)^2), -2 *q(1) * q(4) + 2 * q(2) * q(3), 2 * q(1) * q(3) + 2 * q(2) * q(4)];
                    [2 *q(1) * q(4) + 2 * q(2) * q(3), 1-2*(q(2)^2+q(4)^2), -2 * q(1) * q(2) + 2 * q(3) * q(4)];
                    [-2 *q(1) * q(3) + 2 * q(2) * q(4), 2 * q(1) * q(2) + 2 * q(3) * q(4), 1-2*(q(2)^2+q(3)^2)]
                ];
        end
        function ret = R3(q)
            
            %q = quatnormalize(q); % Hier Fehler?
            % -2 *q(1) * q(3) + 2 * q(2) * q(4), 2 * q(1) * q(2) + 2 * q(3) * q(4), 1-2*(q(2)^2+q(3)^2)
            n   = size(q, 2);
            ret = zeros(3, n);
            ret(1, :) = -2 *q(1, :) .* q(3, :) + 2 * q(2, :) .* q(4, :);
            ret(2, :) = 2 * q(1, :) .* q(2, :) + 2 * q(3, :) .* q(4, :);
            ret(3, :) = 1- 2* (q(2, :).^2+q(3, :).^2);
        end
        
        
    end
    
    methods(Access=private)
        function set_interval(obj, n_int)
            obj.n_int = n_int;
        end
        
    end
    
    methods(Access=protected)
        function emptyResults(obj)
            
            emptyResults@classModell(obj);
            
            obj.T   = [];
            obj.TD  = [];
            obj.TDD = [];
            
            obj.Rges    = [];
            obj.RvD     = [];
            obj.RvDD    = [];
        
            obj.Q       = [];
            obj.QD      = [];
            obj.QDD     = [];
            
            obj.isEmptyM = true;    
            obj.isEmptyMD = true;
            obj.isEmptyMDD = true;
            
            obj.isEmptyTheta= true;    
            obj.isEmptyThetaD= true;
            obj.isEmptyThetaDD= true;
            
            obj.isEmptyT= true;
            obj.isEmptyTD= true;
            obj.isEmptyTDD= true;
            
            obj.isEmptyRges = true;
            obj.isEmptyRvD = true;
            obj.isEmptyRvDD = true;
            
            obj.isEmptyQ = true;
            obj.isEmptyQD = true;
            obj.isEmptyQDD= true;
        end
    end
    
    methods
        function ret = get.I_1(cq)
            ret = [
                [cq.I(1),         0,        0];
                [0,         cq.I(2),        0];
                [0,               0,   cq.I(3)];
            ];
        end
        
        function ret = get.n_state_contr(obj)
            ret = obj.n_contr + obj.n_var;
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
        
        %M Bestimmung der Massenmatrix des Quadrocopters
        function ret = get.M(cq)
            if cq.isEmptyM
                cq.M = zeros(6, 6);
                cq.M(1:3, 1:3) = eye(3) .* cq.m;
                cq.M(4:6, 4:6) = cq.I_1;
                cq.isEmptyM = false;
            end
            ret = cq.M;
        end
        %
        function ret = get.theta(cq)
            
            if cq.isEmptyTheta
                %state = [r, q, v, w]^T in R^13
                
                q   = cq.state(4:7    , :);
                r_g = cq.R3(q) * cq.env.g; %
                v   = cq.state(8:10   , :);
                w   = cq.state(11:13  , :);
                
                
                cq.theta = zeros(6, cq.n_int);
                cq.theta(1, :) = cq.m * (w(2, :) .* v(3, :) - w(3, :) .* v(2, :) + r_g(1, :));
                cq.theta(2, :) = cq.m * (w(3, :) .* v(1, :) - w(1, :) .* v(3, :) + r_g(2, :));
                cq.theta(3, :) = cq.m * (w(1, :) .* v(2, :) - w(2, :) .* v(1, :) + r_g(3, :));
                cq.theta(4, :) = cq.I(3) *  w(2, :) .* v(3, :) - cq.I(2) * w(3, :) .* v(2, :);
                cq.theta(5, :) = cq.I(1) * w(3, :) .* v(1, :)  - cq.I(3) * w(1, :) .* v(3, :);
                cq.theta(6, :) = cq.I(2) *  w(1, :) .* v(2, :) - cq.I(1) * w(2, :) .* v(1, :);
                cq.isEmptyTheta = false;
            end
            ret = cq.theta;
        end
        
        function ret = get.thetaD(cq)
            if cq.isEmptyThetaD
                
                q   = cq.state(4:7    , :);
                %r_g = R3(q) * cq.env.g; %
                v   = cq.state(8:10   , :);
                w   = cq.state(11:13  , :);
                %   state = [r, q, v, w, dot v, dot w]^T in R^19
                
                cq.thetaD = zeros(28, cq.n_int);
                
                %%   Ableitung (groß) Theta 1 nach ...
                cq.thetaD(1, :) = - 2 * cq.m * cq.env.g * q(3, :); % q1
                cq.thetaD(2, :) =   2 * cq.m * cq.env.g * q(4, :); % q2
                cq.thetaD(3, :) = - 2 * cq.m * cq.env.g * q(1, :); % q3
                cq.thetaD(4, :) =   2 * cq.m * cq.env.g * q(2, :); % q4
                cq.thetaD(5, :) =   cq.m * w(3, :);                 % v2
                cq.thetaD(6, :) =   cq.m * w(2, :);                 % v4
                cq.thetaD(7, :) =   cq.m * v(3, :);                 % w2
                cq.thetaD(8, :) =   cq.m * v(2, :);                 % w3
                
                %%   Ableitung (groß) Theta 2 nach ...
                cq.thetaD(9, :)  =   2 * cq.m * cq.env.g * q(2, :); % q1
                cq.thetaD(10, :) =   2 * cq.m * cq.env.g * q(1, :); % q2
                cq.thetaD(11, :) =   2 * cq.m * cq.env.g * q(4, :); % q3
                cq.thetaD(12, :) =   2 * cq.m * cq.env.g * q(3, :); % q4
                cq.thetaD(13, :) =     cq.m * w(3, :);                 % v2
                cq.thetaD(14, :) =   - cq.m * w(1, :);                 % v4
                cq.thetaD(15, :) =   - cq.m * v(3, :);                 % w2
                cq.thetaD(16, :) =     cq.m * v(1, :);                 % w3
                
                %%   Ableitung (groß) Theta 3 nach ...
                cq.thetaD(17, :) =  -4 * cq.m * cq.env.g * q(2, :); % q2
                cq.thetaD(18, :) =  -4 * cq.m * cq.env.g * q(3, :); % q3
                cq.thetaD(19, :) =  -    cq.m * w(2, :);            % v1
                cq.thetaD(20, :) =       cq.m * w(1, :);            % v2
                cq.thetaD(21, :) =       cq.m * v(2, :);              % w1
                cq.thetaD(22, :) =  -    cq.m * v(1, :);              % w2
                
                %%   Ableitung (groß) Theta 3 nach ...
                cq.thetaD(17, :) =  -4 * cq.m * cq.env.g * q(2, :); % q2
                cq.thetaD(18, :) =  -4 * cq.m * cq.env.g * q(3, :); % q3
                cq.thetaD(19, :) =  -    cq.m * w(2, :);            % v1
                cq.thetaD(20, :) =       cq.m * w(1, :);            % v2
                cq.thetaD(21, :) =       cq.m * v(2, :);              % w1
                cq.thetaD(22, :) =  -    cq.m * v(1, :);              % w2
                
                %%   Ableitung (groß) Theta 4 nach ...
                cq.thetaD(23, :) =  (cq.I(3) - cq.I(2)) * w(3, :);  % w2
                cq.thetaD(24, :) =  (cq.I(3) - cq.I(2)) * w(2, :);  % w3
                
                %%   Ableitung (groß) Theta 5 nach ...
                cq.thetaD(25, :) =  (cq.I(1) - cq.I(3)) * w(3, :);  % w1
                cq.thetaD(26, :) =  (cq.I(1) - cq.I(3)) * w(1, :);  % w3
                
                %%   Ableitung (groß) Theta 6 nach ...
                cq.thetaD(27, :) =  (cq.I(2) - cq.I(1)) * w(2, :);  % w1
                cq.thetaD(28, :) =  (cq.I(2) - cq.I(1)) * w(1, :);  % w2
                
                cq.isEmptyThetaD = false;
                
            end
            ret = cq.thetaD;
        end
        
        function ret = get.thetaDD(cq)
            if cq.isEmptyThetaDD
                %q   = cq.state(4:7    , :);
                %r_g = R3(q) * cq.env.g; %
                %v   = cq.state(8:10   , :);
                %w   = cq.state(11:13  , :);
                %   state = [r, q, v, w, dot v, dot w]^T in R^19
                
                cq.thetaDD = zeros(28, cq.n_int);
                
                %%  Hesse - Matrix von Theta[1] 
                cq.thetaDD(1, :) = - 2 * cq.m * cq.env.g; %thetaD(1, :) nach q3 abgeleitet
                cq.thetaDD(2, :) =   2 * cq.m * cq.env.g; %thetaD(2, :) nach q4 abgeleitet
                cq.thetaDD(3, :) = - 2 * cq.m * cq.env.g; %thetaD(3, :) nach q1 abgeleitet
                cq.thetaDD(4, :) =   2 * cq.m * cq.env.g; %thetaD(4, :) nach q2 abgeleitet
                cq.thetaDD(5, :) = -     cq.m;     %thetaD(5, :) nach w3 abgeleitet
                cq.thetaDD(6, :) =       cq.m;     %thetaD(6, :) nach w2 abgeleitet
                cq.thetaDD(7, :) =       cq.m;     %thetaD(7, :) nach v3 abgeleitet
                cq.thetaDD(8, :) = -     cq.m;     %thetaD(8, :) nach v2 abgeleitet
                %% Hesse - Matrix von Theta[2]
                
                cq.thetaDD(9, :)  =  2 * cq.m * cq.env.g; %thetaD(9, :) nach q2 abgeleitet
                cq.thetaDD(10, :) =  2 * cq.m * cq.env.g; %thetaD(10, :) nach q1 abgeleitet
                cq.thetaDD(11, :) =  2 * cq.m * cq.env.g; %thetaD(11, :) nach q4 abgeleitet
                cq.thetaDD(12, :) =  2 * cq.m * cq.env.g; %thetaD(12, :) nach q3 abgeleitet
                cq.thetaDD(13, :) =      cq.m;     %thetaD(13, :) nach w3 abgeleitet
                cq.thetaDD(14, :) = -    cq.m;     %thetaD(14, :) nach w1 abgeleitet
                cq.thetaDD(15, :) = -    cq.m;     %thetaD(15, :) nach v3 abgeleitet
                cq.thetaDD(16, :) =      cq.m;     %thetaD(16, :) nach v1 abgeleitet
                
                %% Hesse - Matrix von Theta[3]
                
                cq.thetaDD(17, :) =  - 4 * cq.m * cq.env.g; %thetaD(17, :) nach q2 abgeleitet
                cq.thetaDD(18, :) =  - 4 * cq.m * cq.env.g; %thetaD(18, :) nach q3 abgeleitet
                cq.thetaDD(19, :) = -      cq.m ;           %thetaD(19, :) nach w2 abgeleitet
                cq.thetaDD(20, :) =        cq.m ;           %thetaD(20, :) nach w1 abgeleitet
                cq.thetaDD(21, :) =        cq.m;            %thetaD(21, :) nach v(2) abgeleitet
                cq.thetaDD(22, :) = -      cq.m;            %thetaD(22, :) nach v1 abgeleitet
                
                %% Hesse - Matrix von Theta[4]
                cq.thetaDD(23, :) =  (cq.I(3) - cq.I(2)) ;  % thetaD(23, :) nach w2 abgeleitet
                cq.thetaDD(24, :) =  (cq.I(3) - cq.I(2)) ;  % thetaD(24, :) nach w3 abgeleitet
                
                %% Hesse - Matrix von Theta[5]
                cq.thetaDD(25, :) =  (cq.I(1) - cq.I(3)) ;  % thetaD(25, :) nach w3 abgeleitet
                cq.thetaDD(26, :) =  (cq.I(1) - cq.I(3)) ;  % thetaD(26, :) nach w1 abgeleitet
                
                %% Hesse - Matrix von Theta[5]
                cq.thetaDD(27, :) =  (cq.I(2) - cq.I(1)) ;  % thetaD(27, :) nach w2 abgeleitet
                cq.thetaDD(28, :) =  (cq.I(2) - cq.I(1)) ;  % thetaD(28, :) nach w1 abgeleitet
                
                cq.isEmptyThetaDD = false;
                
            end
            ret = cq.thetaDD;
        end
        
        function ret = get.T(cq)
            
            if cq.isEmptyT
                %q   = cq.state(4:7    , :);
                %r_g = R3(q) * cq.env.g; %
                %v   = cq.state(8:10   , :);
                u   = cq.contr;
                w   = cq.state(11:13  , :);
                %   state = [r, q, v, w, dot v, dot w]^T in R^19
                
                cq.T = zeros(6, cq.n_int);
                
                cq.T(3, :) = cq.kT * sum(u.^2, 1); %% Testen ob richtig!!
                cq.T(4, :) = cq.kT * cq.d * (u(2, :).^2 - u(4, :).^2) ... 
                    + cq.I_M * w(2, :) .* (-u(1, :) +  u(2, :) - u(3, :) + u(4, :));

                cq.T(5, :) = cq.kT * cq.d * (u(3, :).^2 - u(1, :).^2) ... 
                    + cq.I_M * w(1, :) .* (u(1, :) -  u(2, :) + u(3, :) - u(4, :));

                cq.T(6, :) = cq.kQ * (-u(1, :).^2 +  u(2, :).^2 + u(3, :).^2 + u(4, :).^2);
                cq.isEmptyT = false;
            end
            
            ret = cq.T;
            
        end
        
        function ret = get.TD(cq)
            
            if cq.isEmptyTD
                %q   = cq.state(4:7    , :);
                %r_g = R3(q) * cq.env.g; %
                %v   = cq.state(8:10   , :);
                u   = cq.contr;
                w   = cq.state(11:13  , :);
                %   state = [r, q, v, w, dot v, dot w]^T in R^19
                
                cq.TD = zeros(18, cq.n_int);
                
                %%   Ableitung T[3] nach ...
                cq.TD(1, :) =   2 * cq.kT * u(1, :); % u1
                cq.TD(2, :) =   2 * cq.kT * u(2, :); % u2
                cq.TD(3, :) =   2 * cq.kT * u(3, :); % u3
                cq.TD(4, :) =   2 * cq.kT * u(4, :); % u4
                
                %%   Ableitung T[4] nach ...BB: Ableitungsreihenfolge vertauscht
                cq.TD(5, :) =     cq.I_M * (-u(1, :) + u(2, :) - u(3, :) + u(4, :));%w2
                cq.TD(6, :) =   - cq.I_M * w(2, :); % u1
                cq.TD(7, :) =     cq.I_M * w(2, :) + 2 * cq.kT * cq.d * u(4, :);
                cq.TD(8, :) =   - cq.I_M * w(2, :); % u3
                cq.TD(9, :) =     cq.I_M * w(2, :) - 2 * cq.kT * cq.d * u(2, :); % u4
                
                
                %%   Ableitung T[5] nach ...BB: Ableitungsreihenfolge vertauscht
                cq.TD(10, :) =     cq.I_M * (u(1, :) - u(2, :) + u(3, :) - u(4, :));%w1
                cq.TD(11, :) =     cq.I_M * w(1, :) - 2 * cq.kT * cq.d * u(1, :); %u1
                cq.TD(12, :) =   - cq.I_M * w(1, :); % u2
                cq.TD(13, :) =     cq.I_M * w(1, :) + 2 * cq.kT * cq.d * u(3, :); % u3
                cq.TD(14, :) =   - cq.I_M * w(1, :); % u4
                
                %%   Ableitung T[6] nach ...
                cq.TD(15, :) = - 2 * cq.kQ * u(1, :); %u1
                cq.TD(16, :) =   2 * cq.kQ * u(2, :); %u2
                cq.TD(17, :) = - 2 * cq.kQ * u(3, :); %u3
                cq.TD(18, :) =   2 * cq.kQ * u(4, :); %u4
                
                cq.isEmptyTD = false;
            end
            
            ret = cq.TD;
            
        end
        
        function ret = get.TDD(cq)
            
            if cq.isEmptyTDD
                %q   = cq.state(4:7    , :);
                %r_g = R3(q) * cq.env.g; %
                %v   = cq.state(8:10   , :);
                %u   = cq.control;
                %w   = cq.state(11:12  , :);
                %   state = [r, q, v, w, dot v, dot w]^T in R^19
                
                cq.TDD = zeros(28, cq.n_int);
                
                %%  Hesse - Matrix von T[3] 
                cq.TDD(1, :) =    2 * cq.kT; %TD(1, :) nach u1 abgeleitet
                cq.TDD(2, :) =    2 * cq.kT; %TD(2, :) nach u2 abgeleitet
                cq.TDD(3, :) =    2 * cq.kT; %TD(3, :) nach u3 abgeleitet
                cq.TDD(4, :) =    2 * cq.kT; %TD(4, :) nach u4 abgeleitet
                
                %%  Hesse - Matrix von T[4] 
                cq.TDD(5, :) =  - cq.I_M;           %TD(5, :) nach w2 abgeleitet
                cq.TDD(6, :) =    cq.I_M;           %TD(6, :) nach w2 abgeleitet %BB: Reihenfolge geändert
                cq.TDD(7, :) =    2 * cq.kT * cq.d; %TD(6, :) nach u2 abgeleitet
                cq.TDD(8, :) =  - cq.I_M;           %TD(7, :) nach w2 abgeleitet
                cq.TDD(9, :) =  - 2 * cq.kT * cq.d; %TD(8, :) nach u4 abgeleitet
                
                cq.TDD(10, :) =   cq.I_M;           %TD(8, :) nach w2 abgeleitet
                cq.TDD(11, :) = - cq.I_M;           %TD(9, :) nach u1 abgeleitet
                cq.TDD(12, :) =  	cq.I_M;         %TD(9, :) nach u2 abgeleitet
                cq.TDD(13, :) = - cq.I_M;           %TD(9, :) nach u3 abgeleitet
                cq.TDD(14, :) =   cq.I_M;           %TD(9, :) nach u4 abgeleitet
                
                %%  Hesse - Matrix von T[5] 
                
                cq.TDD(15, :) =    cq.I_M;           %TD(10, :) nach w1 abgeleitet %BB: Reihenfolge geändert
                cq.TDD(16, :) =  - 2 * cq.kT * cq.d; %TD(10, :) nach u1 abgeleitet
                cq.TDD(17, :) =  - cq.I_M;           %TD(11, :) nach w1 abgeleitet
                cq.TDD(18, :) =    cq.I_M;           %TD(12, :) nach w1 abgeleitet %BB: Reihenfolge geändert
                cq.TDD(19, :) =    2 * cq.kT * cq.d; %TD(12, :) nach u3 abgeleitet 
                cq.TDD(20, :) =  - cq.I_M;           %TD(13, :) nach w1 abgeleitet
                cq.TDD(21, :) =    cq.I_M;           %TD(14, :) nach u1 abgeleitet
                cq.TDD(22, :) =  - cq.I_M;           %TD(14, :) nach u2 abgeleitet
                cq.TDD(23, :) =    cq.I_M;           %TD(14, :) nach u3 abgeleitet
                cq.TDD(24, :) =  - cq.I_M;           %TD(14, :) nach u4 abgeleitet
                
                %%  Hesse - Matrix von T[6]
                cq.TDD(25, :) =  - 2 * cq.kQ;           %TD(14, :) nach u1 abgeleitet
                cq.TDD(26, :) =    2 * cq.kQ;           %TD(14, :) nach u2 abgeleitet
                cq.TDD(27, :) =  - 2 * cq.kQ;           %TD(14, :) nach u3 abgeleitet
                cq.TDD(28, :) =    2 * cq.kQ;           %TD(14, :) nach u4 abgeleitet
                
                cq.isEmptyTDD = false;
            end
            
            ret = cq.TDD;
            
        end
        
        function res = get.Rges(obj)
            if obj.isEmptyRges
                n = obj.n_int;
                obj.Rges = cell(1, n);
                
                for i = 1:n
                    q = obj.state(4:7, i);
                    obj.Rges{i} = obj.R(q);
                end
                
                obj.isEmptyRges = false;
            end
            res = obj.Rges;
        end
        
        function res = get.RvD(obj)
            
            
            if obj.isEmptyRvD
                q = obj.state(4:7, :); 
                v = obj.state(8:10, :);
                
                obj.RvD = zeros(12, obj.n_int);
                
                obj.RvD(1, :) = -2*q(4, :).*v(2, :)+2*q(3, :).*v(3, :);
                obj.RvD(2, :) =  2*q(3, :).*v(2, :)+2*q(4, :).*v(3, :);
                obj.RvD(3, :) =  -4*q(3, :).*v(1, :)+2*q(2, :).*v(2, :)+2*q(1, :).*v(3, :);
                obj.RvD(4, :) = -4*q(4, :).*v(1, :)-2*q(1, :).*v(2, :)+2*q(2, :).*v(3, :);
                
                obj.RvD(5, :)   =  2*q(4, :).*v(1, :)-2*q(2, :).*v(3, :);
                obj.RvD(6, :)   =  2*q(3, :).*v(1, :)-4*q(2, :).*v(2, :)-2*q(1, :).*v(3, :);
                obj.RvD(7, :)   =  2*q(2, :).*v(1, :)+2*q(4, :).*v(3, :);
                obj.RvD(8, :)   =  2*q(1, :).*v(1, :)-4*q(4, :).*v(2, :)+2*q(3, :).*v(3, :);
                
                obj.RvD(9, :)   =  -2*q(3, :).*v(1, :)+2*q(2, :).*v(2, :);
                obj.RvD(10, :)  =  2*q(4, :).*v(1, :)+2*q(1, :).*v(2, :)-4*q(2, :).*v(3, :);
                obj.RvD(11, :)  =  -2*q(1, :).*v(1, :)+2*q(4, :).*v(2, :)-4*q(3, :).*v(3, :);
                obj.RvD(12, :)  =  2*q(2, :).*v(1, :)+2*q(3, :).*v(2, :);

                %obj.RvD(1:3, 5:7) = obj.R(q);
                
                obj.isEmptyRvD = false;
            end
            res = obj.RvD;
        end
        
        function res = get.RvDD(obj)
            if obj.isEmptyRvDD
                q = obj.state(4:7, :); 
                v = obj.state(8:10, :);
                
                obj.RvDD = zeros(90, obj.n_int);
                
                %  Hesse - Matrix von Rv[1] 
                obj.RvDD(1, :) = 2*v(3, :);
                obj.RvDD(2, :) =  -2*v(2, :);
                obj.RvDD(3, :) = -2*q(4, :);
                obj.RvDD(4, :) =  2*q(3, :);
                obj.RvDD(5, :) = 2*v(2, :); 
                obj.RvDD(6, :) =  2*v(3, :);
                obj.RvDD(7, :) = 2*q(3, :);
                obj.RvDD(8, :) =  2*q(4, :);
                obj.RvDD(9, :) = 2*v(3, :);
                obj.RvDD(10, :) = 2*v(2, :);
                obj.RvDD(11, :) =  -4*v(1, :);
                obj.RvDD(12, :) =  -4*q(3, :);
                obj.RvDD(13, :) =  2*q(2, :);
                obj.RvDD(14, :) = 2*q(1, :);
                obj.RvDD(15, :) =  -2*v(2, :);
                obj.RvDD(16, :) = 2*v(3, :);
                obj.RvDD(17, :) =  -4*v(2, :);
                obj.RvDD(18, :) =  -4*q(4, :);
                obj.RvDD(19, :) =  -2*q(1, :);
                obj.RvDD(20, :) =  2*q(2, :);
                obj.RvDD(21, :) =  -4*q(3, :);
                obj.RvDD(22, :) =  -4*q(4, :);
                obj.RvDD(23, :) =  -2*q(4, :); 
                obj.RvDD(24, :) =  2*q(3, :);
                obj.RvDD(25, :) = 2*q(2, :);
                obj.RvDD(26, :) =  -2*q(1, :);
                obj.RvDD(27, :) =  2*q(3, :);
                obj.RvDD(28, :) =  2*q(4, :); 
                obj.RvDD(29, :) =  2*q(1, :); 
                obj.RvDD(30, :) =  2*q(2, :);
                
                %  Hesse - Matrix von Rv[2]
                               
                obj.RvDD(31, :) =  -2*v(3, :);%
                obj.RvDD(32, :) =   2*v(1, :);%
                obj.RvDD(33, :) =  2*q(4, :);%
                obj.RvDD(34, :) =   -2*q(2, :);%
                obj.RvDD(35, :) =  -2*v(3, :);%
                obj.RvDD(36, :) =   -4*v(2, :);%
                obj.RvDD(37, :) =  2*v(1, :);%
                obj.RvDD(38, :) =   2*q(3, :);%
                obj.RvDD(39, :) =  -4*q(2, :);%
                obj.RvDD(41, :) =  -2*q(1, :);%
                obj.RvDD(42, :) =  2*v(1, :);%
                obj.RvDD(43, :) =   2*v(3, :);%
                obj.RvDD(44, :) =   2*q(2, :);%
                obj.RvDD(45, :) =  2*q(4, :);%
                obj.RvDD(46, :) =  2*v(1, :);%
                obj.RvDD(47, :) =   2*v(3, :);%
                obj.RvDD(48, :) =   -4*v(2, :);%
                obj.RvDD(49, :) =  2*q(1, :);%
                obj.RvDD(50, :) =   -4*q(4, :);%
                obj.RvDD(51, :) =   2*q(3, :);%
                obj.RvDD(52, :) =   2*q(4, :);%
                obj.RvDD(53, :) =   2*q(3, :);%
                obj.RvDD(54, :) =   2*q(3, :);%
                obj.RvDD(55, :) =   2*q(1, :);%
                obj.RvDD(56, :) =   -4*q(4, :);%
                obj.RvDD(57, :) =   -2*q(2, :);%
                obj.RvDD(58, :) =   -2*q(1, :);% 
                obj.RvDD(59, :) =   2*q(4, :);%
                obj.RvDD(60, :) =   2*q(3, :);%
                
                %  Hesse - Matrix von Rv[3]
                %
                obj.RvDD(61, :) =  2*v(2, :);
                obj.RvDD(62, :) =   -2*v(1, :);
                obj.RvDD(63, :) =   -2*q(3, :);
                obj.RvDD(64, :) =   2*q(2, :);
                obj.RvDD(65, :) =  2*v(2, :);%
                obj.RvDD(66, :) =    -4*v(3, :);
                obj.RvDD(67, :) =   2*v(1, :);
                obj.RvDD(68, :) =   2*q(4, :);%
                obj.RvDD(69, :) =   2*q(1, :);
                obj.RvDD(70, :) =   -4*q(2, :);
                obj.RvDD(71, :) =   -2*v(1, :);
                obj.RvDD(72, :) =   -4*v(3, :);%
                obj.RvDD(73, :) =   2*v(2, :);
                obj.RvDD(74, :) =  -2*q(1, :);
                obj.RvDD(75, :) =  2*q(4, :);%
                
                obj.RvDD(76, :) =   -4*q(3, :);
                obj.RvDD(77, :) =   2*v(1, :);
                obj.RvDD(78, :) =  2*v(2, :);
                obj.RvDD(79, :) =   2*q(2, :);%
                obj.RvDD(80, :) =   2*q(3, :);
                obj.RvDD(81, :) =   -2*q(3, :);
                obj.RvDD(82, :) =   2*q(4, :);
                obj.RvDD(83, :) =  -2*q(1, :); %
                obj.RvDD(84, :) =   2*q(2, :);
                obj.RvDD(85, :) =   2*q(2, :);
                obj.RvDD(86, :) =   2*q(1, :);
                obj.RvDD(87, :) =   2*q(4, :);
                obj.RvDD(88, :) =   2*q(3, :); %
                obj.RvDD(89, :) =  -4*q(2, :); 
                obj.RvDD(90, :) =  -4*q(3, :);
                
                obj.isEmptyRvDD = false;
            end
            res = obj.RvDD;
        end
        
        function res = get.Q(obj)
            
            if obj.isEmptyQ
                q   = obj.state(4:7    , :);
                w   = obj.state(11:13  , :);
                
                obj.Q = zeros(4, obj.n_int);
                obj.Q(1, :) = -q(2, :).*w(1, :)-q(3, :).*w(2, :)-q(4, :).*w(3, :);
                obj.Q(2, :) = q(1, :).*w(1, :)+q(3, :).*w(3, :)-q(4, :).*w(2, :);
                obj.Q(3, :) = q(1, :).*w(2, :)-q(2, :).*w(3, :)+q(4, :).*w(1, :);
                obj.Q(4, :) = q(1, :).*w(3, :)+q(2, :).*w(2, :)-q(3, :).*w(1, :);
                obj.Q = 1/2 * obj.Q;
                obj.isEmptyQ = false;
            end
            res = obj.Q;
        end
        
        function res = get.QD(obj)
            
            if obj.isEmptyQD
                q   = obj.state(4:7    , :);
                w   = obj.state(11:13  , :);
                
                obj.QD = zeros(24, obj.n_int);
                
                obj.QD(1, :) = -(1/2)*w(1, :);
                obj.QD(2, :) = -(1/2)*w(2, :);
                obj.QD(3, :) = -(1/2)*w(3, :);
                obj.QD(4, :) = -(1/2)*q(2, :);
                
                obj.QD(5, :) = -(1/2)*q(3, :);
                obj.QD(6, :) = -(1/2)*q(4, :);
                obj.QD(7, :) = (1/2)*w(1, :);
                obj.QD(8, :) = (1/2)*w(3, :);
                
                obj.QD(9, :) = -(1/2)*w(2, :);
                obj.QD(10, :)= (1/2)*q(1, :);
                obj.QD(11, :)= -(1/2)*q(4, :);
                obj.QD(12, :)= (1/2)*q(3, :);
                obj.QD(13, :)= (1/2)*w(2, :);
                obj.QD(14, :)= -(1/2)*w(3, :);
                obj.QD(15, :)= (1/2)*w(1, :);
                obj.QD(16, :)= (1/2)*q(4, :);
                obj.QD(17, :)= (1/2)*q(1, :);
                obj.QD(18, :)= -(1/2)*q(2, :);
                obj.QD(19, :)= (1/2)*w(3, :);
                obj.QD(20, :)= (1/2)*w(2, :);
                obj.QD(21, :)= -(1/2)*w(1, :);
                obj.QD(22, :)= -(1/2)*q(3, :);
                obj.QD(23, :)= (1/2)*q(2, :);
                obj.QD(24, :)= (1/2)*q(1, :);
                
                obj.isEmptyQD = false;
            end
            
            res = obj.QD;
        end
        
        function res = get.QDD(obj)
            
            if obj.isEmptyQDD
                
                obj.QDD = zeros(24, obj.n_int);
                
                obj.QDD(1, :) = -1/2;
                obj.QDD(2, :) = -1/2;
                obj.QDD(3, :) = -1/2;
                obj.QDD(4, :) = -1/2;
                obj.QDD(5, :) = -1/2;
                obj.QDD(6, :) = -1/2;
                
                obj.QDD(7, :) = 1/2;
                obj.QDD(8, :) = 1/2;
                obj.QDD(9, :) = -1/2;
                obj.QDD(10, :) = 1/2;
                obj.QDD(11, :) = -1/2;
                obj.QDD(12, :) = 1/2;
                
                obj.QDD(13, :) = 1/2;
                obj.QDD(14, :) = -1/2;
                obj.QDD(15, :) = 1/2;
                obj.QDD(16, :) = 1/2;
                obj.QDD(17, :) = 1/2;
                obj.QDD(18, :) = -1/2;
                
                obj.QDD(19, :) = 1/2;
                obj.QDD(20, :) = 1/2;
                obj.QDD(21, :) = -1/2;
                obj.QDD(22, :) = -1/2;
                obj.QDD(23, :) = 1/2;
                obj.QDD(24, :) = 1/2;
                
                obj.isEmptyQDD = false;
            end
            res = obj.QDD;
            
        end
        
        %% 
        
        function res = getM(obj,ind)
            % the actual mass matrix at time instance ind
            res       = obj.M; % Massenmatrix ist konstant
            
            %res         = zeros(2,2);
            %res(1:2,1)  = M_val(1:2,ind);
            %res(1,2)    = res(2,1);
            %res(2,2)    = M_val(3,ind);
        end
        
        function res = getMD(obj,ind)
            % the actual 1st derivative of mass matrix at time instance ind
            res             = cell(1,1);  % Massenmatrix ist konstant         
        end
        
        function res = getMDD(obj,ind)
            % the actual 2nd derivative of mass matrix at time instance ind
            res             = cell(1,1); % Massenmatrix ist konstant
        end       
        
        function res = getTheta(obj, varargin)
            % the actual coriolis terms
            
            if nargin == 1
                res = obj.theta;
            elseif nargin == 2
                res = obj.theta(:, varargin{1});
            else
                error('Anzahl der Parameter ist falsch');
            end
        end
        
        function res = getThetaD(obj,ind)
            % the actual 1st derivative of coriolis terms
            res = zeros(6,10); %Jacobianmatrix von Theta
            
            thetaD_val  = obj.thetaD;
            
            res(1, 1:4)  = thetaD_val(1:4, ind)';
            res(1, 6:7)  = thetaD_val(5:6, ind)';
            res(1, 9:10) = thetaD_val(7:8, ind)';
            
            res(2, 1:5)  = thetaD_val(9:13, ind)';
            res(2, 7:8)  = thetaD_val(14:15, ind)';
            res(2, 10)   = thetaD_val(16, ind)';
            
            res(3, 2:3)  = thetaD_val(17:18, ind)';
            res(3, 5:6)  = thetaD_val(19:20, ind)';
            res(3, 8:9)  = thetaD_val(21:22, ind)';
            
            res(4, 9:10) = thetaD_val(23:24, ind)';
            
            res(5, 8) = thetaD_val(25, ind)';
            res(5, 10) = thetaD_val(26, ind)';
            
            res(6, 8:9) = thetaD_val(27:28, ind)';
        end 
        
        function res = getThetaDD(obj,ind)
            % the actual 2nd derivative of coriolis terms
            res         = cell(10, 10);
            
            iTD1 = 1;
            iTD2 = 1;
            iTD3 = 1;
            iTD4 = 1;
            iTD5 = 1;
            iTD6 = 1;

            thetaDD_val = obj.thetaDD;
            
            res{1, 2}    = zeros(6, 1);
            res{1, 2}(2) = thetaDD_val(obj.startIndexHesseThetaDD2 + iTD2, ind);
            iTD2 = iTD2 + 1; %TD2 = 2
            
            res{1, 3}    = zeros(6, 1);
            res{1, 3}(1) = thetaDD_val(obj.startIndexHesseThetaDD1 + iTD1,  ind);
            iTD1 = iTD1 + 1; %TD1 = 2
            
            res{2, 1}    = zeros(6, 1);
            res{2, 1}(2) = thetaDD_val(obj.startIndexHesseThetaDD2 + iTD2,  ind);
            iTD2 = iTD2 + 1; %TD2 = 3
            
            res{2, 2}    = zeros(6, 1);
            res{2, 2}(3) = thetaDD_val(obj.startIndexHesseThetaDD3 + iTD3,  ind);
            iTD3 = iTD3 + 1; %TD3 = 2
            
            res{2, 4}    = zeros(6, 1);
            res{2, 4}(1) = thetaDD_val(obj.startIndexHesseThetaDD1 + iTD1,  ind);
            iTD1 = iTD1 + 1; %TD1 = 3
            
            res{3, 1}    = zeros(6, 1);
            res{3, 1}(1) = thetaDD_val(obj.startIndexHesseThetaDD1 + iTD1,  ind);
            iTD1 = iTD1 + 1; %TD1 = 4
            
            res{3, 3}    = zeros(6, 1);
            res{3, 3}(3) = thetaDD_val(obj.startIndexHesseThetaDD3 + iTD3,  ind);
            iTD3 = iTD3 + 1; %TD3 = 3
            
            res{3, 4}    = zeros(6, 1);
            res{3, 4}(2) = thetaDD_val(obj.startIndexHesseThetaDD2 + iTD2,  ind);
            iTD2 = iTD2 + 1; %TD2 = 4
            
            res{4, 2}    = zeros(6, 1);
            res{4, 2}(1) = thetaDD_val(obj.startIndexHesseThetaDD1 + iTD1,  ind);
            iTD1 = iTD1 + 1; %TD1 = 5
            
            res{4, 3}    = zeros(6, 1);
            res{4, 3}(2) = thetaDD_val(obj.startIndexHesseThetaDD2 + iTD2,  ind);
            iTD2 = iTD2 + 1; %TD2 = 5
            
            res{5, 9}    = zeros(6, 1);
            res{5, 9}(3) = thetaDD_val(obj.startIndexHesseThetaDD3 + iTD3,  ind);
            iTD3 = iTD3 + 1; %TD3 = 4
            
            res{5, 10}    = zeros(6, 1);
            res{5, 10}(2) = thetaDD_val(obj.startIndexHesseThetaDD2 + iTD2,  ind);
            iTD2 = iTD2 + 1; %TD2 = 6
            
            res{6, 8}    = zeros(6, 1);
            res{6, 8}(3) = thetaDD_val(obj.startIndexHesseThetaDD3 + iTD3,  ind);
            iTD3 = iTD3 + 1; %TD3 = 5
            
            res{6, 10}    = zeros(6, 1);
            res{6, 10}(1) = thetaDD_val(obj.startIndexHesseThetaDD1 + iTD1,  ind);
            iTD1 = iTD1 + 1; %TD1 = 6
            
            res{7, 8}    = zeros(6, 1);
            res{7, 8}(2) = thetaDD_val(obj.startIndexHesseThetaDD2 + iTD2,  ind);
            iTD2 = iTD2 + 1; %TD2 = 7
            
            res{7, 9}    = zeros(6, 1);
            res{7, 9}(1) = thetaDD_val(obj.startIndexHesseThetaDD1 + iTD1,  ind);
            iTD1 = iTD1 + 1; %TD1 = 7
            
            res{8, 6}    = zeros(6, 1);
            res{8, 6}(3) = thetaDD_val(obj.startIndexHesseThetaDD3 + iTD3,  ind);
            iTD3 = iTD3 + 1; %TD3 = 6
            
            res{8, 7}    = zeros(6, 1);
            res{8, 7}(2) = thetaDD_val(obj.startIndexHesseThetaDD2 + iTD2,  ind);
            iTD2 = iTD2 + 1; %TD2 = 8
            
            res{8, 9}    = zeros(6, 1);
            res{8, 9}(6) = thetaDD_val(obj.startIndexHesseThetaDD6 + iTD6,  ind);
            iTD6 = iTD6 + 1; %TD6 = 2
            
            res{8, 10}    = zeros(6, 1);
            res{8, 10}(5) = thetaDD_val(obj.startIndexHesseThetaDD5 + iTD5,  ind);
            iTD5 = iTD5 + 1; %TD5 = 2
            
            res{9, 5}    = zeros(6, 1);
            res{9, 5}(3) = thetaDD_val(obj.startIndexHesseThetaDD3 + iTD3,  ind);
            
            res{9, 7}    = zeros(6, 1);
            res{9, 7}(1) = thetaDD_val(obj.startIndexHesseThetaDD1 + iTD1,  ind);
            iTD1 = iTD1 + 1; %TD1 = 8
            
            res{9, 8}    = zeros(6, 1);
            res{9, 8}(6) = thetaDD_val(obj.startIndexHesseThetaDD6 + iTD6,  ind);
            
            res{9, 10}    = zeros(6, 1);
            res{9, 10}(4) = thetaDD_val(obj.startIndexHesseThetaDD4 + iTD4,  ind);
            iTD4 = iTD4 + 1; %TD4 = 2
            
            res{10, 5}    = zeros(6, 1);
            res{10, 5}(2) = thetaDD_val(obj.startIndexHesseThetaDD2 + iTD2,  ind);
                        
            res{10, 6}    = zeros(6, 1);
            res{10, 6}(1) = thetaDD_val(obj.startIndexHesseThetaDD1 + iTD1,  ind);
            
            res{10, 8}    = zeros(6, 1);
            res{10, 8}(5) = thetaDD_val(obj.startIndexHesseThetaDD5 + iTD5,  ind);
            
            res{10, 9}    = zeros(6, 1);
            res{10, 9}(4) = thetaDD_val(obj.startIndexHesseThetaDD4 + iTD4,  ind);
            
        end
        
        function res = getT(obj, ind)
            res = obj.T(:, ind);
        end
        
        function res = getTD(obj, ind)
            
            res = zeros(6, 7); %Jacobianmatrix von T BB: Anzahl der Ableitungen von 8 auf 7 reduziert
            
            TD_val  = obj.TD;
            %BB: Ableitungsreihenfolge vertauscht
            res(3, 4:7)  = TD_val(1:4, ind)';
            res(4, 2)    = TD_val(5, ind)';
            res(4, 4:7)  = TD_val(6:9, ind)';
            
            res(5, 1)  = TD_val(10, ind)';
            res(5, 4:7)  = TD_val(11:14, ind)';
            res(6, 4:7)  = TD_val(15:18, ind)';
        end
        
       
        function res = getTDD(obj, ind)
            res         = cell(7, 7);
            
            %iTD1 = 1;
            %iTD2 = 1;
            iTD3 = 1;
            iTD4 = 1;
            iTD5 = 1;
            iTD6 = 1;

            TDD_val = obj.TDD;
            % 7 x 7
            % Zeile 1
            res{1, 4}  = zeros(6, 1);
            res{1, 4}(5) = TDD_val(obj.startIndexHesseTDD5 + iTD5, ind);
            iTD5 = iTD5 + 1; %TD5 = 2
            
            res{1, 5}  = zeros(6, 1);
            res{1, 5}(5) = TDD_val(obj.startIndexHesseTDD5 + iTD5, ind);
            iTD5 = iTD5 + 1; %TD5 = 3
            
            res{1, 6}  = zeros(6, 1);
            res{1, 6}(5) = TDD_val(obj.startIndexHesseTDD5 + iTD5, ind);
            iTD5 = iTD5 + 1; %TD5 = 4
            
            res{1, 7}  = zeros(6, 1);
            res{1, 7}(5) = TDD_val(obj.startIndexHesseTDD5 + iTD5, ind);
            iTD5 = iTD5 + 1; %TD5 = 5
            
            % Zeile 2
            
            res{2, 4}  = zeros(6, 1);
            res{2, 4}(4) = TDD_val(obj.startIndexHesseTDD4 + iTD4, ind);
            iTD4 = iTD4 + 1; %TD4 = 2
            
            res{2, 5}  = zeros(6, 1);
            res{2, 5}(4) = TDD_val(obj.startIndexHesseTDD4 + iTD4, ind);
            iTD4 = iTD4 + 1; %TD4 = 3
            
            res{2, 6}  = zeros(6, 1);
            res{2, 6}(4) = TDD_val(obj.startIndexHesseTDD4 + iTD4, ind);
            iTD4 = iTD4 + 1; %TD4 = 4
            
            res{2, 7}  = zeros(6, 1);
            res{2, 7}(4) = TDD_val(obj.startIndexHesseTDD4 + iTD4, ind);
            iTD4 = iTD4 + 1; %TD4 = 5
            
            % Zeile 4
            
            res{4, 1}  = zeros(6, 1);
            res{4, 1}(5) = TDD_val(obj.startIndexHesseTDD5 + iTD5, ind);
            iTD5 = iTD5 + 1; %TD5 = 6
            
            res{4, 2}  = zeros(6, 1);
            res{4, 2}(4) = TDD_val(obj.startIndexHesseTDD4 + iTD4, ind);
            iTD4 = iTD4 + 1; %TD4 = 6
            
            
            res{4, 4}  = zeros(6, 1);
            res{4, 4}(3) = TDD_val(obj.startIndexHesseTDD3 + iTD3, ind);
            iTD3 = iTD3 + 1; %TD3 = 2
            
            res{4, 4}(5) = TDD_val(obj.startIndexHesseTDD5 + iTD5, ind);
            iTD5 = iTD5 + 1; %TD5 = 7
            
            res{4, 4}(6) = TDD_val(obj.startIndexHesseTDD6 + iTD6, ind);
            iTD6 = iTD6 + 1; %TD6 = 2
            
            % Zeile 5
            res{5, 1}  = zeros(6, 1);
            res{5, 1}(5) = TDD_val(obj.startIndexHesseTDD5 + iTD5, ind);
            iTD5 = iTD5 + 1; %TD5 = 8
            
            res{5, 2}  = zeros(6, 1);
            res{5, 2}(4) = TDD_val(obj.startIndexHesseTDD4 + iTD4, ind);
            iTD4 = iTD4 + 1; %TD4 = 7
            
            res{5, 5}  = zeros(6, 1);
            res{5, 5}(3) = TDD_val(obj.startIndexHesseTDD3 + iTD3, ind);
            iTD3 = iTD3 + 1; %TD3 = 3
            
            res{5, 5}(4) = TDD_val(obj.startIndexHesseTDD4 + iTD4, ind);
            iTD4 = iTD4 + 1; %TD4 = 8
            
            res{5, 5}(6) = TDD_val(obj.startIndexHesseTDD6 + iTD6, ind);
            iTD6 = iTD6 + 1; %TD6 = 3
            
            % Zeile 6
            
            res{6, 1}  = zeros(6, 1);
            res{6, 1}(5) = TDD_val(obj.startIndexHesseTDD5 + iTD5, ind);
            iTD5 = iTD5 + 1; %TD5 = 9
            
            res{6, 2}  = zeros(6, 1);
            res{6, 2}(4) = TDD_val(obj.startIndexHesseTDD4 + iTD4, ind);
            iTD4 = iTD4 + 1; %TD4 = 9
            
            
            res{6, 6}  = zeros(6, 1);
            res{6, 6}(3) = TDD_val(obj.startIndexHesseTDD3 + iTD3, ind);
            iTD3 = iTD3 + 1; %TD3 = 4
            
            res{6, 6}(5) = TDD_val(obj.startIndexHesseTDD5 + iTD5, ind);
            iTD5 = iTD5 + 1; %TD5 = 10
            
            res{6, 6}(6) = TDD_val(obj.startIndexHesseTDD6 + iTD6, ind);
            iTD6 = iTD6 + 1; %TD6 = 4
            
            
            % Zeile 7
            
            res{7, 1}  = zeros(6, 1);
            res{7, 1}(5) = TDD_val(obj.startIndexHesseTDD5 + iTD5, ind);
            %iTD5 = iTD5 + 1; %TD5 = 11
            
            res{7, 2}  = zeros(6, 1);
            res{7, 2}(4) = TDD_val(obj.startIndexHesseTDD4 + iTD4, ind);
            iTD4 = iTD4 + 1; %TD4 = 10
            
            
            
            res{7, 7}  = zeros(6, 1);
            res{7, 7}(3) = TDD_val(obj.startIndexHesseTDD3 + iTD3, ind);
            %iTD3 = iTD3 + 1; %TD3 = 5
            
            res{7, 7}(4) = TDD_val(obj.startIndexHesseTDD4 + iTD4, ind);
            %iTD4 = iTD4 + 1; %TD4 = 11
            
            res{7, 7}(6) = TDD_val(obj.startIndexHesseTDD6 + iTD6, ind);
            %iTD6 = iTD6 + 1; %TD6 = 5
            

        end
        
        function res = getRv(obj, ind)
            v   = obj.state(8:10, ind);
            res = obj.Rges{ind} * v;
        end
        
        function res = getRvD(obj, ind)
            
            res = zeros(3, 7);
            
            Rv_val = obj.RvD;
            
            res(1, 1:4) = Rv_val(1:4, ind)';
            res(2, 1:4) = Rv_val(5:8, ind)';
            res(3, 1:4) = Rv_val(9:12, ind)';
            res(1:3, 5:7) = obj.Rges{ind};
            
        end
        
        function res = getRvDD(obj, ind)
            
            res = cell(7, 7);
            %Matrix(7, 7, {(1, 1) = 0, (1, 2) = 0, (1, 3) = 2*v3, (1, 4) = -2*v2, (1, 5) = 0, (1, 6) = -2*q[4], (1, 7) = 2*q[3], (2, 1) = 0, (2, 2) = 0, (2, 3) = 2*v2, (2, 4) = 2*v3, (2, 5) = 0, (2, 6) = 2*q[3], (2, 7) = 2*q[4], (3, 1) = 2*v3, (3, 2) = 2*v2, (3, 3) = -4*v1, (3, 4) = 0, (3, 5) = -4*q[3], (3, 6) = 2*q[2], (3, 7) = 2*q[1], (4, 1) = -2*v2, (4, 2) = 2*v3, (4, 3) = 0, (4, 4) = -4*v1, (4, 5) = -4*q[4], (4, 6) = -2*q[1], (4, 7) = 2*q[2], (5, 1) = 0, (5, 2) = 0, (5, 3) = -4*q[3], (5, 4) = -4*q[4], (5, 5) = 0, (5, 6) = 0, (5, 7) = 0, (6, 1) = -2*q[4], (6, 2) = 2*q[3], (6, 3) = 2*q[2], (6, 4) = -2*q[1], (6, 5) = 0, (6, 6) = 0, (6, 7) = 0, (7, 1) = 2*q[3], (7, 2) = 2*q[4], (7, 3) = 2*q[1], (7, 4) = 2*q[2], (7, 5) = 0, (7, 6) = 0, (7, 7) = 0})
            %startIndexHesseRDD1 = 0;
            %startIndexHesseRDD2 = 30;
            %startIndexHesseRDD3 = 60;
            
            iRvD1 = 1;
            iRvD2 = 1;
            iRvD3 = 1;
            
            RvD_val = obj.RvDD;
            
            res{1, 2}  = zeros(3, 1);
            res{1, 2}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 2
            
            res{1, 2}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 2
            
            %%%
            res{1, 3}  = zeros(3, 1);
            res{1, 3}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 2
            
            res{1, 3}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 3
            
            %%%
            res{1, 4}  = zeros(3, 1);
            res{1, 4}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 3
            
            res{1, 4}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 3
            
            %%%
            res{1, 5}  = zeros(3, 1);
            res{1, 5}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 4
            
            res{1, 5}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 4
            
            %%%
            res{1, 6}  = zeros(3, 1);
            res{1, 6}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 4
            
            res{1, 6}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 5
            
            %%%
            res{1, 7}  = zeros(3, 1);
            res{1, 7}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 5
            
            res{1, 7}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 5
            
            %%%
            res{2, 1}  = zeros(3, 1);
            res{2, 1}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 6
            
            res{2, 1}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 6
            
            %%%
            
            res{2, 2}  = zeros(3, 1);
            res{2, 2}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 7
            
            res{2, 2}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 7
            
            %%%
            res{2, 3}  = zeros(3, 1);
            res{2, 3}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 6
            
            res{2, 3}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 8
            
            %%%
            res{2, 4}  = zeros(3, 1);
            res{2, 4}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 7
            
            res{2, 4}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 8
            
            %%%
            
            res{2, 5}  = zeros(3, 1);
            res{2, 5}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 9
            
            res{2, 5}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 9
            
            %%%
            
            res{2, 6}  = zeros(3, 1);
            res{2, 6}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 8
            
            res{2, 6}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 10
            
            res{2, 6}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 10
            
            %%%
            
            res{2, 7}  = zeros(3, 1);
            res{2, 7}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 9
            
            res{2, 7}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 11
            
            res{2, 7}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 11
            
            %%%
            res{3, 1}  = zeros(3, 1);
            res{3, 1}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 10
            
            res{3, 1}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 12
            
            %%%
            
            res{3, 2}  = zeros(3, 1);
            res{3, 2}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 11
            
            res{3, 2}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 12
            
            %%%
            res{3, 3}  = zeros(3, 1);
            res{3, 3}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 12
            
            res{3, 3}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 13
            
            %%%
            
            res{3, 4}  = zeros(3, 1);
            res{3, 4}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 13
            
            res{3, 4}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 14
            
            %%%
            
            res{3, 5}  = zeros(3, 1);
            res{3, 5}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 13
            
            res{3, 5}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 14
            
            res{3, 5}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 15
            
            %%%
            res{3, 6}  = zeros(3, 1);
            res{3, 6}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 14
            
            res{3, 6}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 16
            
            %%%
            
            res{3, 7}  = zeros(3, 1);
            res{3, 7}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 15
            
            res{3, 7}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 15
            
            res{3, 7}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 17
            
            %%%
            
            res{4, 1}  = zeros(3, 1);
            res{4, 1}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 16
            
            res{4, 1}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 16
            
            %%%
            res{4, 2}  = zeros(3, 1);
            res{4, 2}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 17
            
            res{4, 2}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 18
            
            %%%
            
            res{4, 3}  = zeros(3, 1);
            res{4, 3}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 17
            
            res{4, 3}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 19
            
            %%%
            
            res{4, 4}  = zeros(3, 1);
            res{4, 4}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 18
            
            res{4, 4}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 18
            
            %%%
            
            res{4, 5}  = zeros(3, 1);
            res{4, 5}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 19
            
            res{4, 5}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 19
            
            res{4, 5}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 20
            
            %%%
            
            res{4, 6}  = zeros(3, 1);
            res{4, 6}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 20
            
            res{4, 6}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 20
            
            res{4, 6}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 21
            
            %%%
            
            res{4, 7}  = zeros(3, 1);
            res{4, 7}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 21
            
            res{4, 7}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 21
            
            
            %%%
            res{5, 1}  = zeros(3, 1);
            res{5, 1}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 22
            
            res{5, 1}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 22
            
            %%%
            res{5, 2}  = zeros(3, 1);
            res{5, 2}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 23
            
            res{5, 2}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 23
            
            %%%
            
            res{5, 3}  = zeros(3, 1);
            res{5, 3}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 22
            
            res{5, 3}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 24
            
            res{5, 3}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 24
            
            %%%
            
            res{5, 4}  = zeros(3, 1);
            res{5, 4}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 23
            
            res{5, 4}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 25
            
            res{5, 4}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 25
            
            %%%
            
            res{6, 1}  = zeros(3, 1);
            res{6, 1}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 24
            
            res{6, 1}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 26
            
            %%%
            
            res{6, 2}  = zeros(3, 1);
            res{6, 2}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 25
            
            res{6, 2}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 26
            
            res{6, 2}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 27
            
            %%%
            
            res{6, 3}  = zeros(3, 1);
            res{6, 3}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 26
            
            res{6, 3}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 28
            
            %%%
            
            res{6, 4}  = zeros(3, 1);
            res{6, 4}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 27
            
            res{6, 4}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 27
            
            res{6, 4}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 29
            
            %%%
            
            res{7, 1}  = zeros(3, 1);
            res{7, 1}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 28
            
            res{7, 1}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 28
            
            %%%
            
            res{7, 2}  = zeros(3, 1);
            res{7, 2}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 29
            
            res{7, 2}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 29
            
            res{7, 2}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            iRvD3 = iRvD3 + 1; %iRvD3 = 30
            
            %%%
            
            res{7, 3}  = zeros(3, 1);
            res{7, 3}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            iRvD1 = iRvD1 + 1; %iRvD1 = 30
            
            res{7, 3}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            iRvD2 = iRvD2 + 1; %iRvD2 = 30
           
            res{7, 3}(3) = RvD_val(obj.startIndexHesseRDD3 + iRvD3, ind);
            %iRvD3 = iRvD3 + 1; %iRvD3 = 31
            
            %%%
            
            res{7, 4}  = zeros(3, 1);
            res{7, 4}(1) = RvD_val(obj.startIndexHesseRDD1 + iRvD1, ind);
            %iRvD1 = iRvD1 + 1; %iRvD1 = 31
            
            res{7, 4}(2) = RvD_val(obj.startIndexHesseRDD2 + iRvD2, ind);
            %iRvD2 = iRvD2 + 1; %iRvD2 = 31
        end
        
        function res = getQ(obj, ind)
            res = obj.Q(:, ind);
        end
        
        function res = getQD(obj, ind)
            
            res = zeros(4, 7);
            
            QD_val = obj.QD;
           
            res(1, 2:7) = QD_val(1:6, ind)';
            res(2, 1) = QD_val(7, ind)';
            res(2, 3:7) = QD_val(8:12, ind)';
            res(3, 1:2) = QD_val(13:14, ind)';
            res(3, 4:7) = QD_val(15:18, ind)';
            res(4, 1:3) = QD_val(19:21, ind)';
            res(4, 5:7) = QD_val(22:24, ind)';
                       
        end
        
        function res = getQDD(obj, ind)
            res = cell(7, 7);
            
            iQD1 = 1;
            iQD2 = 1;
            iQD3 = 1;
            iQD4 = 1;
            
            QDD_val = obj.QDD;
            
            res{1, 5}  = zeros(4, 1);
            res{1, 5}(2) = QDD_val(obj.startIndexHesseQDD2 + iQD2, ind);
            iQD2 = iQD2 + 1; %iQD2 = 2
            
            res{1, 6}  = zeros(4, 1);
            res{1, 6}(3) = QDD_val(obj.startIndexHesseQDD3 + iQD3, ind);
            iQD3 = iQD3 + 1; %iQD3 = 2
            
            res{1, 7}  = zeros(4, 1);
            res{1, 7}(4) = QDD_val(obj.startIndexHesseQDD4 + iQD4, ind);
            iQD4 = iQD4 + 1; %iQD4 = 2
            
            res{2, 5}  = zeros(4, 1);
            res{2, 5}(1) = QDD_val(obj.startIndexHesseQDD1 + iQD1, ind);
            iQD1 = iQD1 + 1; %iQD1 = 2
            
            res{2, 6}  = zeros(4, 1);
            res{2, 6}(4) = QDD_val(obj.startIndexHesseQDD4 + iQD4, ind);
            iQD4 = iQD4 + 1; %iQD4 = 3
            
            res{2, 7}  = zeros(4, 1);
            res{2, 7}(3) = QDD_val(obj.startIndexHesseQDD3 + iQD3, ind);
            iQD3 = iQD3 + 1; %iQD3 = 3
            
            res{3, 5}  = zeros(4, 1);
            res{3, 5}(4) = QDD_val(obj.startIndexHesseQDD4 + iQD4, ind);
            iQD4 = iQD4 + 1; %iQD4 = 4   
            
            res{3, 6}  = zeros(4, 1);
            res{3, 6}(1) = QDD_val(obj.startIndexHesseQDD1 + iQD1, ind);
            iQD1 = iQD1 + 1; %iQD1 = 3
            
            res{3, 7}  = zeros(4, 1);
            res{3, 7}(2) = QDD_val(obj.startIndexHesseQDD2 + iQD2, ind);
            iQD2 = iQD2 + 1; %iQD2 = 3
            
            res{4, 5}  = zeros(4, 1);
            res{4, 5}(3) = QDD_val(obj.startIndexHesseQDD3 + iQD3, ind);
            iQD3 = iQD3 + 1; %iQD3 = 4
            
            res{4, 6}  = zeros(4, 1);
            res{4, 6}(2) = QDD_val(obj.startIndexHesseQDD2 + iQD2, ind);
            iQD2 = iQD2 + 1; %iQD2 = 4
            
            res{4, 7}  = zeros(4, 1);
            res{4, 7}(1) = QDD_val(obj.startIndexHesseQDD1 + iQD1, ind);
            iQD1 = iQD1 + 1; %iQD1 = 4
            
            res{5, 1}  = zeros(4, 1);
            res{5, 1}(2) = QDD_val(obj.startIndexHesseQDD2 + iQD2, ind);
            iQD2 = iQD2 + 1; %iQD2 = 5
            
            res{5, 2}  = zeros(4, 1);
            res{5, 2}(1) = QDD_val(obj.startIndexHesseQDD1 + iQD1, ind);
            iQD1 = iQD1 + 1; %iQD1 = 5
            
            res{5, 3}  = zeros(4, 1);
            res{5, 3}(4) = QDD_val(obj.startIndexHesseQDD4 + iQD4, ind);
            iQD4 = iQD4 + 1; %iQD4 = 5 
            
            res{5, 4}  = zeros(4, 1);
            res{5, 4}(3) = QDD_val(obj.startIndexHesseQDD3 + iQD3, ind);
            iQD3 = iQD3 + 1; %iQD3 = 5
            
            res{6, 1}  = zeros(4, 1);
            res{6, 1}(3) = QDD_val(obj.startIndexHesseQDD3 + iQD3, ind);
            iQD3 = iQD3 + 1; %iQD3 = 6
            
            res{6, 2}  = zeros(4, 1);
            res{6, 2}(4) = QDD_val(obj.startIndexHesseQDD4 + iQD4, ind);
            iQD4 = iQD4 + 1; %iQD4 = 6
            
            res{6, 3}  = zeros(4, 1);
            res{6, 3}(1) = QDD_val(obj.startIndexHesseQDD1 + iQD1, ind);
            iQD1 = iQD1 + 1; %iQD1 = 6
            
            res{6, 4}  = zeros(4, 1);
            res{6, 4}(2) = QDD_val(obj.startIndexHesseQDD2 + iQD2, ind);
            iQD2 = iQD2 + 1; %iQD2 = 6
            
            res{7, 1}  = zeros(4, 1);
            res{7, 1}(4) = QDD_val(obj.startIndexHesseQDD4 + iQD4, ind);
            %iQD4 = iQD4 + 1; %iQD4 = 7
            
            res{7, 2}  = zeros(4, 1);
            res{7, 2}(3) = QDD_val(obj.startIndexHesseQDD3 + iQD3, ind);
            %iQD3 = iQD3 + 1; %iQD3 = 7
            
            res{7, 3}  = zeros(4, 1);
            res{7, 3}(2) = QDD_val(obj.startIndexHesseQDD2 + iQD2, ind);
            %iQD2 = iQD2 + 1; %iQD2 = 7
            
            res{7, 4}  = zeros(4, 1);
            res{7, 4}(1) = QDD_val(obj.startIndexHesseQDD1 + iQD1, ind);
            %iQD1 = iQD1 + 1; %iQD1 = 7
        end
    end
    
    methods(Test)
        % Test der Rotationsmatrix
        function testR(obj)
            alpha = pi/2;
            alpha = alpha/2;
            r = [0 0 1]; %Drehachse
            r = r/norm(r);  
            r = sin(alpha) .* r;
            q = [cos(alpha) r(1) r(2) r(3)];
            R = classQuadrocopter.R(q);
            obj.verifyEqual(R * [1; 0; 0], [0; 1; 0], 'AbsTol', 1e-15);
        end
        
        
        function testState(obj)
            n_int_ = 50;
            n_contr_ = obj.n_contr;
            n_state_ = obj.n_var;
            
            obj.state = rand(n_state_, n_int_+1);
            obj.contr = rand(n_contr_, n_int_+1);
            
            
            o = zeros(length(obj.state), 1);
            
            for i = 1:length(obj.state)
                o(i) = norm(obj.state(4:7, i));
            end
            obj.verifyEqual(max(o), 1, 'AbsTol', obj.eps); % Quaternionen normalizisiert ? SO(3) ?
            
            obj.verifyEqual(obj.isEmptyM, true);    
            obj.verifyEqual(obj.isEmptyMD, true);
            obj.verifyEqual(obj.isEmptyMDD, true);
            obj.verifyEqual(obj.isEmptyTheta, true);    
            obj.verifyEqual(obj.isEmptyThetaD, true);
            obj.verifyEqual(obj.isEmptyThetaDD, true);
            obj.verifyEqual(obj.isEmptyT, true);
            obj.verifyEqual(obj.isEmptyTD, true);
            obj.verifyEqual(obj.isEmptyTDD, true);
        end
        
        function testM(obj)
            dummy = 10;
            M_ = obj.getM(dummy);
            obj.verifySize(M_, [6, 6]);
            obj.verifyEqual(obj.isEmptyM, false);
        end
        
        function testMD(obj)
            dummy = 10;
            MD_ = obj.getMD(dummy);
            obj.verifySize(MD_, [1, 1]);
        end
        
        function testMDD(obj)
            dummy = 10;
            MDD_ = obj.getMDD(dummy);
            obj.verifySize(MDD_, [1, 1]);
        end
        
        function testTheta(obj)
            n_int_ = 50;
            n_contr_ = obj.n_contr;
            n_state_ = obj.n_var;
            
            obj.state = rand(n_state_, n_int_+1);
            obj.contr = rand(n_contr_, n_int_+1);
            
            theta_ = obj.theta;
            obj.verifyEqual(obj.isEmptyTheta, false);
            obj.verifySize(theta_, [6, n_int_ + 1]);
             
            obj.emptyResults();
            obj.verifyEqual(obj.isEmptyTheta, true);
             
            theta_ = obj.getTheta(1);
            obj.verifyEqual(obj.isEmptyTheta, false);
             
            obj.verifySize(theta_, [6, 1]);
            
        end
        
        function testThetaD(obj)
            n_int_ = 50;
            n_contr_ = obj.n_contr;
            n_state_ = obj.n_var;
            
            obj.state = rand(n_state_, n_int_+1);
            obj.contr = rand(n_contr_, n_int_+1);
            
            
            thetaD_ = obj.thetaD;
            obj.verifySize(thetaD_, [28, n_int_ + 1]);
             
            obj.emptyResults();
            obj.verifyEqual(obj.isEmptyThetaD, true);
             
            thetaD_ = obj.getThetaD(1);
            obj.verifyEqual(obj.isEmptyThetaD, false);
            
         end
         
        function testThetaDD(obj)
            n_int_ = 50;
            n_contr_ = obj.n_contr;
            n_state_ = obj.n_var;
            
            obj.state = rand(n_state_, n_int_+1);
            obj.contr = rand(n_contr_, n_int_+1);
            
            
            thetaDD_ = obj.thetaDD;
            obj.verifySize(thetaDD_, [28, n_int_ + 1]);
             
            obj.emptyResults();
            obj.verifyEqual(obj.isEmptyThetaDD, true);
             
            thetaDD_ = obj.getThetaDD(1);
            obj.verifyEqual(obj.isEmptyThetaDD, false);
         end
         
        function testT(obj)
            n_int_ = 50;
            n_contr_ = obj.n_contr;
            n_state_ = obj.n_var;
            
            obj.state = rand(n_state_, n_int_+1);
            obj.contr = rand(n_contr_, n_int_+1);
            
            T_ = obj.T;
            obj.verifySize(T_, [6, n_int_ + 1]);
             
            obj.emptyResults();
            obj.verifyEqual(obj.isEmptyT, true);
             
            T_ = obj.getT(1);
            obj.verifyEqual(obj.isEmptyT, false);
         end
         
        function testTD(obj)
            n_int_ = 50;
            n_contr_ = obj.n_contr;
            n_state_ = obj.n_var;
            
            obj.state = rand(n_state_, n_int_+1);
            obj.contr = rand(n_contr_, n_int_+1);
            
            
            TD_ = obj.TD;
            obj.verifySize(TD_, [18, n_int_ + 1]);
             
            obj.emptyResults();
            obj.verifyEqual(obj.isEmptyTD, true);
             
            TD_ = obj.getTD(1);
            obj.verifyEqual(obj.isEmptyTD, false);
         end
         
        function testTDD(obj)
            n_int_ = 50;
            n_contr_ = obj.n_contr;
            n_state_ = obj.n_var;
            
            obj.state = rand(n_state_, n_int_+1);
            obj.contr = rand(n_contr_, n_int_+1);
            
            
            TDD_ = obj.TDD;
            obj.verifySize(TDD_, [28, n_int_ + 1]);
             
            obj.emptyResults();
            obj.verifyEqual(obj.isEmptyTDD, true);
             
            TDD_ = obj.getTDD(1);
            obj.verifyEqual(obj.isEmptyTDD, false);
            obj.verifySize(TDD_, [7, 7]);
         end
        
        function testRges(obj)
            n_int_ = 50;
            n_state_ = obj.n_var;
            
            obj.state = rand(n_state_, n_int_+1);
            
            Rges_ = obj.Rges;
            
            obj.verifySize(Rges_, [1, n_int_ + 1]);
            obj.emptyResults();
            obj.verifyEqual(obj.isEmptyRges, true);
            
            Rv_ = obj.getRv(1);
            obj.verifyEqual(obj.isEmptyRges, false);
            obj.verifySize(Rv_, [3, 1]);
             
        end
         
        function testRD(obj)
            
            n_int_ = 50;
            n_state_ = obj.n_var;
            
            obj.state = rand(n_state_, n_int_+1);
            
            RvD_ = obj.RvD;
            obj.verifySize(RvD_, [12, n_int_ + 1]);
             
            obj.emptyResults();
            obj.verifyEqual(obj.isEmptyRvD, true);
             
            RvD_ = obj.getRvD(1);
            obj.verifyEqual(obj.isEmptyRvD, false);
            
            obj.verifySize(RvD_, [3, 7]);
        end
         
        function testRDD(obj)
             
            n_int_ = 50;
            n_state_ = obj.n_var;
            
            obj.state = rand(n_state_, n_int_+1);
            
            RvDD_ = obj.RvDD;
            obj.verifySize(RvDD_, [90, n_int_ + 1]);
             
            obj.emptyResults();
            obj.verifyEqual(obj.isEmptyRvDD, true);
             
            RvDD_ = obj.getRvDD(1);
            obj.verifyEqual(obj.isEmptyRvDD, false);
             
        end
         
        function testQ(obj)
            n_int_ = 50;
            
            n_state_ = obj.n_var;
            obj.state = rand(n_state_, n_int_+1);
            
            Q_ = obj.Q;
            obj.verifySize(Q_, [4, n_int_ + 1]);
             
            obj.emptyResults();
            obj.verifyEqual(obj.isEmptyQ, true);
             
            Q_ = obj.getQ(1);
            obj.verifyEqual(obj.isEmptyQ, false);
        end
        
        function testQD(obj)
            
            n_int_ = 50;
            n_state_ = obj.n_var;
            
            obj.state = rand(n_state_, n_int_+1);
            
            QD_ = obj.QD;
            obj.verifySize(QD_, [24, n_int_ + 1]);
             
            obj.emptyResults();
            obj.verifyEqual(obj.isEmptyQD, true);
             
            QD_ = obj.getQD(1);
            obj.verifyEqual(obj.isEmptyQD, false);
            
            obj.verifySize(QD_, [4, 7]);
            
        end
                
        function testQDD(obj)
             
            n_int_ = 50;
            n_state_ = obj.n_var;
            
            obj.state = rand(n_state_, n_int_+1);
            
            QDD_ = obj.QDD;
            obj.verifySize(QDD_, [24, n_int_ + 1]);
             
            obj.emptyResults();
            obj.verifyEqual(obj.isEmptyQDD, true);
             
            QDD_ = obj.getQDD(1);
            obj.verifyEqual(obj.isEmptyQDD, false);
            obj.verifySize(QDD_, [7, 7]);
        end

        
    end
    
end

