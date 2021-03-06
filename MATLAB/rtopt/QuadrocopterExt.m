classdef QuadrocopterExt < BasisQDyn
    % QuadrocopterExt Erweitert die Quadrocopterdynamik durch aeu�ere Einfluesse, Wind etc.
    properties
        hForceExt;
        hMomentExt;
        
        ForceExt;
        MomentExt;
        
        WindSav;
        
    end
    
    methods
        function QExt = QuadrocopterExt(varargin)
            QExt@BasisQDyn(varargin);
            if nargin >= 2
                QExt.WindSav = zeros(13, QExt.environment.n_timepoints);
            end
        end
        function F = wind(obj, t, st, ctr)
            %wind Erweitert das Modell durch ein externes Moment und Kraft
            solver_ = obj.solver;
            
            solver_.vec = [st; ctr];
            [old_intervals1] = solver_.preToDo();
            
            v = zeros(3, 1);%st(8:10);
            
            obj.ForceExt = obj.hForceExt(v);
            obj.MomentExt = obj.hMomentExt();
            
            
            % Matlab Bug.. Ahh
            n1 = logical(size(obj.ForceExt, 1) ~=3);
            n2 = logical(size(obj.ForceExt, 1) ~=3);
            if(n1 || n2)
                error('Force and Moment vector should by 3 x 1');
            end
            
            [F] =solver_.ode(1);
            if ~isempty(obj.WindSav)
                obj.WindSav(:, t) = F - st;
            end
            solver_.postToDo(old_intervals1);
        end
        
        function res = dot(obj,ind) 
            %dot Ueberschreibt die dot Funktion der BasisQDyn - Klasse zur
            %Erweiterung externe Momente und Kraefte
            res = obj.dot@BasisQDyn(ind);
            
            MomentExt_ = obj.MomentExt;
            ForceExt_ = obj.ForceExt;
            
            
            res(8:10) = res(8:10) + ForceExt_;
            res(11:13) = res(11:13) + MomentExt_;
        end
        
        function setupTest(obj,n_intervals)
            % Quadrocopter soll 5 Meter hoch fliegen
            xbc = [         ... Variablenname L�nge   Name
                ... Anfangsbedingung
                0, 0, 0,    ...     r           3      Ortsvektor
                1, 0, 0, 0, ...     q           4      Quaternion (Einheitsquaternion)
                0, 0, 0,    ...     v           3      Translatorische Geschwindigkeit
                0, 0, 0;    ...     w           3      Winkelgeschwindigkeit
                ... Endbedingung
                0, 0, 5,    ...
                1, 0, 0, 0, ...
                0, 0, 0,    ...
                0, 0, 0     ...
                ];
            
            env = Environment();
            env.xbc = xbc;
            env.setUniformMesh(uint8(n_intervals));
            
            model = Quadrocopter();
            opts_ = odeset('RelTol',1e-4,'AbsTol',1e-5);
            integrator = ode15sM(opts_);
            
            %load('TestData', 'data');
            
            %cBQD = QuadrocopterExt(model, env, integrator);
            obj.robot = model;
            obj.environment = env;
            obj.solver = integrator;
            obj.solver.dyn = obj;
            steadyPoint = obj.steadyPoint;
            obj.vec = repmat(steadyPoint, (n_intervals +1), 1); % Setup a initial estimation
            
            obj.hMomentExt = @obj.testMoment;
            obj.hForceExt = @obj.testForce;
            obj.WindSav = zeros(13, env.n_timepoints);
        end
        
        function res =  testForce(obj, v)
            res = rand(3, 1);
        end
        
        function res =  testMoment(obj)
            res = rand(3, 1);
        end
        
    end
    methods(Test)
        
        function testWind(testCase)
            % testWind Testet die Windfunktion
            n_intervals = 20;
            timepoint = 5;
            testCase.setupTest(n_intervals);
            
            st = rand(13, 1);
            ctr = rand(4, 1);
            res = testCase.wind(timepoint, st, ctr);
            steadyPoint1 = testCase.steadyPoint;
            abs(testCase.steadyPoint(1:13) - res)
            
            
        end
            
    end
    
    
end