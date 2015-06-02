classdef(Abstract) Dyn < handle  & TestEnv
    % Dyn Diese Klasse repräsentiert die Dynamik
    
    properties(SetAccess = public, GetAccess= public)
        environment; % in this object all outside parameters are stored, like gravity, wind, time mesh,...
        robot; % handle for a classRobot element providing the mass matrix and the coriolis terms
    end
    
    properties(Abstract)
        state; % the state matrix [r, q, v, w] in R^(13xn) for all time instances
        contr; % the control matrix [u1, u2, u3, u4] in R^(4xn) for all time instances
        backdoor_vec;
    end
    
    methods(Abstract)
        dot(obj,ind);
        dotD(obj,ind);
        dotDD(obj,ind);
    end
    
    methods(Test)
        
        function testdotD(obj)
            % TESTDOTD This method nummerially derives dot and compares it
            % with dotD
            n_intervals = 50;
            obj.setupTest(n_intervals);
            
            for timepoint = 1:(n_intervals+1)
                
                func = @() obj.dot(timepoint);
                numDiff = obj.numDiff_nD(timepoint,func);
                
                %Calculate analytic solution
                ana_dotD = obj.dotD(timepoint);
                
                %Assert that the result has the correct form
                obj.assertSize(ana_dotD, [13,17]);
                
                %Compare the results
                obj.assertLessThan(max(abs(ana_dotD - numDiff)), obj.tol);
            end
        end
        
        function testdotDD(obj)
            % TESTDOTDD This method nummericaly derives dotD and compares
            % it with dotDD
            n_intervals = 50;
            obj.setupTest(n_intervals);
            
            
            for timepoint = 1:(n_intervals+1)
                
                func = @() obj.dotD(timepoint);
                numDiff = obj.numDiff_nxnD(timepoint,func);
                
                ana_dotDD = obj.dotDD(timepoint);
                for j = 1:length(ana_dotDD)
                    num_dotDD = reshape(numDiff(j,:,:), [17 17] );
                    obj.assertLessThan(max(abs(ana_dotDD{j} - num_dotDD)),obj.tol);
                end
            end
        end
    end
    
    methods
        
        %Overwrite function in TestEnv, because we don't want normed
        %quaternios here.
        function func_p = plusEpsShift(obj,i,t,vec_old,func)
            vec_p = vec_old;
            vec_p((t-1)* obj.robot.n_var + i) = vec_p((t-1)* obj.robot.n_var + i) + obj.eps;
            obj.backdoor_vec = vec_p;
            func_p = func();
            obj.vec = vec_old;
        end
        
        %Overwrite function in TestEnv, because we don't want normed
        %quaternios here.
        function func_n = minusEpsShift(obj,i,t,vec_old,func)
            vec_n = vec_old;
            vec_n((t-1)* obj.robot.n_var + i) = vec_n((t-1)* obj.robot.n_var + i) - obj.eps;
            obj.backdoor_vec = vec_n;
            func_n = func();
            obj.vec = vec_old;
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
            obj.environment = env;
            obj.robot = model;
            %obj.state=rand(13,n_intervals+1);
            %obj.contr=rand(4,n_intervals+1);
            
            obj.vec = rand(17* (n_intervals+1), 1);
            
        end
    end
end