classdef(Abstract) Dyn < handle  & TestEnv
    % DYN Zentrale Klasse der Dynamik
    
    properties(SetAccess = public, GetAccess= public)
        environment; % in this object all outside parameters are stored, like gravity, wind, time mesh,...
        robot; % handle for a classRobot element providing the mass matrix and the coriolis terms
        solver; % handle of a Solver
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
                for j = 1:size(ana_dotDD, 1)
                    num_dotDD = reshape(numDiff(j,:,:), [17 17] );
                    anadotDDM = reshape(ana_dotDD(j, :, :), [17, 17]);
                    obj.assertLessThan(max(abs(anadotDDM - num_dotDD)),obj.tol);
                end
            end
        end
    end
    
    methods
       
        function vec = getVecFromCells(o,s,q)
            % GETVECFROMCELLS Takes as input cells s and q, where s_{i} is
            % the state at timepoint i and q{i} is the control at timepoint
            % i. This is rearanged, such that the data fits to the format
            % used in o.vec = [ s_1; q_1; ... ]. The last control value can
            % be set to zero, because we stop the calculation at that point.
            length_s = length(s);
            vec = zeros(length_s * o.robot.n_var,1);
            for i = 1: length_s-1
                vec( (i-1) * o.robot.n_var +1 : i * o.robot.n_var ) = ...
                    [ s{i} ; q{i}] ;
            end
            i = length_s;
            vec( (i-1) * o.robot.n_var +1 : end - o.robot.n_contr ) = s{i};
        end
        
        function func_p = plusEpsShift(obj,i,t,vec_old,func,n,dyn)
            % PLUSEPSSHIFT This funktion overwrites the one in TestEnv, as
            % we need unnormed quaternions here.
            vec_p = vec_old;
            vec_p((t-1)* obj.robot.n_var + i) = vec_p((t-1)* obj.robot.n_var + i) + obj.eps;
            obj.backdoor_vec = vec_p;
            func_p = func();
            obj.vec = vec_old;
        end
        
        function func_n = minusEpsShift(obj,i,t,vec_old,func,n,dyn)
            % MINUSEPSSHIFT Same function as PLUSEPSSHFIT only with a
            % minus.
            vec_n = vec_old;
            vec_n((t-1)* obj.robot.n_var + i) = vec_n((t-1)* obj.robot.n_var + i) - obj.eps;
            obj.backdoor_vec = vec_n;
            func_n = func();
            obj.vec = vec_old;
        end
        
    end
    
    methods
        function setupTest(obj,n_intervals)
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