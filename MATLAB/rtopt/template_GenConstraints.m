classdef(Abstract) GenConstraints < handle
    %BASISGENQDYN wird von MAPLE/PYTHON generiert und enthaelt die
    %Berechnung der Ableitungen
    
    properties
        dode;   % handle for the ForwEuler element providing the discretization of the ode
        dyn;    % Dynamik
    end
    
    properties
        EqCon;
        EqConD; 
        EqConDD;
    end
    
    properties(Dependent)
        vec;
    end
    
    %setter;
    methods
        function set.vec(obj, vec_)
            obj.dyn.backdoor_vec = vec_;
        end
    end
    
    properties(GetAccess=private, SetAccess = protected)
        %isEmptyF;
        %isEmptyJ;
        %isEmptyH;
    end
    
    methods
        function cGC = GenConstraints(dode)
            % constructor based on two input values
            % a classForwEuler element and a classOCPparam element
            cGC.dode = dode;
            cGC.dyn = dode.dyn;
        end
        
        function res = get.EqCon(obj)
            %if obj.isEmptyF
                [q,v,omega,u,Iges,IM,m,kT,kQ,d,g, n_int] = getParams(obj);
                
                $4$ %CountConstraints
                res = zeros((n_int+ 1) * CountConstraints, 1);
                
                for timestep = 1:n_int+1
                    $0$
                    res(CountConstraints * (timestep-1)+1:CountConstraints * timestep) = t1;
                end
        end
        
        function res = get.EqConD(obj)
            [q,v,omega,u,Iges,IM,m,kT,kQ,d,g, n_int, n_state, n_contr] = getParams(obj);
            
            $4$ %CountConstraints
            $3$ %CountJacobi
            
            srow = zeros(CountJacobi * (n_int+1), 1);
            scol = zeros(CountJacobi * (n_int+1), 1);
            sval = zeros(CountJacobi * (n_int+1), 1);
            
            for timestep=1:(n_int+1)
                $1$
                srow(CountJacobi* (timestep-1)+1:CountJacobi * timestep)= uint16(t1(:, 1) + (timestep-1)*CountConstraints);
                scol(CountJacobi* (timestep-1)+1:CountJacobi * timestep) = uint16(t1(:, 2) + (timestep-1)*(n_state+n_contr));
                sval(CountJacobi* (timestep-1)+1:CountJacobi * timestep) = t1(:, 3);
            end
            
            res = sparse(srow, scol, sval, CountConstraints * (n_int + 1), (n_int+1)*(n_state+n_contr));
        end
        
        function res = get.EqConDD(obj)
            [q,v,omega,u,Iges,IM,m,kT,kQ,d,g, n_int, n_state, n_contr] = getParams(obj);
            
            $4$ %CountConstraints
            $2$ %HesseMatrix
            
            res = cell(uint16((n_int+1)*CountConstraints), 1);
            
            n_var = n_state+n_contr;
            ind = 1;
            for timestep=1:(n_int +1)
                for ConstraintStep = 1:CountConstraints
                    
                    tmp = t1(uint16(t1(:, 1)) == ConstraintStep, :);
                    
                    srow = (tmp(:, 2) + (timestep-1)*n_var);
                    scol = (tmp(:, 3) + (timestep-1)*n_var);
                    sval = tmp(:, 4);
                    res{ind} = sparse(srow, scol, sval,(n_int+1) * n_var,(n_int+1)*n_var);
                    ind = ind + 1;
                end
            end
        end
        
        function [q,v,omega,u,Iges,IM,m,kT,kQ,d,g, n_int, n_state, n_contr] = getParams(obj)
            
            q   = obj.dyn.state(4:7    , :);
            v   = obj.dyn.state(8:10   , :);
            omega   = obj.dyn.state(11:13  , :);
            u   = obj.dyn.contr;
            n_int = obj.dyn.environment.n_intervals; 
            n_state = obj.dyn.robot.n_state;
            n_contr = obj.dyn.robot.n_contr;
            
            Iges = obj.dyn.robot.I;
            IM = obj.dyn.robot.I_M;
            m = obj.dyn.robot.m;
            
            kT  = obj.dyn.robot.kT;
            kQ  = obj.dyn.robot.kQ;
            d   = obj.dyn.robot.d;
            g   = obj.dyn.environment.g;
        end
        
    end
    
end