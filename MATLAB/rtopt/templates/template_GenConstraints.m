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
        InEqCon;
        InEqConD; 
        InEqConDD;
    end
    
    properties(Dependent)
        vec;
    end
    
    %setter;
    methods
        function set.vec(obj, vec_)
            obj.dyn.backdoor_vec = vec_;
            obj.dode.flag_h = false;
            obj.dode.flag_hDD = false;
        end
    end
    
    methods
        function cGC = GenConstraints(dode)
            % constructor based on two input values
            % a classForwEuler element and a classOCPparam element
            cGC.dode = dode;
            if ~isempty(cGC.dode)
                cGC.dyn = dode.dyn;
            end
        end
        
        function res = get.EqCon(obj)
            %if obj.isEmptyF
            [r, q,v,omega,u,Iges,IM,m,kT,kQ,d,g, n_int] = getParams(obj);
                
            $4$ %eqCountConstraints
            res = zeros((n_int+ 1) * eqCountConstraints, 1);
                
            for timestep = 1:n_int+1
                $0$
                res(eqCountConstraints * (timestep-1)+1:eqCountConstraints * timestep) = t1;
            end
        end
        
        function res = get.InEqCon(obj)
            %if obj.isEmptyF
            [r, q,v,omega,u,Iges,IM,m,kT,kQ,d,g, n_int, n_state, n_contr, n_var, umin, umax] = getParams(obj);
                
            $10$ %inEqCountConstraints
            res = zeros((n_int+ 1) * inEqCountConstraints, 1);
                
            for timestep = 1:n_int+1
                $6$ %inEqCountConstraintsFunction
                res(inEqCountConstraints * (timestep-1)+1:inEqCountConstraints * timestep) = t1;
            end
        end
        
        function res = get.EqConD(obj)
            [r, q,v,omega,u,Iges,IM,m,kT,kQ,d,g, n_int, n_state, n_contr] = getParams(obj);
            
            $4$ %eqCountConstraints
            $3$ %eqCountJacobi
            
            srow = zeros(eqCountJacobi * (n_int+1), 1);
            scol = zeros(eqCountJacobi * (n_int+1), 1);
            sval = zeros(eqCountJacobi * (n_int+1), 1);
            
            for timestep=1:(n_int+1)
                $1$
                srow(eqCountJacobi* (timestep-1)+1:eqCountJacobi * timestep)= uint16(t1(:, 1) + (timestep-1)*eqCountConstraints);
                scol(eqCountJacobi* (timestep-1)+1:eqCountJacobi * timestep) = uint16(t1(:, 2) + (timestep-1)*(n_state+n_contr));
                sval(eqCountJacobi* (timestep-1)+1:eqCountJacobi * timestep) = t1(:, 3);
            end
            
            res = sparse(srow, scol, sval, eqCountConstraints * (n_int + 1), (n_int+1)*(n_state+n_contr));
        end
        
        function res = get.InEqConD(obj)
            [r, q,v,omega,u,Iges,IM,m,kT,kQ,d,g, n_int, n_state, n_contr] = getParams(obj);
            
            $10$ %inEqCountConstraints
            $9$ %inEqCountJacobi
            
            srow = zeros(inEqCountJacobi * (n_int+1), 1);
            scol = zeros(inEqCountJacobi * (n_int+1), 1);
            sval = zeros(inEqCountJacobi * (n_int+1), 1);
            
            for timestep=1:(n_int+1)
                $7$ %inEqConstraintsJacobi
                srow(inEqCountJacobi* (timestep-1)+1:inEqCountJacobi * timestep)= uint16(t1(:, 1) + (timestep-1)*inEqCountConstraints);
                scol(inEqCountJacobi* (timestep-1)+1:inEqCountJacobi * timestep) = uint16(t1(:, 2) + (timestep-1)*(n_state+n_contr));
                sval(inEqCountJacobi* (timestep-1)+1:inEqCountJacobi * timestep) = t1(:, 3);
            end
            
            res = sparse(srow, scol, sval, inEqCountConstraints * (n_int + 1), (n_int+1)*(n_state+n_contr));
        end
        
        function res = get.EqConDD(obj)
            [r, q,v,omega,u,Iges,IM,m,kT,kQ,d,g, n_int, n_state, n_contr] = getParams(obj);
            
            $4$ %eqCountConstraints
            $2$ %eqHesseMatrix
            
            res = cell(uint16((n_int+1)*eqCountConstraints), 1);
            
            n_var = n_state+n_contr;
            ind = 1;
            for timestep=1:(n_int +1)
                for ConstraintStep = 1:eqCountConstraints
                    
                    tmp = t1(uint16(t1(:, 1)) == ConstraintStep, :);
                    
                    srow = (tmp(:, 2) + (timestep-1)*n_var);
                    scol = (tmp(:, 3) + (timestep-1)*n_var);
                    sval = tmp(:, 4);
                    res{ind} = sparse(srow, scol, sval,(n_int+1) * n_var,(n_int+1)*n_var);
                    ind = ind + 1;
                end
            end
        end
        
        function res = get.InEqConDD(obj)
            [r, q,v,omega,u,Iges,IM,m,kT,kQ,d,g, n_int, n_state, n_contr] = getParams(obj);
            
            $10$ %inEqCountConstraints
            $8$ %inEqHesseMatrix
            
            res = cell(uint16((n_int+1)*inEqCountConstraints), 1);
            
            n_var = n_state+n_contr;
            ind = 1;
            for timestep=1:(n_int +1)
                for ConstraintStep = 1:inEqCountConstraints
                    
                    tmp = t1(uint16(t1(:, 1)) == ConstraintStep, :);
                    
                    srow = (tmp(:, 2) + (timestep-1)*n_var);
                    scol = (tmp(:, 3) + (timestep-1)*n_var);
                    sval = tmp(:, 4);
                    res{ind} = sparse(srow, scol, sval,(n_int+1) * n_var,(n_int+1)*n_var);
                    ind = ind + 1;
                end
            end
        end
        
        function [r, q,v,omega,u,Iges,IM,m,kT,kQ,d,g, n_int, n_state, n_contr, n_var, umin, umax] = getParams(obj)
            
            r   = obj.dyn.state(1:3    , :);
            q   = obj.dyn.state(4:7    , :);
            v   = obj.dyn.state(8:10   , :);
            omega   = obj.dyn.state(11:13  , :);
            u   = obj.dyn.contr;
            n_int = obj.dyn.environment.n_intervals; 
            n_state = obj.dyn.robot.n_state;
            n_contr = obj.dyn.robot.n_contr;
            n_var = obj.dyn.robot.n_var;
            
            Iges = obj.dyn.robot.I;
            IM = obj.dyn.robot.I_M;
            m = obj.dyn.robot.m;
            
            kT  = obj.dyn.robot.kT;
            kQ  = obj.dyn.robot.kQ;
            d   = obj.dyn.robot.d;
            g   = obj.dyn.environment.g;
            
            umin = obj.dyn.robot.u_min;
            umax = obj.dyn.robot.u_max;
        end
        
    end
    
end