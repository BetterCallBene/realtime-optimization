classdef(Abstract) BasisGenQDyn < Dyn
    %BASISGENQDYN wird von MAPLE/PYTHON generiert und enthaelt die
    %Berechnung der Ableitungen
    
    properties
        F;
        J;
        H;
        
    end
    
    properties(GetAccess=private, SetAccess = protected)
        isEmptyF;
        isEmptyJ;
        isEmptyH;
    end
    
    methods
        function cBGQD = BasisGenQDyn()
            cBGQD@Dyn();
            cBGQD.emptyResults();
        end
        
        function emptyResults(obj)
            
            obj.F   = [];
            obj.J  = [];
            obj.H = [];
            
            obj.isEmptyF = true;
            obj.isEmptyJ = true;
            obj.isEmptyH = true;
        end
        
        function res = get.F(obj)
            if obj.isEmptyF
                [q,v,omega,u,Iges,IM,m,kT,kQ,d,g] = getParams(obj);
                
                $0$
                
                obj.F = t1;
                obj.isEmptyF = false;
            end
            res = obj.F;
        end
        
        function res = get.J(obj)
            if obj.isEmptyJ
                [q,v,omega,u,Iges,IM,m,kT,kQ,d,g] = getParams(obj);
                
                $1$
                
                obj.J = t1;
                obj.isEmptyJ = false;
            end
            res = obj.J;
        end
        
        function res = get.H(obj)
            if obj.isEmptyH
                [q,v,omega,u,Iges,IM,m,kT,kQ,d,g] = getParams(obj);
                
                onesV = ones(1, obj.environment.n_timepoints); %Korregiere die Matrix Dimension
                
                $2$
                
                obj.H = t1;
                obj.isEmptyH = false;
            end
            res = obj.H;
        end
        
        function res = FTilde(obj, y, u)
            [q,v,omega,u,Iges,IM,m,kT,kQ,d,g] = getParams(obj, y, u);
            $3$
            res = t1';
        end
        
        function res = JTilde(obj, y, u)
            [q,v,omega,u,Iges,IM,m,kT,kQ,d,g] = getParams(obj, y, u);
            $4$
            res = t1;
        end
        
        function res = getJTilde(obj, y, u)
            J_ = obj.JTilde(y, u);
            res = sparse(J_(:, 1), J_(:, 2), J_(:, 3), 13, 17);
        end
        
        function res = HTilde(obj, y, u)
            [q,v,omega,u,Iges,IM,m,kT,kQ,d,g] = getParams(obj, y, u);
            ones = 1;
            $5$
            res = t1;
        end
        
        function res = getHTilde(obj, y, u)
            
            res = cell(1, 13);
            H_ = obj.HTilde(y, u);
            for i = 1:size(H_, 1)
                if isempty(res{H_(i, 1)})
                    res{H_(i, 1)} = sparse(17, 17);
                end
                res{H_(i, 1)}(H_(i, 2), H_(i, 3)) = H_(i, 4);
            end
            
        end
        
        function [q,v,omega,u,Iges,IM,m,kT,kQ,d,g] = getParams(obj, varargin)
            
            if nargin == 3            
                y = varargin{1};
                u = varargin{2};
            else
                y = obj.state;
                u   = obj.contr;
            end
            q   = y(4:7    , :);
            v   = y(8:10   , :);
            omega   = y(11:13  , :);
            
            
            
            Iges = obj.robot.I;
            IM = obj.robot.I_M;
            m = obj.robot.m;
            
            kT  = obj.robot.kT;
            kQ  = obj.robot.kQ;
            d   = obj.robot.d;
            g   = obj.environment.g;
        end
        
    end
    
end