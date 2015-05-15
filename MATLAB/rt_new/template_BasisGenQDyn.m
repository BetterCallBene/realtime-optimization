classdef(Abstract) BasisGenQDyn < Dyn
    %BASISGENQDYN wird von MAPLE/PYTHON generiert und enthält die
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
                
                {0}
                
                obj.F = cg;
                obj.isEmptyF = false;
            end
            res = obj.F;
        end
        
        function res = get.J(obj)
            if obj.isEmptyJ
                [q,v,omega,u,Iges,IM,m,kT,kQ,d,g] = getParams(obj);
                
                {1}
                
                obj.J = t1;
                obj.isEmptyJ = false;
            end
            res = obj.J;
        end
        
        function res = get.H(obj)
            if obj.isEmptyH
                [q,v,omega,u,Iges,IM,m,kT,kQ,d,g] = getParams(obj);
                
                onesV = ones(1, obj.environment.n_intervals); %Korregiere die Matrix Dimension
                
                {2}
                
                obj.H = t1;
                obj.isEmptyH = false;
            end
            res = obj.H;
        end
        
        function [q,v,omega,u,Iges,IM,m,kT,kQ,d,g] = getParams(obj)
            q   = obj.state(4:7    , :);
            v   = obj.state(8:10   , :);
            omega   = obj.state(11:13  , :);
            u   = obj.contr;
            
            
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