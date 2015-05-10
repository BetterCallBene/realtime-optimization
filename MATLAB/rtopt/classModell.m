classdef(Abstract) classModell < handle & classTestEnv
    %classModell Dies ist die Interface/Abstrakte Klasse nachfolgende 
    
   properties(Abstract, ...
              GetAccess = public, ...
              SetAccess = public ...
             )
         state;      % the state matrix [r, q, v, w, dot v, dot w] in R^19 for all time instances
   end
   
   properties( Abstract, ... 
                SetAccess = protected, ... 
                GetAccess = public ...
    )
        M;          % matrix of mass matrix entries for all time instances
        MD;         % matrix of Jacobian entries of mass matrix for all time instances
        MDD;        % matrix of Hessian entries of mass matrix for all time instances
        
        theta;      % matrix of theta entries for all time instances
        thetaD;     % matrix of Jacobian entries of theta for all time instances
        thetaDD;    % matrix of Hessian entries of theta for all time instances
   end
   
   
   methods(Access=protected)
       % other methods
        function emptyResults(obj)
            
            % if state or param are change all these values
            % have to be computed again.
            obj.M       = [];
            obj.MD      = [];
            obj.MDD     = [];
            
            obj.theta   = [];
            obj.thetaD  = [];
            obj.thetaDD = [];
            
        end
   end
   
   methods
       %constructor
        function cR = classModell(varargin)
            cR.emptyResults();
        end
   end
   
    
end

