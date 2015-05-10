classdef classTestEnv < matlab.unittest.TestCase
    %CLASSTESTENV Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        eps = 1e-2;
    end
    
    methods
        %numDiffObj 
        function numdObj = numDiffObj(obj, func)   
%             %ToDo: Implementierung der numerischen Differentation
%                         
%               n   = length(obj.vec);
%               vec = obj.vec;
%  
%               for i=1:n
%                  % plus epsilon shift
%                   vecn    = vec;
%                   vecn(i) = vecn(i) + obj.eps;
%                   obj.vec = vecn;
%                   hp      = func();
% %     
% %                 % minus epsilon shift
%                   vecn    = vec;
%                   vecn(i) = vecn(i) - obj.eps;
%                   obj.vec = vecn;
%                   hm      = func();
% %     
% %                 % central difference
%                   numdObj    = (hp - hm)/2/obj.eps;
%                   max(numdObj)
%                   disp(numdObj);
%                   
%                   
% %                 % comparison numerical to analytical derivative
% %                 
%              end
        end
            % display result
    end


% %second order
% H       = obj.get_hess();
% m       = length(H);
% g_num   = cell(m,1);
% 
% for i=1:m
%     g_num{i} = zeros(n,n);
% end
% 
% error_val = 0;
% 
% for i=1:n
%     % plus epsilon shift
%     vecn    = vec;
%     vecn(i) = vecn(i) +eps;
%     obj.vec = vecn;
%     hp      = obj.get_jac();
%     
%     % minus epsilon shift
%     vecn    = vec;
%     vecn(i) = vecn(i) -eps;
%     obj.vec = vecn;
%     hm      = obj.get_jac();
%     
%     % central difference
%     num = (hp - hm)/2/eps;
%     for j =1:m
%         g_num{j}(:,i) = num(:,j);
%     end
% end
% 
% 
% for j=1:m
%     diff = H{j}- g_num{j};
%     diff_val = max(max(abs(diff)));
%     if (diff_val > error_val)
%         %disp(diff);
%         error_val = diff_val;
%     end
% end
% % display result
% disp(error_val)
%    end
    
end

