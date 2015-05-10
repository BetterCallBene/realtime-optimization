function numDiffObj(obj)

% numerical differentiation for objects of classes
% using functions get_func, get_jac, get_hess

vec = obj.vec;
eps = 1e-2;
n   = length(vec);

% first order
g   = obj.get_jac();
error_val = 0;

for i=1:n
    % plus epsilon shift
    vecn    = vec;
    vecn(i) = vecn(i) +eps;
    obj.vec = vecn;
    hp      = obj.get_func();
    
    % minus epsilon shift
    vecn    = vec;
    vecn(i) = vecn(i) -eps;
    obj.vec = vecn;
    hm      = obj.get_func();
    
    % central difference
    gnum    = (hp - hm)/2/eps;
    
    % comparison numerical to analytical derivative
    diff = g(i,:)' - gnum;
    diff_val = max(max(abs(diff)));
    if (diff_val > error_val)
        error_val = diff_val;
    end
end
% display result
disp(error_val)

%second order
H       = obj.get_hess();
m       = length(H);
g_num   = cell(m,1);

for i=1:m
    g_num{i} = zeros(n,n);
end

error_val = 0;

for i=1:n
    % plus epsilon shift
    vecn    = vec;
    vecn(i) = vecn(i) +eps;
    obj.vec = vecn;
    hp      = obj.get_jac();
    
    % minus epsilon shift
    vecn    = vec;
    vecn(i) = vecn(i) -eps;
    obj.vec = vecn;
    hm      = obj.get_jac();
    
    % central difference
    num = (hp - hm)/2/eps;
    for j =1:m
        g_num{j}(:,i) = num(:,j);
    end
end


for j=1:m
    diff = H{j}- g_num{j};
    diff_val = max(max(abs(diff)));
    if (diff_val > error_val)
        %disp(diff);
        error_val = diff_val;
    end
end
% display result
disp(error_val)