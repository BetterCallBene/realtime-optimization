function [ deltay ] = riccati2(LD, LDD )
% RICCATI2

[n,m] = size(LDD);

if n~=m
    error('LDD not quadratic');
end

cRic = RiccatiManager(m);
for i = n:-1:1
    
    if ( LDD(i,i) == 0)
        error('LDD(i,i) == 0');
    end
    
    %TODO: find better right limit than n
    cRic.doStep(i,LDD(i,:), LD(i),n);
    
end

   deltay = cRic.resolveRecursion();

end

