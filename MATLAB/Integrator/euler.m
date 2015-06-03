function [t, y] = euler(number_of_steps, stepsize, f, y0, t0)

% Wie bekomme ich f?

n = number_of_steps;
y = zeros(17,n+1);
t = zeros(1,n+1);

y(1) = y0; 
t(1) = t0;

for j = 1:n;
    y(j+1) = y(j) + stepsize .* feval(f,y(j));
    t(j+1) = t(j) + stepsize;
end
