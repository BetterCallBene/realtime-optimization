function [y] = diag_beispiel(x0)

y = @(x) diag([exp(x),exp(2.*x),exp(3.*x)])*x0';


