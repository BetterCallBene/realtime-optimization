clear
load('ode15i', 'F', 'J')

ode15iF = F;
ode15iJ = J;

load('ForwEuler', 'F', 'J')

ForwEulerF = F;
ForwEulerJ = J;


ode15iF - ForwEulerF
spy(abs(ForwEulerJ - ode15iJ))