function visPositionAttitude(t, position, attitude )
%VISPOSITIONATTITUDE position and attitude over time
% t             is a row vector containing every time step
% position      the ith row consists of [X,Y,Y] coordinates of the position
%               at time t(i)
% attitude      the ith row consists of [phi, theta, psi] angles at time t(i) 
% 
% An error is thrown, when the window is closed during a running
% visualization. To avoid that, press the STOP button first.
%
% At the moment only + configuration can be displayed.
global visualization

n = length(t);

visualization.yout = [ zeros(n,3) , attitude, zeros(n,3), position];
visualization.tout = t;

QuadAnim4;

end
