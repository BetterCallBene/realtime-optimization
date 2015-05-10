function visualize( t, y )
%VISUALIZE   produces a plot of the quadrocopter over time
% t     is a row vector containing every time step
% y     the ith row contains angular velocity, attitude, position and
%       velocity. Each of these consists of 3 values.
%
% An error is thrown, when the window is closed during a running
% visualization. To avoid that, press the STOP button first.
%
% At the moment only + configuration can be displayed.
global visualization

visualization.tout=t;
visualization.yout=y;

QuadAnim4;
end

