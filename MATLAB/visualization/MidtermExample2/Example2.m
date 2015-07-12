clear

currentpath = cd('..');
addpath(pwd)
cd(currentpath)

load Circle.mat

global visualization
            
Q = zeros(length(t), 12);
array = reshape(v, [length(v)/length(t), length(t)])'; 

Q(:, 1:3)  = array(:, 11:13); %omega
[Q(:, 4), Q(:, 5), Q(:, 6)] = quat2angle((array(:, 4:7))); % Winkel
Q(:, 7:9) =  array(:, 8:10); %v
Q(:, 10:12) = array(:, 1:3); %position


visualization.yout = Q;
visualization.tout = t;

QuadAnim4;
