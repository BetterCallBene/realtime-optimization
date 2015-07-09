clear 
close all

currentpath = cd('..');
addpath(pwd)
cd(currentpath)


global visualization;

steadyPoint = [2, 0, 5];

load RData.mat
ProcessTime1 = ProcessTime;

n_intervals = length(intervals);


steadyPoints = repmat(steadyPoint, n_intervals, 1); 
resR = v;
           
QR = zeros(n_intervals, 16);
array = reshape(resR, [length(resR)/n_intervals, n_intervals])'; 
% 
QR(:, 1:3)  = array(:, 11:13); %omega
[QR(:, 4), QR(:, 5), QR(:, 6)] = quat2angle((array(:, 4:7))); % Winkel
QR(:, 7:9) =  array(:, 8:10); %v
QR(:, 10:12) = array(:, 1:3); %position
%QR(:, 13:16) = array(:, 14:17); %u
QR(:, 13:15) =  Wind(1:3, :)';
QR(:, 16) = sqrt(sum((array(:, 1:3) - steadyPoints).^2, 2));










visualization.yout = QR;
visualization.tout = intervals;
% 
QuadAnim41;
