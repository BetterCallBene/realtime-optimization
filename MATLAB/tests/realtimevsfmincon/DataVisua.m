clear 
close all

currentpath = cd('..');
addpath(pwd)
cd(currentpath)

steadyPoint = [2, 0, 7];

load RData.mat
 
resR = v;
 
global visualization
%           
QR = zeros(length(intervals), 16);
array = reshape(resR, [length(resR)/length(intervals), length(intervals)])'; 
% 
QR(:, 1:3)  = array(:, 11:13); %omega
[QR(:, 4), QR(:, 5), QR(:, 6)] = quat2angle((array(:, 4:7))); % Winkel
QR(:, 7:9) =  array(:, 8:10); %v
QR(:, 10:12) = array(:, 1:3); %position
QR(:, 13:16) = array(:, 14:17); %u

% 
load FData.mat

resF = v;

QF = zeros(length(intervals), 16);
array = reshape(resF, [length(resF)/length(intervals), length(intervals)])'; 

QF(:, 1:3)  = array(:, 11:13); %omega
[QF(:, 4), QF(:, 5), QF(:, 6)] = quat2angle((array(:, 4:7))); % Winkel
QF(:, 7:9) =  array(:, 8:10); %v
QF(:, 10:12) = array(:, 1:3); %position
QF(:, 13:16) = array(:, 14:17); %u
% 
pair3 = ones(length(intervals), 1);
steadyPointM = pair3 * steadyPoint;
% 
%aF = sqrt(sum((QF(:, 10:12) - steadyPointM).^2, 2));

%plot(QF(:, 13:16))
%aR = sqrt(sum((QR(:, 10:12) - steadyPointM).^2, 2));
aF = sqrt(sum((QF(:, 10:12) - steadyPointM).^2, 2));

plot(intervals, aF)
%QF(end, 10:12) - steadyPoint
%QR(end, 10:12) - steadyPoint
% % 
%figure
%plot(intervals, aF, intervals, aR)
% figure
% for timepoint = 1:length(intervals)
%     plot(timepoint, QR(timepoint, 12), 'o') 
%     xlim([0, length(intervals)]);
%     ylim([min(QR(:, 12))-0.1, max(QR(:, 12))+0.1]);
%     hold on
%     pause(0.5)
% end
% 
% 
% figure 
% subplot(2, 1, 1)
% plot(intervals, QF(:, 13))
% plot(intervals, QF(:, 14))
% 
% subplot(2, 1, 2)
% plot(intervals, QF(:, 15))
% plot(intervals, QF(:, 16))
% 
% figure 
% subplot(2, 2, 1)
% plot(intervals, QR(:, 13))
% plot(intervals, QR(:, 14))
% % 
% subplot(2, 2, 2)
% plot(intervals, QR(:, 15))
% plot(intervals, QR(:, 16))
% 

% visualization.yout = QR;
% visualization.tout = intervals;
% 
% QuadAnim4;
