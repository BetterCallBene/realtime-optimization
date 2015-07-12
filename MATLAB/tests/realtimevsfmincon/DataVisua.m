clear 
close all

currentpath = cd('..');
addpath(pwd)
cd(currentpath)

steadyPoint = [2, 0, 10];

load RData.mat
ProcessTime1 = ProcessTime;
 
resR = v;
           
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
ProcessTime2 = ProcessTime;

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

aR = sqrt(sum((QR(:, 10:12) - steadyPointM).^2, 2));
aF = sqrt(sum((QF(:, 10:12) - steadyPointM).^2, 2));

mFigure = figure();
%title('Jump from [2 0 5] to [2 0 10]');
set(mFigure,'name','Jump from [2 0 5] to [2 0 10]','numbertitle','off')
set(mFigure , 'Position', [0 0, 1920 , 1080 ]);

mTime = uicontrol('style', 'text');

strProcessTime = ['Process Time of fmincon: ', num2str(100), '%  Realtime: ', num2str((ProcessTime1/ProcessTime2)*100, 2), '%'];
set(mTime, 'String', strProcessTime);
set(mTime,'Units','characters')
set(mTime, 'Position', [-1, 48, 20, 5]);
set(mTime, 'FontSize', 14);

mSolver = uicontrol('style', 'text');

set(mSolver, 'String', 'Solver: ode15s with RelTol = 1e-2, AbsTol = 1e-3');
set(mSolver,'Units','characters')
set(mSolver, 'Position', [-1, 41, 20, 5]);
set(mSolver, 'FontSize', 14);

for timepoint = 1:length(intervals)
    subplot(2, 2, 1)
    
    plot(1:timepoint, QR(1:timepoint, 12), '-ro', 'markersize', 7);
    title('Current height of the quadrocopter', 'FontSize', 14);
    
    hold on
    plot(1:timepoint, QF(1:timepoint, 12), '-bo', 'markersize', 7);
    
    set(gca, 'FontSize', 12)
    xlim([0, length(intervals)]);
    ylim([min(QR(:, 12))-0.1, max(QR(:, 12))+0.1]);
    xlabel('t in s', 'FontSize', 14);
    ylabel('h in m', 'FontSize', 14);
    
    
    hold off
    
    subplot(2, 2, 2)
    
    plot(1:timepoint, aR(1:timepoint), '-ro', 'markersize', 7);
    title('Distance between the point [2 0 5] to [2 0 10]', 'FontSize', 14)
    hold on
    plot(1:timepoint, aF(1:timepoint), '-bo', 'markersize', 7);
    
    set(gca, 'FontSize', 12)
    xlim([0, length(intervals)]);
    ylim([min(aF(:))-0.1, max(aR())+0.1]);
    xlabel('t in s', 'FontSize', 14);
    ylabel('d in m', 'FontSize', 14);
    
    hold off
    
    subplot(2, 2, 3)
    
    plot(1:timepoint, QR(1:timepoint, 14), '-r');
    title('Motor speed of the realtime ansatz', 'FontSize', 14)
    
    set(gca, 'FontSize', 12)
    xlim([0, length(intervals)]);
    ylim([min(QR(:, 14))-0.1, max(QR(:, 14))+0.1]);
    xlabel('t in s', 'FontSize', 14);
    ylabel('u in 10^{-4} rpm', 'FontSize', 14);
    
    hold off
    
    subplot(2, 2, 4)
    
    plot(1:timepoint, QF(1:timepoint, 14), '-b');
    title('Motor speed of the fmincon ansatz', 'FontSize', 14)
    
    set(gca, 'FontSize', 12)
    xlim([0, length(intervals)]);
    ylim([min(QR(:, 14))-0.1, max(QR(:, 14))+0.1]);
    xlabel('t in s', 'FontSize', 14);
    ylabel('u in 10^{-4} rpm', 'FontSize', 14);
    
    hold off
    
    pause(0.2)
end
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

%visualization.yout = QR;
%visualization.tout = intervals;
% 
%QuadAnim4;
