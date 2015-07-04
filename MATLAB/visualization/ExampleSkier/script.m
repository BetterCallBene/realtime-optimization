close all
clear all
load('skier2.mat');

costF = 0;
pos = zeros(3,n_timepoints);
for i = 1:n_timepoints
    costF(i) = res{i,5};
    tmp = res{i,1};
    pos(:,i) = tmp(1:3);
end

% Plot cost function
plot(costF)
title('Cost function value');

% Plot CamPos and Pos
figure
plot3(cam_pos_t(1,:), cam_pos_t(2,:), cam_pos_t(3,:),'r')
hold on
plot3(pos(1,:), pos(2,:), pos(3,:),'b')
axis([-10 10  0 1400 0 800])
view(-116,16);
title(' CamPosition (red) and computed position (blue)');


%Plot
norm_t = zeros(1,1200);
for i=1:1200
norm_t(i) = norm(cam_pos_t(:,i) - pos(:,i));
end

figure
plot(norm_t);
title('Distance between camPosition and computed position');


%For Screencast
% figure;
% pause(10);
% for i = 1:1400
% plot3(pos(1,1:i), pos(2,1:i), pos(3,1:i));
% axis([-10 10  0 1400 0 800])
% view(-116,16);
% pause(0.1);
% end