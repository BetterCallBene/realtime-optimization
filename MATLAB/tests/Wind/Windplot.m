load('WindRand.mat');

figure

for timepoint = 1:size(Wind, 2)

    xyz = 0.3 * Wind(1:3, timepoint);

    quiver3(0, 0, 0, xyz(1), xyz(2), xyz(3), 'LineWidth', 3)
    title('Change in position caused by wind', 'FontSize', 14)
    view(-30, 45)
    
    xlim([0 0.3])
    ylim([0 0.3])
    zlim([0 0.3])    
   
    hold off
    pause(0.3)
end
