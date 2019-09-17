%% Translation Autoencoders

load('/home/tem/Documents/MATLAB/BMI/KD/6-28-19_Tem/MAT_Data/jaco_task1_translation.mat');

offset = 0.3;
xmin = min(jaco_tasks_translation(:, 1)); xmax = max(jaco_tasks_translation(:, 1)); xrange = [xmin-offset, xmax+offset];
ymin = min(jaco_tasks_translation(:, 2)); ymax = max(jaco_tasks_translation(:, 2)); yrange = [ymin-offset, ymax+offset];
zmin = min(jaco_tasks_translation(:, 3)); zmax = max(jaco_tasks_translation(:, 3)); zrange = [zmin-offset, zmax+offset];


figure(1); hold on; grid on;

% Scatterplot
s1 = scatter3(jaco_tasks_translation(:, 1), jaco_tasks_translation(:, 2), jaco_tasks_translation(:, 3), 'k')
xlabel('X'); ylabel('Y'); zlabel('Z');
axis([xrange, yrange, zrange]);
axis square;
view([142,31]);

% Autoencoders
translation_autoenc = trainAutoencoder(jaco_tasks_translation);
jaco_tasks_translation_decoded = predict(translation_autoenc, jaco_tasks_translation);
% biplot(coeff);

% Reconstructed Scatterplot
s2 = scatter3(jaco_tasks_translation_decoded(:, 1), jaco_tasks_translation_decoded(:, 2), jaco_tasks_translation_decoded(:, 3), 'b')
xlabel('X'); ylabel('Y'); zlabel('Z');
axis([xrange, yrange, zrange]);
axis square;
view([142,31]);

% Axes
% start_point = mu;
% end_point = coeff + mu';
% l1 = line([start_point(1) end_point(1, 1)],[start_point(2) end_point(2, 1)],[start_point(3) end_point(3, 1)],'Color', 'r', 'LineWidth', 1.5);
% l2 = line([start_point(1) end_point(1, 2)],[start_point(2) end_point(2, 2)],[start_point(3) end_point(3, 2)],'Color', 'g', 'LineWidth', 1.5);
% l3 = line([start_point(1) end_point(1, 3)],[start_point(2) end_point(2, 3)],[start_point(3) end_point(3, 3)],'Color', 'b', 'LineWidth', 1.5); 
%  
% l4 = line(xrange, [0,0], [0,0], 'Color', 'r', 'LineWidth', 1.5); %draw x and y axes.
% l5 = line([0,0], yrange, [0,0], 'Color', 'g','LineWidth', 1.5);
% l6 = line([0,0], [0,0], zrange, 'Color', 'b','LineWidth', 1.5);

title('Task 1 Translation PCA');
legend([s1 s2], {'Raw', 'Decoded'});
hold off