clear all; clc; close all
%%
global start_time bagfile;
%%
tasks = {'Task1', 'Task2', 'Task3', 'Task4', 'Task5'};
bag_folder = '/home/tem/Documents/MATLAB/BMI/KD/6-28-19_Tem/Task4/'; %path to appropriate folder containing the .bag files. 
fnames = dir(bag_folder);
numfids = length(fnames);

source_frame = 'j2s7s300_link_base';
target_frame = 'j2s7s300_end_effector';
jaco_task4_translation = []; %matrix for storing translation data
jaco_task4_rotation = []; %matrix for storing quaternion data
for i=3:numfids
    filename = fnames(i).name;  
    disp(filename);
    
    bagfile = ros.Bag.load(strcat(bag_folder, filename));
    bagfile.info()
    
    %% get start_time;
    get_start_time();
    
    %% get tf tree;
    jaco_poses = get_tf(source_frame, target_frame);
    jaco_translation = zeros(length(jaco_poses), 3); 
    jaco_rotation = zeros(length(jaco_poses), 4);
    for j=1:length(jaco_poses)
        jaco_translation(j, :) = jaco_poses(j).translation';
        jaco_rotation(j, :) = jaco_poses(j).rotation';
    end
    
    %append translation data from each trial
    jaco_task4_translation = [jaco_task4_translation; jaco_translation];
    %append orientation data from each trial
    jaco_task4_rotation = [jaco_task4_rotation; jaco_rotation];
end

%%

save('/home/tem/Documents/MATLAB/BMI/KD/6-28-19_Tem/MAT_Data/jaco_task4_translation.mat', 'jaco_task4_translation');
save('/home/tem/Documents/MATLAB/BMI/KD/6-28-19_Tem/MAT_Data/jaco_task4_rotation.mat', 'jaco_task4_rotation');

%% 
clear all; clc; close all;

%% Translation PCA

load('/home/tem/Documents/MATLAB/BMI/KD/6-28-19_Tem/MAT_Data/jaco_task4_translation.mat');
offset = 0.2;
xmin = min(jaco_task4_translation(:, 1)); xmax = max(jaco_task4_translation(:, 1)); xrange = [xmin-offset, xmax+offset];
ymin = min(jaco_task4_translation(:, 2)); ymax = max(jaco_task4_translation(:, 2)); yrange = [ymin-offset, ymax+offset];
zmin = min(jaco_task4_translation(:, 3)); zmax = max(jaco_task4_translation(:, 3)); zrange = [zmin-offset, zmax+offset];
offset = 0.3;
xrange = [-0.8, 0.8]; yrange = [-0.8, 0.8]; zrange = [-0.8, 0.8];


figure(1); hold on; grid on;

% Scatterplot
s1 = scatter3(jaco_task4_translation(:, 1), jaco_task4_translation(:, 2), jaco_task4_translation(:, 3), 'k')
xlabel('X'); ylabel('Y'); zlabel('Z');
axis([xrange, yrange, zrange]);
axis square;
view([142,31]);

% PCA
[coeff, score, latent, tsquared, explained, mu] = pca(jaco_task4_translation);
jaco_tasks_translation_centered = score*coeff;
% biplot(coeff);

% Centered Scatterplot
s2 = scatter3(jaco_tasks_translation_centered(:, 1), jaco_tasks_translation_centered(:, 2), jaco_tasks_translation_centered(:, 3), 'b')
xlabel('X'); ylabel('Y'); zlabel('Z');
axis([xrange, yrange, zrange]);
axis square;
view([142,31]);

% Axes
start_point = mu;
end_point = coeff + mu';
l1 = line([start_point(1) end_point(1, 1)],[start_point(2) end_point(2, 1)],[start_point(3) end_point(3, 1)],'Color', 'r', 'LineWidth', 1.5);
l2 = line([start_point(1) end_point(1, 2)],[start_point(2) end_point(2, 2)],[start_point(3) end_point(3, 2)],'Color', 'g', 'LineWidth', 1.5);
l3 = line([start_point(1) end_point(1, 3)],[start_point(2) end_point(2, 3)],[start_point(3) end_point(3, 3)],'Color', 'b', 'LineWidth', 1.5); 
 
l4 = line(xrange, [0,0], [0,0], 'Color', 'r', 'LineWidth', 1.5); %draw x and y axes.
l5 = line([0,0], yrange, [0,0], 'Color', 'g','LineWidth', 1.5);
l6 = line([0,0], [0,0], zrange, 'Color', 'b','LineWidth', 1.5);

title('Task 4 Translation PCA');
legend([s1 s2], {'Raw', 'Centered'});
hold off

%% Rotation PCA

load('/home/tem/Documents/MATLAB/BMI/KD/6-28-19_Tem/MAT_Data/jaco_task4_rotation.mat');
offset = 0.2;
xmin = min(jaco_task4_rotation(:, 1)); xmax = max(jaco_task4_rotation(:, 1)); xrange = [xmin-offset, xmax+offset];
ymin = min(jaco_task4_rotation(:, 2)); ymax = max(jaco_task4_rotation(:, 2)); yrange = [ymin-offset, ymax+offset];
zmin = min(jaco_task4_rotation(:, 3)); zmax = max(jaco_task4_rotation(:, 3)); zrange = [zmin-offset, zmax+offset];
offset = 0.3;
xrange = [-0.8, 0.8]; yrange = [-0.8, 0.8]; zrange = [-0.8, 0.8];

figure(2); hold on; grid on;

% Scatterplot
s1 = scatter3(jaco_task4_rotation(:, 1), jaco_task4_rotation(:, 2), jaco_task4_rotation(:, 3), 'k')
xlabel('X'); ylabel('Y'); zlabel('Z');
axis([xrange, yrange, zrange]);
axis square;
view([142,31]);

% PCA
[coeff, score, latent, tsquared, explained, mu] = pca(jaco_task4_rotation);
jaco_tasks_rotation_centered = score*coeff;
% biplot(coeff);

% Centered Scatterplot
s2 = scatter3(jaco_tasks_rotation_centered(:, 1), jaco_tasks_rotation_centered(:, 2), jaco_tasks_rotation_centered(:, 3), 'b')
xlabel('X'); ylabel('Y'); zlabel('Z');
axis([xrange, yrange, zrange]);
axis square;
view([142,31]);

% Axes
start_point = mu;
end_point = coeff + mu';
l1 = line([start_point(1) end_point(1, 1)],[start_point(2) end_point(2, 1)],[start_point(3) end_point(3, 1)],'Color', 'r', 'LineWidth', 1.5);
l2 = line([start_point(1) end_point(1, 2)],[start_point(2) end_point(2, 2)],[start_point(3) end_point(3, 2)],'Color', 'g', 'LineWidth', 1.5);
l3 = line([start_point(1) end_point(1, 3)],[start_point(2) end_point(2, 3)],[start_point(3) end_point(3, 3)],'Color', 'b', 'LineWidth', 1.5); 
 
l4 = line(xrange, [0,0], [0,0], 'Color', 'r', 'LineWidth', 1.5); %draw x and y axes.
l5 = line([0,0], yrange, [0,0], 'Color', 'g','LineWidth', 1.5);
l6 = line([0,0], [0,0], zrange, 'Color', 'b','LineWidth', 1.5);

title('Task 4 Rotation PCA');
legend([s1 s2], {'Raw', 'Centered'});
hold off

%% UTIL FUNCTIONS FOR BAG STUFF
function  get_start_time()
    global start_time;
    topic = '/tf';
    [~, meta] = get_bag_data(topic);
    start_time = meta{1}.time.time;
end

function [msgs, meta] = get_bag_data(topic)
    global bagfile;
    [msgs, meta] = bagfile.readAll({topic});
end


function tf_pose = get_tf(source_frame, target_frame)
    global bagfile;
    tree = ros.TFTree(bagfile);
    tree.allFrames();
    topic = '/tf';
    [~, meta] = get_bag_data(topic);
    ts = zeros(length(meta), 1);
    for i=1:length(meta)
        ts(i) = meta{i}.time.time;
    end
    ts(ts > (tree.time_end - 1)) = tree.time_end - 1;
    ts(ts < (tree.time_begin + 1)) = tree.time_begin + 1; 
    tf_pose = tree.lookup(source_frame, target_frame, ts);
       
end