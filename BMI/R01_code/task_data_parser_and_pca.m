clear all; clc; close all
%%
global start_time bagfile;
%%
tasks = {'Task1', 'Task2', 'Task3', 'Task4', 'Task5'};
bag_folder = '/home/tem/Documents/MATLAB/BMI/KD/6-28-19_Tem/Task1/'; %path to appropriate folder containing the .bag files. 
fnames = dir(bag_folder);
numfids = length(fnames);

source_frame = 'j2s7s300_link_base';
target_frame = 'j2s7s300_end_effector';
jaco_tasks_translation = []; %matrix for storing translation data
jaco_tasks_rotation = []; %matrix for storing quaternion data
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
    jaco_tasks_translation = [jaco_tasks_translation; jaco_translation];
    %append orientation data from each trial
    jaco_tasks_rotation = [jaco_tasks_rotation; jaco_rotation];
end

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