function [ee_cell, joint_cell] = bag_parser(path)
    % bag_parser -> Load and parse ROS bag files for JACO (j2s7s300)
    %
    % PURPOSE:
    %
    % Extract pose data in the end-effector control space and joint control
    % space from ROS bag files of the Kinova JACO (j2s7s300) robot arm 
    % executing a single Activities of Daily Living (ADL) task.
    %
    % The user specifies the path to the bag files. The user should have 
    % already placed the bag files for this task into their own folder.
    % Only bag files should be in placed in this folder. Each bag file 
    % represents a trial (i.e. repetition of the task).
    %
    % INPUTS:
    %
    % path - the search path (string) MATLAB needs to find the bag files
    % 
    % OUTPUTS:
    % 
    % ee_cell - A 1xn cell array where each cell is a 1x2 cell array 
    %           containing the pose data in end-effector space as a 
    %           translation matrix and rotation matrix for a single bag 
    %           file (n - total number of bag files/trials).
    % 
    % joint_cell - A 1xn cell array where each cell is a matrix
    %              representing the velocity data in joint space for a 
    %              single bag file
    
    files = dir(path); % return list of files in path directory
    num_files = length(files);
    
    source = 'j2s7s300_link_base'; % coordinate frames for EE space
    target = 'j2s7s300_end_effector';
    
    ee_cell = {};
    joint_cell = {};
    
    for i = 3:num_files
        
        T = []; % matrix to store translation (pos) data (m)
        R = []; % matrix to store rotation data (ori) (quaternions)
        
        J = []; % matrix to store joint state data (rad/sec)
        
        folder = strcat(files(i).folder, '/');
        name = files(i).name;
         
%%%%%%%%%%%%%%%%%%%%%% Robotic Systems Toolbox %%%%%%%%%%%%%%%%%%%%%%%%%%%

%         bag = rosbag(strcat(folder, name); % load ROS bag file
%         
%         % End-Effector Space
%         tf_bag = select(bag, 'Topic', '/tf');
%         tf_msgs = readMessages(tf_bag);
%         
%         
%         
%         % Joint Space
%         js_bag = select(bag, 'Topic', '/j2s7s300_driver/out/joint_state');
%         js_msgs = readMessages(js_bag);
%         for j = 1:length(js_msgs)
%             J = [J; js_msgs{j}.Velocity(1:7)'];
%         end
%           
%         joint_cell{end+1} = J;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%% +ros package %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        bag = ros.Bag.load(strcat(folder, name)); % +ros package
        
        % End-Effector Space (borrowed from Deepak's task_data script)
        tree = ros.TFTree(bag); % transformation tree
        [tf_msgs, tf_meta] = bag.readAll({'/tf'});
        
        ts = zeros(length(tf_meta),1); % timestamps for /tf messages
        for j = 1:length(tf_meta)
            ts(j) = tf_meta{j}.time.time;
        end
        
        ts(ts > (tree.time_end - 1)) = tree.time_end - 1; % may hold key to issues with Robotics toolbox
        ts(ts < (tree.time_begin + 1)) = tree.time_begin + 1; % why does this work??
        
        ee_poses = tree.lookup(source, target, ts)'; % Ex1 struct
        
        for j = 1:length(ee_poses)
            T = [T; ee_poses(j).translation']; % Ex3 matrix
            R = [R; ee_poses(j).rotation']; % Ex4 matrix (quaternions)
        end
        
        ee_cell{end+1} = {T,R};
        
        % Joint Space
        js_msgs = bag.readAll({'/j2s7s300_driver/out/joint_state'});
        for j = 1:length(js_msgs)
            J = [J; js_msgs{j}.velocity(1:7)']; % Jx7 matrix
        end
        
        joint_cell{end+1} = J;
        
    end
    
    ee_cell = ee_cell';
    joint_cell = joint_cell';
end
        
        
        
        
        
        
        
    