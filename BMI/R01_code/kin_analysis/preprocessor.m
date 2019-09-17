function [ee_data, joint_data] = preprocessor(ee_cell, joint_cell, path, name)
    % preprocessor -> Preprocess bag files' data from bag_parser and save 
    %                 CSV files for dimensionality reduction
    %
    % PURPOSE:
    %
    % Convert rotation data in end-effector space from quaternions to Euler
    % angles. Obtain velocity data from end-effector space data by taking 
    % the differences between consecutive data points and dividing by the 
    % time step. Then apply a low-pass filter to reduce noise in the 
    % dataset. Joint and EE data are the outputs and they are saved as CSV 
    % files.
    % 
    % The user specifies the directory where the CSV files should be saved.
    % The user also determines the name of the CSV files, but the string,
    % "joint" or "ee" will be appended to the file name so the user knows
    % which control space is represented in the dataset. 
    %
    % INPUTS:
    % 
    % ee_cell - A 1xn cell array where each cell is a 1x2 cell array 
    %           containing the pose data in end-effector space as a 
    %           translation matrix and rotation matrix for a single bag 
    %           file (n - total number of bag files/trials).
    % 
    % joint_cell - A 1xn cell array where each cell is a matrix
    %              representing the velocity data in joint space for a 
    %              single bag file
    %
    % path - The location of the directory where the CSV files will be
    %        saved.
    %
    % name - The name of the CSV files. This name should reflect which
    %        robot and which task are represented by the data in the files.
    %
    % OUTPUTS:
    %
    % ee_data - A 1xn cell array where each cell is a 1x2 cell array
    %           containing the change in pose data in end-effector space 
    %           as a translation matrix and rotation matrix for a single 
    %           bag file
    %
    % joint_data - A 1xn cell array where each cell is a matrix
    %              representing the velocity data in joint space for a 
    %              single bag file
    
    ee_data = {};
    joint_data = {};
    
    % Calculate time (s) between /tf messages
    bag = rosbag('/home/tem/Documents/MATLAB/BMI/KD/6-28-19_Tem/Task1/T1_T1.bag');
    tf_bag = select(bag, 'Topic', '/tf');
    dt = (tf_bag.EndTime - tf_bag.StartTime)/(tf_bag.NumMessages-1);
    
    % Moving Average Filter Box Size
    box = 10;
    b = (1/box)*ones(box,1);
    
    max_ee = [];
    max_js = [];
         
    for i = 1:length(ee_cell) % loop through each bag file (trial)
        
%%%%%%%%%%%%%%%%%%%%%%%%%% End-Effector Space %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Convert rotation data from quaternions to Euler Angles (rad)
        R = quat2eul(ee_cell{i}{2}, 'XYZ'); % Ex3 matrix
        
        % Store translation and rotation data together
        ee_pose = [ee_cell{i}{1} R]; % Ex6 matrix
        
        V = [];
        
        for j = 2:length(ee_pose)
            
            % Calculate change in pose
            v = ee_pose(j,:) - ee_pose(j-1,:);
            
            % If all values less than 1e-4, then robot is static (ignore)
            if max(abs(v)) >= 1e-4
                V = [V; v];
            end
            
        end
        
        % Divide by time to get average velocity (m/s, rad/s)
        V = V/dt;
        
        % Find the maxes across the dataset
        max_ee = [max_ee; max(abs(V),[],1)];
        max_ee = max(max_ee,[],1);
        
        % Apply moving average filter
        V = mavefilter(V,10,4);
                
        % Append to output
        ee_data{end+1} = V;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Joint Space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        J = [];
        joint_speed = joint_cell{i};
        
        for j = 1:length(joint_speed)
            
            w = joint_speed(j,:);
            
            % If all values less than 1e-4, then robot is static (ignore)
            if max(abs(w)) >= 1e-4
                J = [J; w];
            end
        
        % Find the maxes across the dataset
        max_js = [max_js; max(abs(J),[],1)];
        max_js = max(max_js,[],1);
            
        end
        
        % Append to output
        joint_data{end+1} = J;
    
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normalize the Datasets %%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:length(ee_data)
        
        ee_data{i} = ee_data{i}./max_ee;
        joint_data{i} = joint_data{i}./max_js;
    
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save CSV files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    ee_data = ee_data';
    ee_file = strcat(path,'/',name,'_ee.csv');
    csvwrite(ee_file, cell2mat(ee_data));
    
    joint_data = joint_data';
    joint_file = strcat(path,'/',name,'_js.csv');
    csvwrite(joint_file, cell2mat(joint_data));
    
end
            