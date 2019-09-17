% Load and parse files
route1 = '/home/tem/Documents/MATLAB/BMI/KD/6-28-19_Tem/Task1';
route2 = '/home/tem/Documents/MATLAB/BMI/KD/6-28-19_Tem/Task2';
route3 = '/home/tem/Documents/MATLAB/BMI/KD/6-28-19_Tem/Task3';
route4 = '/home/tem/Documents/MATLAB/BMI/KD/6-28-19_Tem/Task4';
route5 = '/home/tem/Documents/MATLAB/BMI/KD/6-28-19_Tem/Task5';

[ee_cell1, joint_cell1] = bag_parser(route1);
[ee_cell2, joint_cell2] = bag_parser(route2);
[ee_cell3, joint_cell3] = bag_parser(route3);
[ee_cell4, joint_cell4] = bag_parser(route4);
[ee_cell5, joint_cell5] = bag_parser(route5);

% Preprocess and normalize data
name1 = 'task1';
name2 = 'task2';
name3 = 'task3';
name4 = 'task4';
name5 = 'task5';
route = '/home/tem/Documents/MATLAB/BMI/KD/CSV';

[ee_data1, joint_data1] = preprocessor(ee_cell1, joint_cell1, route, name1);
[ee_data2, joint_data2] = preprocessor(ee_cell2, joint_cell2, route, name2);
[ee_data3, joint_data3] = preprocessor(ee_cell3, joint_cell3, route, name3);
[ee_data4, joint_data4] = preprocessor(ee_cell4, joint_cell4, route, name4);
[ee_data5, joint_data5] = preprocessor(ee_cell5, joint_cell5, route, name5);

% Boxplots of velocity distributions
ee_boxplot(ee_data1, ee_data2, ee_data3, ee_data4, ee_data5);
joint_boxplot(joint_data1, joint_data2, joint_data3, joint_data4, joint_data5);

