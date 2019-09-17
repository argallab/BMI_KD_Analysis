function ee_boxplot(varargin)
    % ee_boxplots -> Takes in velocity data in the end-effector (EE) space 
    %                for multiple tasks and produces boxplots for each task
    %                The boxplots are compared across axes and tasks
    %
    % PURPOSE:
    % 
    % Using cell arrays of velocity data in the EE space for a variable 
    % number of tasks, convert each cell array into a data matrix.Produce 
    % boxplots that compare the velocity distribution of an axis across 
    % different tasks. Produce boxplots that compare the velocity 
    % distributions of axes within a task. One-Way Analysis of Variance is 
    % performed to compare all distributions within a graph.
    %
    % INPUTS:
    % 
    % varargin - A variable number of cell arrays where each cell array
    %            contains the velocity data for a single task. Each cell in
    %            an array contains the velocity data for a single trial of 
    %            that task.
    % OUTPUTS:
    % 
    % Figure 1: N graphs where a single graph displays the velocity boxplot
    % of each axis for a single task (assuming N tasks)
    %
    % Figure 2: 6 graphs where a single graph displays the velocity boxplot
    % of the same axis for each task.
    
    axis_cell = {[];[];[];[];[];[]}; % varargin to be rearranged by axis
    
    n = length(varargin); 
    g = ceil(n/3);
    
    figure(1)
    sgtitle('EE Velocity vs Axis for each Task')
    
    for i = 1:n
        
        % Convert each cell array to a matrix
        varargin{i} = cell2mat(varargin{i});
        
        % Rearrange data by axis
        for j = 1:length(axis_cell)
            
            vec = varargin{i}(:,j);
            if i > 1
                diff = length(vec) - size(axis_cell{j},1);
                if diff > 0
                    axis_cell{j} = [axis_cell{j}; nan(diff, size(axis_cell{j},2))];
                else
                    vec = [vec; nan(abs(diff),1)];
                end
            end
            axis_cell{j} = [axis_cell{j} vec];
            
        end
        
        
        % Conduct one-way ANOVA tests on all pairs of axes
        [task_group, task_pval] = anova_pair(varargin{i});
        
        % Boxplots of axes for a task
        subplot(g,3,i)
        boxplot(varargin{i})
        xticklabels({'TX','TY','TZ','RX','RY','RZ'})
        xlabel('Axis')
        ylabel('Velocity')
        title(strcat('Task ',int2str(i)))
        hold on
        sigstar(task_group, task_pval);
        hold off
            
    end
    
    figure(2)
    sgtitle('EE Velocity vs Task for each Axis')
    titles = {'Trans X (TX)','Trans Y (TY)','Trans Z (TZ)','Rot X (RX)','Rot Y (RY)','Rot Z (RZ)'};
    
    xtl = {};
        
    % Cell array for xticklabels
    for j = 1:size(axis_cell{i},2)
        s = int2str(j);
        xtl{end+1} = s;
    end
    
    for i = 1:length(axis_cell)  
        
        % Conduct one-way ANOVA tests on all pairs of tasks
        [axis_group, axis_pval] = anova_pair(axis_cell{i});
        
        % Boxplots of tasks for an axis
        subplot(2,3,i)
        boxplot(axis_cell{i})
        xticklabels(xtl);
        xlabel('Task')
        ylabel('Velocity')
        title(titles{i})
        hold on
        sigstar(axis_group, axis_pval);
        hold off
        
    end
end
    
        