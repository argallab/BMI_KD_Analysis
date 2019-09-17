function joint_boxplot(varargin)
    % joint_boxplot -> Takes in velocity data in the joint space for 
    %                  multiple tasks and produces boxplots for each task. 
    %                  The boxplots are compared across axes and tasks
    %
    % PURPOSE:
    %
    % Using cell arrays of velocity data in the joint space for a variable 
    % number of tasks, convert each cell array into a data matrix. Produce 
    % boxplots that compare the velocity distribution of a joint across 
    % different tasks. Produce boxplots that compare the velocity 
    % distributions of joints within a task. One-Way Analysis of Variance 
    % is performed to compare all distributions within a graph.
    %
    % INPUTS:
    %
    % varargin - A variable number of cell arrays where each cell array
    %            contains the velocity data for a single task. Each cell in
    %            an array contains the velocity data for a single trial of 
    %            that task.
    %
    % OUTPUTS:
    %
    % Figure 1: N graphs where a single graph displays the velocity boxplot
    % of each joint for a single task (assuming N tasks)
    %
    % Figure 2: 6 graphs where a single graph displays the velocity boxplot
    % of the same joint for each task.
    
    j_cell = {[];[];[];[];[];[];[]}; % varargin to be rearranged by joint
    
    n = length(varargin); 
    g = ceil(n/3);
    
    figure(1)
    sgtitle('Joint Velocity vs Axis for each Task')
    
    for i = 1:n
        
        % Convert each cell array into a matrix
        varargin{i} = cell2mat(varargin{i});
        
        % Rearrange data by axis
        for j = 1:length(j_cell)
            
            vec = varargin{i}(:,j);
            if i > 1
                diff = length(vec) - size(j_cell{j},1);
                if diff > 0
                    j_cell{j} = [j_cell{j}; nan(diff, size(j_cell{j},2))];
                else
                    vec = [vec; nan(abs(diff),1)];
                end
            end
            j_cell{j} = [j_cell{j} vec];
            
        end
        
        % Conduct one-way ANOVA tests on all pairs of axes
        [task_group, task_pval] = anova_pair(varargin{i});
        
        xtl = {};
        
        % Cell array for xticklabels
        for j = 1:size(varargin{i},2)
            s = int2str(j);
            xtl{end+1} = s;
        end
        
        % Boxplots of axes for a task
        subplot(g,3,i)
        boxplot(varargin{i})
        xticklabels(xtl)
        xlabel('Joint')
        ylabel('Velocity')
        title(strcat('Task ',int2str(i)))
        hold on
        sigstar(task_group, task_pval);
        hold off
            
    end
    
    figure(2)
    sgtitle('Joint Velocity vs Task for each Axis')
   
    
    for i = 1:length(j_cell)
        
        xtl = {};
        
        % Cell array for xticklabels
        for j = 1:size(j_cell{i},2)
            s = int2str(j);
            xtl{end+1} = s;
        end
        
        % Conduct one-way ANOVA tests on all pairs of tasks
        [j_group, j_pval] = anova_pair(j_cell{i});
        
        % Boxplots of tasks for an axis
        subplot(4,2,i)
        boxplot(j_cell{i})
        xticklabels(xtl);
        xlabel('Task')
        ylabel('Velocity')
        title(strcat('Joint ',int2str(i)))
        hold on
        sigstar(j_group, j_pval);
        hold off
        
    end
end
    