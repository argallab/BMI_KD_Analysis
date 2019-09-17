function [group, p_val] = anova_pair(data)
       % anova_pair -> Take in a data matrix and conduct a one-way analysis
       %               of variance (ANOVA) on all possible pairs of column 
       %               data
       %
       % PURPOSE:
       %
       % Given a data matrix, where the columns are features and the rows
       % are samples, conduct a one-way ANOVA test on all possible pairs of
       % column data to determine if there are significant difference in
       % mean between columns. Differences in mean are considered to be
       % significant when a one-way ANOVA test yields a p-value less than
       % 0.05.
       %
       % Those pairs that have a significant difference in mean will be
       % saved into a cell array as a 1x2 vector and the corresponding
       % p-value will be saved into the vector p_val.
       %
       % The outputs will then be used as inputs to the sigstar function,
       % so that these results can be displayed on graphs.
       %
       % INPUTS:
       %
       % data - A data matrix where the columns are the features and the
       %        rows are the samples (data points). The distributions of
       %        the features will be compared using One-Way ANOVA.
       %
       % OUTPUTS:
       %
       % group - A cell array where each cell is a 1x2 vector representing
       %         the features that had a significant difference in mean. 
       %         The position of the vector in the cell array is the same
       %         as the corresponding p-value in the vector p-val.
       %
       % p_val - A vector containing all p-values that resulted from
       %         One-Way ANOVA tests where differences in means were 
       %         significant
       
       n = size(data,2); % number of features
       
       pairs = nchoosek(1:n,2); % all possible pairs
       
           
       
       p_val = [];
       group = {};
       
       for i = 1:size(pairs,1)
           
           g = pairs(i,:);
           test = [data(:,g(1)) data(:,g(2))];
           p = anova1(test,g,'off'); % One-way ANOVA
           
           % If ANOVA test is significant, then add p-value and pair
           if p <= 0.05
               p_val = [p_val; p];
               group{end+1} = g;
           end
       end
end