function M = mavefilter(data, n, z)
    % mavefilter -> Filters data by applying a moving average (MAV) and
    %               removing data that is considered erroneous (E).
    %
    % PURPOSE:
    %
    % Calculate the moving average and standard deviation of the last n 
    % data points. Determine which points (if any) are more than z standard
    % deviations away from the mean. Those points will be considered noise
    % and deleted. If the data has multiple features, then a point will
    % only be removed from the dataset if it can be considered noise with
    % respect to all features. Thus, the data is filtered.
    %
    % The moving average for the first n-1 data points will have a lookback
    % window equal to the position of the current data point (i.e. if n =
    % 10, but we are at the ith data point -where i < n-  then n = i until
    % i >= 10).
    %
    % INPUTS:
    %
    % data - Data to be filtered. It can be a vector or matrix, however
    %        columns are always considered features and rows are points.
    %
    % n - The number of data points to be averaged (also known as the box
    %     or lookback window)
    % 
    % z - The number of standard deviations required for a data point to be
    %     considered erroneous.
    %
    % OUTPUTS:
    % 
    % M - The data after it has been filtered. 
    
    M = data(1,:);
    
    for i = 2:length(data) % Ignore first data point (trivial)
        
        % Lookback Window
        if i < n
            w = 1;
        else
            w = i-n+1;
        end
        
        D = data(w:i,:);
        
        % Calculate mean and standard deviation
        means = mean(D);
        stds = std(D);
        
        % Determine if data is erroneous
        d = (data(i,:)-means)./stds;
        
        if min(abs(d)) <= z
            M = [M; means]; % append moving average
        end
    end
end
        
        
        