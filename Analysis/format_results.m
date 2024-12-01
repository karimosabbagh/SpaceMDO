function resultTable = format_results(output, PB, filename)
    % Function to create a result table from the output vector and PB struct
    
    % Get the variable names, lower bounds, and upper bounds from PB
    numVars = numel(PB.var);
    varNames = cell(numVars, 1);
    lowerBounds = zeros(numVars, 1);
    upperBounds = zeros(numVars, 1);
    
    for i = 1:numVars
        varNames{i} = PB.var{i}{1}; % Get variable name
        lowerBounds(i) = PB.var{i}{6}; % Get lower bound
        upperBounds(i) = PB.var{i}{7}; % Get upper bound
    end
    
    % Ensure the output vector has the same number of elements as variables
    if numel(output.x) ~= numVars
        error('Mismatch between the number of variables and the output vector length.');
    end

    % Determine if the result is active (at bounds)
    isActive = (output.x(:) == lowerBounds) | (output.x(:) == upperBounds);


    % Create a copy of the results to adjust for date and FOV conversion
    resultsAdjusted = output.x(:); % Start with a copy of the original results
    
    % for i = 1:numVars
    %     % Check if the variable name contains 'date' (for date conversion)
    %     if contains(varNames{i}, 'date', 'IgnoreCase', true)
    %         % Convert from datenum to datetime
    %         resultsAdjusted(i) = datetime(fix(resultsAdjusted(i)), 'ConvertFrom', 'datenum');
    %         lowerBounds(i) = datetime(lowerBounds(i), 'ConvertFrom', 'datenum'); % Convert lower bound to a date
    %         upperBounds(i) = datetime(upperBounds(i), 'ConvertFrom', 'datenum'); % Convert upper bound to a date
    %     end
    % 
    %     % Check if the variable name contains 'FOV' (for FOV conversion)
    %     if contains(varNames{i}, 'FOV', 'IgnoreCase', true)
    %         % Convert from radians to degrees
    %         resultsAdjusted(i) = rad2deg(resultsAdjusted(i));
    %         lowerBounds(i) = rad2deg(lowerBounds(i)); % Convert lower bound from radians to degrees
    %         upperBounds(i) = rad2deg(upperBounds(i)); % Convert upper bound from radians to degrees
    %     end
    % end

    % Create the table
    resultTable = table(varNames, resultsAdjusted(:), lowerBounds, upperBounds, isActive, ...
                        'VariableNames', {'VariableName', 'Result', 'LowerBound', 'UpperBound', 'isActive'});
    % 
    % % Create the table
    % resultTable = table(varNames, output.x(:), lowerBounds, upperBounds, isActive, ...
    %                     'VariableNames', {'VariableName', 'Result', 'LowerBound', 'UpperBound', 'isActive'});
    
    % Display the table
    disp(resultTable);

    save(filename, 'resultTable');
end
