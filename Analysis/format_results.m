function resultTable = format_results(output, PB) %, filename)

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
    
    % Create the table
    resultTable = table(varNames, output.x(:), lowerBounds, upperBounds, isActive, ...
                        'VariableNames', {'VariableName', 'Result', 'LowerBound', 'UpperBound', 'isActive'});

    
    % Display the table
    disp(resultTable);

end
