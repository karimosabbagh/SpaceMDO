function plot_evaluation_results(folder_names, x_values, x_label, desired_evaluations)
    
     % Paths for the Analysis and Experiments folders
    analysis_folder = pwd; % Current directory for Analysis
    experiments_folder = fullfile(analysis_folder, '..', 'Experiments'); % Experiments folder is assumed to be one level up
    
    % Paths for the Analysis and Experiments folders
    analysis_folder = pwd; % Current directory for Analysis
    experiments_folder = fullfile(analysis_folder, '..', 'Experiments'); % Experiments folder is assumed to be one level up

    % Check if string_names is a cell array of strings
    if ~iscell(desired_evaluations) || ~all(cellfun(@ischar, desired_evaluations))
        error('string_names must be a cell array of strings.');
    end
    
    % Check if x_values is a cell array of strings
    if iscell(x_values) && all(cellfun(@ischar, x_values))
        x_is_strings = true;
        x_categories = categorical(x_values); % Convert strings to categorical for plotting
    elseif isnumeric(x_values)
        x_is_strings = false;
    else
        error('x_values must be a cell array of strings or a numeric array.');
    end
    
    % Number of folders
    numFolders = length(folder_names);
    
    % Initialize results storage
    y_values = cell(numFolders, length(desired_evaluations));
    
    % Loop through each folder
    for i = 1:numFolders
        folder = fullfile(experiments_folder,folder_names{i});
        
        % Load output.mat file
        outputFilePath = fullfile(folder, 'output.mat');
        if ~isfile(outputFilePath)
            error('File output.mat not found in folder: %s', folder);
        end
        load(outputFilePath, 'output', 'PB'); % Load variables output and PB
        
        % Run the evaluation function
        [results, evaluations] = evaluate_subsystems(output.x, PB);
        
        % Extract required values
        for j = 1:length(desired_evaluations)
            try
                y_values{i, j} = eval(desired_evaluations{j});
            catch ME
                error('Field "%s" not found.', desired_evaluations{j})
            end
            % if isfield(evaluations, desired_evaluations{j})
            %     y_values{i, j} = evaluations.(desired_evaluations{j});
            %     %y_values{i, j} = desired_evaluations(j);
            % else
            %     error('Field "%s" not found in evaluations.', desired_evaluations{j});
            % end
        end
    end
    
    % Plot results
    figure;
    hold on;

    % Process labels to remove anything before the dot
    processed_labels = cellfun(@(s) regexp(s, '\.(.*)', 'tokens'), desired_evaluations, 'UniformOutput', false);
    processed_labels = cellfun(@(c) c{1}{1}, processed_labels, 'UniformOutput', false); % Extract the token
    
    % Replace underscores with spaces in labels
    processed_labels = cellfun(@(s) strrep(s, '_', ' '), processed_labels, 'UniformOutput', false);

    if length(desired_evaluations) == 1
        % Single y-value plot
        y_combined = abs(cell2mat(y_values(:, 1)'));
        if x_is_strings
            plot(x_categories, y_combined, '-o');
        else
            plot(x_values, y_combined, '-o');
        end
        
        % ylabel(strrep(desired_evaluations{1}, '_', ' '));
        ylabel(processed_labels{1});
    elseif length(desired_evaluations) == 2
        % Two y-value plot with two y-axes
        yyaxis left;
        y_left = abs(cell2mat(y_values(:, 1)'));
        if x_is_strings
            plot(x_categories, y_left, '-o');
        else
            plot(x_values, y_left, '-o');
        end
        % ylabel(strrep(desired_evaluations{1}, '_', ' '));
        ylabel(processed_labels{1});
        yyaxis right;
        y_right = abs(cell2mat(y_values(:, 2)'));
        if x_is_strings
            plot(x_categories, y_right, '-s');
        else
            plot(x_values, y_right, '-s');
        end
        % ylabel(strrep(desired_evaluations{2}, '_', ' '));
        ylabel(processed_labels{2});
    else
        error('Only one or two string_names are supported.');
    end
    
    % Common x-axis labeling
    xlabel(x_label);
    %legend(desired_evaluations, 'Location', 'best');
    % legend(cellfun(@(s) strrep(s, '_', ' '), desired_evaluations, 'UniformOutput', false), 'Location', 'best');
    legend(processed_labels, 'Location', 'best');
    title('Evaluation Results');
    grid on;
    hold off;
end
