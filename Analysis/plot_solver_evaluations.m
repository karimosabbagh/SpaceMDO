function plot_solver_evaluations()
    % Define solver folders and types
    solver_folders = {'Solver_Types/sqp', 'Solver_Types/mads', 'Solver_Types/ip'};
    solver_types = {'SQP', 'MADS', 'Interior-Point'};

    % Paths for the Analysis and Experiments folders
    analysis_folder = pwd; % Current directory for Analysis
    experiments_folder = fullfile(analysis_folder, '..', 'Experiments'); % Experiments folder is assumed to be one level up

    % Number of solvers
    numSolvers = length(solver_folders);

    % Initialize results storage
    cost_values = zeros(numSolvers, 1);
    percent_coverage_values = zeros(numSolvers, 1);

    % Loop through each solver folder
    for i = 1:numSolvers
        folder = fullfile(experiments_folder, solver_folders{i});

        % Load output.mat file
        outputFilePath = fullfile(folder, 'output.mat');
        if ~isfile(outputFilePath)
            error('File output.mat not found in folder: %s', folder);
        end
        load(outputFilePath, 'output', 'PB'); % Load variables output and PB

        % Run the evaluation function
        [results, evaluations] = evaluate_subsystems(output.x, PB);

        % Extract required values
        cost_values(i) = evaluations.cost;
        percent_coverage_values(i) = abs(evaluations.percent_coverage); % Ensure conversion to percentage
    end

    % Scale Adjustment for Percent Coverage
    scale_factor = max(cost_values) / max(percent_coverage_values); % Adjust scale to match cost
    scaled_percent_coverage = percent_coverage_values * scale_factor;

    % Create the grouped bar chart
    figure;
    hold on;

    % Plot bars for cost and scaled percent coverage
    bar_width = 0.4;
    x_positions = 1:numSolvers;

    % Bar for cost
    bar(x_positions - bar_width/2, cost_values, bar_width, 'FaceColor', [0 0.4470 0.7410], 'DisplayName', 'Cost');

    % Bar for percent coverage (scaled)
    bar(x_positions + bar_width/2, scaled_percent_coverage, bar_width, 'FaceColor', [0.8500 0.3250 0.0980], 'DisplayName', 'Percent Coverage');

    % Add secondary y-axis
    yyaxis left;
    ylabel('Cost');
    yyaxis right;
    ylabel('Percent Coverage');

    % X-axis labels
    set(gca, 'XTick', x_positions, 'XTickLabel', solver_types);
    xlabel('Solver Types');

    % Title and legend
    title('Evaluations for Different Solver Types');
    legend({'Cost', 'Percent Coverage (%)'}, 'Location', 'southoutside', 'Orientation', 'horizontal');

    % Add grid for better readability
    grid on;
    hold off;
end

% function plot_solver_evaluations()
%     % Define solver folders and types
%     solver_folders = {'Solver_Types/sqp', 'Solver_Types/mads', 'Solver_Types/ip'};
%     solver_types = {'SQP', 'MADS', 'Interior-Point'};
% 
%     % Paths for the Analysis and Experiments folders
%     analysis_folder = pwd; % Current directory for Analysis
%     experiments_folder = fullfile(analysis_folder, '..', 'Experiments'); % Experiments folder is assumed to be one level up
% 
%     % Number of solvers
%     numSolvers = length(solver_folders);
% 
%     % Initialize results storage
%     cost_values = zeros(numSolvers, 1);
%     percent_coverage_values = zeros(numSolvers, 1);
% 
%     % Loop through each solver folder
%     for i = 1:numSolvers
%         folder = fullfile(experiments_folder, solver_folders{i});
% 
%         % Load output.mat file
%         outputFilePath = fullfile(folder, 'output.mat');
%         if ~isfile(outputFilePath)
%             error('File output.mat not found in folder: %s', folder);
%         end
%         load(outputFilePath, 'output', 'PB'); % Load variables output and PB
% 
%         % Run the evaluation function
%         [results, evaluations] = evaluate_subsystems(output.x, PB);
% 
%         % Extract required values
%         cost_values(i) = evaluations.cost;
%         percent_coverage_values(i) = abs(evaluations.percent_coverage) * 100; % Convert to positive percentage
%     end
% 
%     % Scale Adjustment for Percent Coverage
%     scale_factor = max(cost_values) / max(percent_coverage_values); % Adjust scale to match cost
%     scaled_percent_coverage = percent_coverage_values * scale_factor;
% 
%     % Create the grouped bar chart
%     figure;
%     hold on;
% 
%     % Plot bars for cost and scaled percent coverage
%     bar_width = 0.4;
%     x_positions = 1:numSolvers;
% 
%     % Bar for cost
%     bar(x_positions - bar_width/2, cost_values, bar_width, 'FaceColor', [0 0.4470 0.7410], 'DisplayName', 'Cost');
% 
%     % Bar for percent coverage (scaled)
%     bar(x_positions + bar_width/2, scaled_percent_coverage, bar_width, 'FaceColor', [0.8500 0.3250 0.0980], 'DisplayName', 'Percent Coverage');
% 
%     % Add secondary y-axis
%     yyaxis left;
%     ylabel('Cost');
%     yyaxis right;
%     ylabel('Percent Coverage (%)');
% 
%     % X-axis labels
%     set(gca, 'XTick', x_positions, 'XTickLabel', solver_types);
%     xlabel('Solver Types');
% 
%     % Title and legend
%     title('Evaluations for Different Solver Types');
%     legend({'Cost', 'Percent Coverage (%)'}, 'Location', 'best');
% 
%     % Add grid for better readability
%     grid on;
%     hold off;
% end
