function plot_cost_component_distribution()
    % Define cost folders and values
    cost_folders = {'Cost/structuremassdominant', 'Cost/propellantmassdominant', ...
                    'Cost/ispdominant', 'Cost/FOVdominant'};
    cost_values = {'Structure', 'Propellant', 'ISP', 'FOV'};

    % Paths for the Analysis and Experiments folders
    analysis_folder = pwd; % Current directory for Analysis
    experiments_folder = fullfile(analysis_folder, '..', 'Experiments'); % Experiments folder is assumed to be one level up

    % Number of folders
    numFolders = length(cost_folders);

    % Initialize results storage
    percentages = zeros(numFolders, length(cost_values)); % Preallocate for speed

    % Loop through each folder
    for i = 1:numFolders
        folder = fullfile(experiments_folder, cost_folders{i});

        % Load output.mat file
        outputFilePath = fullfile(folder, 'output.mat');
        if ~isfile(outputFilePath)
            error('File output.mat not found in folder: %s', folder);
        end
        load(outputFilePath, 'output', 'PB'); % Load variables output and PB

        % Run the evaluation function
        [results, evaluations] = evaluate_subsystems(output.x, PB);

        % Extract percentage contributions for each cost component
        percentages(i, 1) = evaluations.cost_split(1);
        percentages(i, 2) = evaluations.cost_split(2);
        percentages(i, 3) = evaluations.cost_split(3);
        percentages(i, 4) = evaluations.cost_split(4);
    end

    % Plot stacked bar chart
    figure;
    bar(percentages, 'stacked');
    
    % Set x-axis labels to the folder names without the path
    xticks(1:numFolders);
    xticklabels({'Structure Mass  Dominant', 'Propellant Mass Dominant', ...
                 'ISP Dominant', 'FOV Dominant'});
    
    % Add legend, labels, and title
    legend(cost_values, 'Location', 'bestoutside');
    xlabel('Dominant Cost Factor');
    ylabel('Percentage Contribution (%)');
    title('Cost Component Distribution Across Exponent Dominance');
    grid on;
end


% cost_folders = {'Cost/structuremassdominant', 'Cost/propellantmassdominant','Cost/ispdominant','Cost/FOVdominant'};
% cost_values = {' Structure ' , ' Propellant', 'ISP' , 'FOV'};
% 
%  % Paths for the Analysis and Experiments folders
% analysis_folder = pwd; % Current directory for Analysis
% experiments_folder = fullfile(analysis_folder, '..', 'Experiments'); % Experiments folder is assumed to be one level up
% 
% 
% % Number of folders
% numFolders = length(cost_folders);
% 
% % Initialize results storage
% y_values = cell(numFolders, length(desired_evaluations));
% 
% % Loop through each folder
% for i = 1:numFolders
%     folder = fullfile(experiments_folder,cost_folders{i});
% 
%     % Load output.mat file
%     outputFilePath = fullfile(folder, 'output.mat');
%     if ~isfile(outputFilePath)
%         error('File output.mat not found in folder: %s', folder);
%     end
%     load(outputFilePath, 'output', 'PB'); % Load variables output and PB
% 
%     % Run the evaluation function
%     [results, evaluations] = evaluate_subsystems(output.x, PB);
% 
%     [Structure_Perc, Propellant_Perc, ISP_perc, FOV_perc] = evaluations.cost_split;
% 
% 
% 