%% main_plotting.m
% This script serves as the main file for generating all SVD projection plots.

%% 1. Initialization and Setup
clear;
clc;
addpath(genpath(pwd));

%% 2. Load Data
% Load the simulation index and the calculated SVD features.
try
    load('SimIndex.mat', 'SimIndex');
    load('Sim_SVD_Features.mat', 'all_sim_svd_features');
catch
    error('Required data files (SimIndex.mat or Sim_SVD_Features.mat) not found. Please run analyze_simulation_features.m first.');
end

%% 3. Generate Systematic Sensitivity Plots
% This section runs the systematic plotting analysis, which creates plots
% for all combinations of filtering by one parameter and coloring by another.
fprintf('Running systematic plotting analysis...\n');

% --- Configuration ---
connect_lines = true; % Set to true to plot with connecting lines

output_base_directory = 'SVD_plots';
if connect_lines
    output_base_directory = [output_base_directory, '_with_lines'];
end
if ~exist(output_base_directory, 'dir')
    mkdir(output_base_directory);
end

% --- Load or Initialize Plotting Index ---
plot_index_file = 'PlotIndex.mat';
if exist(plot_index_file, 'file')
    load(plot_index_file, 'PlotIndex');
else
    PlotIndex = struct('run_id', {}, 'filter', {}, 'color_filter', {}, 'output_directory', {});
end

% --- Define Parameters ---
param_names = {'m1', 'm2', 'm3', 'm4', 'm5', 'm6'};
base_values = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]; % Define a base case for constants

% --- Generate All Plotting Combinations ---
all_params_indices = 1:length(param_names);

for i = 1:length(all_params_indices)
    for j = 1:length(all_params_indices)
        if i == j
            continue;
        end

        x_axis_param_idx = i;
        color_param_idx = j;
    
    % --- Define Filter and Color Filter ---
    current_filter = struct();
    current_color_filter = struct();
    
    % Set all other parameters to their base value
    constant_param_indices = setdiff(all_params_indices, [x_axis_param_idx, color_param_idx]);
    for cp_idx = constant_param_indices
        current_filter.(param_names{cp_idx}) = base_values(cp_idx);
    end
    
    % The color filter will be the second parameter in the combination
    current_color_filter.(param_names{color_param_idx}) = unique(SimIndex.inputs(color_param_idx, :));
    
    % --- Check if this combination has been run ---
    % (For simplicity, we'll run all and overwrite for now. A more robust
    % implementation would check for existing identical runs.)
    
    run_id = length(PlotIndex) + 1;
    
    % --- Create Output Directory ---
    output_dir_name = sprintf('plot_run_%03d_x-%s_color-%s', run_id, param_names{x_axis_param_idx}, param_names{color_param_idx});
    output_directory = fullfile(output_base_directory, output_dir_name);
    if ~exist(output_directory, 'dir')
        mkdir(output_directory);
    end
    
    fprintf('--- Running Plotting Analysis ID: %d ---\n', run_id);
    fprintf('X-axis: %s, Color: %s\n', param_names{x_axis_param_idx}, param_names{color_param_idx});
    
    % --- Generate Plots ---
    generate_sensitivity_plots(current_filter, current_color_filter, output_directory, SimIndex, all_sim_svd_features, connect_lines);
    
    % --- Update and Save Plot Index ---
    PlotIndex(run_id).run_id = run_id;
    PlotIndex(run_id).filter = current_filter;
    PlotIndex(run_id).color_filter = current_color_filter;
    PlotIndex(run_id).output_directory = output_directory;
    
    save(plot_index_file, 'PlotIndex');
    
    fprintf('Finished Plotting Analysis ID: %d. Plots saved to: %s\n', run_id, output_directory);
    
    end
end

disp('All plotting analyses complete.');
fprintf('Systematic plotting analysis complete.\n\n');

%% 4. Generate Default SVD Sensitivity Plots
% This section runs the default SVD sensitivity plotting, which creates
% scatter plots of SVD projection coefficients against each parameter multiplier.
fprintf('Running default SVD sensitivity plotting...\n');

% --- Configuration ---
output_base_directory = 'plots_svd_sensitivity';
if ~exist(output_base_directory, 'dir')
    mkdir(output_base_directory);
end

% --- Validate Data Dimensions ---
if isempty(SimIndex.inputs) || isempty(all_sim_svd_features)
    error('SimIndex.inputs or all_sim_svd_features is empty. No data to plot.');
end

if size(SimIndex.inputs, 2) ~= size(all_sim_svd_features, 2)
    error('Mismatch in number of simulations between SimIndex.inputs and all_sim_svd_features.');
end

% --- Define Parameters ---
param_names = {'m1', 'm2', 'm3', 'm4', 'm5', 'm6'};
num_params = length(param_names);

% --- Define SVD Component Names ---
num_svd_components = size(all_sim_svd_features, 1);
svd_component_names = cell(num_svd_components, 1);
for k = 1:num_svd_components
    svd_component_names{k} = sprintf('SVD_Component_%d', k);
end

% --- Define Filters ---
filters = {};

% Filter 1: everything 
filter1 = struct();
filters{end+1} = filter1;

% Filter 2: Vary m1 and m2, others constant (color combinations)
filter2 = struct();
filter2.m3 = 1.0;
filter2.m4 = 1.0;
filter2.m5 = 1.0;
filter2.m6 = 1.0;
filters{end+1} = filter2;

% Filter 3: Highlight m1 = 1.0 with a color filter
filter3 = struct();
filter3.m3 = 1.0;
filter3.m4 = 1.0;
filter3.m5 = 1.0;
filter3.m6 = 1.0;
filter3.color_filter.m1 = [0.7,1.0, 1.3];
filters{end+1} = filter3;

% Filter 4: Vary m1, others constant
filter4 = struct();
filter4.m2 = 1.0;
filter4.m3 = 1.0;
filter4.m4 = 1.0;
filter4.m5 = 1.0;
filter4.m6 = 1.0;
filters{end+1} = filter4;


% --- Iterate through each filter ---
for f_idx = 1:length(filters)
    current_filter = filters{f_idx};
    
    % --- Create Filter-Specific Directory ---
    filter_dir_name = sprintf('filter_%d', f_idx);
    filter_dir_path = fullfile(output_base_directory, filter_dir_name);
    if ~exist(filter_dir_path, 'dir')
        mkdir(filter_dir_path);
    end
    
    % --- Write Filter to Text File ---
    filter_file_path = fullfile(filter_dir_path, 'filter.txt');
    fid = fopen(filter_file_path, 'w');
    if fid == -1
        error('Could not open filter.txt for writing.');
    end
    fprintf(fid, 'Filter %d:\n', f_idx);
    filter_fields = fieldnames(current_filter);
    for i = 1:length(filter_fields)
        field = filter_fields{i};
        if strcmp(field, 'color_filter')
            fprintf(fid, '  color_filter:\n');
            color_filter_fields = fieldnames(current_filter.color_filter);
            for j = 1:length(color_filter_fields)
                color_field = color_filter_fields{j};
                values = current_filter.color_filter.(color_field);
                fprintf(fid, '    %s: [%s]\n', color_field, num2str(values));
            end
        else
            values = current_filter.(field);
            fprintf(fid, '  %s: [%s]\n', field, num2str(values));
        end
    end
    fclose(fid);
    
    % --- Apply Filter ---
    if isfield(current_filter, 'color_filter')
        main_filter = rmfield(current_filter, 'color_filter');
    else
        main_filter = current_filter;
    end
    [filtered_inputs, filtered_svd_features] = apply_filter(SimIndex, all_sim_svd_features, main_filter, param_names);
    
    % --- Determine which parameters are varying ---
    varying_param_indices = get_varying_parameters(filtered_inputs, num_params);
    
    % --- Plotting Logic ---
    if isfield(current_filter, 'color_filter')
        for p_idx = varying_param_indices
            plot_with_color_filter(filtered_inputs, filtered_svd_features, p_idx, current_filter.color_filter, param_names, svd_component_names, filter_dir_path, connect_lines);
        end
    else
        switch length(varying_param_indices)
            case 1
                plot_single_parameter_sensitivity(filtered_inputs, filtered_svd_features, varying_param_indices, param_names, svd_component_names, filter_dir_path);
            case 2
                plot_two_parameter_color_combinations(filtered_inputs, filtered_svd_features, varying_param_indices, param_names, svd_component_names, filter_dir_path);
            otherwise
                plot_default_sensitivity(filtered_inputs, filtered_svd_features, varying_param_indices, param_names, svd_component_names, filter_dir_path);
        end
    end
end

disp('Plotting complete.');
fprintf('Default SVD sensitivity plotting complete.\n\n');

%% 5. Generate Violin Plots
fprintf('Generating violin plots...\n');
output_violin_directory = 'violin_plots';
if ~exist(output_violin_directory, 'dir')
    mkdir(output_violin_directory);
end

param_names = {'m1', 'm2', 'm3', 'm4', 'm5', 'm6'};
num_svd_components = size(all_sim_svd_features, 1);

for p_idx = 1:length(param_names)
    param_name = param_names{p_idx};
    
    % Create subdirectory for the parameter
    param_dir = fullfile(output_violin_directory, param_name);
    if ~exist(param_dir, 'dir')
        mkdir(param_dir);
    end
    
    % Use all simulation data for the plots
    all_inputs = SimIndex.inputs;
    all_svd_features = all_sim_svd_features;
    
    multipliers = unique(all_inputs(p_idx, :));
    
    for svd_idx = 1:num_svd_components
        hFig = figure('Name', ['Violin plot for SVD ' num2str(svd_idx) ' vs ' param_name]);
       
        % Collect data for violin plot
        violin_data = [];
        groups = [];
        for m_idx = 1:length(multipliers)
            multiplier = multipliers(m_idx);
            indices = all_inputs(p_idx, :) == multiplier;
            data = all_svd_features(svd_idx, indices);
            violin_data = [violin_data, data];
            groups = [groups, repmat(multiplier, 1, length(data))];
        end
        
        violinplot(groups, violin_data);
        
        %set(gca, 'XTickLabel', cellstr(num2str(multipliers(:))));
        
        title(['SVD ' num2str(svd_idx) ' vs ' param_name]);
        xlabel(['Multiplier for ' param_name]);
        ylabel('SVD Feature Value');
        
        saveas(hFig, fullfile(param_dir, ['SVD' num2str(svd_idx) '_violin.png']));
        close(hFig);
    end
end
fprintf('Violin plots generated.\n\n');


