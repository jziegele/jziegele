function plot_summary_metric_by_parameters(filter_params, summary_metric_name, x_axis_param, group_by_param, show_constants_in_title)
% plot_summary_metric_by_parameters - Plots a summary metric from SimIndex.
%
% Syntax: plot_summary_metric_by_parameters(filter_params, summary_metric_name, x_axis_param, group_by_param, show_constants_in_title)
%
% Inputs:
%   filter_params         - A struct where each field is a parameter name (e.g., 'm1')
%                           and the value is the desired value (or array of values).
%   summary_metric_name   - String, the name of the summary metric to plot (e.g., 'QPA_mean').
%   x_axis_param          - String, the name of the parameter to use for the X-axis.
%                           Must be a parameter that varies (has multiple values) in filter_params.
%   group_by_param        - Optional string, the name of another parameter to group by.
%                           Must also be a varying parameter in filter_params.
%   show_constants_in_title - Optional boolean (default: true). If true,
%                             the constant parameters are listed in the plot title.
%
% Example:
%   filter.m1 = 1.0;
%   filter.m2 = [0.8, 1.0, 1.2];
%   filter.m3 = [0.8, 1.0, 1.2];
%   plot_summary_metric_by_parameters(filter, 'QPA_mean', 'm2', 'm3');

%% 0. Argument Handling
if nargin < 5
    show_constants_in_title = true; % Default to showing the title info
end
if nargin < 4
    show_constants_in_title = true;
    group_by_param = ''; % No grouping by default
end

% Validate required inputs
if nargin < 3
    error('Not enough input arguments. filter_params, summary_metric_name, and x_axis_param are required.');
end

%% 1. Load and Validate Index File
sim_index_file = 'SimIndex.mat';
output_directory = 'simulation_results'; % Not directly used here, but good to keep context

if ~exist(sim_index_file, 'file')
    error('SimIndex.mat not found. Please run main.m to generate simulations.');
end
load(sim_index_file, 'SimIndex');

if isempty(SimIndex)
    error('SimIndex is empty. No simulations to plot.');
end

%% 2. Filter the Simulation Index
param_names = {'m1', 'm2', 'm3', 'm4', 'm5', 'm6'};
filter_fields = fieldnames(filter_params);

% Initialize a mask for all simulations
combined_mask = true(1, size(SimIndex.inputs, 2));

% Apply filters
for i = 1:length(filter_fields)
    field = filter_fields{i};
    value = filter_params.(field);
    param_index = find(strcmp(param_names, field));
    
    if isempty(param_index)
        warning('Parameter %s not found. Skipping.', field);
        continue;
    end
    
    % Update the mask
    combined_mask = combined_mask & ismember(SimIndex.inputs(param_index, :), value);
end

% Get the indices of the simulations that passed the filter
subset_indices = find(combined_mask);

if isempty(subset_indices)
    disp('No simulations found matching the specified filter criteria.');
    return;
end

%% 3. Prepare Data for Plotting
% Validate x_axis_param
x_axis_param_index = find(strcmp(param_names, x_axis_param));
if isscalar(filter_params.(x_axis_param))
    error('x_axis_param (%s) must be a varying parameter (i.e., have multiple values in filter_params).', x_axis_param);
end

% Extract all relevant data points from the subset
all_x_values = SimIndex.inputs(x_axis_param_index, subset_indices);
all_y_values = cellfun(@(s) s.(summary_metric_name), SimIndex.summaryMetrics(subset_indices));

unique_x_values = unique(all_x_values);

if isempty(unique_x_values)
    error('No unique X-axis values found for parameter %s. Check filter_params and SimIndex data.', x_axis_param);
end

if isempty(group_by_param)
    % Simple plot (no grouping)
    plot_data = zeros(1, length(unique_x_values));
    for i = 1:length(unique_x_values)
        idx = (all_x_values == unique_x_values(i));
        plot_data(i) = mean(all_y_values(idx));
    end
    plot_type = 'simple';
else
    % Grouped plot
    group_by_param_index = find(strcmp(param_names, group_by_param));
    if isscalar(filter_params.(group_by_param))
        error('group_by_param (%s) must be a varying parameter (i.e., have multiple values in filter_params).', group_by_param);
    end
    all_group_values = SimIndex.inputs(group_by_param_index, subset_indices);
    unique_group_values = unique(all_group_values);

    if isempty(unique_group_values)
        error('No unique grouping values found for parameter %s. Check filter_params and SimIndex data.', group_by_param);
    end

    plot_data = NaN(length(unique_x_values), length(unique_group_values));

    for i = 1:length(subset_indices)
        sim_index = subset_indices(i);
        current_params = SimIndex.inputs(:, sim_index);
        [~, x_idx] = ismember(current_params(x_axis_param_index), unique_x_values);
        [~, group_idx] = ismember(current_params(group_by_param_index), unique_group_values);

        if x_idx > 0 && group_idx > 0
            plot_data(x_idx, group_idx) = SimIndex.summaryMetrics{sim_index}.(summary_metric_name);
        else
            warning('Could not find valid indices for x_axis_param or group_by_param for simulation %d. Skipping.', sim_index);
        end
    end
    plot_type = 'grouped';
end

%% 4. Plotting
figure;
hold on;
grid on;

if strcmp(plot_type, 'simple')
    bar(unique_x_values, plot_data);
    xticks(unique_x_values);
    % plot(unique_x_values, plot_data, 'o-', 'LineWidth', 1.5); % Alternative line plot
else % grouped
    b = bar(unique_x_values, plot_data, 'grouped');
    xticks(unique_x_values);
    % Add text labels to bars (optional, similar to plot_sensitivity_by_output)
    for k = 1:length(b)
        xtips = b(k).XEndPoints;
        ytips = b(k).YEndPoints;
        labels = string(round(b(k).YData, 2));
        text(xtips, ytips, labels, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom','FontSize', 8);
    end
    legend_labels = arrayfun(@(val) sprintf('%s: %.2f', group_by_param, val), unique_group_values, 'UniformOutput', false);
    legend(legend_labels, 'Location', 'bestoutside');
    xticklabels(unique_x_values);
end

%% 5. Finalize Plot
% Create a title string that includes the constant parameters
title_str = [summary_metric_name, ' vs. ', x_axis_param];
if show_constants_in_title
    constants_str = '';
    for i = 1:length(filter_fields)
        field = filter_fields{i};
        value = filter_params.(field);
        % Only add to title if it's a single, constant value AND not the x or group param
        if isscalar(value) && ~strcmp(field, x_axis_param) && ~strcmp(field, group_by_param)
            constants_str = [constants_str, sprintf('%s=%.2f, ', field, value)];
        end
    end

    if ~isempty(constants_str)
        title_str = [title_str, ' (Constants: ' constants_str(1:end-2) ')'];
    end
end

hold off;
title(title_str);
xlabel(x_axis_param);
y_label = strrep(summary_metric_name, '_', ' '); % Replace underscore for better display
ylabel(y_label);

end