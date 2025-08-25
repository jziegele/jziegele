
function plot_pressure_volume_by_parameters(filter_params, pressure_to_plot, volume_to_plot, show_constants_in_title)
% plot_pressure_volume_by_parameters - Filters simulations and plots a pressure and volume waveform on the same graph.
%
% Syntax: plot_pressure_volume_by_parameters(filter_params, pressure_to_plot, volume_to_plot, show_constants_in_title)
%
% Inputs:
%   filter_params - A struct where each field is a parameter name (e.g., 'm1')
%                   and the value is the desired value (or array of values).
%   pressure_to_plot - String name of the pressure to plot (e.g., 'P11').
%   volume_to_plot - String name of the volume to plot (e.g., 'V11').
%   show_constants_in_title - Optional boolean (default: true). If true,
%                             the constant parameters are listed in the plot title.

%% 0. Argument Handling
if nargin < 4
    show_constants_in_title = true; % Default to showing the title info
end

%% 1. Load and Validate Index File
sim_index_file = 'SimIndex.mat';
output_directory = 'simulation_results';

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

% Identify varying parameters for the legend
legend_vars = {};
for i = 1:length(filter_fields)
    field = filter_fields{i};
    value = filter_params.(field);
    if ~isscalar(value)
        legend_vars{end+1} = field;
    end
end

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

%% 3. Map Signal Names to Indices
pressure_map = containers.Map({'P11', 'P21', 'P12', 'P22', 'P13', 'P23'}, [3, 4, 5, 6, 7, 8]);
volume_map = containers.Map({'V11', 'V12', 'V13', 'V21', 'V22', 'V23'}, [10, 11, 12, 13, 14, 15]);

if ~isKey(pressure_map, pressure_to_plot)
    error('Invalid pressure_to_plot name. Use one of: P11, P21, P12, P22, P13, P23');
end
if ~isKey(volume_map, volume_to_plot)
    error('Invalid volume_to_plot name. Use one of: V11, V12, V13, V21, V22, V23');
end

pressure_idx = pressure_map(pressure_to_plot);
volume_idx = volume_map(volume_to_plot);

%% 4. Plot the Waveforms
figure;

for i = 1:length(subset_indices)
    sim_index = subset_indices(i);
    
    % --- Load Data ---
    filename = sprintf('simulation%05d.mat', sim_index);
    file_to_load = fullfile(output_directory, filename);
    if ~exist(file_to_load, 'file')
        warning('Data file %s not found. Skipping.', filename);
        continue;
    end
    data = load(file_to_load);
    
    % --- Extract Variables ---
    t = data.results.t;
    Pressure = data.results.X(:,pressure_idx);
    Volume = data.results.o(:,volume_idx);
    tper = data.results.Outputs.tper;
    
    % --- Define Time Window (e.g., 5th to 7th cardiac cycle) ---
    x_start = 5 * tper;
    x_end = 7 * tper;
    t_idx = t >= x_start & t < x_end;
    
    % --- Create Legend Entry ---
    legend_str = '';
    current_params = SimIndex.inputs(:, sim_index);
    for j = 1:length(legend_vars)
        field = legend_vars{j};
        param_index = find(strcmp(param_names, field));
        if ~isempty(param_index)
            legend_str = [legend_str, sprintf('%s: %.2f, ', field, current_params(param_index))];
        else
            warning('Variable %s not found in simulation inputs for legend.', field);
        end
    end
    if isempty(legend_str)
        legend_str = sprintf('Run %d', sim_index);
    else
        legend_str = legend_str(1:end-2);
    end

    % --- Plot ---
    yyaxis left
    plot(t(t_idx) - t(find(t_idx,1)), Pressure(t_idx), '-', 'LineWidth', 1.5, 'DisplayName', [pressure_to_plot ' ' legend_str]);
    hold on;
    ylabel('Pressure (mmHg)');

    yyaxis right
    plot(t(t_idx) - t(find(t_idx,1)), Volume(t_idx), '--', 'LineWidth', 1.5, 'DisplayName', [volume_to_plot ' ' legend_str]);
    ylabel('Volume (mL)');
end

%% 5. Finalize Plot
grid on;
hold off;
legend('-DynamicLegend', 'Location', 'best');
xlabel('Time (s)');

% Create a title string that includes the constant parameters if requested
title_str = sprintf('%s and %s Waveform Comparison', pressure_to_plot, volume_to_plot);
if show_constants_in_title
    constants_str = '';
    filter_fields = fieldnames(filter_params);
    for i = 1:length(filter_fields)
        field = filter_fields{i};
        value = filter_params.(field);
        % Only add to title if it's a single, constant value
        if isscalar(value)
            constants_str = [constants_str, sprintf('%s=%.2f, ', field, value)];
        end
    end

    if ~isempty(constants_str)
        title_str = [title_str, ' (Constants: ' constants_str(1:end-2) ')'];
    end
end

title(title_str);

end
