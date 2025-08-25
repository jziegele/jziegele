
function plot_qpa_by_parameters(filter_params, show_constants_in_title)
% plot_qpa_by_parameters - Filters simulations and plots QPA waveforms.
%
% Syntax: plot_qpa_by_parameters(filter_params, show_constants_in_title)
%
% Inputs:
%   filter_params - A struct where each field is a parameter name (e.g., 'm1')
%                   and the value is the desired value (or array of values).
%   show_constants_in_title - Optional boolean (default: true). If true,
%                             the constant parameters are listed in the plot title.
%
% Example:
%   filter.m1 = 1.0;
%   filter.m2 = [1.0, 1.3];
%   plot_qpa_by_parameters(filter); % Legend automatically shows m2, title shows m1 as constant.

%% 0. Argument Handling
if nargin < 2
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

%% 3. Plot the QPA Waveforms
figure;
hold on;
grid on;

legend_entries = {};

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
    Q_PA = data.results.X(:,2); % QPA is the 2nd column of X
    tper = data.results.Outputs.tper;
    
    % --- Define Time Window (e.g., 5th to 7th cardiac cycle) ---
    x_start = 5 * tper;
    x_end = 7 * tper;
    t_idx = t >= x_start & t < x_end;
    
    % --- Plot ---
    plot(t(t_idx) - t(find(t_idx,1)), Q_PA(t_idx), 'LineWidth', 1.5);
    
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
        legend_entries{end+1} = sprintf('Run %d', sim_index); % Fallback
    else
        legend_entries{end+1} = legend_str(1:end-2); % Remove trailing comma
    end
end

%% 4. Finalize Plot
% Create a title string that includes the constant parameters if requested
title_str = 'QPA Waveform Comparison';
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

hold off;
title(title_str, 'fontsize',8);
xlabel('Time (s)');
ylabel('Penetrating Artery Flow (Q_{PA}) [mL/s]');
legend(legend_entries, 'Location', 'best');

end
