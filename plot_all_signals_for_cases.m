function plot_all_signals_for_cases(filter_params)
% plot_all_signals_for_cases - For each simulation matching the filter, 
% plots arterial/venous pressures/volumes on separate graphs.
%
% Syntax: plot_all_signals_for_cases(filter_params)
%
% Inputs:
%   filter_params - A struct where each field is a parameter name (e.g., 'm1')
%                   and the value is the desired value (or array of values).

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

%% 3. Plot the Waveforms for Each Case
art_pressure_names = {'P11', 'P12', 'P13'};
art_pressure_indices = [3, 5, 7];
ven_pressure_names = {'P21', 'P22', 'P23'};
ven_pressure_indices = [4, 6, 8];

art_volume_names = {'V11', 'V12', 'V13'};
art_volume_indices = [10, 11, 12];
ven_volume_names = {'V21', 'V22', 'V23'};
ven_volume_indices = [13, 14, 15];

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
    X = data.results.X;
    o = data.results.o;
    tper = data.results.Outputs.tper;
    
    % --- Define Time Window (e.g., 5th to 7th cardiac cycle) ---
    x_start = 5 * tper;
    x_end = 7 * tper;
    t_idx = t >= x_start & t < x_end;
    
    % --- Get current parameters for title ---
    current_params = SimIndex.inputs(:, sim_index);
    title_str_params = strjoin(arrayfun(@(j) sprintf('%s=%.2f', param_names{j}, current_params(j)), 1:6, 'UniformOutput', false), ', ');

    % --- Plot Arterial Pressures ---
    figure;
    hold on;
    grid on;
    for p_idx = 1:3
        plot(t(t_idx) - t(find(t_idx,1)), X(t_idx, art_pressure_indices(p_idx)), 'LineWidth', 1.5, 'DisplayName', art_pressure_names{p_idx});
    end
    hold off;
    xlabel('Time (s)');
    ylabel('Pressure (mmHg)');
    title(sprintf('Arterial Pressures for Case: %s', title_str_params));
    legend('Location', 'best');

    % --- Plot Venous Pressures ---
    figure;
    hold on;
    grid on;
    for p_idx = 1:3
        plot(t(t_idx) - t(find(t_idx,1)), X(t_idx, ven_pressure_indices(p_idx)), 'LineWidth', 1.5, 'DisplayName', ven_pressure_names{p_idx});
    end
    hold off;
    xlabel('Time (s)');
    ylabel('Pressure (mmHg)');
    title(sprintf('Venous Pressures for Case: %s', title_str_params));
    legend('Location', 'best');

    % --- Plot Arterial Volumes ---
    figure;
    hold on;
    grid on;
    for v_idx = 1:3
        plot(t(t_idx) - t(find(t_idx,1)), o(t_idx, art_volume_indices(v_idx)), 'LineWidth', 1.5, 'DisplayName', art_volume_names{v_idx});
    end
    hold off;
    xlabel('Time (s)');
    ylabel('Volume (mL)');
    title(sprintf('Arterial Volumes for Case: %s', title_str_params));
    legend('Location', 'best');

    % --- Plot Venous Volumes ---
    figure;
    hold on;
    grid on;
    for v_idx = 1:3
        plot(t(t_idx) - t(find(t_idx,1)), o(t_idx, ven_volume_indices(v_idx)), 'LineWidth', 1.5, 'DisplayName', ven_volume_names{v_idx});
    end
    hold off;
    xlabel('Time (s)');
    ylabel('Volume (mL)');
    title(sprintf('Venous Volumes for Case: %s', title_str_params));
    legend('Location', 'best');
end

end