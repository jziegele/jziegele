%% main.m
% This script runs a sensitivity analysis for the Myocardial Perfusion Model.

%% 1. Initialization and Setup
clear;
clc;
addpath(genpath(pwd));

% --- Configuration ---
overwrite_existing = false; % Set to true to re-run and overwrite existing simulations
simulation_input = 0; % ExerciseLevel
extend_patient_data = false;
output_directory = 'simulation_results';
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

% --- Load Master Index File ---
sim_index_file = 'SimIndex.mat';
if exist(sim_index_file, 'file')
    load(sim_index_file, 'SimIndex');
else
    % Initialize SimIndex with the new structure
    SimIndex = struct('inputs', [], 'summaryMetrics', {{}});
end

%% 2. Define Parameter Space
% Define the values for each parameter you want to test
m1_values = [0.7, 1.0, 1.3];
m2_values = [0.9, 1.0, 1.2, 1.4];
m3_values = [0.9, 1.0, 1.2, 1.4];
m4_values = [0.6, 0.8, 1.0, 1.2, 1.4];
m5_values = [1.0, 1.25, 1.5];
m6_values = [0.5, 1.0, 2.0];

% Create a grid of all parameter combinations
[m1_grid, m2_grid, m3_grid, m4_grid, m5_grid, m6_grid] = ndgrid(m1_values, m2_values, m3_values, m4_values, m5_values, m6_values);
parameter_combinations = [m1_grid(:), m2_grid(:), m3_grid(:), m4_grid(:), m5_grid(:), m6_grid(:)];

%% 3. Run Simulations
for i = 1:size(parameter_combinations, 1)
    % --- Get Current Parameter Combination ---
    m_params = parameter_combinations(i, :);

    % --- Check if this simulation has already been run ---
    simulation_exists = false;
    existing_run_index = -1;
    if ~isempty(SimIndex.inputs)
        % Compare the current m_params with all columns in SimIndex.inputs
        [simulation_exists, existing_run_index] = ismember(m_params, SimIndex.inputs', 'rows');
    end

    % --- Run Simulation if it doesn't exist or if overwriting is enabled ---
    if simulation_exists && ~overwrite_existing
        fprintf('Skipping simulation for parameters [%s] (already exists).\n', num2str(m_params));
        continue; % Skip to the next iteration
    else
        if simulation_exists
            fprintf('Overwriting simulation for parameters [%s].\n', num2str(m_params));
            runID = existing_run_index;
        else
            fprintf('Running simulation for parameters [%s].\n', num2str(m_params));
            runID = size(SimIndex.inputs, 2) + 1;
        end
        
        % --- Load Input Data and Set Parameters ---
        inputData = Load_input_data(simulation_input, extend_patient_data);
        m_myo = ones(1, 10);
        m_myo(1:6) = m_params; % Assuming the first 6 parameters are m1-m6
        
        % --- Run Simulation ---
        [t, X, o] = run_myocardial_ODE(inputData, m_myo, simulation_input);
        
        % --- Compute Outputs ---
        [Outputs, ~] = compute_outputs(inputData, t, X, o);
        
        % --- Save Results ---
        results.t = t;
        results.X = X;
        results.o = o;
        results.Outputs = Outputs;
        
        filename = sprintf('simulation%05d.mat', runID);
        save(fullfile(output_directory, filename), 'results');
        
        % --- Update SimIndex ---
        SimIndex.inputs(:, runID) = m_params;
        SimIndex.summaryMetrics{runID} = Outputs;
        
        % --- Save SimIndex ---
        save(sim_index_file, 'SimIndex');
    end
end

disp('Sensitivity analysis complete.');

%% 4. Example: Plotting Specific Results
% After running the analysis, you can use the following code in the command window
% to plot the QPA waveforms for a specific subset of your simulations.

% Example 1: Plotting with a single constant and automatically determined legend variables
% filter1.m1 = 1.0;
% filter1.m2 = [1.0, 1.3]; % m2 will automatically be in the legend
% plot_qpa_by_parameters(filter1);

% Example 2: Hiding constants from the title
% filter2.m1 = 1.0;
% filter2.m2 = [1.0, 1.3];
% plot_qpa_by_parameters(filter2, false);

% Sanbox
filter.m1 = [1.0];
filter.m2 = [1.0];
filter.m3 = [1.0];
filter.m4 = [1.0];
filter.m5 = [1.0];
filter.m6 = [0.5, 1.0, 2.0];
plot_qpa_by_parameters(filter);
%% 5. Example: Plotting Summary Metrics
% After running the analysis, you can use the following code in the command window
% to plot summary metrics.

% Example 1: Plot QPA_mean vs. m2, holding m1 constant
filter_metric.m1 = [0.7, 1.0, 1.3];
filter_metric.m2 = [1.0];% m2 will be on the X-axis
filter_metric.m3 = [0.7, 1.0, 1.3];
filter_metric.m4 = [1.0];
filter_metric.m5 = [1.0];
filter_metric.m6 = [1.0];
plot_summary_metric_by_parameters(filter_metric, 'EM_ratio', 'm1','m3', true);

% Example 2: Plot QPA_mean vs. m2, grouped by m3, holding m1 constant
% filter_metric.m1 = 1.0;
% filter_metric.m2 = [0.8, 1.0, 1.2];
% filter_metric.m3 = [0.9, 1.1]; % m3 will be used for grouping
% plot_summary_metric_by_parameters(filter_metric, 'QPA_mean', 'm2', 'm3');

% Example 3: Plot without constants in the title
% filter_metric.m1 = 1.0;
% filter_metric.m2 = [0.8, 1.0, 1.2];
% plot_summary_metric_by_parameters(filter_metric, 'QPA_mean', 'm2', '', false);

%% 6. Example: Plotting Sensitivity by Output
% This section demonstrates how to use plot_sensitivity_by_output.
% It assumes a 'one-at-a-time' sensitivity analysis where only one parameter
% is varied from its baseline (1.0) while others are held constant at 1.0.
% The SimIndex is filtered to extract these specific variations.

% Define the parameters and multipliers used in the sensitivity analysis
param_names = {'m1', 'm2', 'm3', 'm4', 'm5', 'm6'};
multipliers = {
    [0.7, 1.0, 1.3],    % m1 multipliers
    [0.9, 1.0, 1.2, 1.4], % m2 multipliers
    [0.9, 1.0, 1.2, 1.4], % m3 multipliers
    [0.6, 0.8, 1.0, 1.2, 1.4], % m4 multipliers
    [1.0, 1.25, 1.5],   % m5 multipliers
    [0.5, 1.0, 2.0]     % m6 multipliers
};

% Define the base case parameters from which variations will be made
% You can change these values to specify a different base case.
base_case_params = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]; 

% Define all target field names to plot
target_fieldNames = {
    'QPA_mean', 'EE_ratio', 'EM_ratio', ...
    'avgVolEpi', 'avgVolMid', 'avgVolEndo', 'IMBV_change', ...
    'R11avg', 'R11min', 'R11max', ...
    'R12avg', 'R12min', 'R12max', ...
    'R13avg', 'R13min', 'R13max', ...
    'Rm1avg', 'Rm1min', 'Rm1max', ...
    'Rm2avg', 'Rm2min', 'Rm2max', ...
    'Rm3avg', 'Rm3min', 'Rm3max', ...
    'R21avg', 'R21min', 'R21max', ...
    'R22avg', 'R22min', 'R22max', ...
    'R23avg', 'R23min', 'R23max'
};
target_structType = 'summaryMetrics'; % This matches the field in SimIndex

for field_idx = 1:length(target_fieldNames)
    target_fieldName = target_fieldNames{field_idx};

    fprintf('\nPreparing data for plot_sensitivity_by_output for %s...\n', target_fieldName);

    % Initialize SensitivityResults structure for each fieldName
    SensitivityResults = struct();

    % Populate SensitivityResults for one-at-a-time sensitivity
    % Iterate through each parameter that was varied (m1 to m6)
    for p_idx = 1:length(param_names)
        current_param_name = param_names{p_idx};
        current_multipliers = multipliers{p_idx}; % Get multipliers for the current parameter
        
        % Iterate through each multiplier
        for m_val_multiplier = current_multipliers
            % Convert multiplier to string for field name (e.g., 'mult_0p7')
            % Use the multiplier value for the field name, not the absolute m_val
            mult_str = strrep(sprintf('mult_%.2f', m_val_multiplier), '.', 'p');
            
            % Create a temporary parameter vector with only the current parameter varied from the base case
            temp_params = base_case_params;
            temp_params(p_idx) = m_val_multiplier; % Use the multiplier directly as the absolute value for the varied parameter
            
            % Search SimIndex for this specific parameter combination
            found_index = -1;
            for k = 1:size(SimIndex.inputs, 2)
                if isequal(SimIndex.inputs(:, k)', temp_params)
                    found_index = k;
                    break;
                end
            end

            if found_index ~= -1
                % Populate SensitivityResults in the format expected by plot_sensitivity_by_output
                SensitivityResults.m_myo_params.(sprintf('index%d', p_idx)).(mult_str).(target_structType).(target_fieldName) = ...
                    SimIndex.summaryMetrics{found_index}.(target_fieldName);
            else
                fprintf('Warning: No simulation found for %s = %.2f (others at base case [%s]). This combination might not have been run in the full factorial analysis.\n', current_param_name, m_val_multiplier, num2str(base_case_params));
            end 
        end 
    end 
    % Call the plotting function
    fprintf('Calling plot_sensitivity_by_output for %s...\n', target_fieldName);

    % Find the absolute baseline value (all m parameters = 1.0)
    absolute_baseline_value = NaN;
    baseline_params_all_one = ones(1, 6);
    for k = 1:size(SimIndex.inputs, 2)
        if isequal(SimIndex.inputs(:, k)', baseline_params_all_one)
            absolute_baseline_value = SimIndex.summaryMetrics{k}.(target_fieldName);
            break;
        end
    end

    if isnan(absolute_baseline_value)
        fprintf('Warning: Absolute baseline (all m parameters = 1.0) not found in SimIndex. Cannot plot baseline bar.\n');
    end

    % Define the baseline multipliers (all 1.0 for the absolute baseline)
    baselineMultipliers = [1.3, 1.3, 1.3, 1.3, 1.3, 1.3];

    % Note: The call to plot_sensitivity_by_output might need to be adjusted
    % depending on how it handles the new 'multipliers' cell array.
    % This example passes the entire cell array.
    plot_sensitivity_by_output(SensitivityResults, 'm_myo_params', target_structType, target_fieldName, 0, multipliers, baselineMultipliers, base_case_params);

    fprintf('plot_sensitivity_by_output example complete for %s.\n', target_fieldName);
end % End of loop for target_fieldNames



%% 7. Example: Plotting Pressure and Volume Waveforms
% After running the analysis, you can use the following code in the command window
% to plot the pressure and volume waveforms for a specific subset of your simulations.

filter.m1 = 1.0;
filter.m2 = 1.0;
filter.m3 = 1.0;
filter.m4 = 1.0;
filter.m5 = 1.0;
filter.m6 = m6_values;
% --- Pressure Plot Example ---
plot_pressures_by_parameters(filter);

% --- Volume Plot Example ---
plot_volumes_by_parameters(filter);


%% 8. Example: Plotting All Signals for Specific Cases
% After running the analysis, you can use the following code in the command window
% to generate plots for specific parameter combinations.

filter.m1 = m1_values;
filter.m2 = 1.0;
filter.m3 = 1.0;
filter.m4 = 1.0;
filter.m5 = 1.0;
filter.m6 = 1.0; % Will generate plots for m6=0.5 and m6=1.0
plot_all_signals_for_cases(filter);