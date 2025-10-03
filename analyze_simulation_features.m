%% Simplified analysis script with feature error vector calculation

clc;
clear;

fid = fopen('analysis_log.txt', 'w');
fprintf(fid, 'Log file created.\n');

% --- Load Data ---
svd_basis_file = 'SVD_Basis.mat';
if ~exist(svd_basis_file, 'file')
    error('SVD_Basis.mat not found.');
end
load(svd_basis_file, 'leftSV_raw');
fprintf(fid, 'SVD_Basis.mat loaded.\n');

sim_index_file = 'SimIndex.mat';
if ~exist(sim_index_file, 'file')
    error('SimIndex.mat not found.');
end
load(sim_index_file, 'SimIndex');
fprintf(fid, 'SimIndex.mat loaded.\n');

output_directory = 'simulation_results';
num_simulations = size(SimIndex.inputs, 2);
fprintf(fid, 'Number of simulations: %d\n', num_simulations);

% Load the Pdata and PLV for the rest case
InputData_rest = Load_input_data(0);
Pdata_rest = InputData_rest.Pdata;
Plv_rest = InputData_rest.Plv; % Get PLV data
fprintf(fid, 'Pdata and Plv data loaded.\n');

% --- Initialize Results ---
num_features = size(leftSV_raw, 1);
all_sim_svd_features = NaN(size(leftSV_raw, 2), num_simulations);
all_sim_feature_errors = NaN(num_features, num_simulations); % To store the feature error vectors
all_sim_fourier_features = NaN(num_features, num_simulations); % To store the Fourier feature vectors (N_harmonics*2+1)

% --- Main Loop ---
% Clear the function from memory to ensure the latest version is used
clear calculate_svd_features;

for i = 1:num_simulations
    filename = sprintf('simulation%05d.mat', i);
    file_to_load = fullfile(output_directory, filename);

    if ~exist(file_to_load, 'file')
        fprintf(fid, 'Simulation file not found: %s\n', file_to_load);
        continue;
    end

    load(file_to_load, 'results');

    try
        % Call the updated function
        [svd_features, feature_errors, fourier_features] = calculate_svd_features(results, leftSV_raw, Pdata_rest, Plv_rest);

        % Store results
        all_sim_svd_features(:, i) = svd_features;
        all_sim_feature_errors(:, i) = feature_errors;
        all_sim_fourier_features(:, i) = fourier_features; % Store Fourier features
        
        if any(isnan(svd_features))
            fprintf(fid, 'NaNs found in simulation %d\n', i);
        else
            fprintf(fid, 'Processed sim %d.\n', i);
        end
        
    catch e
        fprintf(fid, 'Error in simulation %d: %s\n', i, e.message);
        % Do nothing on error, just leave as NaN
    end
end

% --- Save Results ---
save('Sim_SVD_Features.mat', 'all_sim_svd_features');
save('Sim_Feature_Error.mat', 'all_sim_feature_errors');
save('Sim_Fourier_Features.mat', 'all_sim_fourier_features'); % Save Fourier features
fprintf(fid, '\nAnalysis complete.\n');
fprintf(fid, 'SVD features saved to Sim_SVD_Features.mat\n');
fprintf(fid, 'Feature error vectors saved to Sim_Feature_Error.mat\n');
fprintf(fid, 'Fourier feature vectors saved to Sim_Fourier_Features.mat\n'); % Log Fourier features save
fclose(fid);
