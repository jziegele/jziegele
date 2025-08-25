%Simplified analysis script

clc;
clear;

fid = fopen('analysis_log.txt', 'w');
fprintf(fid, 'Log file created.\n');

svd_basis_file = 'SVD_Basis.mat';
if ~exist(svd_basis_file, 'file')
    fprintf(fid, 'SVD_Basis.mat not found.\n');
    error('SVD_Basis.mat not found.');
end
load(svd_basis_file, 'leftSV_raw');
fprintf(fid, 'SVD_Basis.mat loaded.\n');

sim_index_file = 'SimIndex.mat';
if ~exist(sim_index_file, 'file')
    fprintf(fid, 'SimIndex.mat not found.\n');
    error('SimIndex.mat not found.');
end
load(sim_index_file, 'SimIndex');
fprintf(fid, 'SimIndex.mat loaded.\n');

output_directory = 'simulation_results';
num_simulations = size(SimIndex.inputs, 2);
fprintf(fid, 'Number of simulations: %d\n', num_simulations);

num_svd_features = size(leftSV_raw, 2);
all_sim_svd_features = NaN(num_svd_features, num_simulations);

for i = 1:1%1num_simulations
    filename = sprintf('simulation%05d.mat', i);
    file_to_load = fullfile(output_directory, filename);

    if ~exist(file_to_load, 'file')
        fprintf(fid, 'Simulation file not found: %s\n', file_to_load);
        continue;
    end

    load(file_to_load, 'results');

    try
        svd_features = calculate_svd_features(results, leftSV_raw);
        if any(isnan(svd_features))
            fprintf(fid, 'NaNs found in simulation %d\n', i);
        else
            fprintf(fid, 'Successfully calculated features for simulation %d\n', i);
        end
        all_sim_svd_features(:, i) = svd_features;
    catch e
        fprintf(fid, 'Error in simulation %d: %s\n', i, e.message);
        % Do nothing on error, just leave as NaN
    end
end

save('Sim_SVD_Features.mat', 'all_sim_svd_features');
fprintf(fid, 'Analysis complete.\n');
fclose(fid);