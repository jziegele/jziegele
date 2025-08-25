%% convert_sim_index.m
% This script converts the old, struct-based SimIndex.mat to the new, streamlined format.

clear;
clc;

sim_index_file = 'SimIndex.mat';
backup_file = 'SimIndex_backup.mat';

% --- 1. Load the data using a try-catch block for robustness ---
try
    fprintf('Attempting to load %s...\n', sim_index_file);
    load(sim_index_file, 'SimIndex');
    fprintf('%s loaded successfully.\n', sim_index_file);
catch
    fprintf('%s not found or failed to load. Attempting to load backup...\n', sim_index_file);
    try
        load(backup_file, 'SimIndex');
        fprintf('%s loaded successfully.\n', backup_file);
    catch e
        error('Failed to load both %s and %s. Error: %s', sim_index_file, backup_file, e.message);
    end
end

% --- 2. Check if conversion is needed ---
if ~isfield(SimIndex, 'filename')
    fprintf('The loaded SimIndex appears to be in the NEW format already. No conversion needed.\n');
    return; % Exit the script
end

% --- 3. Backup the old file (if it doesn't already exist) ---
if ~exist(backup_file, 'file')
    fprintf('Backing up original file to %s...\n', backup_file);
    copyfile(sim_index_file, backup_file);
end

% --- 4. Initialize the new structure ---
fprintf('Initializing new SimIndex structure...\n');
num_simulations = length(SimIndex);
NewSimIndex = struct('inputs', zeros(6, num_simulations), 'summaryMetrics', {cell(1, num_simulations)});

% --- 5. Loop through the old data and populate the new structure ---
fprintf('Converting %d simulation records...\n', num_simulations);
for i = 1:num_simulations
    params_vector = [
        SimIndex(i).inputs.m1;
        SimIndex(i).inputs.m2;
        SimIndex(i).inputs.m3;
        SimIndex(i).inputs.m4;
        SimIndex(i).inputs.m5;
        SimIndex(i).inputs.m6
    ];
    NewSimIndex.inputs(:, i) = params_vector;
    NewSimIndex.summaryMetrics{i} = SimIndex(i).summaryMetrics;
end

% --- 6. Overwrite the old SimIndex variable and save ---
fprintf('Saving converted SimIndex back to %s...\n', sim_index_file);
SimIndex = NewSimIndex;
save(sim_index_file, 'SimIndex');

fprintf('Conversion complete!\n');