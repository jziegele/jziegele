function [sim_svd_features] = calculate_svd_features(simulation_result, svd_basis)
% calculate_svd_features - Calculates SVD features for a simulated waveform.
%
% This function takes the output from a single simulation run, processes its
% flow waveform, performs a Fourier analysis, and then projects the resulting
% Fourier features onto a pre-calculated SVD basis from real experimental data.
%
% Inputs:
%   simulation_result - A struct containing the simulation results for one run.
%                       This should be the 'results' struct from a 'simulationXXXXX.mat' file.
%   svd_basis         - The left singular vectors (leftSV) from the SVD analysis
%                       of the real pig data. This should be a matrix (e.g., 30x17).
%
% Outputs:
%   sim_svd_features  - A vector of the calculated SVD coefficients for the
%                       simulated data.

%% 1. Extract and Preprocess Simulated Data

% Extract flow and time from the results struct
Q_sim_raw = simulation_result.X(:,2); % Q_PA is the 2nd column of X
t_sim = simulation_result.t;
tper = simulation_result.Outputs.tper;

% Repeat the cardiac cycle to get at least 30 beats
num_beats_existing = floor((t_sim(end) - t_sim(1)) / tper);
if num_beats_existing == 0
    num_beats_existing = 1;
end
num_repeats = ceil(30 / num_beats_existing);

Q_repeated = repmat(Q_sim_raw, num_repeats, 1);
t_repeated = t_sim;
for i = 2:num_repeats
    t_repeated = [t_repeated; t_sim(2:end) + t_repeated(end)];
end


% The original analysis used a 30-beat segment. We will replicate that here.
% We select a segment from the end of the simulation to ensure it has reached a steady state.
num_beats_to_capture = 30;
end_time = t_repeated(end);
start_time = end_time - (num_beats_to_capture * tper);
if start_time < t_repeated(1)
    start_time = t_repeated(1);
end

% Create a logical index for the desired time window
t_idx = t_repeated >= start_time & t_repeated <= end_time;
t_segment = t_repeated(t_idx);
Q_segment = Q_repeated(t_idx);

% Normalize the flow waveform in the same way as the real data
Q_mean = mean(Q_segment);
Q_new = Q_segment - Q_mean;
Q_new_max = max(Q_new);

if Q_new_max == 0
    sim_svd_features = zeros(size(svd_basis, 2), 1);
    return;
end

Q_processed = Q_new / Q_new_max;

%% 2. Perform Fourier Analysis

% Ensure the signal is long enough and has an even number of points for FFT
if mod(length(Q_processed), 2) ~= 0
    Q_processed = Q_processed(1:end-1);
    t_segment = t_segment(1:end-1);
end

[~, F_Q, fshift] = CalcFourierModes(Q_processed, t_segment);

% The Fourier analysis on the real data used a specific frequency resolution
% that resulted in harmonics being spaced 30 points apart. We will mimic
% that indexing here.
N_harmonics = 15;

% Find the indices of the desired harmonics (1Hz, 2Hz, ..., 15Hz)
harmonic_indices = zeros(1, N_harmonics);
for s = 1:N_harmonics
    [~, idx] = min(abs(fshift - s));
    harmonic_indices(s) = idx;
end

% Add logging for debugging
global analysis_log_fid; % Declare global variable for file ID
if isempty(analysis_log_fid) || analysis_log_fid == -1
    analysis_log_fid = fopen('analysis_log.txt', 'a'); % Open in append mode
end
fprintf(analysis_log_fid, 'Debug: length(F_Q) = %d, N_harmonics = %d\n', length(F_Q), N_harmonics);
fprintf(analysis_log_fid, 'Debug: fshift min = %f, max = %f\n', min(fshift), max(fshift));
fprintf(analysis_log_fid, 'Debug: Harmonic indices: %s\n', mat2str(harmonic_indices));


% Check if the FFT output is long enough to extract all harmonics
if max(harmonic_indices) > length(F_Q)
    error('The simulated signal is too short or sampled too slowly to extract the required harmonics. Increase simulation time or check sampling rate.');
end

%% 3. Create Fourier Feature Vector

% Initialize the 31-component feature vector (15 real, 15 imaginary, 1 raw flow)
sim_fourier_features = zeros(N_harmonics * 2 + 1, 1);

idx = 1;
for s = 1:N_harmonics
    j = harmonic_indices(s);

    % Get the complex fourier coefficient for the j-th harmonic
    complex_coeff = F_Q(j);

    % Store the real and imaginary parts
    sim_fourier_features(idx) = real(complex_coeff);
    sim_fourier_features(idx + N_harmonics) = imag(complex_coeff);
    idx = idx + 1;
end

% Add the raw flow as the 31st feature
sim_fourier_features(N_harmonics * 2 + 1) = Q_mean;

%% 4. Project onto SVD Basis

% Project the simulation's Fourier features onto the SVD basis from the real data.
% The SVD basis (leftSV) should have dimensions of [31 x num_components].
if size(svd_basis, 1) ~= 31
    error('The provided svd_basis must have 31 rows to match the number of Fourier features (30 Fourier features + 1 raw flow).');
end

sim_svd_features = svd_basis' * sim_fourier_features;

fprintf('Successfully calculated SVD features for the simulation.\n');

end