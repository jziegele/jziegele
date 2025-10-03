function [sim_svd_features, feature_error_vector, sim_fourier_features] = calculate_svd_features(simulation_result, svd_basis, Pdata, Plv)
% This function calculates the SVD features for a single simulation result.
% It now also returns the vector of errors for each feature.

try
    %% 1. Find the start of the cardiac cycle from Pdata
    t_pdata = Pdata(1,:);
    dplvdt_pdata = Pdata(4,:);

    % Find the time of the first peak of dPLV/dt
    mean_dplvdt_pdata = mean(dplvdt_pdata);
    std_dplvdt_pdata = std(dplvdt_pdata);
    minpeakheight = mean_dplvdt_pdata + 0.8*std_dplvdt_pdata;
    [~, locs] = findpeaks(dplvdt_pdata, 'MinPeakHeight', minpeakheight);
    if isempty(locs)
        [~, max_idx] = max(dplvdt_pdata);
        t_start_of_cycle = t_pdata(max_idx);
    else
        % Find the first peak after a short stabilization period (e.g., 0.1s)
        first_peak_loc = locs(find(t_pdata(locs) > 0.1, 1));
        if isempty(first_peak_loc)
            t_start_of_cycle = t_pdata(locs(1));
        else
            t_start_of_cycle = t_pdata(first_peak_loc);
        end
    end

    %% 2. Extract and Preprocess Simulated Data
    Q_sim_raw = simulation_result.X(:,2); % Q_PA is the 2nd column of X
    Q_sim_raw = Q_sim_raw * 15; % change units from ml/s to ml/min
    t_sim = simulation_result.t;
    tper = simulation_result.Outputs.tper;

    % --- New DYNAMIC Segmentation Logic ---
    t_start_cycle = t_start_of_cycle; % From previous section
    available_duration = t_sim(end) - t_start_cycle;
    if available_duration <= 0, error('No usable signal duration found after first peak.'); end
    num_cycles_available = floor(available_duration / tper);

    if num_cycles_available < 1, error('Signal is less than one full cycle long.'); end

    t_end_cycle = t_start_cycle + (num_cycles_available * tper);
    t1_sim_idx = find(t_sim >= t_start_cycle, 1, 'first');
    t2_sim_idx = find(t_sim < t_end_cycle, 1, 'last');

    if isempty(t1_sim_idx) || isempty(t2_sim_idx), error('Could not find valid start/end indices for segmentation.'); end

    Q_segment = Q_sim_raw(t1_sim_idx:t2_sim_idx);
    t_segment = t_sim(t1_sim_idx:t2_sim_idx);

    %% 2. Perform Fourier Transform
    % Normalize the signal before resampling
    Q_mean = mean(Q_segment);
    % Q_new = Q_segment - Q_mean;
    % Q_new_max = max(Q_new);
    % Q_processed = Q_new / Q_new_max;

    dt_mean = mean(diff(t_segment));
    t_resampled = (t_segment(1):dt_mean:t_segment(end))';
    Q_resampled = interp1(t_segment, Q_segment, t_resampled);
    if mod(length(Q_resampled), 2) ~= 0
        Q_resampled = Q_resampled(1:end-1);
        t_resampled = t_resampled(1:end-1);
    end

    [~, F_Q, fshift] = CalcFourierModes(Q_resampled, t_resampled);

    %% 3. Create Fourier Feature Vector
    N_harmonics = 15;
    sim_fourier_features = zeros(N_harmonics * 2 + 1, 1);
    sim_fourier_features(1) = Q_mean; % Q_mean is the first feature

    f_beat = 1/ tper;
    N_harmonics = 15;
    F_Q_positive = F_Q(fshift >= 0);
    fshift_positive = fshift(fshift >= 0);

    for s = 1:N_harmonics
        f_target = s * f_beat;
        [~, j] = min(abs(fshift_positive - f_target));
        complex_coeff = F_Q_positive(j);
        % Real parts are features 2-16
        sim_fourier_features(s + 1) = real(complex_coeff);
        % Imaginary parts are features 17-31
        sim_fourier_features(s + 1 + N_harmonics) = imag(complex_coeff);
    end

    %% 4. Project onto SVD Basis
    sim_svd_features = svd_basis' * sim_fourier_features;

    %% 5. Calculate Feature-Space Error Vector
    reconstructed_fourier_features = svd_basis * sim_svd_features;
    feature_error_vector = sim_fourier_features - reconstructed_fourier_features;

catch ME
    % If any error occurs, return NaNs and print the error message.
    fprintf('Error in calculate_svd_features: %s\n', ME.message);
    sim_svd_features = NaN(size(svd_basis, 2), 1);
    feature_error_vector = NaN(size(svd_basis, 1), 1);
    sim_fourier_features = NaN(size(svd_basis, 1), 1);
end
end