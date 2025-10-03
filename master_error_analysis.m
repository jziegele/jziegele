%% master_error_analysis.m
% This script consolidates the main error analysis workflows and several
% diagnostic tools into a single file.
% Use the configuration flags in Section 1 to choose which analysis to run.


%% 1. CONFIGURATION
clc;
clear;
close all;

% --- Main Analysis Flags ---
% This analysis answers: "Why did we choose 15 modes?"
% It checks how error changes as you increase the number of Fourier and SVD modes.
run_mode_sensitivity_analysis = true;

% This analysis answers: "Is the error low for all my simulations?"
% It calculates the final reconstruction error for the chosen model design
% across all simulations and generates a histogram.
run_global_reconstruction_error_analysis = true;

% --- Diagnostic Tool Flags ---
run_fft_scaling_check = false;      % Checks your MATLAB FFT function. Doesn't use project data.
run_simple_validation = false;      % Quick check that error decreases from 14 to 15 modes.
run_deep_dive_verification = false; % In-depth visual debug for a SINGLE simulation.

% --- Configuration for Deep Dive ---
simulation_to_verify = 1010; % Simulation index for the deep-dive tool.


%% 2. INITIALIZATION AND COMMON DATA LOADING
fprintf('--- Master Error Analysis & Diagnostic Script ---\n');

output_folder = 'ErrorResults';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
    fprintf('Created output directory: %s\n', output_folder);
end

if run_mode_sensitivity_analysis || run_global_reconstruction_error_analysis
    fprintf('\nLoading data for main analysis...\n');
    svd_basis_file = 'SVD_Basis.mat';
    sim_index_file = 'SimIndex.mat';
    output_directory = 'simulation_results';
    if ~exist(svd_basis_file, 'file')
        error('SVD_Basis.mat not found.');
    end
    load(svd_basis_file, 'leftSV_raw');
    svd_basis = leftSV_raw;
    if ~exist(sim_index_file, 'file')
        error('SimIndex.mat not found.');
    end
    load(sim_index_file, 'SimIndex');
    num_simulations = size(SimIndex.inputs, 2);
    InputData_rest = Load_input_data(0);
    Pdata_rest = InputData_rest.Pdata;
    fprintf('Data loaded.\n');
end


%% 3. ANALYSIS 1: MODE SENSITIVITY
if run_mode_sensitivity_analysis
    fprintf('\n--- Running Analysis 1: Mode Sensitivity ---\n');
    num_fourier_modes_to_check = 15;
    num_svd_modes_to_check = 15;
    all_mse_fourier_recon = zeros(num_fourier_modes_to_check, num_simulations);
    all_nmse_fourier_recon = zeros(num_fourier_modes_to_check, num_simulations);
    all_mse_svd_recon = zeros(num_svd_modes_to_check, num_simulations);
    all_nmse_svd_recon = zeros(num_svd_modes_to_check, num_simulations);

    for sim_idx = 1:num_simulations
        fprintf('Processing simulation %d of %d...\n', sim_idx, num_simulations);
        filename = sprintf('simulation%%05d.mat', sim_idx);
        file_to_load = fullfile(output_directory, filename);
        if ~exist(file_to_load, 'file')
            continue;
        end
        load(file_to_load, 'results');
        try
            Q_sim_raw = results.X(:,2);
            t_sim = results.t;
            tper = results.Outputs.tper;
            PLV = InputData_rest.Plv;
            t_pdata = Pdata_rest(1,:);
            dPLVdt = PLV(2:end) - PLV(1:end-1);
            dPLVdt_smooth = smoothdata(dPLVdt, "movmean", 50);
            [~, locs] = findpeaks(dPLVdt_smooth, 'MinPeakHeight', mean(dPLVdt_smooth) + 0.8*std(dPLVdt_smooth), 'MinPeakDistance', 300);
            t_start_cycle = t_pdata(locs(1));
            num_cycles_available = floor((t_sim(end) - t_start_cycle) / tper);
            t_end_cycle = t_start_cycle + (num_cycles_available * tper);
            t1_sim_idx = find(t_sim >= t_start_cycle, 1, 'first');
            t2_sim_idx = find(t_sim < t_end_cycle, 1, 'last');
            Q_segment = Q_sim_raw(t1_sim_idx:t2_sim_idx);
            t_segment = t_sim(t1_sim_idx:t2_sim_idx);
            if mod(length(t_segment), 2) ~= 0
                t_segment = t_segment(1:end-1);
                Q_segment = Q_segment(1:end-1);
            end
            Q_mean = mean(Q_segment);
            Q_new = Q_segment - Q_mean;
            Q_new_max = max(Q_new);
            Q_processed = Q_new / Q_new_max;
            dt_mean = mean(diff(t_segment));
            t_resampled = (t_segment(1):dt_mean:t_segment(end))';
            Q_resampled = interp1(t_segment, Q_processed, t_resampled);
            if mod(length(Q_resampled), 2) ~= 0
                Q_resampled = Q_resampled(1:end-1);
                t_resampled = t_resampled(1:end-1);
            end
            [~, F_Q, fshift] = CalcFourierModes(Q_resampled, t_resampled);
            f_beat = 1 / tper;
            for n_modes = 1:num_fourier_modes_to_check
                F_Q_reconstructed = zeros(size(F_Q));
                for s = 1:n_modes
                    f_target = s * f_beat;
                    [~, p_idx] = min(abs(fshift - f_target));
                    [~, n_idx] = min(abs(fshift - (-f_target)));
                    F_Q_reconstructed(p_idx) = F_Q(p_idx);
                    F_Q_reconstructed(n_idx) = F_Q(n_idx);
                end
                Q_recon_norm = ifft(ifftshift(F_Q_reconstructed)) * length(t_resampled);
                all_mse_fourier_recon(n_modes, sim_idx) = mean((Q_resampled(:) - Q_recon_norm(:)).^2);
                all_nmse_fourier_recon(n_modes, sim_idx) = sum((Q_resampled(:) - Q_recon_norm(:)).^2) / sum(Q_resampled(:).^2);
            end
            sim_fourier_features = zeros(31, 1);
            sim_fourier_features(1) = Q_mean;
            N_harmonics = 15;
            F_Q_positive = F_Q(fshift >= 0);
            fshift_positive = fshift(fshift >= 0);
            for s = 1:N_harmonics
                f_target = s * f_beat;
                [~, j] = min(abs(fshift_positive - f_target));
                complex_coeff = F_Q_positive(j);
                sim_fourier_features(s + 1) = real(complex_coeff);
                sim_fourier_features(s + 1 + N_harmonics) = imag(complex_coeff);
            end
            sim_svd_features = svd_basis' * sim_fourier_features;
            Q_segment_interp = interp1(t_segment, Q_segment, t_resampled, 'linear', 'extrap');
            for n_svd_modes = 1:num_svd_modes_to_check
                recon_fourier_features = svd_basis(:, 1:n_svd_modes) * sim_svd_features(1:n_svd_modes);
                F_Q_recon_svd = zeros(size(F_Q));
                recon_Q_mean = recon_fourier_features(end);
                for s = 1:15
                    complex_coeff = recon_fourier_features(s) + 1i * recon_fourier_features(s + 15);
                    f_target = s * f_beat;
                    [~, p_idx] = min(abs(fshift - f_target));
                    [~, n_idx] = min(abs(fshift - (-f_target)));
                    F_Q_recon_svd(p_idx) = complex_coeff;
                    F_Q_recon_svd(n_idx) = conj(complex_coeff);
                end
                Q_recon_norm_svd = ifft(ifftshift(F_Q_recon_svd)) * length(t_resampled);
                Q_recon_svd = (Q_recon_norm_svd * Q_new_max) + recon_Q_mean;
                all_mse_svd_recon(n_svd_modes, sim_idx) = mean((Q_segment_interp(:) - Q_recon_svd(:)).^2);
                all_nmse_svd_recon(n_svd_modes, sim_idx) = sum((Q_segment_interp(:) - Q_recon_svd(:)).^2) / sum(Q_segment_interp(:).^2);
            end
        catch e
            fprintf('Error in Analysis 1 for sim %d: %s\n', sim_idx, e.message);
        end
    end
    save(fullfile(output_folder, 'Mode_Sensitivity_Error.mat'), 'all_mse_fourier_recon', 'all_nmse_fourier_recon', 'all_mse_svd_recon', 'all_nmse_svd_recon');
    figure;
    hold on;
    for sim_idx = 1:num_simulations
        plot(1:num_fourier_modes_to_check, all_mse_fourier_recon(:, sim_idx), 'Color', [0.7 0.7 0.7 0.3]);
    end
    plot(1:num_fourier_modes_to_check, mean(all_mse_fourier_recon, 2), 'k-', 'LineWidth', 2);
    hold off;
    title('Fourier MSE vs. Modes');
    xlabel('Modes');
    ylabel('MSE');
    grid on;
    saveas(gcf, fullfile(output_folder, 'Mode_Sensitivity_Fourier_MSE.png'));
    figure;
    hold on;
    for sim_idx = 1:num_simulations
        plot(1:num_fourier_modes_to_check, all_nmse_fourier_recon(:, sim_idx), 'Color', [0.7 0.7 0.7 0.3]);
    end
    plot(1:num_fourier_modes_to_check, mean(all_nmse_fourier_recon, 2), 'k-', 'LineWidth', 2);
    hold off;
    title('Fourier NMSE vs. Modes');
    xlabel('Modes');
    ylabel('NMSE');
    grid on;
    saveas(gcf, fullfile(output_folder, 'Mode_Sensitivity_Fourier_NMSE.png'));
    figure;
    hold on;
    for sim_idx = 1:num_simulations
        plot(1:num_svd_modes_to_check, all_mse_svd_recon(:, sim_idx), 'Color', [0.7 0.7 0.7 0.3]);
    end
    plot(1:num_svd_modes_to_check, mean(all_mse_svd_recon, 2), 'k-', 'LineWidth', 2);
    hold off;
    title('SVD MSE vs. Modes');
    xlabel('Modes');
    ylabel('MSE');
    grid on;
    saveas(gcf, fullfile(output_folder, 'Mode_Sensitivity_SVD_MSE.png'));
    figure;
    hold on;
    for sim_idx = 1:num_simulations
        plot(1:num_svd_modes_to_check, all_nmse_svd_recon(:, sim_idx), 'Color', [0.7 0.7 0.7 0.3]);
    end
    plot(1:num_svd_modes_to_check, mean(all_nmse_svd_recon, 2), 'k-', 'LineWidth', 2);
    hold off;
    title('SVD NMSE vs. Modes');
    xlabel('Modes');
    ylabel('NMSE');
    grid on;
    saveas(gcf, fullfile(output_folder, 'Mode_Sensitivity_SVD_NMSE.png'));
    fprintf('--- Analysis 1: Mode Sensitivity Complete ---\n');
end


%% 4. ANALYSIS 2: GLOBAL RECONSTRUCTION ERROR
if run_global_reconstruction_error_analysis
    fprintf('\n--- Running Analysis 2: Global Reconstruction Error ---\n');
    all_sim_svd_features_global = NaN(size(svd_basis, 2), num_simulations);
    all_sim_errors_global = NaN(2, num_simulations);
    for i = 1:num_simulations
        fprintf('Processing simulation %d of %d...\n', i, num_simulations);
        filename = sprintf('simulation%%05d.mat', i);
        file_to_load = fullfile(output_directory, filename);
        if ~exist(file_to_load, 'file')
            continue;
        end
        load(file_to_load, 'results');
        try
            Q_sim_raw = results.X(:,2);
            t_sim = results.t;
            tper = results.Outputs.tper;
            PLV = InputData_rest.Plv;
            t_pdata = Pdata_rest(1,:);
            dPLVdt = PLV(2:end) - PLV(1:end-1);
            dPLVdt_smooth = smoothdata(dPLVdt, "movmean", 50);
            [~, locs] = findpeaks(dPLVdt_smooth, 'MinPeakHeight', mean(dPLVdt_smooth) + 0.8*std(dPLVdt_smooth), 'MinPeakDistance', 300);
            t_start_cycle = t_pdata(locs(1));
            num_cycles_available = floor((t_sim(end) - t_start_cycle) / tper);
            t_end_cycle = t_start_cycle + (num_cycles_available * tper);
            t1_sim_idx = find(t_sim >= t_start_cycle, 1, 'first');
            t2_sim_idx = find(t_sim < t_end_cycle, 1, 'last');
            Q_segment = Q_sim_raw(t1_sim_idx:t2_sim_idx);
            t_segment = t_sim(t1_sim_idx:t2_sim_idx);
            if mod(length(t_segment), 2) ~= 0
                t_segment = t_segment(1:end-1);
                Q_segment = Q_segment(1:end-1);
            end
            Q_mean = mean(Q_segment);
            Q_new = Q_segment - Q_mean;
            Q_new_max = max(Q_new);
            Q_processed = Q_new / Q_new_max;
            dt_mean = mean(diff(t_segment));
            t_resampled = (t_segment(1):dt_mean:t_segment(end))';
            Q_resampled = interp1(t_segment, Q_processed, t_resampled);
            if mod(length(Q_resampled), 2) ~= 0
                Q_resampled = Q_resampled(1:end-1);
                t_resampled = t_resampled(1:end-1);
            end
            [~, F_Q, fshift] = CalcFourierModes(Q_resampled, t_resampled);
            N_harmonics = 15;
            f_beat = 1/tper;
            sim_fourier_features = zeros(N_harmonics * 2 + 1, 1);
            sim_fourier_features(1) = Q_mean;
            F_Q_positive = F_Q(fshift >= 0);
            fshift_positive = fshift(fshift >= 0);
            for s = 1:N_harmonics
                f_target = s * f_beat;
                [~,j] = min(abs(fshift_positive - f_target));
                complex_coeff = F_Q_positive(j);
                sim_fourier_features(s + 1) = real(complex_coeff);
                sim_fourier_features(s + 1 + N_harmonics) = imag(complex_coeff);
            end
            sim_svd_features = svd_basis' * sim_fourier_features;
            recon_fourier_features = svd_basis * sim_svd_features;
            F_Q_recon_svd = zeros(size(F_Q));
            recon_Q_mean = recon_fourier_features(end);
            for s = 1:N_harmonics
                complex_coeff = recon_fourier_features(s) + 1i*recon_fourier_features(s+N_harmonics);
                f_target = s * f_beat;
                [~,p_idx] = min(abs(fshift - f_target));
                [~,n_idx] = min(abs(fshift - (-(s*f_beat))));
                F_Q_recon_svd(p_idx) = complex_coeff;
                F_Q_recon_svd(n_idx) = conj(complex_coeff);
            end
            Q_recon_norm_svd = ifft(ifftshift(F_Q_recon_svd)) * length(t_resampled);
            Q_recon_svd = (Q_recon_norm_svd * Q_new_max) + recon_Q_mean;
            Q_segment_interp = interp1(t_segment, Q_segment, t_resampled);
            nrmse_svd = sqrt(sum((Q_segment_interp - Q_recon_svd).^2) / sum(Q_segment_interp.^2));
            F_Q_recon_fourier = zeros(size(F_Q));
            for s = 1:N_harmonics
                f_target = s * f_beat;
                [~,p_idx] = min(abs(fshift - f_target));
                [~,n_idx] = min(abs(fshift - (-(s*f_beat))));
                F_Q_recon_fourier(p_idx) = F_Q(p_idx);
                F_Q_recon_fourier(n_idx) = F_Q(n_idx);
            end
            Q_recon_norm_fourier = ifft(ifftshift(F_Q_recon_fourier)) * length(t_resampled);
            Q_recon_fourier = (Q_recon_norm_fourier * Q_new_max) + Q_mean;
            nrmse_fourier = sqrt(sum((Q_segment_interp - Q_recon_fourier).^2) / sum(Q_segment_interp.^2));
            all_sim_svd_features_global(:, i) = sim_svd_features;
            all_sim_errors_global(1, i) = nrmse_svd;
            all_sim_errors_global(2, i) = nrmse_fourier;
        catch e
            fprintf('Error in Analysis 2 for sim %d: %s\n', i, e.message);
        end
    end
    save(fullfile(output_folder, 'Global_Error_and_Features.mat'), 'all_sim_svd_features_global', 'all_sim_errors_global');
    svd_errors = all_sim_errors_global(1, :);
    fourier_errors = all_sim_errors_global(2, :);
    svd_errors = svd_errors(~isnan(svd_errors));
    fourier_errors = fourier_errors(~isnan(fourier_errors));
    figure('Position', [100, 100, 1200, 500]);
    max_err = max([svd_errors, fourier_errors]) * 100;
    if isempty(max_err) || max_err == 0
        max_err = 10;
    end
    subplot(1, 2, 1);
    histogram(fourier_errors * 100, 'BinWidth', 0.5);
    title('Global Fourier Error');
    xlabel('NRMSE %');
    ylabel('Count');
    grid on;
    xlim([0, max_err + 1]);
    subplot(1, 2, 2);
    histogram(svd_errors * 100, 'BinWidth', 0.5);
    title('Global SVD Error');
    xlabel('NRMSE %');
    ylabel('Count');
    grid on;
    xlim([0, max_err + 1]);
    sgtitle('Distribution of Reconstruction Errors');
    saveas(gcf, fullfile(output_folder, 'Global_Reconstruction_Error_Histogram.png'));
    fprintf('--- Analysis 2: Global Reconstruction Error Complete ---\n');
end


%% 5. DIAGNOSTIC TOOL: FFT SCALING VERIFICATION
if run_fft_scaling_check
    fprintf('\n--- Running Diagnostic: FFT Scaling Verification ---\n');
    close all;
    fs = 1000;
    t = 0:1/fs:1-1/fs;
    f1 = 10;
    f2 = 30;
    x_original = sin(2*pi*f1*t) + 0.5*cos(2*pi*f2*t);
    x_reconstructed = ifft(fft(x_original));
    max_abs_error = max(abs(x_original - x_reconstructed));
    fprintf('Max error between original and FFT-reconstructed signal: %e\n', max_abs_error);
    if max_abs_error < 1e-12
        fprintf('Conclusion: FFT/IFFT scaling is correct.\n');
    else
        fprintf('Conclusion: FFT/IFFT scaling issue detected.\n');
    end
    figure;
    plot(t, x_original, 'b-', 'LineWidth', 2);
    hold on;
    plot(t, real(x_reconstructed), 'r--');
    title('Diagnostic: FFT Scaling Verification');
    legend('Original', 'Reconstructed');
    grid on;
    saveas(gcf, fullfile(output_folder, 'Diagnostic_FFT_Scaling_Check.png'));
    fprintf('--- Diagnostic Complete ---\n');
end


%% 6. DIAGNOSTIC TOOL: SIMPLE SVD VALIDATION
if run_simple_validation
    fprintf('\n--- Running Diagnostic: Simple SVD Validation ---\n');
    close all;
    try
        load('SVD_Basis.mat', 'leftSV_raw');
        load(fullfile('simulation_results', sprintf('simulation%%05d.mat', 1)), 'results');
        InputData_rest = Load_input_data(0);
        Pdata_rest = InputData_rest.Pdata;
        Q_sim_raw = results.X(:,2);
        t_sim = results.t;
        tper = results.Outputs.tper;
        PLV = InputData_rest.Plv;
        t_pdata = Pdata_rest(1,:);
        dPLVdt = PLV(2:end) - PLV(1:end-1);
        dPLVdt_smooth = smoothdata(dPLVdt, "movmean", 50);
        [~, locs] = findpeaks(dPLVdt_smooth, 'MinPeakHeight', mean(dPLVdt_smooth) + 0.8*std(dPLVdt_smooth), 'MinPeakDistance', 300);
        t_start_cycle = t_pdata(locs(1));
        num_cycles_available = floor((t_sim(end) - t_start_cycle) / tper);
        t_end_cycle = t_start_cycle + (num_cycles_available * tper);
        t1_sim_idx = find(t_sim >= t_start_cycle, 1, 'first');
        t2_sim_idx = find(t_sim < t_end_cycle, 1, 'last');
        Q_segment = Q_sim_raw(t1_sim_idx:t2_sim_idx);
        t_segment = t_sim(t1_sim_idx:t2_sim_idx);
        if mod(length(t_segment),2)~=0
            t_segment=t_segment(1:end-1);
            Q_segment=Q_segment(1:end-1);
        end
        Q_mean = mean(Q_segment);
        Q_new = Q_segment - Q_mean;
        Q_new_max = max(Q_new);
        Q_processed = Q_new / Q_new_max;
        dt_mean = mean(diff(t_segment));
        t_resampled = (t_segment(1):dt_mean:t_segment(end))';
        Q_resampled = interp1(t_segment, Q_processed, t_resampled);
        if mod(length(Q_resampled),2)~=0
            Q_resampled=Q_resampled(1:end-1);
            t_resampled=t_resampled(1:end-1);
        end
        [~, F_Q, fshift] = CalcFourierModes(Q_resampled, t_resampled);
        f_beat = 1/tper;
        sim_fourier_features = zeros(31,1);
        sim_fourier_features(1) = Q_mean;
        N_harmonics = 15;
        F_Q_positive = F_Q(fshift>=0);
        fshift_positive = fshift(fshift>=0);
        for s=1:N_harmonics
            [~,j] = min(abs(fshift_positive - (s*f_beat)));
            complex_coeff = F_Q_positive(j);
            sim_fourier_features(s + 1) = real(complex_coeff);
            sim_fourier_features(s + 1 + N_harmonics) = imag(complex_coeff);
        end
        sim_svd_features = leftSV_raw' * sim_fourier_features;
        recon_fourier_A = leftSV_raw(:,1:14) * sim_svd_features(1:14);
        F_Q_recon_A = zeros(size(F_Q));
        for s=1:15
            complex_coeff = recon_fourier_A(s) + 1i*recon_fourier_A(s+15);
            [~,p_idx]=min(abs(fshift-(s*f_beat)));
            [~,n_idx]=min(abs(fshift-(-(s*f_beat))));
            F_Q_recon_A(p_idx)=complex_coeff;
            F_Q_recon_A(n_idx)=conj(complex_coeff);
        end
        Q_recon_norm_A = ifft(ifftshift(F_Q_recon_A)) * length(t_resampled);
        error_14 = mean((Q_resampled(:) - Q_recon_norm_A(:)).^2);
        recon_fourier_B = leftSV_raw(:,1:15) * sim_svd_features(1:15);
        F_Q_recon_B = zeros(size(F_Q));
        for s=1:15
            complex_coeff = recon_fourier_B(s) + 1i*recon_fourier_B(s+15);
            [~,p_idx]=min(abs(fshift-(s*f_beat)));
            [~,n_idx]=min(abs(fshift-(-(s*f_beat))));
            F_Q_recon_B(p_idx)=complex_coeff;
            F_Q_recon_B(n_idx)=conj(complex_coeff);
        end
        Q_recon_norm_B = ifft(ifftshift(F_Q_recon_B)) * length(t_resampled);
        error_15 = mean((Q_resampled(:) - Q_recon_norm_B(:)).^2);
        fprintf('MSE with 14 modes: %e\n', error_14);
        fprintf('MSE with 15 modes: %e\n', error_15);
        if error_15 < error_14
            fprintf('SUCCESS: Error decreased with more modes.\n');
        else
            fprintf('FAILURE: Error did NOT decrease.\n');
        end
    catch e
        fprintf('Error during Simple SVD Validation: %s\n', e.message);
    end
    fprintf('--- Diagnostic Complete ---\n');
end


%% 7. DIAGNOSTIC TOOL: DEEP DIVE SVD VERIFICATION
if run_deep_dive_verification
    fprintf('\n--- Running Diagnostic: Deep Dive for Sim %d ---\n', simulation_to_verify);
    close all;
    try
        load('SVD_Basis.mat', 'leftSV_raw');
        sim_filename = sprintf('simulation%%05d.mat', simulation_to_verify);
        load(fullfile('simulation_results', sim_filename), 'results');
        InputData_rest = Load_input_data(0);
        Pdata_rest = InputData_rest.Pdata;
        Q_sim_raw = results.X(:,2);
        t_sim = results.t;
        tper = results.Outputs.tper;
        PLV = InputData_rest.Plv;
        t_pdata = Pdata_rest(1,:);
        dPLVdt = PLV(2:end) - PLV(1:end-1);
        dPLVdt_smooth = smoothdata(dPLVdt, "movmean", 50);
        [~, locs] = findpeaks(dPLVdt_smooth, 'MinPeakHeight', mean(dPLVdt_smooth) + 0.8*std(dPLVdt_smooth), 'MinPeakDistance', 300);
        t_start_cycle = t_pdata(locs(1));
        t_end_cycle = t_start_cycle + (10 * tper);
        t1_sim_idx = find(t_sim >= t_start_cycle, 1, 'first');
        t2_sim_idx = find(t_sim < t_end_cycle, 1, 'last');
        Q_segment = Q_sim_raw(t1_sim_idx:t2_sim_idx);
        t_segment = t_sim(t1_sim_idx:t2_sim_idx);
        if mod(length(t_segment),2)~=0
            t_segment=t_segment(1:end-1);
            Q_segment=Q_segment(1:end-1);
        end
        Q_mean = mean(Q_segment);
        Q_new = Q_segment - Q_mean;
        Q_new_max = max(Q_new);
        Q_processed = Q_new / Q_new_max;
        dt_mean = mean(diff(t_segment));
        t_resampled = (t_segment(1):dt_mean:t_segment(end))';
        Q_resampled = interp1(t_segment, Q_processed, t_resampled);
        if mod(length(Q_resampled),2)~=0
            Q_resampled=Q_resampled(1:end-1);
            t_resampled=t_resampled(1:end-1);
        end
        [~, F_Q, fshift] = CalcFourierModes(Q_resampled, t_resampled);
        N_harmonics = 15;
        f_beat = 1/tper;
        sim_fourier_features = zeros(N_harmonics*2+1,1);
        sim_fourier_features(1) = Q_mean;
        F_Q_positive = F_Q(fshift>=0);
        fshift_positive = fshift(fshift>=0);
        for s=1:N_harmonics
            [~,j]=min(abs(fshift_positive-(s*f_beat)));
            complex_coeff=F_Q_positive(j);
            sim_fourier_features(s+1)=real(complex_coeff);
            sim_fourier_features(s+1+N_harmonics)=imag(complex_coeff);
        end
        sim_svd_features = leftSV_raw' * sim_fourier_features;
        F_Q_recon_15 = zeros(size(F_Q));
        for s=1:N_harmonics
            [~,p_idx]=min(abs(fshift-(s*f_beat)));
            [~,n_idx]=min(abs(fshift-(-(s*f_beat))));
            complex_coeff=F_Q(p_idx);
            F_Q_recon_15(p_idx)=complex_coeff;
            F_Q_recon_15(n_idx)=conj(complex_coeff);
        end
        Q_recon_norm_15 = ifft(ifftshift(F_Q_recon_15)) * length(t_resampled);
        Q_recon_from_15 = (Q_recon_norm_15 * Q_new_max) + Q_mean;
        recon_fourier_features = leftSV_raw * sim_svd_features;
        F_Q_recon_svd = zeros(size(F_Q));
        recon_Q_mean = recon_fourier_features(end);
        for s=1:N_harmonics
            complex_coeff = recon_fourier_features(s) + 1i*recon_fourier_features(s+N_harmonics);
            [~,p_idx]=min(abs(fshift-(s*f_beat)));
            [~,n_idx]=min(abs(fshift-(-(s*f_beat))));
            F_Q_recon_svd(p_idx)=complex_coeff;
            F_Q_recon_svd(n_idx)=conj(complex_coeff);
        end
        Q_recon_norm_svd = ifft(ifftshift(F_Q_recon_svd)) * length(t_resampled);
        Q_recon_from_svd = (Q_recon_norm_svd * Q_new_max) + recon_Q_mean;
        figure('Position',[100,100,1200,600]);
        hold on;
        plot(t_segment, Q_segment, 'b-', 'LineWidth', 2.5);
        plot(t_resampled, Q_resampled, 'k.', 'MarkerSize', 8);
        plot(t_resampled, Q_recon_from_15, '--', 'Color', '#FFA500', 'LineWidth', 2);
        plot(t_resampled, Q_recon_from_svd, 'r:', 'LineWidth', 2);
        legend('Original', 'Resampled', 'Recon (15 Harmonics)', 'Recon (SVD)');
        grid on;
        title(sprintf('Deep Dive for Sim %d', simulation_to_verify));
        saveas(gcf, fullfile(output_folder, sprintf('Diagnostic_Deep_Dive_Sim%%05d.png', simulation_to_verify)));
    catch e
        fprintf('Error during Deep Dive Verification: %s\n', e.message);
    end
    fprintf('--- Diagnostic Complete ---\n');
end


if ~run_mode_sensitivity_analysis && ~run_global_reconstruction_error_analysis && ~run_fft_scaling_check && ~run_simple_validation && ~run_deep_dive_verification
    fprintf('\nWarning: All analysis and diagnostic flags were set to false. Nothing to run.\n');
end

fprintf('--- Master Script Finished ---\n');