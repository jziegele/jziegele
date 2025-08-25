function [Outputs,OxygenInput] = compute_outputs(inputData, t, X, o)
%%
% COMPUTE_OUTPUTS
%   Calculates various physiological output metrics from the raw simulation results
%   (time-series data from the ODE solver).
%
% Inputs:
%   inputData - Structure containing input parameters (HR, tper, etc.)
%   t         - Time vector from ODE solver (s)
%   X         - State variables from ODE solver
%   o         - Detailed intermediate outputs from dXdT_myocardium (O array)
%
% Outputs:
%   Outputs     - Structure containing calculated output metrics (e.g., mean flows, ratios, volumes)
%   OxygenInput - Structure containing selected time-series data for potential oxygen transport calculations

HR = inputData.HR; % Heart Rate (bpm)
tper = inputData.tper; % Cardiac Period (s)

% Extracting relevant time-series data from X and o arrays
% Note: The indexing for 'o' corresponds to the order defined in the 'O' output of dXdT_myocardium.m
Q_PA = X(:,2); % Flow at the start of the circuit (Penetrating Artery Flow) (mL/s)

% Arterial Flows (Q1x)
Q_11 = o(:,1); % Arterial flow - Epicardial layer (mL/s)
Q_12 = o(:,2); % Arterial flow - Mid-wall layer (mL/s)
Q_13 = o(:,3); % Arterial flow - Endocardial layer (mL/s)

% Microvascular Flows (Qm x)
Q_m1 = o(:,4); % Microvascular flow - Epicardial layer (mL/s)
Q_m2 = o(:,5); % Microvascular flow - Mid-wall layer (mL/s)
Q_m3 = o(:,6); % Microvascular flow - Endocardial layer (mL/s)

% Venous Flows (Q2x)
Q_21 = o(:,7); % Venous flow - Epicardial layer (mL/s)
Q_22 = o(:,8); % Venous flow - Mid-wall layer (mL/s)
Q_23 = o(:,9); % Venous flow - Endocardial layer (mL/s)

% Arterial Volumes (V1x)
V_11 = o(:,10); % Arterial volume - Epicardial layer (mL)
V_12 = o(:,11); % Arterial volume - Mid-wall layer (mL)
V_13 = o(:,12); % Arterial volume - Endocardial layer (mL)

% Venous Volumes (V2x)
V_21 = o(:,13); % Venous volume - Epicardial layer (mL)
V_22 = o(:,14); % Venous volume - Mid-wall layer (mL)
V_23 = o(:,15); % Venous volume - Endocardial layer (mL)

% Resistances (R1x, R2x, Rmx) - Note: These are instantaneous resistances
% R1x (Arterial Resistances)
R_11 = o(:,19); % Arterial resistance - Epicardial layer (mmHg/(mL/s))
R_12 = o(:,20); % Arterial resistance - Mid-wall layer (mmHg/(mL/s))
R_13 = o(:,21); % Arterial resistance - Endocardial layer (mmHg/(mL/s))

% Rmx (Microvascular Resistances)
R_m1 = o(:,25); % Microvascular resistance - Epicardial layer (mmHg/(mL/s))
R_m2 = o(:,26); % Microvascular resistance - Mid-wall layer (mmHg/(mL/s))
R_m3 = o(:,27); % Microvascular resistance - Endocardial layer (mmHg/(mL/s))

% R2x (Venous Resistances)
R_21 = o(:,22); % Venous resistance - Epicardial layer (mmHg/(mL/s))
R_22 = o(:,23); % Venous resistance - Mid-wall layer (mmHg/(mL/s))
R_23 = o(:,24); % Venous resistance - Endocardial layer (mmHg/(mL/s))

% Other Volumes and Flows
V_PA = o(:,28); % Penetrating Artery Volume (mL)
V_PV = o(:,29); % Penetrating Vein Volume (mL)
Q_PV = o(:,30); % Flow at the end of the circuit (Penetrating Vein Flow) (mL/s)

% Intramyocardial Pressures (P_im) - from dXdT_myocardium O array
P_im1 = o(:,16); % Epicardial Intramyocardial Pressure (mmHg)
P_im3 = o(:,18); % Endocardial Intramyocardial Pressure (mmHg)

% State variables from X array (Pressures)
P11 = X(:,3); % Epicardial Arterial Pressure (mmHg)
P21 = X(:,4); % Epicardial Venous Pressure (mmHg)
P12 = X(:,5); % Mid-wall Arterial Pressure (mmHg)
P22 = X(:,6); % Mid-wall Venous Pressure (mmHg)
P13 = X(:,7); % Endocardial Arterial Pressure (mmHg)
P23 = X(:,8); % Endocardial Venous Pressure (mmHg)

% Calculate total volumes for each layer
volume_epi = V_11 + V_21; % Total Epicardial Volume (mL)
volume_mid = V_12 + V_22; % Total Mid-wall Volume (mL)
volume_endo = V_13 + V_23; % Total Endocardial Volume (mL)

% Define the desired time window for the 5th to 7th cardiac cycle
x_start_desired = 5 * tper;
x_end_desired = 7 * tper;

% Determine the actual time window to plot
if t(end) < x_start_desired
    % Not even 5 cardiac cycles available, plot all available data
    x_start = t(1);
    x_end = t(end);
else
    % At least 5 cardiac cycles available
    x_start = x_start_desired;
    if t(end) < x_end_desired
        % Less than 7 cardiac cycles available, plot up to the end of data
        x_end = t(end);
    else
        % 7 or more cardiac cycles available, plot the desired range
        x_end = x_end_desired;
    end
end

% Create a logical index for the determined time window
t_idx = t >= x_start & t < x_end;
idx_first = find(t_idx,1,'first');
idx_last = find(t_idx,1,'last');
t_range = t(idx_first:idx_last);

% --- Helper function for correct time-averaged mean ---
function avg = time_avg(time_vec, data_vec)
    if length(time_vec) < 2
        avg = mean(data_vec); % Fallback for single point
        return;
    end
    avg = trapz(time_vec, data_vec) / (time_vec(end) - time_vec(1));
end

% --- Calculate Output Metrics using correct time-averaged mean ---
epsilon = 1e-9; % A small number to prevent division by zero

Outputs.QPA_mean = time_avg(t_range, Q_PA(idx_first:idx_last));
Outputs.EE_ratio = time_avg(t_range, Q_13(idx_first:idx_last)) / (time_avg(t_range, Q_11(idx_first:idx_last)) + epsilon);
Outputs.EM_ratio = time_avg(t_range, Q_13(idx_first:idx_last)) / (time_avg(t_range, Q_12(idx_first:idx_last)) + epsilon);
Outputs.avgVolEpi = time_avg(t_range, volume_epi(idx_first:idx_last));
Outputs.avgVolMid = time_avg(t_range, volume_mid(idx_first:idx_last));
Outputs.avgVolEndo = time_avg(t_range, volume_endo(idx_first:idx_last));

Vol_1per = volume_epi(t_idx) + volume_mid(t_idx) + volume_endo(t_idx);
Outputs.IMBV_change = (max(Vol_1per) - min(Vol_1per)) / (max(Vol_1per) + epsilon); % Intramyocardial Blood Volume Change
Outputs.tper = tper; % Add tper to Outputs struct

% Average, Min, and Max Resistances
Outputs.R11avg = time_avg(t_range, R_11(idx_first:idx_last)); Outputs.R11min = min(R_11(idx_first:idx_last)); Outputs.R11max = max(R_11(idx_first:idx_last));
Outputs.R12avg = time_avg(t_range, R_12(idx_first:idx_last)); Outputs.R12min = min(R_12(idx_first:idx_last)); Outputs.R12max = max(R_12(idx_first:idx_last));
Outputs.R13avg = time_avg(t_range, R_13(idx_first:idx_last)); Outputs.R13min = min(R_13(idx_first:idx_last)); Outputs.R13max = max(R_13(idx_first:idx_last));
Outputs.Rm1avg = time_avg(t_range, R_m1(idx_first:idx_last)); Outputs.Rm1min = min(R_m1(idx_first:idx_last)); Outputs.Rm1max = max(R_m1(idx_first:idx_last));
Outputs.Rm2avg = time_avg(t_range, R_m2(idx_first:idx_last)); Outputs.Rm2min = min(R_m2(idx_first:idx_last)); Outputs.Rm2max = max(R_m2(idx_first:idx_last));
Outputs.Rm3avg = time_avg(t_range, R_m3(idx_first:idx_last)); Outputs.Rm3min = min(R_m3(idx_first:idx_last)); Outputs.Rm3max = max(R_m3(idx_first:idx_last));
Outputs.R21avg = time_avg(t_range, R_21(idx_first:idx_last)); Outputs.R21min = min(R_21(idx_first:idx_last)); Outputs.R21max = max(R_21(idx_first:idx_last));
Outputs.R22avg = time_avg(t_range, R_22(idx_first:idx_last)); Outputs.R22min = min(R_22(idx_first:idx_last)); Outputs.R22max = max(R_22(idx_first:idx_last));
Outputs.R23avg = time_avg(t_range, R_23(idx_first:idx_last)); Outputs.R23min = min(R_23(idx_first:idx_last)); Outputs.R23max = max(R_23(idx_first:idx_last));

% Package selected time-series data into OxygenInput struct
OxygenInput = struct('T', t, 'V_PA', V_PA, 'V_11', V_11, 'V_12', V_12, 'V_13', V_13, ...
    'V_21', V_21, 'V_22', V_22, 'V_23', V_23, 'V_PV', V_PV, 'Q_PA', Q_PA, ...
    'Q_11', Q_11, 'Q_m1', Q_m1, 'Q_21', Q_21, 'Q_12', Q_12, 'Q_m2', Q_m2, ...
    'Q_22', Q_22, 'Q_13', Q_13, 'Q_m3', Q_m3, 'Q_23', Q_23, 'Q_PV', Q_PV);
end