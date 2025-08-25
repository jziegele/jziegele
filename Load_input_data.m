function [InputData] = Load_input_data(simulation_input, extend_data)
%
% Load_input_data
%   Loads input data based on the type of input provided.
%   It can handle a numeric ExerciseLevel or a string patientFileName.
%
% Inputs:
%   simulation_input - Can be either:
%                      - A numeric scalar (0-4) for ExerciseLevel
%                      - A string for the patient .mat file name
%   extend_data      - Optional boolean (true/false). If true and using a 
%                      patient file, the single-cycle data will be repeated
%                      to create a longer time-series.
%
% Outputs:
%   InputData     - Structure containing processed signal characteristics

% Set default for extend_data if not provided
if nargin < 2
    extend_data = false;
end

dt = 1/1000; % Time step for interpolation (s)

if isnumeric(simulation_input)
    % Handle numeric input (ExerciseLevel)
    ExerciseLevel = simulation_input;
    load("AllResults.mat"); % Loads Myo_Inputs_Rest, Myo_Inputs_MildE1, etc.
    
    switch ExerciseLevel
        case 0
            Myo = Myo_Inputs_Rest;
            HR = 64;
        case 1
            Myo = Myo_Inputs_MildE1;
            HR = 90;
        case 2
            Myo = Myo_Inputs_MildE2;
            HR = 120;
        case 3
            Myo = Myo_Inputs_MildE3;
            HR = 150;
        case 4
            Myo = Myo_Inputs_MaxE;
            HR = 180;
        otherwise
            error('Invalid ExerciseLevel. Must be between 0 and 4.');
    end
    
    t_orig = Myo.t; 
    Pao_orig = Myo.P_SA;
    Plv_orig = Myo.P_LV;
    SIP_orig = Myo.Y; 
    Pra_orig = Myo.P_RA;
    
else
    % Handle string input (patientFileName)
    patientFileName = simulation_input;
    patientDataPath = fullfile('StartingPatients4Vic', patientFileName);
    
    if ~exist(patientDataPath, 'file')
        error('Patient data file not found: %s', patientDataPath);
    end
    
    loaded_data = load(patientDataPath);
    
    if isfield(loaded_data, 'timecourse') && isfield(loaded_data, 'inputs')
        Myo = loaded_data.timecourse;
        HR = loaded_data.inputs.HR;
    else
        error('The patient file %s does not contain the required ''timecourse'' and ''inputs'' structs.', patientFileName);
    end
    
    t_orig = Myo.T; 
    Pao_orig = Myo.P_SA;
    Plv_orig = Myo.P_LV;
    SIP_orig = Myo.Y;
    Pra_orig = Myo.P_RA;

end

% --- Optional Data Extension ---
if extend_data
    num_cycles_to_create = 30; % Increased to 30 to be safely longer than the 25-cycle sim
    tper_single = t_orig(end) - t_orig(1);
    
    % Create long vectors by repeating the single cycle
    t_long = [];
    Pao_long = [];
    Plv_long = [];
    SIP_long = [];
    Pra_long = [];

    for i = 0:num_cycles_to_create-1
        % Append all but the last point of the cycle, shifted in time
        t_long = [t_long; t_orig(1:end-1) + (i * tper_single)];
        Pao_long = [Pao_long; Pao_orig(1:end-1)];
        Plv_long = [Plv_long; Plv_orig(1:end-1)];
        SIP_long = [SIP_long; SIP_orig(1:end-1)];
        Pra_long = [Pra_long; Pra_orig(1:end-1)];
    end
    
    % Add the very last point to close the final cycle
    last_time_point = t_orig(end) + ((num_cycles_to_create-1) * tper_single);
    t_long = [t_long; last_time_point];
    Pao_long = [Pao_long; Pao_orig(end)];
    Plv_long = [Plv_long; Plv_orig(end)];
    SIP_long = [SIP_long; SIP_orig(end)];
    Pra_long = [Pra_long; Pra_orig(end)];

    % Overwrite original vectors with the new long ones
    t_orig = t_long;
    Pao_orig = Pao_long;
    Plv_orig = Plv_long;
    SIP_orig = SIP_long;
    Pra_orig = Pra_long;
end

% --- Common Calculations ---
tstart = t_orig(1);
tend = floor(t_orig(end));
tdata = linspace(tstart,tend,(tend-tstart)/dt);

Pao = interp1(t_orig, Pao_orig, tdata);
Plv = interp1(t_orig, Plv_orig, tdata);
SIP = interp1(t_orig, SIP_orig, tdata);
Pra = interp1(t_orig, Pra_orig, tdata);

% Scaling SIP based on Plv
SIP = SIP * (0.2 * max(Plv) / max(SIP)); 

tper = 60 / HR; % Cardiac Period (s)

dPlvdt = gradient(Plv, dt);
dSIPdt = gradient(SIP, dt);

% Combine all pressure data into a single matrix for ODE solver
Pdata = [tdata; Pao; Plv; dPlvdt; SIP; dSIPdt; Pra];

% Package all input data into a struct
InputData = struct('Pdata', Pdata, 'HR', HR, 'tdata', tdata, 'dt', dt, 'tper', tper, 'Plv', Plv, 'Pao', Pao, 'SIP', SIP, 'Pra', Pra);

end