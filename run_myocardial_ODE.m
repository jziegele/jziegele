function [t, X, o] = run_myocardial_ODE(inputData, m_myo, ExerciseLevel)
% Purpose 
    % Runs ODE model based on input parameters and outputs 
% inputs 
    % Exerccise level  - manually selected in driver
% outputs
    % InputData  - signal characteristics used in ODE modeling 

switch ExerciseLevel
    case 0
        f=[1,1,1,1];
    case 1
        f=[2.5,3.3,5.0,2.5];
    case 2
        f=[3.7,5.1, 8.7,3.7];
    case 3
        f=[4.9,7.0, 13.0,4.9];
    case 4
        f=[6.3,9.2,19.2,6.3];
end

Xo_myo = [60 1 50 50 85 85 120 120 5]';
opts=odeset('MaxStep',1/200, 'RelTol', 1e-5);

%End time should be set such that it is at the end of a cardiac cycle
%Set to either 30 or 60 seconds or a set number of cardiac cycles 

% Determine the total duration of available input data
total_data_duration = inputData.Pdata(1, end);
max_cardiac_cycles = floor(total_data_duration / inputData.tper)

% Determine the duration for the first ODE run
if total_data_duration >= 60
    first_run_duration = 60;
elseif total_data_duration >= 30
    first_run_duration = 30;
else
    first_run_duration = max_cardiac_cycles * inputData.tper;
end

% Determine the duration for the second ODE run
second_run_duration = first_run_duration;

%Run the ode once to reach steady state 
%and a second time for input data to the oxygen transport model
[~,X_temp] = ode15s(@dXdT_myo,[0.0, first_run_duration],Xo_myo,opts,inputData.Pdata, m_myo,0);
[t,X] = ode15s(@dXdT_myo,[0.0, second_run_duration],X_temp(end,:),opts,inputData.Pdata, m_myo,0);

[~, o] = dXdT_myo(t, X, inputData.Pdata, m_myo, 1);

end