
% run_single_simulation.m
%
% This script runs a single simulation of the myocardial perfusion model
% with default parameters and saves the results to simulationTEST.mat.

clear all; close all; clc;

% Add the path to the project and its subfolders
addpath(genpath(pwd));

% Load the baseline model parameters
m_myo = ones(1, 10);

simulation_input = 0; % ExerciseLevel
extend_patient_data = false;
% Load the input data
inputData = Load_input_data(simulation_input, extend_patient_data);

% Run the simulation
[t, X, o] = run_myocardial_ODE(inputData, m_myo, 0);

% Save the results
save('simulationTEST.mat', 't', 'X', 'o', 'm_myo');

disp('Single simulation completed and results saved to simulationTEST.mat');

%% 
Q_sim_raw = X(:,2)*60;
figure(); hold on;
plot(t,Q_sim_raw)
xlim([0, 5])
