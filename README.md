# Myocardial Perfusion Model Sensitivity Analysis

## Description

This project performs a comprehensive sensitivity analysis for a lumped-parameter model of myocardial perfusion. It allows for the systematic variation of model parameters, execution of simulations, and analysis of the results to understand the influence of each parameter on various physiological outputs.

The core of the analysis involves:
1.  Running a series of simulations over a multi-dimensional parameter space.
2.  Calculating key physiological metrics from the simulation results (e.g., mean flows, pressures, volumes, and endo-epicardial ratios).
3.  Utilizing Singular Value Decomposition (SVD) to extract dominant features from the time-series data.
4.  Generating a variety of plots to visualize the sensitivity of the outputs to the input parameters.

## Features

-   **Automated Simulation Runner:** The `main.m` script automates the process of running simulations for all combinations of specified parameter multipliers.
-   **Robust Simulation Core:** Uses the `ode15s` solver for stiff ODEs to simulate the myocardial perfusion model.
-   **Comprehensive Output Calculation:** The `compute_outputs.m` script calculates a wide range of physiological metrics from the raw simulation data.
-   **Advanced Feature Extraction:** The `calculate_svd_features.m` script uses SVD to reduce the dimensionality of the simulation results and extract key features.
-   **Systematic Plotting:** The `main_plotting.m` script provides a framework for systematically generating a large number of plots to explore the data.
-   **Variety of Plotting Functions:** Includes functions to generate various types of plots, such as:
    -   Time-series plots of pressures, flows, and volumes.
    -   Summary metric plots.
    -   Sensitivity plots.
    -   Violin plots for visualizing distributions.

## Getting Started

### Prerequisites

-   MATLAB (version R2018a or newer recommended).
-   No special toolboxes are required for the core functionality, but some plotting functions might have specific dependencies (e.g., `violinplot`). Please ensure you have the necessary files for all plotting functions.

### Installation

1.  Clone or download this repository to your local machine.
2.  Open MATLAB and navigate to the project directory.
3.  Add the project directory and its subfolders to the MATLAB path:

```matlab
addpath(genpath(pwd));
```

## Usage

The analysis is typically performed in three steps:

### 1. Running the Simulations (`main.m`)

This script is the primary driver for the sensitivity analysis. It runs the simulations for all parameter combinations and saves the results. It also contains commented-out examples for plotting the direct results of the simulations (e.g., pressure and flow waveforms).

To run the simulations, execute the `main.m` script in the MATLAB command window:

```matlab
main
```

### 2. Performing the SVD Projection (`analyze_simulation_features.m`)

After running the simulations, this script should be run to perform the SVD projection on the simulation results. It calculates the SVD features for all simulations and saves them to `Sim_SVD_Features.mat`.

To run the SVD projection analysis, execute the `analyze_simulation_features.m` script:

```matlab
analyze_simulation_features
```

### 3. Plotting the SVD Projection Analysis (`main_plotting.m`)

This script is used to generate plots specifically for the SVD projection analysis. It loads the SVD features and creates a variety of plots to visualize the results of the sensitivity analysis on the SVD components.

To generate the SVD analysis plots, execute the `main_plotting.m` script:

```matlab
main_plotting
```

## File Structure

Here is a brief overview of the project's file structure:

```
SVDProjection/
|-- main.m                      # Main script to run the sensitivity analysis
|-- main_plotting.m             # Main script to generate plots
|-- compute_outputs.m           # Function to calculate physiological outputs
|-- run_myocardial_ODE.m        # Function to run the ODE simulation
|-- dXdT_myo.m                  # The core ODE model of myocardial perfusion
|-- calculate_svd_features.m    # Function to calculate SVD features
|-- CalcFourierModes.m          # Function to calculate Fourier modes
|-- plot_...                    # Various plotting functions
|-- SimIndex.mat                # Master index of all simulations
|-- SVD_Basis.mat               # SVD basis vectors
|-- AllResults.mat              # All simulation results
|-- simulation_results/         # Directory containing the raw results of each simulation
|-- plots_svd_sensitivity/      # Directory for sensitivity plots
|-- violin_plots/               # Directory for violin plots
|-- README.md                   # This file
```

## Contributing

Contributions to this project are welcome. If you find any issues or have suggestions for improvements, please feel free to open an issue or submit a pull request.

For questions or collaborations, please contact the project owner.
