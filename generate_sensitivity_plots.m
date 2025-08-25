function generate_sensitivity_plots(filter_config, color_filter_config, output_dir, SimIndex, all_sim_svd_features, connect_lines)
% generate_sensitivity_plots - Generates and saves sensitivity plots based on filter configurations.

% --- Define Parameters ---
param_names = {'m1', 'm2', 'm3', 'm4', 'm5', 'm6'};
num_params = length(param_names);

% --- Define SVD Component Names ---
num_svd_components = size(all_sim_svd_features, 1);
svd_component_names = cell(num_svd_components, 1);
for k = 1:num_svd_components
    svd_component_names{k} = sprintf('SVD_Component_%d', k);
end

% --- Apply Main Filter ---
[filtered_inputs, filtered_svd_features] = apply_filter(SimIndex, all_sim_svd_features, filter_config, param_names);

% --- Determine which parameter is on the x-axis ---
varying_param_indices = get_varying_parameters(filtered_inputs, num_params);
% We expect two varying parameters, one for x-axis and one for color.
% The one NOT in the color filter is the x-axis parameter.
color_filter_param_name = fieldnames(color_filter_config);
x_axis_param_idx = -1;
for idx = varying_param_indices
    if ~strcmp(param_names{idx}, color_filter_param_name{1})
        x_axis_param_idx = idx;
        break;
    end
end

if x_axis_param_idx == -1
    warning('Could not determine x-axis parameter. Skipping plot generation for this run.');
    return;
end

% --- Call the plotting function ---
% The plot_with_color_filter function is now the main plotting engine.
plot_with_color_filter(filtered_inputs, filtered_svd_features, x_axis_param_idx, color_filter_config, param_names, svd_component_names, output_dir, connect_lines);

end