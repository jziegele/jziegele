function plot_single_parameter_sensitivity(filtered_inputs, filtered_svd_features, varying_param_indices, param_names, svd_component_names, filter_dir_path)
% PLOT_SINGLE_PARAMETER_SENSITIVITY - Plots sensitivity for a single varying parameter.

varying_param_idx = varying_param_indices(1);
varying_param_name = param_names{varying_param_idx};
unique_values = unique(filtered_inputs(varying_param_idx, :));
colors = lines(length(unique_values));

num_svd_components = size(filtered_svd_features, 1);

for comp_idx = 1:num_svd_components
    current_component_name = svd_component_names{comp_idx};
    hFig = figure('Visible', 'off');
    hold on;
    legend_entries = {};
    for i = 1:length(unique_values)
        val = unique_values(i);
        indices = filtered_inputs(varying_param_idx, :) == val;
        scatter(filtered_inputs(varying_param_idx, indices), filtered_svd_features(comp_idx, indices), [], colors(i,:), 'filled');
        legend_entries{end+1} = sprintf('%s = %.2f', varying_param_name, val);
    end
    hold off;
    
    xlabel(varying_param_name);
    ylabel(current_component_name);
    title(sprintf('Sensitivity of %s to %s', current_component_name, varying_param_name));
    legend(legend_entries);
    grid on;
    
    plot_filename = fullfile(filter_dir_path, sprintf('%s_vs_%s.png', current_component_name, varying_param_name));
    saveas(hFig, plot_filename);
    close(hFig);
end

end