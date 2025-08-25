function plot_default_sensitivity(filtered_inputs, filtered_svd_features, varying_param_indices, param_names, svd_component_names, filter_dir_path)
% PLOT_DEFAULT_SENSITIVITY - Plots sensitivity for more than two varying parameters.

num_svd_components = size(filtered_svd_features, 1);

for p_idx = varying_param_indices
    current_param_name = param_names{p_idx};
    param_dir = fullfile(filter_dir_path, current_param_name);
    if ~exist(param_dir, 'dir')
        mkdir(param_dir);
    end
    
    for comp_idx = 1:num_svd_components
        current_component_name = svd_component_names{comp_idx};
        x_values = filtered_inputs(p_idx, :);
        y_values = filtered_svd_features(comp_idx, :);
        
        hFig = figure('Visible', 'off');
        scatter(x_values, y_values, 'filled');
        
        xlabel_str = sprintf('%s Multiplier', current_param_name);
        ylabel_str = sprintf('SVD Component Coefficient: %s', current_component_name);
        title_str = sprintf('Sensitivity of %s to %s', current_component_name, current_param_name);
        
        xlabel(xlabel_str);
        ylabel(ylabel_str);
        title(title_str);
        grid on;
        
        plot_filename = fullfile(param_dir, sprintf('%s_vs_%s.png', current_component_name, current_param_name));
        saveas(hFig, plot_filename);
        close(hFig);
    end
end

end