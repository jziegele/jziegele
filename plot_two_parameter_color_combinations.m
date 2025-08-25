function plot_two_parameter_color_combinations(filtered_inputs, filtered_svd_features, varying_param_indices, param_names, svd_component_names, filter_dir_path)
% PLOT_TWO_PARAMETER_COLOR_COMBINATIONS - Plots the interaction between two varying parameters with unique colors for each combination.

param1_idx = varying_param_indices(1);
param2_idx = varying_param_indices(2);
param1_name = param_names{param1_idx};
param2_name = param_names{param2_idx};

% Get unique values of the second parameter for coloring
param2_values = unique(filtered_inputs(param2_idx, :));
num_colors = length(param2_values);
colors = lines(num_colors);

num_svd_components = size(filtered_svd_features, 1);

for comp_idx = 1:num_svd_components
    current_component_name = svd_component_names{comp_idx};
    hFig = figure('Visible', 'off');
    hold on;
    
    for i = 1:length(param2_values)
        val = param2_values(i);
        indices = filtered_inputs(param2_idx, :) == val;
        
        scatter(filtered_inputs(param1_idx, indices), filtered_svd_features(comp_idx, indices), [], colors(i,:), 'filled');
    end
    
    hold off;
    
    xlabel(param1_name);
    ylabel(current_component_name);
    title(sprintf('Interaction of %s and %s on %s', param1_name, param2_name, current_component_name));
    
    % Create legend
    legend_entries = arrayfun(@(x) sprintf('%s = %.2f', param2_name, x), param2_values, 'UniformOutput', false);
    legend(legend_entries);
    grid on;
    
    plot_filename = fullfile(filter_dir_path, sprintf('%s_vs_%s_and_%s_combinations.png', current_component_name, param1_name, param2_name));
    saveas(hFig, plot_filename);
    close(hFig);
end

end