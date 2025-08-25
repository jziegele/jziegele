function plot_with_color_filter(filtered_inputs, filtered_svd_features, p_idx, color_filter, param_names, svd_component_names, filter_dir_path, connect_lines)
% PLOT_WITH_COLOR_FILTER - Plots sensitivity with a color filter, with an option to connect the points with lines.

current_param_name = param_names{p_idx};

num_svd_components = size(filtered_svd_features, 1);

% --- Get the parameter and values from the color_filter ---
color_filter_fields = fieldnames(color_filter);
if isempty(color_filter_fields)
    warning('Color filter is empty. No colored points will be plotted.');
    return;
end
% Assuming only one field is used for color filtering
color_param_name = color_filter_fields{1}; 
color_param_idx = find(strcmp(param_names, color_param_name));
if isempty(color_param_idx)
    warning('Color filter parameter ''%s'' not found.', color_param_name);
    return;
end
color_values = color_filter.(color_param_name);
colors = lines(length(color_values)); % Generate a set of colors

for comp_idx = 1:num_svd_components
    current_component_name = svd_component_names{comp_idx};
    x_values = filtered_inputs(p_idx, :);
    y_values = filtered_svd_features(comp_idx, :);
    
    hFig = figure('Visible', 'off');
    hold on;
    
    % Plot all points in gray first as a background
    scatter(x_values, y_values, 'filled', 'MarkerFaceColor', [0.8 0.8 0.8], 'DisplayName', 'All Data');
    
    % --- Loop through each multiplier value and plot with a unique color ---
    for i = 1:length(color_values)
        val = color_values(i);
        
        % Find indices for the current multiplier value
        color_indices = (filtered_inputs(color_param_idx, :) == val);
        
        if any(color_indices)
            % Get the x and y data for the current color group
            group_x = x_values(color_indices);
            group_y = y_values(color_indices);
            
            legend_text = sprintf('%s = %.2f', color_param_name, val);

            if connect_lines
                % Sort the data by x-value to connect them correctly
                [sorted_x, sort_order] = sort(group_x);
                sorted_y = group_y(sort_order);
                % Plot points and connecting line in one command
                plot(sorted_x, sorted_y, '-o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', legend_text);
            else
                % Plot only the scatter points
                scatter(group_x, group_y, 'filled', 'MarkerFaceColor', colors(i,:), 'DisplayName', legend_text);
            end
        end
    end
    
    hold off;
    
    xlabel_str = sprintf('%s Multiplier', current_param_name);
    ylabel_str = sprintf('SVD Component Coefficient: %s', current_component_name);
    title_str = sprintf('Sensitivity of %s to %s (Colored by %s)', current_component_name, current_param_name, color_param_name);
    
    xlabel(xlabel_str);
    ylabel(ylabel_str);
    title(title_str);
    legend('Location', 'best');
    grid on;
    
    plot_filename = fullfile(filter_dir_path, sprintf('%s_vs_%s_color_filtered_by_%s.png', current_component_name, current_param_name, color_param_name));
    saveas(hFig, plot_filename);
    close(hFig);
end

end