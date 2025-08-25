function [filtered_inputs, filtered_svd_features] = apply_filter(SimIndex, all_sim_svd_features, current_filter, param_names)
% APPLY_FILTER - Applies a filter to the simulation data.

sim_indices_to_plot = true(1, size(SimIndex.inputs, 2));
filter_fields = fieldnames(current_filter);
for i = 1:length(filter_fields)
    field = filter_fields{i};
    param_idx = find(strcmp(param_names, field));
    if ~isempty(param_idx)
        values = current_filter.(field);
        sim_indices_to_plot = sim_indices_to_plot & ismember(SimIndex.inputs(param_idx, :), values);
    end
end

filtered_inputs = SimIndex.inputs(:, sim_indices_to_plot);
filtered_svd_features = all_sim_svd_features(:, sim_indices_to_plot);

end