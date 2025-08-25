function [varying_param_indices] = get_varying_parameters(filtered_inputs, num_params)
% GET_VARYING_PARAMETERS - Determines which parameters are varying in a dataset.

varying_param_indices = [];
for i = 1:num_params
    if length(unique(filtered_inputs(i, :))) > 1
        varying_param_indices(end+1) = i;
    end
end

end