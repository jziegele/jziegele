load('Sim_SVD_Features.mat');
fid = fopen('nan_check_result.txt', 'w');
if any(isnan(all_sim_svd_features(:)))
    fprintf(fid, 'File contains NaNs');
else
    fprintf(fid, 'File does not contain NaNs');
end
fclose(fid);