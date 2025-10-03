%check the flow values for the training data and the rsults of the
%simulations 

% check the training data 
load('SVD_basis.mat','leftSV_raw')
load('matrixA_raw.mat')
matrixA_raw(:,end) = [];
figure(); 
subplot(1,2,1);
hold on; 
histogram(leftSV_raw(1,:))
xlabel('LSV')
subplot(1,2,2)
histogram(matrixA_raw(1,:))
xlabel('Qpa mean ml/min')
hold off;

figure();
scatter(leftSV_raw(1,:), matrixA_raw(1,:))


%% 
fourier_base_directory = "fourier_verification";
if ~exist(fourier_base_directory, 'dir')
    mkdir(fourier_base_directory)
end 
load('matrixA_raw.mat')
matrixA_raw(:,end) = [];
load("Sim_Fourier_Features.mat")

for i = 1:size(all_sim_fourier_features,1)

    sim_features = all_sim_fourier_features(i,:);
    pig_features = matrixA_raw(i,:);
    data = [sim_features, pig_features];
    group = [repmat({'sim'}, length(sim_features), 1); repmat({'pig'}, length(pig_features), 1)];

  
    figure();
    hold on; 
    boxplot(data, group);
    xlabel('Simulated Features');
    ylabel('Pig Features');
    title(['Feature Comparison for Index ' num2str(i)]);
    hold off;

end 

%% 

Q_sim_raw = results.X(:,2);
figure(); hold on;
plot(results.t,Q_sim_raw)