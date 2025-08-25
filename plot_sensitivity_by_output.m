function plot_sensitivity_by_output(SensResults, var, structType, fieldName, displayOption, multipliers, baselineMultipliers, base_case_params)
% PLOT_SENSITIVITY_BY_OUTPUT plots how a single output is affected by all inputs

    mainDir = fullfile(pwd, 'SensitivityResults');
    if ~exist(mainDir, 'dir')
        mkdir(mainDir)
    end
    
    nParams = numel(fieldnames(SensResults.(var)));
    max_multipliers = max(cellfun(@length, multipliers));

    valuesMatrix = NaN(max_multipliers, nParams);

    for paramIndex = 1:nParams
        index_str = sprintf("index%d", paramIndex);
        current_multipliers = multipliers{paramIndex};

        for i = 1:length(current_multipliers)
            m_val_multiplier = current_multipliers(i);
            mult_str = strrep(sprintf('mult_%.2f', m_val_multiplier), '.', 'p');
            if isfield(SensResults.(var).(index_str), mult_str) && ...
               isfield(SensResults.(var).(index_str).(mult_str), structType) && ...
               isfield(SensResults.(var).(index_str).(mult_str), structType) && ...
               isfield(SensResults.(var).(index_str).(mult_str).(structType), fieldName)

                valuesMatrix(i, paramIndex) = SensResults.(var).(index_str).(mult_str).(structType).(fieldName);
            end
        end
    end

    if displayOption == 0
        fig = figure('Position', [100, 100, 1200, 800]);
        ax = axes('Parent', fig);
    elseif displayOption == 1
        fig = figure('Visible', 'off', 'Position', [100, 100, 1200, 800]);
        ax = axes('Parent', fig);
    else
        disp('DisplayOption value is not 0 or 1')
        return
    end

    b = bar(ax, valuesMatrix', 'grouped');
    colors = parula(max_multipliers);

    for i = 1:length(b)
        b(i).FaceColor = colors(i, :);
    end
    xticklabels(arrayfun(@(x) sprintf('index %d', x), 1:nParams, 'UniformOutput', false));
    xlabel(ax, sprintf('%s Index', var));
    ylabel(ax, fieldName);
    
    legend_labels = arrayfun(@(x) sprintf('Multiplier %d', x), 1:max_multipliers, 'UniformOutput', false);
    legend(ax, legend_labels, 'Location', 'bestoutside');

    title(ax, ['Grouped Bar: Sensitivity of ', fieldName, ' across all m_myo indices (Base Case: ', num2str(base_case_params), ')']);

    for k = 1:length(b) % k iterates through the multiplier index (1 to max_multipliers)
        xtips = b(k).XEndPoints;
        
        % Labels for Multiplier (at half height)
        ytips_multiplier = b(k).YData / 2; % Half the height of the bar
        
        current_multiplier_labels = cell(1, nParams);
        for j = 1:nParams % j iterates through the parameter index (1 to nParams)
            if k <= length(multipliers{j}) % Check if this multiplier exists for this parameter
                current_multiplier_labels{j} = sprintf('%.2f', multipliers{j}(k));
            else
                current_multiplier_labels{j} = ''; % No label for NaN bars
            end
        end
        
        text(xtips, ytips_multiplier, current_multiplier_labels, ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', ...
             'Color', 'white', ...
             'FontSize', 8, ...
             'Rotation', 90); % Added rotation

        % Labels for Bar Value (at top)
        ytips_value = b(k).YEndPoints; % Top of the bar
        labels_value = string(round(b(k).YData, 2));
        
        text(xtips, ytips_value, labels_value, ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', ...
             'Color', 'black', ...
             'FontSize', 8);
    end

    if displayOption == 0
        % Create a subdirectory based on the fieldName
        saveDir = fullfile(mainDir, fieldName);
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
        end

        % Create a unique filename with a timestamp
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        baseSaveName = ['sensitivity_plot_', timestamp];
        
        % Save the plot
        plotSaveName = [baseSaveName, '.png'];
        ax.Toolbar.Visible = 'off';
        exportgraphics(fig, fullfile(saveDir, plotSaveName), 'Resolution', 300);
        close(fig);

        % Save the multiplier information to a text file
        txtSaveName = [baseSaveName, '_multipliers.txt'];
        fid = fopen(fullfile(saveDir, txtSaveName), 'w');
        if fid ~= -1
            fprintf(fid, 'Multipliers used for the sensitivity analysis plot generated at %s:\n\n', timestamp);
            for i = 1:length(multipliers)
                fprintf(fid, 'Parameter %d: %s\n', i, mat2str(multipliers{i}));
            end
            fclose(fid);
        end

        fprintf('Saved sensitivity analysis plot and multiplier info for %s in %s\n', fieldName, saveDir);
    end

end 

