clc; clear; close all;

%% Directories
dataDir = './DELIVERABLE1_RES';
outputDir = './DELIVERABLE1_IMAGES';

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

files = dir(fullfile(dataDir, '*_chunked.csv'));

if isempty(files)
    error('No *_chunked.csv files found in %s', dataDir);
end

%% Process each CSV independently
for k = 1:length(files)
    filename = fullfile(dataDir, files(k).name);
    [~, baseName, ~] = fileparts(files(k).name);
    matrixName = erase(baseName, '_chunked');
    fprintf('\n📄 Processing: %s\n', files(k).name);

    try
        T = readtable(filename);
    catch ME
        warning('⚠️ Could not read %s: %s', files(k).name, ME.message);
        continue;
    end

    % Ensure types are consistent
    if iscell(T.Executable)
        T.Executable = string(T.Executable);
    end

    executables = unique(string(T.Executable));
    threads = unique(T.Threads);
    chunkSizes = unique(T.ChunkSize);

    %% Compute 90th percentile times
    summaryCells = {};
    for e = 1:length(executables)
        for c = 1:length(chunkSizes)
            for t = 1:length(threads)
                filt = (string(T.Executable) == executables(e)) & ...
                       (T.Threads == threads(t)) & ...
                       (T.ChunkSize == chunkSizes(c));

                times = T.Time(filt);
                if ~isempty(times)
                    summaryCells = [summaryCells; {executables(e), threads(t), chunkSizes(c), prctile(times, 90)}];
                end
            end
        end
    end

    if isempty(summaryCells)
        warning('No valid data for %s', files(k).name);
        continue;
    end

    summaryData = cell2table(summaryCells, ...
        'VariableNames', {'Executable','Threads','ChunkSize','Time90'});

    %% ===== LINE PLOT =====
    fig = figure('Name', sprintf('Line Plot - %s', matrixName), ...
        'Color', 'w', 'Position', [100, 100, 900, 600]);
    set(fig, 'Renderer', 'opengl');
    hold on; grid on; box on;

    colors = lines(length(chunkSizes));
    markers = {'o', 's', 'd', '^', 'v', 'p', 'h'};

    xPositions = 1:length(threads);

    for e = 1:length(executables)
        for c = 1:length(chunkSizes)
            dataC = summaryData(summaryData.ChunkSize == chunkSizes(c) & ...
                                summaryData.Executable == executables(e), :);
            if isempty(dataC)
                continue;
            end

            [~, idx] = ismember(dataC.Threads, threads);
            
            % Create clean display name for legend
            execStr = char(executables(e));
            chunkStr = sprintf('Chunk %d', chunkSizes(c));
            displayName = sprintf('%s - %s', execStr, chunkStr);
            
            plot(xPositions(idx), dataC.Time90, '-', ...
                'Color', colors(c,:), ...
                'DisplayName', displayName, ...
                'LineWidth', 1.5, ...
                'Marker', markers{mod(e-1, length(markers))+1}, ...
                'MarkerSize', 8, ...
                'MarkerFaceColor', colors(c,:));

            % Add smart value labels
            for j = 1:height(dataC)
                yVal = dataC.Time90(j);
                if yVal >= 0.01
                    labelText = sprintf('%.3f', yVal);
                elseif yVal >= 0.001
                    labelText = sprintf('%.4f', yVal);
                else
                    labelText = sprintf('%.2e', yVal);
                end

                text(xPositions(idx(j)), yVal, labelText, ...
                    'VerticalAlignment', 'bottom', ...
                    'HorizontalAlignment', 'center', ...
                    'FontSize', 7, 'Color', colors(c,:), ...
                    'FontWeight', 'bold', ...
                    'BackgroundColor', [1 1 1 0.7], ...
                    'Margin', 1);
            end
        end
    end

    % Clean matrix name for title (replace underscores properly)
    cleanMatrixName = strrep(matrixName, '_', ' ');
    
    xlabel('Number of Threads', 'FontSize', 12);
    ylabel('Execution Time (s)', 'FontSize', 12);
    title(sprintf('90th Percentile Time vs Threads - %s', cleanMatrixName), ...
        'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    set(gca, 'XTick', xPositions, 'XTickLabel', string(threads), ...
        'FontSize', 11, 'YScale', 'log');
    axis tight;
    hold off;

    % Save line plot
    outFile = fullfile(outputDir, sprintf('LinePlot_%s.png', matrixName));
    exportgraphics(fig, outFile, 'Resolution', 300);
    fprintf('✅ Saved: %s\n', outFile);
    close(fig);

    %% ===== HEATMAP =====
    for e = 1:length(executables)
        heatData = NaN(length(chunkSizes), length(threads));
        for c = 1:length(chunkSizes)
            for t = 1:length(threads)
                filt = (summaryData.Executable == executables(e)) & ...
                       (summaryData.ChunkSize == chunkSizes(c)) & ...
                       (summaryData.Threads == threads(t));
                if any(filt)
                    heatData(c,t) = summaryData.Time90(filt);
                end
            end
        end

        if any(~isnan(heatData(:)))
            figH = figure('Name', sprintf('Heatmap - %s (%s)', matrixName, executables(e)), ...
                'Color','w', 'Position', [150, 150, 900, 700]);
            set(figH, 'Renderer', 'opengl');

            % Clean executable name and matrix name for title
            execStr = char(executables(e));
            cleanMatrixName = strrep(matrixName, '_', ' ');
            
            % Create modern styled heatmap with cell borders
            hold on;
            
            % Plot each cell as a rectangle
            for c = 1:length(chunkSizes)
                for t = 1:length(threads)
                    val = heatData(c,t);
                    if ~isnan(val)
                        % Normalize value for color mapping
                        minVal = min(heatData(:), [], 'omitnan');
                        maxVal = max(heatData(:), [], 'omitnan');
                        normVal = (val - minVal) / (maxVal - minVal);
                        
                        % Get color from parula colormap
                        cmap = parula(256);
                        colorIdx = max(1, min(256, round(normVal * 255) + 1));
                        cellColor = cmap(colorIdx, :);
                        
                        % Draw rectangle (full size, no gaps)
                        rectangle('Position', [t-0.5, c-0.5, 1, 1], ...
                            'FaceColor', cellColor, ...
                            'EdgeColor', [0.3 0.3 0.3], ...
                            'LineWidth', 1.5);
                        
                        % Add value label inside cell
                        if val >= 0.01
                            lbl = sprintf('%.3f', val);
                        elseif val >= 0.001
                            lbl = sprintf('%.4f', val);
                        else
                            lbl = sprintf('%.2e', val);
                        end
                        
                        text(t, c, lbl, ...
                            'HorizontalAlignment', 'center', ...
                            'VerticalAlignment', 'middle', ...
                            'FontSize', 10, 'Color', 'w', ...
                            'FontWeight', 'bold');
                    end
                end
            end
            hold off;
            
            % Set axis limits and labels
            xlim([0.5, length(threads)+0.5]);
            ylim([0.5, length(chunkSizes)+0.5]);
            
            xlabel('Number of Threads', 'FontSize', 13, 'FontWeight', 'bold');
            ylabel('Chunk Size', 'FontSize', 13, 'FontWeight', 'bold');
            title(sprintf('90th Percentile Time Heatmap - %s (%s)', ...
                cleanMatrixName, execStr), ...
                'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
            
            % Format axis
            set(gca, 'XTick', 1:length(threads), 'XTickLabel', threads, ...
                'YTick', 1:length(chunkSizes), 'YTickLabel', chunkSizes, ...
                'FontSize', 11, 'Box', 'on', 'LineWidth', 1.5);
            
            % Add colorbar
            colormap('parula');
            caxis([min(heatData(:), [], 'omitnan'), max(heatData(:), [], 'omitnan')]);
            cb = colorbar;
            cb.Label.String = 'Execution Time (s)';
            cb.Label.FontSize = 12;
            cb.Label.FontWeight = 'bold';
            cb.FontSize = 11;
            
            % Save heatmap
            heatFile = fullfile(outputDir, sprintf('Heatmap_%s_%s.png', matrixName, execStr));
            exportgraphics(figH, heatFile, 'Resolution', 300);
            fprintf('✅ Saved: %s\n', heatFile);
            close(figH);
        end
    end
end

fprintf('\n✅ All line plots and heatmaps generated successfully!\n');