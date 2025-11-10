%% plot_all_spmv_results.m
% Processes all CSV files in /results directory
% Creates separate execution time and speedup plots for each matrix type
clear; close all; clc;

%% Configuration
resultsDir = './results'; % Change to your results directory path
outputDir = './DELIVERABLE1_IMAGES'; % Folder to save figures

% Create output folder if it doesn't exist
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Find all CSV files
csvFiles = dir(fullfile(resultsDir, '*.csv'));
nFiles = numel(csvFiles);

if nFiles == 0
    error('No CSV files found in %s', resultsDir);
end

fprintf('Found %d CSV file(s) in %s\n\n', nFiles, resultsDir);

%% Process each CSV file separately
for f = 1:nFiles
    filename = fullfile(csvFiles(f).folder, csvFiles(f).name);
    [~, matrixName, ~] = fileparts(csvFiles(f).name);
    
    fprintf('=== Processing: %s ===\n', csvFiles(f).name);
    
    try
        data = readtable(filename);
        data.Executable = categorical(data.Executable);
        fprintf('  Loaded %d rows\n', height(data));
        
        %% Compute 90th percentile statistics
        executables = categories(data.Executable);
        nExec = numel(executables);
        
        p90Times = struct();
        stdTimes = struct();
        threadCounts = struct();
        
        for i = 1:nExec
            exec = executables{i};
            subset = data(data.Executable == exec, :);
            
            threads = unique(subset.Threads);
            pTimes = zeros(size(threads));
            sTimes = zeros(size(threads));
            
            for j = 1:numel(threads)
                t = threads(j);
                runs = subset.Time(subset.Threads == t);
                pTimes(j) = prctile(runs, 90); % 90th percentile
                sTimes(j) = std(runs);
            end
            
            p90Times.(exec) = pTimes;
            stdTimes.(exec) = sTimes;
            threadCounts.(exec) = threads;
        end
        
        %% Plot settings
        colors = lines(nExec);
        markers = {'o', 's', 'd', '^', 'v', 'p', 'h', '*', 'x', '+'};
        
        %% ===== EXECUTION TIME PLOT =====
        figure('Position', [100, 100, 900, 600]);
        hold on; grid on; box on;
        
        allThreads = unique(data.Threads);
        xPositions = 1:numel(allThreads);
        
        for i = 1:nExec
            exec = executables{i};
            threads = threadCounts.(exec);
            pTimes = p90Times.(exec);
            sTimes = stdTimes.(exec);
            [~, idx] = ismember(threads, allThreads);
            
            errorbar(xPositions(idx), pTimes, sTimes, ...
                'LineWidth', 1.5, ...
                'Marker', markers{mod(i-1, numel(markers))+1}, ...
                'Color', colors(i,:), ...
                'MarkerFaceColor', colors(i,:), ...
                'MarkerSize', 8);
        end
        
        xlabel('Number of Threads', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel('Execution Time (s)', 'Interpreter', 'latex', 'FontSize', 12);
        title(sprintf('SpMV Performance - %s (90th Percentile)', strrep(matrixName, '_', '\_')), ...
            'Interpreter', 'tex', 'FontSize', 13);
        legend(strrep(executables, 'spmv_', ''), 'Location', 'northwest', 'Interpreter', 'latex');
        set(gca, 'XTick', xPositions, 'XTickLabel', string(allThreads), 'FontSize', 11, 'YScale', 'log');
        axis tight;
        grid on;
        
        %% Smart label positioning for execution time
        allLabels = [];
        for i = 1:nExec
            exec = executables{i};
            threads = threadCounts.(exec);
            pTimes = p90Times.(exec);
            [~, idx] = ismember(threads, allThreads);
            
            for j = 1:numel(pTimes)
                allLabels = [allLabels; struct('x', xPositions(idx(j)), 'y', pTimes(j), 'exec', i)];
            end
        end
        
        [~, sortIdx] = sortrows([[allLabels.x]', [allLabels.y]'], [1, 2]);
        allLabels = allLabels(sortIdx);
        
        for k = 1:length(allLabels)
            lbl = allLabels(k);
            pTime = lbl.y;
            
            if pTime >= 0.01
                labelText = sprintf('%.3f', pTime);
            elseif pTime >= 0.001
                labelText = sprintf('%.4f', pTime);
            else
                labelText = sprintf('%.2e', pTime);
            end
            
            sameX = [allLabels.x] == lbl.x;
            yValues = [allLabels(sameX).y];
            if length(yValues) > 1
                yRank = find(sort(yValues) == pTime, 1);
                if mod(yRank, 2) == 0
                    hAlign = 'left'; vAlign = 'middle'; xOffset = 0.15;
                else
                    hAlign = 'center'; vAlign = 'bottom'; xOffset = 0;
                end
            else
                hAlign = 'center'; vAlign = 'bottom'; xOffset = 0;
            end
            
            text(lbl.x + xOffset, pTime, labelText, ...
                'VerticalAlignment', vAlign, 'HorizontalAlignment', hAlign, ...
                'FontSize', 7, 'Color', colors(lbl.exec,:), 'FontWeight', 'bold', ...
                'BackgroundColor', [1 1 1 0.7], 'EdgeColor', 'none', 'Margin', 1);
        end
        
        % Save execution time figure
        execTimeFilename = fullfile(outputDir, sprintf('%s_exectime.png', matrixName));
        exportgraphics(gcf, execTimeFilename, 'Resolution', 300);
        fprintf('  ✅ Saved: %s\n', execTimeFilename);
        
        %% ===== SPEEDUP PLOT =====
        seqExec = executables(contains(string(executables), 'sequential', 'IgnoreCase', true));
        
        if ~isempty(seqExec)
            seqExec = seqExec{1};
            seqTime = mean(p90Times.(seqExec));
            fprintf('  Using %s as baseline (90th percentile: %.6f s)\n', seqExec, seqTime);
            
            figure('Position', [150, 150, 900, 600]);
            hold on; grid on; box on;
            
            for i = 1:nExec
                exec = executables{i};
                if strcmp(exec, seqExec), continue; end
                
                threads = threadCounts.(exec);
                pTimes = p90Times.(exec);
                speedup = seqTime ./ pTimes;
                [~, idx] = ismember(threads, allThreads);
                
                plot(xPositions(idx), speedup, ...
                    'LineWidth', 1.5, ...
                    'Marker', markers{mod(i-1, numel(markers))+1}, ...
                    'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 8);
            end
            
            % Ideal linear speedup line
            if any(allThreads > 1)
                parallelThreads = allThreads(allThreads > 1);
                [~, idx] = ismember(parallelThreads, allThreads);
                plot(xPositions(idx), parallelThreads, '--k', 'LineWidth', 1.2, 'DisplayName', 'Ideal Linear');
            end
            
            xlabel('Number of Threads', 'Interpreter', 'latex', 'FontSize', 12);
            ylabel('Speedup (×)', 'Interpreter', 'latex', 'FontSize', 12);
            title(sprintf('SpMV Speedup vs Sequential - %s', strrep(matrixName, '_', '\_')), ...
                'Interpreter', 'tex', 'FontSize', 13);
            
            legendEntries = executables(~strcmp(executables, seqExec));
            legend([strrep(legendEntries, 'spmv_', ''); {'Ideal Linear'}], ...
                'Location', 'northwest', 'Interpreter', 'latex');
            set(gca, 'XTick', xPositions, 'XTickLabel', string(allThreads), 'FontSize', 11);
            axis tight; grid on;
            
            %% Smart label positioning for speedup
            allSpeedupLabels = [];
            for i = 1:nExec
                exec = executables{i};
                if strcmp(exec, seqExec), continue; end
                
                threads = threadCounts.(exec);
                pTimes = p90Times.(exec);
                speedup = seqTime ./ pTimes;
                [~, idx] = ismember(threads, allThreads);
                
                for j = 1:numel(speedup)
                    allSpeedupLabels = [allSpeedupLabels; struct('x', xPositions(idx(j)), 'y', speedup(j), 'exec', i)];
                end
            end
            
            [~, sortIdx] = sortrows([[allSpeedupLabels.x]', [allSpeedupLabels.y]'], [1, 2]);
            allSpeedupLabels = allSpeedupLabels(sortIdx);
            
            for k = 1:length(allSpeedupLabels)
                lbl = allSpeedupLabels(k);
                labelText = sprintf('%.2fx', lbl.y);
                
                sameX = [allSpeedupLabels.x] == lbl.x;
                yValues = [allSpeedupLabels(sameX).y];
                if length(yValues) > 1
                    yRank = find(sort(yValues) == lbl.y, 1);
                    if mod(yRank, 2) == 0
                        hAlign = 'left'; vAlign = 'middle'; xOffset = 0.15;
                    else
                        hAlign = 'center'; vAlign = 'bottom'; xOffset = 0;
                    end
                else
                    hAlign = 'center'; vAlign = 'bottom'; xOffset = 0;
                end
                
                text(lbl.x + xOffset, lbl.y, labelText, ...
                    'VerticalAlignment', vAlign, 'HorizontalAlignment', hAlign, ...
                    'FontSize', 7, 'Color', colors(lbl.exec,:), 'FontWeight', 'bold', ...
                    'BackgroundColor', [1 1 1 0.7], 'EdgeColor', 'none', 'Margin', 1);
            end
            
            % Save speedup figure
            speedupFilename = fullfile(outputDir, sprintf('%s_speedup.png', matrixName));
            exportgraphics(gcf, speedupFilename, 'Resolution', 300);
            fprintf('  Saved: %s\n', speedupFilename);
        else
            warning('  No sequential implementation found for %s', matrixName);
        end
        
        fprintf('\n');
        
    catch ME
        warning('Could not process %s: %s', csvFiles(f).name, ME.message);
        fprintf('\n');
    end
end

fprintf('All processing complete!\n');
