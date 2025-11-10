clc; clear; close all;

%% Directory and Files
dataDir = './results';
files = dir(fullfile(dataDir, '*_chunked.csv'));

if isempty(files)
    error('No CSV files found in the specified directory');
end

allData = table();

% Read all CSVs
for k = 1:length(files)
    filename = fullfile(dataDir, files(k).name);
    try
        T = readtable(filename);
        % Add consistent matrix name
        T.MatrixName = repmat({erase(files(k).name,'_chunked.csv')}, height(T), 1);
        
        % Convert Executable to string if needed
        if iscell(T.Executable)
            T.Executable = string(T.Executable);
        end
        
        allData = [allData; T];
    catch ME
        warning('Failed to read file %s: %s', filename, ME.message);
    end
end

if isempty(allData)
    error('No data was successfully loaded');
end

%% Unique identifiers
executables = unique(string(allData.Executable));
threads = unique(allData.Threads);          
chunkSizes = unique(allData.ChunkSize);    

numCombinations = length(executables) * length(threads) * length(chunkSizes);
summaryCells = cell(numCombinations, 4);
idx = 1;

%% Compute 90th percentile times
for e = 1:length(executables)
    for c = 1:length(chunkSizes)
        for t = 1:length(threads)
            filt = (string(allData.Executable) == executables(e)) & ...
                   (allData.Threads == threads(t)) & ...
                   (allData.ChunkSize == chunkSizes(c));
            
            times = allData.Time(filt);
            if ~isempty(times)
                summaryCells(idx,:) = {executables(e), threads(t), chunkSizes(c), prctile(times,90)};
                idx = idx + 1;
            end
        end
    end
end

if idx == 1
    error('No valid combinations found');
end

summaryCells = summaryCells(1:idx-1,:);
summaryData = cell2table(summaryCells, ...
    'VariableNames', {'Executable','Threads','ChunkSize','Time90'});

%% Line Plots
figure('Name','Line Plots of 90th Percentile Time','Color','w'); 
set(gcf, 'Renderer', 'opengl');
set(gca, 'FontName', 'Arial', 'FontSize', 12);
hold on;

colors = lines(length(chunkSizes));
for c = 1:length(chunkSizes)
    dataC = summaryData(summaryData.ChunkSize == chunkSizes(c), :);
    for e = 1:length(executables)
        dataE = dataC(dataC.Executable == executables(e), :);
        if ~isempty(dataE)
            plot(dataE.Threads, dataE.Time90, '-o', ...
                'Color', colors(c,:), ...
                'DisplayName', sprintf('%s - Chunk %d', executables(e), chunkSizes(c)), ...
                'LineWidth', 2, 'MarkerSize', 8);
        end
    end
end

xlabel('Threads'); ylabel('Time (s)');
title('90th Percentile Time vs Threads');
legend('Location','best'); grid on; hold off;

%% Heatmaps
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
        figure('Name', sprintf('Heatmap for %s', executables(e)),'Color','w');
        set(gcf, 'Renderer', 'opengl');
        imagesc(threads, chunkSizes, heatData);
        set(gca,'YDir','normal');                
        colorbar;
        set(gca, 'FontName', 'Arial', 'FontSize', 12);
        
        % Log color scale if supported
        if verLessThan('matlab','9.1') == 0 % 2016b+ supports log scale
            set(gca, 'ColorScale', 'log');           
        end
        
        xlabel('Threads'); ylabel('Chunk Size');
        title(sprintf('90th Percentile Time Heatmap - %s', executables(e)));
    else
        warning('No data for executable %s', executables(e));
    end
end
