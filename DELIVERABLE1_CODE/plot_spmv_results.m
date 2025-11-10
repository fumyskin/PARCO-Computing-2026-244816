%% plot_spmv_results_balanced.m
% Clear workspace
clear; close all; clc;

%% Load data
filename = 'results_utm.csv'; % <-- change to your CSV filename
data = readtable(filename);

% Convert Executable column to categorical for easy grouping
data.Executable = categorical(data.Executable);

%% Compute statistics
executables = categories(data.Executable);
nExec = numel(executables);

meanTimes = struct();
stdTimes = struct();
threadCounts = struct();

for i = 1:nExec
    exec = executables{i};
    subset = data(data.Executable == exec, :);
    
    threads = unique(subset.Threads);
    mTimes = zeros(size(threads));
    sTimes = zeros(size(threads));
    
    for j = 1:numel(threads)
        t = threads(j);
        runs = subset.Time(subset.Threads == t);
        mTimes(j) = mean(runs);
        sTimes(j) = std(runs);
    end
    
    meanTimes.(exec) = mTimes;
    stdTimes.(exec) = sTimes;
    threadCounts.(exec) = threads;
end

%% Plot settings
colors = lines(nExec);
markers = {'o', 's', 'd', '^', 'v'};
figure;
hold on; grid on; box on;

% Create evenly spaced x positions
allThreads = unique(data.Threads);
xPositions = 1:numel(allThreads); % evenly spaced indices

for i = 1:nExec
    exec = executables{i};
    threads = threadCounts.(exec);
    mTimes = meanTimes.(exec);
    sTimes = stdTimes.(exec);
    
    % Map thread numbers to evenly spaced x positions
    [~, idx] = ismember(threads, allThreads);
    
    errorbar(xPositions(idx), mTimes, sTimes, ...
        'LineWidth', 1.5, ...
        'Marker', markers{mod(i-1, numel(markers))+1}, ...
        'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:));
end

%% Aesthetics
xlabel('Number of Threads', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Execution Time (s)', 'Interpreter', 'latex', 'FontSize', 12);
title('SpMV Performance with Different Scheduling Policies', ...
    'Interpreter', 'latex', 'FontSize', 13);
legend(strrep(executables, 'spmv_', ''), 'Location', 'northwest', 'Interpreter', 'latex');

% X-axis evenly spaced with real thread labels
set(gca, 'XTick', xPositions);
set(gca, 'XTickLabel', string(allThreads));
set(gca, 'FontSize', 11);

% Use a log scale for the y-axis to make small sequential runs visible
set(gca, 'YScale', 'log');
ylim([1e-5, 1]); % adjust upper bound if you have higher times

axis tight;
grid on;

%% Save figure (optional)
exportgraphics(gcf, 'spmv_performance_balanced.png', 'Resolution', 300);
