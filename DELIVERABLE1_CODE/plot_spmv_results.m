%% plot_spmv_results.m
% Clear workspace
clear; close all; clc;

%% Load data
filename = 'results.csv'; % <-- change to your CSV filename
data = readtable(filename);

% Convert Executable column to categorical for easy grouping
data.Executable = categorical(data.Executable);

%% Compute statistics
executables = categories(data.Executable);
nExec = numel(executables);

meanTimes = struct();
stdTimes = struct();

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
markers = {'o', 's', 'd'};
figure;
hold on; grid on; box on;

for i = 1:nExec
    exec = executables{i};
    threads = threadCounts.(exec);
    mTimes = meanTimes.(exec);
    sTimes = stdTimes.(exec);
    
    errorbar(threads, mTimes, sTimes, ...
        'LineWidth', 1.5, ...
        'Marker', markers{i}, ...
        'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:));
end

%% Aesthetics
xlabel('Number of Threads', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Execution Time (s)', 'Interpreter', 'latex', 'FontSize', 12);
title('SpMV Performance with Different Scheduling Policies', ...
    'Interpreter', 'latex', 'FontSize', 13);
legend(strrep(executables, 'spmv_', ''), 'Location', 'northwest', 'Interpreter', 'latex');
set(gca, 'YScale', 'log'); % log scale often helps with time data
set(gca, 'FontSize', 11);
xticks(unique(data.Threads));
axis tight;

%% Save figure (optional)
exportgraphics(gcf, 'spmv_performance.png', 'Resolution', 300);
