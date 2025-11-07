%% ========== DECAY WEIGHT CALCULATOR AND VISUALIZER ==========
% Calculates and visualizes the exponential decay weights used in 
% connectivity analysis. Can use actual trajectory data or simulate.
clear; clc;

%% Configuration
YEAR = 2019;
QUARTER = 4;
DECAY_HALFLIFE_DAYS = 7;   % Half-life for exponential decay
TIMESTEP_MINUTES = 15;      % Time resolution (from trajectory data)

% Platform selection
USE_MAC = true;  % Set to true for Mac, false for Windows

fprintf('=== DECAY WEIGHT CALCULATOR ===\n');
fprintf('Decay half-life: %.1f days\n', DECAY_HALFLIFE_DAYS);
fprintf('Timestep resolution: %d minutes\n', TIMESTEP_MINUTES);
if USE_MAC
    fprintf('Platform: Mac\n');
else
    fprintf('Platform: Windows\n');
end

%% Setup paths (cross-platform)
projectPath = matlab.project.rootProject().RootFolder;
dataPath = fullfile(projectPath, 'data');
quarter_name = sprintf('Q%d_%d', QUARTER, YEAR);
outputPath = fullfile(projectPath, 'output');

if USE_MAC
    % Mac paths
    tempPath = fullfile('/Volumes', 'Farmer_diss', 'Dissertation', 'CMS_traj', quarter_name);
else
    % Windows paths
    tempPath = fullfile('D:', 'Dissertation', 'CMS_traj', quarter_name);
end

fprintf('Quarter: %s\n', quarter_name);
fprintf('Trajectory path: %s\n', tempPath);
fprintf('Output path: %s\n', outputPath);

%% Find and load a trajectory file
fprintf('\nSearching for trajectory files...\n');
trajlist = dir(fullfile(tempPath, 'traj*.nc'));

if isempty(trajlist)
    error('No trajectory files found in %s', tempPath);
end

fprintf('Found %d trajectory files\n', length(trajlist));

% Use first trajectory file
traj_file = fullfile(tempPath, trajlist(1).name);
fprintf('Loading trajectory file: %s\n', trajlist(1).name);

try
    % Read trajectory data to get timestep count
    lon = ncread(traj_file, 'lon');
    n_timesteps = size(lon, 1);
    fprintf('Found %d timesteps in trajectory file\n', n_timesteps);
catch ME
    error('Could not read trajectory file: %s', ME.message);
end

%% Calculate decay weights
fprintf('\nCalculating decay weights...\n');

% Decay constant (lambda for exponential decay)
decay_constant = log(2) / DECAY_HALFLIFE_DAYS;
fprintf('Decay constant (λ): %.6f day^-1\n', decay_constant);

% Time vector
timesteps_per_day = (24 * 60) / TIMESTEP_MINUTES;
timestep_indices = (0:n_timesteps-1)';  % Start at 0 (release time)
time_in_days = timestep_indices / timesteps_per_day;
time_in_hours = timestep_indices * (TIMESTEP_MINUTES / 60);

% Calculate decay weights
decay_weights = exp(-decay_constant * time_in_days);

fprintf('Number of timesteps: %d\n', n_timesteps);
fprintf('Total duration: %.2f days (%.1f hours)\n', max(time_in_days), max(time_in_hours));
fprintf('Initial weight (t=0): %.4f\n', decay_weights(1));
fprintf('Final weight (t=%.1f days): %.6f\n', max(time_in_days), decay_weights(end));
fprintf('Weight at half-life (t=%.1f days): %.4f\n', DECAY_HALFLIFE_DAYS, ...
    exp(-decay_constant * DECAY_HALFLIFE_DAYS));

%% Create summary table for key timepoints
key_days = [0, 1, 3, 7, 14, 21, 30];
key_days = key_days(key_days <= max(time_in_days));

fprintf('\n=== DECAY WEIGHTS AT KEY TIMEPOINTS ===\n');
fprintf('Day\tWeight\t\tRelative to t=0\n');
fprintf('---\t------\t\t---------------\n');
for d = key_days
    weight = exp(-decay_constant * d);
    relative = weight / decay_weights(1);
    fprintf('%.1f\t%.6f\t%.1f%%\n', d, weight, relative * 100);
end

%% Save results to file
output_file = fullfile(outputPath, sprintf('decay_weights_halflife_%.1fdays.mat', DECAY_HALFLIFE_DAYS));
save(output_file, 'decay_weights', 'time_in_days', 'time_in_hours', ...
    'timestep_indices', 'decay_constant', 'DECAY_HALFLIFE_DAYS', ...
    'TIMESTEP_MINUTES', 'n_timesteps');
fprintf('\nDecay weights saved to: %s\n', output_file);

% Also save as CSV for easy inspection
csv_file = fullfile(outputPath, sprintf('decay_weights_halflife_%.1fdays.csv', DECAY_HALFLIFE_DAYS));
decay_table = table(timestep_indices, time_in_hours, time_in_days, decay_weights, ...
    'VariableNames', {'Timestep', 'Hours', 'Days', 'DecayWeight'});
writetable(decay_table, csv_file);
fprintf('Decay table saved to: %s\n', csv_file);

%% Visualization 1: Decay curve over time
figure('Position', [100 100 1200 400]);

% Plot in days
subplot(1,3,1);
plot(time_in_days, decay_weights, 'b-', 'LineWidth', 2);
hold on;
% Mark half-life point
if max(time_in_days) >= DECAY_HALFLIFE_DAYS
    plot(DECAY_HALFLIFE_DAYS, 0.5, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    text(DECAY_HALFLIFE_DAYS, 0.5, sprintf('  Half-life\n  (%.1f days)', DECAY_HALFLIFE_DAYS), ...
        'VerticalAlignment', 'bottom', 'FontSize', 10);
end
xlabel('Time since release (days)', 'FontSize', 11);
ylabel('Decay Weight', 'FontSize', 11);
title('Exponential Decay Curve', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
ylim([0 1.05]);

% Plot in hours (first 7 days)
subplot(1,3,2);
max_hours_to_show = min(7 * 24, max(time_in_hours));
mask = time_in_hours <= max_hours_to_show;
plot(time_in_hours(mask), decay_weights(mask), 'b-', 'LineWidth', 2);
xlabel('Time since release (hours)', 'FontSize', 11);
ylabel('Decay Weight', 'FontSize', 11);
title('First Week Detail', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
ylim([0 1.05]);

% Log scale plot
subplot(1,3,3);
semilogy(time_in_days, decay_weights, 'b-', 'LineWidth', 2);
hold on;
if max(time_in_days) >= DECAY_HALFLIFE_DAYS
    plot(DECAY_HALFLIFE_DAYS, 0.5, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
end
xlabel('Time since release (days)', 'FontSize', 11);
ylabel('Decay Weight (log scale)', 'FontSize', 11);
title('Log-Scale View', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
ylim([min(decay_weights)/2, 1.5]);

sgtitle(sprintf('Decay Weights (Half-life = %.1f days, %d-min timesteps)', ...
    DECAY_HALFLIFE_DAYS, TIMESTEP_MINUTES), 'FontSize', 14, 'FontWeight', 'bold');

% Save figure
fig_file = fullfile(outputPath, sprintf('decay_visualization_halflife_%.1fdays.png', DECAY_HALFLIFE_DAYS));
saveas(gcf, fig_file);
fprintf('Figure saved to: %s\n', fig_file);

%% Visualization 2: Comparison of different half-lives
figure('Position', [100 100 900 600]);
hold on;

halflife_values = [3, 7, 14, 21, 30];
colors = lines(length(halflife_values));

for i = 1:length(halflife_values)
    hl = halflife_values(i);
    lambda = log(2) / hl;
    weights = exp(-lambda * time_in_days);
    plot(time_in_days, weights, 'LineWidth', 2, 'Color', colors(i,:), ...
        'DisplayName', sprintf('%.0f-day half-life', hl));
end

xlabel('Time since release (days)', 'FontSize', 12);
ylabel('Decay Weight', 'FontSize', 12);
title('Comparison of Different Half-Lives', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 11);
grid on;
ylim([0 1.05]);
xlim([0 min(30, max(time_in_days))]);

fig_file2 = fullfile(outputPath, 'decay_comparison_multiple_halflifes.png');
saveas(gcf, fig_file2);
fprintf('Comparison figure saved to: %s\n', fig_file2);

%% Summary statistics
fprintf('\n=== SUMMARY STATISTICS ===\n');
fprintf('Sum of all decay weights: %.2f\n', sum(decay_weights));
fprintf('Mean decay weight: %.4f\n', mean(decay_weights));
fprintf('Median decay weight: %.4f\n', median(decay_weights));
fprintf('Std deviation: %.4f\n', std(decay_weights));

% Calculate effective duration (when weight drops to 1% of initial)
threshold_idx = find(decay_weights <= 0.01, 1);
if ~isempty(threshold_idx)
    fprintf('Time to 1%% weight: %.1f days (timestep %d)\n', ...
        time_in_days(threshold_idx), threshold_idx-1);
end

fprintf('\n=== COMPLETE ===\n');