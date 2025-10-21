clear;clc

%% Initialize paths
projectPath = matlab.project.rootProject().RootFolder;
tempPath = fullfile(projectPath, 'temp');
outputPath = fullfile(projectPath, 'output');

%% =================================================================
%% EFFICIENT VERSION with NA percentage monitoring
%% =================================================================

trajlist = dir(fullfile(tempPath,'traj*.nc'));

%% Check NA percentages across all files
fprintf('\n=== Checking NA percentages in trajectory files ===\n');
fprintf('This helps determine simulation completion status\n\n');

% Sample files to check (check all if < 50 files, otherwise sample)
files_to_check = min(length(trajlist), 50);
if length(trajlist) > 50
    check_indices = round(linspace(1, length(trajlist), files_to_check));
else
    check_indices = 1:length(trajlist);
end

% Initialize tracking
na_stats = struct();

for i = 1:length(check_indices)
    file_idx = check_indices(i);
    filename = fullfile(tempPath, trajlist(file_idx).name);
    
    fprintf('File %d/%d: %s\n', i, length(check_indices), trajlist(file_idx).name);
    
    try
        % Get file info
        file_info = ncinfo(filename);
        
        % Check each variable
        for j = 1:length(file_info.Variables)
            var_name = file_info.Variables(j).Name;
            
            try
                var_data = ncread(filename, var_name);
                
                % Calculate NA percentage
                if isnumeric(var_data)
                    total_elements = numel(var_data);
                    na_elements = sum(isnan(var_data(:)));
                    na_percent = (na_elements / total_elements) * 100;
                    
                    % Store in structure (accumulate across files)
                    if ~isfield(na_stats, var_name)
                        na_stats.(var_name) = [];
                    end
                    na_stats.(var_name)(end+1) = na_percent;
                    
                    fprintf('  %s: %.1f%% NA (%d/%d elements)\n', ...
                        var_name, na_percent, na_elements, total_elements);
                end
            catch
                fprintf('  %s: Could not read\n', var_name);
            end
        end
        fprintf('\n');
        
    catch ME
        fprintf('  Error reading file: %s\n\n', ME.message);
    end
end

%% Summary statistics
fprintf('=== SUMMARY: Average NA percentages across checked files ===\n');
var_names = fieldnames(na_stats);
for i = 1:length(var_names)
    avg_na = mean(na_stats.(var_names{i}));
    min_na = min(na_stats.(var_names{i}));
    max_na = max(na_stats.(var_names{i}));
    fprintf('%s: %.1f%% avg (range: %.1f%% - %.1f%%)\n', ...
        var_names{i}, avg_na, min_na, max_na);
end
fprintf('\n');

%% Plot trajectories
target_trajectories = 200;
rng(randi(10000));

% Pre-allocate
selected_data = struct('lon', cell(target_trajectories, 1), 'lat', cell(target_trajectories, 1));
trajectories_collected = 0;

% Sample from random files
selected_files = randsample(length(trajlist), min(10, length(trajlist)), false);

for i = selected_files'
    if trajectories_collected >= target_trajectories, break; end
    
    filename = fullfile(tempPath, trajlist(i).name);
    try
        location = ncread(filename, 'location');
        lon = ncread(filename, 'lon');
        lat = ncread(filename, 'lat');
        
        unique_locs = unique(location);
        n_take = min(length(unique_locs), target_trajectories - trajectories_collected);
        sampled_locs = randsample(unique_locs, n_take, false);
        
        for loc = sampled_locs'
            if trajectories_collected >= target_trajectories, break; end
            loc_idx = find(location == loc, 1);
            trajectories_collected = trajectories_collected + 1;
            selected_data(trajectories_collected).lon = lon(:, loc_idx);
            selected_data(trajectories_collected).lat = lat(:, loc_idx);
        end
    catch
        continue;
    end
end

% Trim unused entries
selected_data = selected_data(1:trajectories_collected);

%% Plot trajectories with starting location stars
figure('Position', [100 100 800 600]);

% Load landmask
try
    mapshow(shaperead(fullfile(outputPath, 'landmask_dissolved.shp')));
catch
end
hold on;

% Plot trajectories and collect start points
colors = lines(10);
start_lons = zeros(trajectories_collected, 1);
start_lats = zeros(trajectories_collected, 1);

for i = 1:trajectories_collected
    % Plot trajectory
    plot(selected_data(i).lon, selected_data(i).lat, 'Color', colors(mod(i-1,10)+1,:), 'LineWidth', 0.5);
    
    % Store starting position
    valid_idx = find(~isnan(selected_data(i).lon) & ~isnan(selected_data(i).lat), 1);
    if ~isempty(valid_idx)
        start_lons(i) = selected_data(i).lon(valid_idx);
        start_lats(i) = selected_data(i).lat(valid_idx);
    end
end

% Plot starting locations as red stars
valid_starts = ~isnan(start_lons);
plot(start_lons(valid_starts), start_lats(valid_starts), 'r*', 'MarkerSize', 8, 'LineWidth', 1.5);

% Calculate trajectory duration for title
trajectory_duration_hours = 0;
trajectory_duration_days = 0;
if trajectories_collected > 0
    % Get time data from first trajectory to calculate duration
    sample_file = fullfile(tempPath, trajlist(selected_files(1)).name);
    try
        time_data = ncread(sample_file, 'time');
        if length(time_data) > 1
            total_duration_seconds = time_data(end) - time_data(1);
            trajectory_duration_hours = total_duration_seconds / 3600;
            trajectory_duration_days = trajectory_duration_hours / 24;
        end
    catch
        % If we can't read time data, skip duration calculation
    end
end

axis equal; axis([294, 296, 17.6, 19.2]);

% Create title with trajectory duration
if trajectory_duration_days >= 1
    title(sprintf('%d trajectories with start points (%.1f days)', trajectories_collected, trajectory_duration_days));
else
    title(sprintf('%d trajectories with start points (%.1f hours)', trajectories_collected, trajectory_duration_hours));
end

fprintf('Plotted %d trajectories\n', trajectories_collected);
if trajectory_duration_hours > 0
    fprintf('Trajectory duration: %.1f hours (%.1f days)\n', trajectory_duration_hours, trajectory_duration_days);
end