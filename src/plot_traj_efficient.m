clear;clc

%% Initialize paths
projectPath = matlab.project.rootProject().RootFolder;
tempPath = fullfile(projectPath, 'temp');
outputPath = fullfile(projectPath, 'output');

%% =================================================================
%% OPTION A: DETAILED VERSION (with debugging and verbose output)
%% =================================================================
%{
trajlist = dir(fullfile(tempPath,'traj*.nc'));
fprintf('Found %d trajectory files\n', length(trajlist));

if isempty(trajlist)
    error('No trajectory files found in %s', tempPath);
end

% Quick metadata scan to understand data structure
sample_file = fullfile(tempPath, trajlist(1).name);
try
    sample_location = ncread(sample_file, 'location');
    sample_lon = ncread(sample_file, 'lon');
    fprintf('Sample file structure: %d locations, %dx%d lon matrix\n', ...
        length(sample_location), size(sample_lon,1), size(sample_lon,2));
catch ME
    error('Cannot read sample file %s: %s', sample_file, ME.message);
end

% Smart sampling parameters
target_trajectories = 200; % Fixed reasonable number for testing
fprintf('Target trajectories to plot: %d\n', target_trajectories);

% Random seed
seed = randi(10000);
rng(seed);
fprintf('Using random seed: %d\n', seed);

% Pre-allocate trajectory storage for efficiency
selected_data = struct('lon', cell(target_trajectories, 1), 'lat', cell(target_trajectories, 1));
trajectories_collected = 0;

% Calculate how many files to sample (at least 25% of files, up to all)
files_to_sample = max(ceil(length(trajlist) * 0.25), min(10, length(trajlist)));
selected_file_indices = randsample(length(trajlist), files_to_sample, false);

fprintf('Sampling from %d files...\n', files_to_sample);

tic;
for i = 1:length(selected_file_indices)
    if trajectories_collected >= target_trajectories
        break;
    end
    
    file_idx = selected_file_indices(i);
    filename = fullfile(tempPath, trajlist(file_idx).name);
    
    fprintf('Processing file %d/%d: %s\n', i, length(selected_file_indices), trajlist(file_idx).name);
    
    try
        % Read file data
        location = ncread(filename, 'location');
        lon = ncread(filename, 'lon');
        lat = ncread(filename, 'lat');
        
        fprintf('  File contains %d location entries\n', length(location));
        
        % Get unique locations in this file
        unique_locs = unique(location);
        fprintf('  Unique locations: %d\n', length(unique_locs));
        
        % Calculate how many trajectories to take from this file
        remaining_needed = target_trajectories - trajectories_collected;
        remaining_files = length(selected_file_indices) - i + 1;
        trajectories_from_this_file = min(length(unique_locs), ...
            ceil(remaining_needed / remaining_files));
        
        fprintf('  Taking %d trajectories from this file\n', trajectories_from_this_file);
        
        if trajectories_from_this_file > 0
            % Randomly select locations from this file
            if length(unique_locs) <= trajectories_from_this_file
                sampled_locs = unique_locs;
            else
                sampled_locs = randsample(unique_locs, trajectories_from_this_file, false);
            end
            
            % Extract data for each selected location
            for j = 1:length(sampled_locs)
                if trajectories_collected >= target_trajectories
                    break;
                end
                
                loc = sampled_locs(j);
                loc_indices = find(location == loc);
                
                if ~isempty(loc_indices)
                    % Take all particles for this location (or sample if too many)
                    if length(loc_indices) > 10 % Limit particles per location for visualization
                        loc_indices = randsample(loc_indices, 10, false);
                    end
                    
                    for k = 1:length(loc_indices)
                        if trajectories_collected >= target_trajectories
                            break;
                        end
                        trajectories_collected = trajectories_collected + 1;
                        idx = loc_indices(k);
                        selected_data(trajectories_collected).lon = lon(:, idx);
                        selected_data(trajectories_collected).lat = lat(:, idx);
                    end
                end
            end
        end
        
    catch ME
        warning('ME:FileError', 'Error reading file %s: %s', filename, ME.message);
        continue;
    end
end

load_time = toc;

% Trim unused pre-allocated entries
if trajectories_collected < target_trajectories
    selected_data = selected_data(1:trajectories_collected);
end

fprintf('\nSuccessfully loaded %d trajectories in %.2f seconds\n', trajectories_collected, load_time);

if trajectories_collected == 0
    error('No trajectories were successfully loaded. Check your data files.');
end

%% Plot the sampled data (DETAILED VERSION)
fprintf('Creating plot...\n');
figure('Position', [100 100 800 600]);

% Load and plot landmask
try
    landmask = shaperead(fullfile(outputPath, 'landmask_dissolved.shp'));
    mapshow(landmask);
    fprintf('Landmask loaded successfully\n');
catch ME
    warning('ME:LandmaskError', 'Could not load landmask: %s', ME.message);
end

hold on;

% Plot trajectories and collect start points
tic;
colors = lines(min(10, trajectories_collected));
start_lons = zeros(trajectories_collected, 1);
start_lats = zeros(trajectories_collected, 1);

for i = 1:trajectories_collected
    color_idx = mod(i-1, size(colors,1)) + 1;
    
    % Plot trajectory
    plot(selected_data(i).lon, selected_data(i).lat, 'Color', colors(color_idx,:), 'LineWidth', 0.5);
    
    % Store starting position
    valid_indices = ~isnan(selected_data(i).lon) & ~isnan(selected_data(i).lat);
    if any(valid_indices)
        first_valid = find(valid_indices, 1);
        start_lons(i) = selected_data(i).lon(first_valid);
        start_lats(i) = selected_data(i).lat(first_valid);
    else
        start_lons(i) = NaN;
        start_lats(i) = NaN;
    end
end

% Plot starting locations as red stars
valid_starts = ~isnan(start_lons) & ~isnan(start_lats);
plot(start_lons(valid_starts), start_lats(valid_starts), 'r*', 'MarkerSize', 8, 'LineWidth', 1.5);

plot_time = toc;

% Finalize plot
axis equal; 
try
    axis([294, 296, 17.6, 19.2]); % Your specified bounds
catch
    fprintf('Using automatic axis bounds\n');
end

title(sprintf('%d trajectories (seed:%d, load:%.1fs, plot:%.2fs)', ...
    trajectories_collected, seed, load_time, plot_time));

fprintf('Plotting completed in %.2f seconds\n', plot_time);
fprintf('Total time: %.2f seconds\n', load_time + plot_time);

% Show some statistics
if trajectories_collected > 0
    all_lons = [selected_data.lon];
    all_lats = [selected_data.lat];
    fprintf('\nTrajectory bounds:\n');
    fprintf('  Longitude: %.3f to %.3f\n', min(all_lons(:)), max(all_lons(:)));
    fprintf('  Latitude: %.3f to %.3f\n', min(all_lats(:)), max(all_lats(:)));
end
%}

%% =================================================================
%% OPTION B: EFFICIENT VERSION (minimal output, maximum speed)
%% =================================================================

trajlist = dir(fullfile(tempPath,'traj*.nc'));
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





%% =================================================================
%% TIME-STEP FREQUENCY ANALYSIS
%% =================================================================

fprintf('\n=== Analyzing time-step frequency ===\n');

% Use the already loaded trajectory data to analyze time steps
if trajectories_collected > 0
    % Analyze first few trajectories for time stepping
    sample_file = fullfile(tempPath, trajlist(selected_files(1)).name);
    
    try
        time_data = ncread(sample_file, 'time');
        time_info = ncinfo(sample_file, 'time');
        
        % Get units
        units_attr = time_info.Attributes(strcmp({time_info.Attributes.Name}, 'units'));
        if ~isempty(units_attr)
            fprintf('Time units: %s\n', units_attr.Value);
        end
        
        if length(time_data) > 1
            time_diffs = diff(time_data);
            mean_diff = mean(time_diffs);
            
            fprintf('Time steps: %d total points\n', length(time_data));
            fprintf('Mean time step: %.3f time units\n', mean_diff);
            
            % Convert to hours if units suggest days
            if ~isempty(units_attr) && contains(lower(units_attr.Value), 'day')
                hours = mean_diff * 24;
                fprintf('Mean time step: %.1f hours\n', hours);
                if abs(hours - 24) < 1
                    fprintf('*** CONFIRMED: ~24 hour time steps ***\n');
                end
            end
            
            % Show first few time points
            fprintf('First 5 time values: ');
            fprintf('%.2f ', time_data(1:min(5,length(time_data))));
            fprintf('\n');
        end
        
    catch ME
        fprintf('Could not analyze time steps: %s\n', ME.message);
    end
else
    fprintf('No trajectories loaded for analysis\n');
end