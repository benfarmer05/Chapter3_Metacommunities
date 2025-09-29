clear;clc

%% Initialize paths
projectPath = matlab.project.rootProject().RootFolder;
tempPath = fullfile(projectPath, 'temp');
outputPath = fullfile(projectPath, 'output');

%% =================================================================
%% SEED CONFIGURATION
%% =================================================================
USE_FIXED_SEED = true;  % Set to true to use fixed seed, false for random
FIXED_SEED = 5648;      % Change this to your desired seed value

if USE_FIXED_SEED
    seed = FIXED_SEED;
    rng(seed);
    fprintf('Using fixed seed: %d\n', seed);
else
    seed = randi(10000);
    rng(seed);
    fprintf('Using random seed: %d\n', seed);
end

%% =================================================================
%% TRULY RANDOM PARTICLE SAMPLING VERSION
%% =================================================================

trajlist = dir(fullfile(tempPath,'traj*.nc'));
target_trajectories = 200;

fprintf('Found %d trajectory files\n', length(trajlist));

% Pre-allocate
selected_data = struct('lon', cell(target_trajectories, 1), 'lat', cell(target_trajectories, 1));
trajectories_collected = 0;

% Sample from random files (use more files for better randomness)
files_to_sample = min(length(trajlist), max(5, ceil(length(trajlist) * 0.5)));
selected_files = randsample(length(trajlist), files_to_sample, false);

fprintf('Sampling from %d files for maximum randomness...\n', files_to_sample);

for file_idx = selected_files'
    if trajectories_collected >= target_trajectories, break; end
    
    filename = fullfile(tempPath, trajlist(file_idx).name);
    fprintf('Processing file: %s\n', trajlist(file_idx).name);
    
    try
        % Read all data from file
        lon = ncread(filename, 'lon');
        lat = ncread(filename, 'lat');
        
        % Get total number of particles in this file
        num_particles = size(lon, 2);
        fprintf('  File contains %d particles\n', num_particles);
        
        % Calculate how many particles to sample from this file
        remaining_needed = target_trajectories - trajectories_collected;
        remaining_files = sum(selected_files > file_idx); % files left to process
        
        % Take a reasonable number from this file, but not all from one file
        max_from_this_file = min(num_particles, ceil(remaining_needed / max(1, remaining_files + 1)));
        particles_to_take = min(max_from_this_file, remaining_needed);
        
        fprintf('  Taking %d random particles from this file\n', particles_to_take);
        
        % RANDOMLY sample particle indices (this is the key change!)
        if num_particles <= particles_to_take
            % Take all particles if file is small
            sampled_particle_indices = 1:num_particles;
        else
            % Randomly sample particle indices
            sampled_particle_indices = randsample(num_particles, particles_to_take, false);
        end
        
        % Extract trajectories for randomly selected particles
        for i = 1:length(sampled_particle_indices)
            if trajectories_collected >= target_trajectories, break; end
            
            particle_idx = sampled_particle_indices(i);
            trajectories_collected = trajectories_collected + 1;
            
            % Store the entire trajectory for this randomly selected particle
            selected_data(trajectories_collected).lon = lon(:, particle_idx);
            selected_data(trajectories_collected).lat = lat(:, particle_idx);
        end
        
    catch ME
        warning('Error reading file %s: %s', filename, ME.message);
        continue;
    end
end

% Trim unused entries
selected_data = selected_data(1:trajectories_collected);

fprintf('Successfully loaded %d random trajectories\n', trajectories_collected);

%% Plot trajectories with starting location stars
figure('Position', [100 100 800 600]);

% Load landmask
try
    mapshow(shaperead(fullfile(outputPath, 'landmask_dissolved.shp')));
    fprintf('Landmask loaded successfully\n');
catch
    fprintf('Could not load landmask\n');
end
hold on;

% Plot trajectories and collect start points
colors = lines(10);
start_lons = zeros(trajectories_collected, 1);
start_lats = zeros(trajectories_collected, 1);

fprintf('Plotting %d trajectories...\n', trajectories_collected);

for i = 1:trajectories_collected
    % Plot trajectory
    plot(selected_data(i).lon, selected_data(i).lat, 'Color', colors(mod(i-1,10)+1,:), 'LineWidth', 0.5);
    
    % Store starting position (first valid point)
    valid_idx = find(~isnan(selected_data(i).lon) & ~isnan(selected_data(i).lat), 1);
    if ~isempty(valid_idx)
        start_lons(i) = selected_data(i).lon(valid_idx);
        start_lats(i) = selected_data(i).lat(valid_idx);
    else
        start_lons(i) = NaN;
        start_lats(i) = NaN;
    end
end

% Plot starting locations as red stars
valid_starts = ~isnan(start_lons);
plot(start_lons(valid_starts), start_lats(valid_starts), 'r*', 'MarkerSize', 8, 'LineWidth', 1.5);

fprintf('Plotted %d starting points\n', sum(valid_starts));

% Calculate trajectory duration for title
trajectory_duration_hours = 0;
trajectory_duration_days = 0;
if trajectories_collected > 0
    % Get time data from first file to calculate duration
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
    title(sprintf('%d random trajectories with start points (%.1f days, seed:%d)', trajectories_collected, trajectory_duration_days, seed));
else
    title(sprintf('%d random trajectories with start points (%.1f hours, seed:%d)', trajectories_collected, trajectory_duration_hours, seed));
end

fprintf('Plot completed!\n');
if trajectory_duration_hours > 0
    fprintf('Trajectory duration: %.1f hours (%.1f days)\n', trajectory_duration_hours, trajectory_duration_days);
end

% Show spatial distribution of start points
if sum(valid_starts) > 0
    fprintf('\nStart point distribution:\n');
    fprintf('  Longitude range: %.3f to %.3f\n', min(start_lons(valid_starts)), max(start_lons(valid_starts)));
    fprintf('  Latitude range: %.3f to %.3f\n', min(start_lats(valid_starts)), max(start_lats(valid_starts)));
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