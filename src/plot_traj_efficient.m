clear;clc

%% Initialize paths
projectPath = matlab.project.rootProject().RootFolder;
tempPath = fullfile(projectPath, 'temp');
outputPath = fullfile(projectPath, 'output');

%% =================================================================
%% EFFICIENT VERSION with TRUE randomness
%% =================================================================
    
trajlist = dir(fullfile(tempPath,'traj*.nc'));
target_trajectories = 200;

% TRUE random seed based on current time
rng('shuffle');  % Initialize with current time
seed = randi(10000);  % Now this will be truly random
rng(seed);  % Set the seed (optional: for reproducibility if you save the seed)

fprintf('Using random seed: %d\n', seed);

% Pre-allocate
selected_data = struct('lon', cell(target_trajectories, 1), 'lat', cell(target_trajectories, 1));
trajectories_collected = 0;

% Sample from random files (now truly random)
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

% Create title with trajectory duration and seed
if trajectory_duration_days >= 1
    title(sprintf('%d trajectories (seed:%d, %.1f days)', trajectories_collected, seed, trajectory_duration_days));
else
    title(sprintf('%d trajectories (seed:%d, %.1f hours)', trajectories_collected, seed, trajectory_duration_hours));
end

fprintf('Plotted %d trajectories\n', trajectories_collected);
if trajectory_duration_hours > 0
    fprintf('Trajectory duration: %.1f hours (%.1f days)\n', trajectory_duration_hours, trajectory_duration_days);
end






% %% version just plotting locations
% 
% 
% clear;clc
% 
% %% Initialize paths
% projectPath = matlab.project.rootProject().RootFolder;
% tempPath = fullfile(projectPath, 'temp');
% outputPath = fullfile(projectPath, 'output');
% 
% %% =================================================================
% %% Plot all unique locations (no trajectories)
% %% =================================================================
% 
% trajlist = dir(fullfile(tempPath,'traj*.nc'));
% 
% % Collect all unique starting locations
% all_start_lons = [];
% all_start_lats = [];
% 
% fprintf('Reading trajectory files...\n');
% 
% for i = 1:length(trajlist)
%     filename = fullfile(tempPath, trajlist(i).name);
%     try
%         location = ncread(filename, 'location');
%         lon = ncread(filename, 'lon');
%         lat = ncread(filename, 'lat');
% 
%         unique_locs = unique(location);
% 
%         % Get starting position for each unique location
%         for loc = unique_locs'
%             loc_idx = find(location == loc, 1);
%             valid_idx = find(~isnan(lon(:, loc_idx)) & ~isnan(lat(:, loc_idx)), 1);
%             if ~isempty(valid_idx)
%                 all_start_lons(end+1) = lon(valid_idx, loc_idx);
%                 all_start_lats(end+1) = lat(valid_idx, loc_idx);
%             end
%         end
%     catch
%         continue;
%     end
% end
% 
% %% Plot unique locations only
% figure('Position', [100 100 800 600]);
% 
% % Load landmask
% try
%     mapshow(shaperead(fullfile(outputPath, 'landmask_dissolved.shp')));
% catch
% end
% hold on;
% 
% % Plot all starting locations as red stars
% plot(all_start_lons, all_start_lats, 'r*', 'MarkerSize', 8, 'LineWidth', 1.5);
% 
% axis equal; axis([294, 296, 17.6, 19.2]);
% 
% title(sprintf('%d unique trajectory starting locations', length(all_start_lons)));
% 
% fprintf('Plotted %d unique starting locations\n', length(all_start_lons));





% %% check traj files to make totally certain their IDs are good
% 
% %% Simple Trajectory to Polygon Check - Multiple Locations
% clear; clc;
% 
% %% Load data
% projectPath = matlab.project.rootProject().RootFolder;
% dataPath = fullfile(projectPath, 'data');
% tempPath = fullfile('D:\Dissertation\CMS_traj', 'Q1_2019');
% 
% % Load polygons
% centroids = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
% unique_IDs = centroids(:,1);
% Xs = [centroids(:,8) centroids(:,10) centroids(:,12) centroids(:,14) centroids(:,8)];
% Ys = [centroids(:,9) centroids(:,11) centroids(:,13) centroids(:,15) centroids(:,9)];
% 
% % Load release file
% releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019_Q1.txt'));
% release_IDs = releasefile(:,1);
% 
% %% Get random trajectory file
% trajlist = dir(fullfile(tempPath, 'traj*.nc'));
% rand_file_idx = randi(length(trajlist));
% filename = fullfile(tempPath, trajlist(rand_file_idx).name);
% 
% fprintf('Selected file: %s\n', trajlist(rand_file_idx).name);
% 
% %% Read trajectory data
% location = ncread(filename, 'location');  % Zero-based line numbers
% lon = ncread(filename, 'lon');
% lat = ncread(filename, 'lat');
% 
% % Convert lon if needed
% if min(lon(:)) > 180
%     lon = lon - 360;
% end
% 
% %% Get unique locations in this file
% unique_locations = unique(location);
% n_unique = length(unique_locations);
% fprintf('Found %d unique locations in this file\n', n_unique);
% fprintf('Location values range: %d to %d\n', min(unique_locations), max(unique_locations));
% fprintf('Release file has %d entries (indices 1 to %d)\n', length(release_IDs), length(release_IDs));
% 
% %% Check a few random ones (pick 5)
% n_check = min(5, n_unique);
% check_locs = unique_locations(randperm(n_unique, n_check));
% 
% figure('Position', [100 100 1400 900]);
% 
% for i = 1:n_check
%     subplot(2, 3, i);
%     hold on;
% 
%     loc_zero = check_locs(i);
%     fprintf('\n--- Checking location %d (zero-based) ---\n', loc_zero);
% 
%     % CRITICAL: Check if this is zero-based or one-based
%     % If location values start at 0, add 1 for MATLAB indexing
%     % If location values start at 1, use as-is
%     if min(unique_locations) == 0
%         loc_one = loc_zero + 1;  % Zero-based, add 1
%         fprintf('  Converting to one-based: %d\n', loc_one);
%     else
%         loc_one = loc_zero;  % Already one-based
%         fprintf('  Already one-based: %d\n', loc_one);
%     end
% 
%     % Get reef ID from release file
%     reef_id = release_IDs(loc_one);
%     fprintf('  Reef ID from release file line %d: %d\n', loc_one, reef_id);
% 
%     % Find polygon
%     poly_idx = find(unique_IDs == reef_id);
%     poly_x = Xs(poly_idx, :);
%     poly_y = Ys(poly_idx, :);
% 
%     % Get particles from this location
%     mask = location == loc_zero;
%     loc_lon = lon(:, mask);
%     loc_lat = lat(:, mask);
% 
%     % Plot polygon
%     plot(poly_x, poly_y, 'r-', 'LineWidth', 2);
% 
%     % Plot trajectories
%     plot(loc_lon, loc_lat, 'b-', 'LineWidth', 0.5);
% 
%     % Mark start points
%     start_lon = loc_lon(1, :);
%     start_lat = loc_lat(1, :);
%     plot(start_lon, start_lat, 'r*', 'MarkerSize', 8);
% 
%     % Check if inside
%     inside = inpolygon(start_lon, start_lat, poly_x, poly_y);
% 
%     axis equal tight;
%     grid on;
%     title(sprintf('ID %d: %d/%d inside', reef_id, sum(inside), length(inside)));
% end
% 
% sgtitle(sprintf('File: %s', trajlist(rand_file_idx).name), 'Interpreter', 'none');