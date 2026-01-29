%% Plot 1000 random particle trajectories
clear; clc;

%% Setup paths
projectPath = matlab.project.rootProject().RootFolder;
dataPath = fullfile(projectPath, 'data');
tempPath = 'D:\Dissertation\CMS_traj\Q1_2019';

%% Configuration
target_trajectories = 1000;
lon_min = -65.2;
lon_max = -64.6;
lat_min = 18.0;
lat_max = 18.5;

%% Plotting parameters
grid_alpha = 0.3;        % Grid cell transparency (0 = invisible, 1 = solid)
trajectory_alpha = 0.75;  % Trajectory transparency (0 = invisible, 1 = solid)
star_size = 10;          % Red star marker size

rng('shuffle');
seed = randi(10000);
% seed = 939; %for showing complexity in open ocean
rng(seed);
fprintf('Using random seed: %d\n', seed);

%% Load reef polygons and filter to bounds
fprintf('Loading reef polygons...\n');
centroids = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
Xs = [centroids(:,8) centroids(:,10) centroids(:,12) centroids(:,14) centroids(:,8)];
Ys = [centroids(:,9) centroids(:,11) centroids(:,13) centroids(:,15) centroids(:,9)];

% Filter to only polygons in bounded region
poly_in_bounds = false(size(Xs, 1), 1);
for i = 1:size(Xs, 1)
    poly_center_x = mean(Xs(i, 1:4));
    poly_center_y = mean(Ys(i, 1:4));
    if poly_center_x >= lon_min && poly_center_x <= lon_max && ...
       poly_center_y >= lat_min && poly_center_y <= lat_max
        poly_in_bounds(i) = true;
    end
end

Xs_bounded = Xs(poly_in_bounds, :);
Ys_bounded = Ys(poly_in_bounds, :);
fprintf('Filtered to %d polygons in bounded region\n', sum(poly_in_bounds));

%% Collect trajectories within spatial bounds
fprintf('Collecting up to %d trajectories from bounded region...\n', target_trajectories);
trajlist = dir(fullfile(tempPath, 'traj*.nc'));

selected_data = struct('lon', cell(target_trajectories, 1), 'lat', cell(target_trajectories, 1));
trajectories_collected = 0;

selected_files = randsample(length(trajlist), min(50, length(trajlist)), false);

for i = selected_files'
    if trajectories_collected >= target_trajectories, break; end
    
    filename = fullfile(tempPath, trajlist(i).name);
    try
        lon = ncread(filename, 'lon');
        lat = ncread(filename, 'lat');
        
        if min(lon(:)) > 180
            lon = lon - 360;
        end
        
        n_particles = size(lon, 2);
        for p = 1:n_particles
            if trajectories_collected >= target_trajectories, break; end
            
            start_lon = lon(1, p);
            start_lat = lat(1, p);
            
            if start_lon < lon_min || start_lon > lon_max || ...
               start_lat < lat_min || start_lat > lat_max || ...
               isnan(start_lon) || isnan(start_lat)
                continue;
            end
            
            trajectories_collected = trajectories_collected + 1;
            selected_data(trajectories_collected).lon = lon(:, p);
            selected_data(trajectories_collected).lat = lat(:, p);
        end
        
    catch
        continue;
    end
end

selected_data = selected_data(1:trajectories_collected);
fprintf('Collected %d trajectories from bounded region\n', trajectories_collected);

%% Create plot
figure('Position', [100 100 1200 900]);
hold on;

% Plot bounded reef polygons in light gray with transparency
fprintf('Plotting reef polygons...\n');
for i = 1:size(Xs_bounded, 1)
    fill(Xs_bounded(i,1:4), Ys_bounded(i,1:4), [0.85 0.85 0.85], 'EdgeColor', [0.85 0.85 0.85], 'FaceAlpha', grid_alpha, 'LineWidth', 0.5);
end

% Plot trajectories and collect start points
fprintf('Plotting trajectories...\n');
colors = lines(10);
start_lons = zeros(trajectories_collected, 1);
start_lats = zeros(trajectories_collected, 1);

for i = 1:trajectories_collected
    valid = ~isnan(selected_data(i).lon) & ~isnan(selected_data(i).lat);
    color_with_alpha = [colors(mod(i-1,10)+1,:), trajectory_alpha];
    plot(selected_data(i).lon(valid), selected_data(i).lat(valid), ...
        'Color', color_with_alpha, 'LineWidth', 0.3);
    
    valid_idx = find(valid, 1);
    if ~isempty(valid_idx)
        start_lons(i) = selected_data(i).lon(valid_idx);
        start_lats(i) = selected_data(i).lat(valid_idx);
    end
end

% Plot starting locations as red stars
valid_starts = ~isnan(start_lons);
plot(start_lons(valid_starts), start_lats(valid_starts), 'r*', 'MarkerSize', star_size, 'LineWidth', 1);

xlabel('Longitude');
ylabel('Latitude');
xlim([min(Xs(:)), max(Xs(:))]);
ylim([min(Ys(:)), max(Ys(:))]);
axis equal;

fprintf('Plot complete!\n');