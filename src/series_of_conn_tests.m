clear; clc;
fprintf('=== SIMPLE POLYGON TEST ===\n');

%% Setup paths
projectPath = matlab.project.rootProject().RootFolder;
dataPath = fullfile(projectPath, 'data');
tempPath = fullfile(projectPath, 'temp');

%% Load reef centroids
centroids = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
unique_IDs = centroids(:,1);

% Build polygon coordinate arrays
Xs = [centroids(:,8) centroids(:,10) centroids(:,12) centroids(:,14) centroids(:,8)];
Ys = [centroids(:,9) centroids(:,11) centroids(:,13) centroids(:,15) centroids(:,9)];

fprintf('Loaded %d reef polygons.\n', size(centroids,1));

%% Pick one trajectory file
trajlist = dir(fullfile(tempPath, 'traj*.nc'));
if isempty(trajlist)
    error('No trajectory files found!');
end

% Randomly select one file
file_idx = randi(numel(trajlist));
file = fullfile(tempPath, trajlist(file_idx).name);
fprintf('Selected file: %s\n', trajlist(file_idx).name);

%% Load release file to get reef IDs
releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019_Q1.txt'));
release_reef_IDs = releasefile(:,1); % Column 1 has the unique reef IDs

%% Load data from this file
lon = ncread(file, 'lon');
lat = ncread(file, 'lat');
loc = ncread(file, 'location') + 1; % Line numbers in release file (1-based)

% Convert longitude from 0-360 to -180 to 180
lon(lon > 180) = lon(lon > 180) - 360;

% Get release locations (first time step)
release_lon = lon(1,:)';
release_lat = lat(1,:)';
release_line_nums = loc;

% Map: line number -> reef ID -> row index in centroids
reef_IDs_for_particles = release_reef_IDs(release_line_nums);
[~, release_loc] = ismember(reef_IDs_for_particles, unique_IDs);

fprintf('Found %d particles.\n', numel(release_lon));
fprintf('Unique source locations: %d\n', numel(unique(release_loc)));

%% Plot
figure('Position', [100 100 1000 800]);
hold on;

% Get active locations
active_locs = unique(release_loc);

% Plot only the active polygons that have releases (in blue)
for i = 1:numel(active_locs)
    idx = active_locs(i);
    if idx >= 1 && idx <= size(Xs,1)
        plot(Xs(idx,:), Ys(idx,:), 'b-', 'LineWidth', 2);
    end
end

% Plot release points (red dots)
scatter(release_lon, release_lat, 30, 'r', 'filled', 'MarkerEdgeColor', 'k');

xlabel('Longitude');
ylabel('Latitude');
title(sprintf('Release Locations from %s', trajlist(file_idx).name), 'Interpreter', 'none');
legend({'Source polygons', 'Release points'}, 'Location', 'best');
grid on;

% Set axis limits with buffer BEFORE axis equal
lon_range = [min(release_lon), max(release_lon)];
lat_range = [min(release_lat), max(release_lat)];
lon_buffer = 0.2 * diff(lon_range);
lat_buffer = 0.2 * diff(lat_range);
xlim([lon_range(1) - lon_buffer, lon_range(2) + lon_buffer]);
ylim([lat_range(1) - lat_buffer, lat_range(2) + lat_buffer]);

axis equal;

%% Diagnostic info
fprintf('\n=== DIAGNOSTIC INFO ===\n');
fprintf('Polygon longitude range: [%.2f, %.2f]\n', min(Xs(:)), max(Xs(:)));
fprintf('Release longitude range: [%.2f, %.2f]\n', min(release_lon), max(release_lon));
fprintf('Polygon latitude range: [%.2f, %.2f]\n', min(Ys(:)), max(Ys(:)));
fprintf('Release latitude range: [%.2f, %.2f]\n', min(release_lat), max(release_lat));

fprintf('\nSample polygon coords (first location):\n');
fprintf('  X: [%.2f, %.2f, %.2f, %.2f]\n', Xs(1,1), Xs(1,2), Xs(1,3), Xs(1,4));
fprintf('  Y: [%.2f, %.2f, %.2f, %.2f]\n', Ys(1,1), Ys(1,2), Ys(1,3), Ys(1,4));

fprintf('\nSample release coords (first 3 particles):\n');
for i = 1:min(3, numel(release_lon))
    fprintf('  Particle %d: lon=%.2f, lat=%.2f, source_loc=%d\n', ...
        i, release_lon(i), release_lat(i), release_loc(i));
end

fprintf('\n=== TEST COMPLETE ===\n');





















clear; clc;
fprintf('=== SIMPLE CONNECTIVITY TEST ===\n');

%% Setup paths
projectPath = matlab.project.rootProject().RootFolder;
dataPath = fullfile(projectPath, 'data');
tempPath = fullfile(projectPath, 'temp');

%% Load reef centroids and build polygons
centroids = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
unique_IDs = centroids(:,1);
Xs = [centroids(:,8) centroids(:,10) centroids(:,12) centroids(:,14) centroids(:,8)];
Ys = [centroids(:,9) centroids(:,11) centroids(:,13) centroids(:,15) centroids(:,9)];
n_locations = size(centroids,1);

fprintf('Loaded %d reef polygons.\n', n_locations);

%% Load release file
releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019_Q1.txt'));
release_reef_IDs = releasefile(:,1);

%% Pick one trajectory file
trajlist = dir(fullfile(tempPath, 'traj*.nc'));
file_idx = randi(numel(trajlist));
file = fullfile(tempPath, trajlist(file_idx).name);
fprintf('Using file: %s\n', trajlist(file_idx).name);

%% Load trajectory data
lon = ncread(file, 'lon');
lat = ncread(file, 'lat');
loc = ncread(file, 'location') + 1;

% Convert longitude
lon(lon > 180) = lon(lon > 180) - 360;

% Get endpoints (last time step)
end_lon = lon(end,:)';
end_lat = lat(end,:)';

% Map release locations
reef_IDs_for_particles = release_reef_IDs(loc);
[~, source_idx] = ismember(reef_IDs_for_particles, unique_IDs);

fprintf('Processing %d particles.\n', numel(end_lon));

%% Find destination polygon for each endpoint
dest_idx = zeros(numel(end_lon), 1);

for j = 1:n_locations
    in = inpolygon(end_lon, end_lat, Xs(j,:), Ys(j,:));
    dest_idx(in) = j;
end

%% Build connectivity matrix
valid = dest_idx > 0 & source_idx > 0;
src = source_idx(valid);
dst = dest_idx(valid);

fprintf('Valid connections: %d / %d particles (%.1f%%)\n', ...
    sum(valid), numel(valid), 100*sum(valid)/numel(valid));

% Simple connectivity matrix (counts)
ConnMatrix = zeros(n_locations, n_locations);
for i = 1:numel(src)
    ConnMatrix(src(i), dst(i)) = ConnMatrix(src(i), dst(i)) + 1;
end

% Get active sources and destinations
active_sources = unique(src);
active_dests = unique(dst);

fprintf('\nActive sources: %d\n', numel(active_sources));
fprintf('Active destinations: %d\n', numel(active_dests));
fprintf('Self-recruitment: %d particles\n', sum(src == dst));

%% Simple visualization
figure('Position', [100 100 800 600]);

% Extract submatrix for active locations only
active_locs = unique([active_sources; active_dests]);
sub_conn = ConnMatrix(active_locs, active_locs);

imagesc(sub_conn);
colorbar;
xlabel('Destination Polygon');
ylabel('Source Polygon');
title(sprintf('Connectivity Matrix (%s)', trajlist(file_idx).name), 'Interpreter', 'none');
axis square;

fprintf('\n=== TEST COMPLETE ===\n');

%% Calculate network metrics
fprintf('\n=== NETWORK METRICS ===\n');

% Out-degree: how many particles each source sent out
out_degree = sum(ConnMatrix, 2);

% In-degree: how many particles each destination received
in_degree = sum(ConnMatrix, 1)';

% Self-recruitment rate
self_recruitment = diag(ConnMatrix) ./ max(out_degree, 1);

% Total connectivity (out + in)
total_connectivity = out_degree + in_degree;

% Report for active locations
for i = 1:min(5, numel(active_locs))
    loc_idx = active_locs(i);
    fprintf('Location %d (ID=%d): Out=%d, In=%d, Self=%.2f%%\n', ...
        loc_idx, unique_IDs(loc_idx), out_degree(loc_idx), ...
        in_degree(loc_idx), 100*self_recruitment(loc_idx));
end

%% Spatial plot of connectivity
figure('Position', [100 100 1000 800]);
hold on;

% Plot only active polygons
for i = 1:numel(active_locs)
    idx = active_locs(i);
    plot(Xs(idx,:), Ys(idx,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
end

% Color polygons by total connectivity (using centroids)
centroids_x = mean(Xs(:,1:4), 2);
centroids_y = mean(Ys(:,1:4), 2);

% Normalize connectivity for sizing
conn_normalized = total_connectivity(active_locs);
conn_normalized = conn_normalized / max(conn_normalized) * 200 + 20;

scatter(centroids_x(active_locs), centroids_y(active_locs), ...
    conn_normalized, total_connectivity(active_locs), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

colorbar;
xlabel('Longitude');
ylabel('Latitude');
title('Reef Connectivity (size & color = total in+out degree)');
axis equal tight;
grid on;

fprintf('\n=== ANALYSIS COMPLETE ===\n');























clear; clc;
fprintf('=== SINGLE EVENT CONNECTIVITY TEST ===\n');

%% Setup paths
projectPath = matlab.project.rootProject().RootFolder;
dataPath = fullfile(projectPath, 'data');
tempPath = fullfile(projectPath, 'temp');

%% Load reef centroids and build polygons
centroids = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
unique_IDs = centroids(:,1);
Xs = [centroids(:,8) centroids(:,10) centroids(:,12) centroids(:,14) centroids(:,8)];
Ys = [centroids(:,9) centroids(:,11) centroids(:,13) centroids(:,15) centroids(:,9)];
n_locations = size(centroids,1);

fprintf('Loaded %d reef polygons.\n', n_locations);

%% Load release file
releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019_Q1.txt'));
release_reef_IDs = releasefile(:,1);

%% Identify all release events across trajectory files
fprintf('\nScanning trajectory files...\n');
trajlist = dir(fullfile(tempPath, 'traj*.nc'));
fprintf('Found %d trajectory files.\n', numel(trajlist));

all_release_dates = [];
for i = 1:numel(trajlist)
    filename = fullfile(tempPath, trajlist(i).name);
    try
        releasedate = ncread(filename, 'releasedate');
        all_release_dates = [all_release_dates; releasedate(:)];
    catch
        warning('Could not read %s', trajlist(i).name);
    end
end

unique_release_dates = unique(all_release_dates);
n_events = length(unique_release_dates);

fprintf('Found %d unique release events.\n', n_events);
for i = 1:min(5, n_events)
    fprintf('  Event %d: Julian date %.2f\n', i, unique_release_dates(i));
end

%% Pick one event to analyze
event_idx = randi(n_events);
target_date = unique_release_dates(event_idx);
fprintf('\nSelected Event %d (Julian date %.2f)\n', event_idx, target_date);

%% Collect all particles from this event across all files
fprintf('\nCollecting particles from all files for this event...\n');

all_end_lon = [];
all_end_lat = [];
all_source_idx = [];

for i = 1:numel(trajlist)
    if mod(i, 10) == 0
        fprintf('  Processing file %d/%d...\n', i, numel(trajlist));
    end
    
    filename = fullfile(tempPath, trajlist(i).name);
    
    try
        % Read data
        lon = ncread(filename, 'lon');
        lat = ncread(filename, 'lat');
        location = ncread(filename, 'location') + 1;
        releasedate = ncread(filename, 'releasedate');
        
        % Filter to target event
        event_mask = abs(releasedate - target_date) < 0.01;
        if ~any(event_mask)
            continue;
        end
        
        % Convert longitude
        lon(lon > 180) = lon(lon > 180) - 360;
        
        % Get endpoints for this event's particles
        end_lon = lon(end, event_mask);
        end_lat = lat(end, event_mask);
        loc = location(event_mask);
        
        % Map to source indices
        reef_IDs = release_reef_IDs(loc);
        [~, src_idx] = ismember(reef_IDs, unique_IDs);
        
        % Accumulate
        all_end_lon = [all_end_lon; end_lon(:)];
        all_end_lat = [all_end_lat; end_lat(:)];
        all_source_idx = [all_source_idx; src_idx(:)];
        
    catch ME
        warning('Error processing %s: %s', trajlist(i).name, ME.message);
    end
end

fprintf('Collected %d total particles for this event.\n', numel(all_end_lon));

%% Find destination polygons for all endpoints
fprintf('\nFinding destinations...\n');

dest_idx = zeros(numel(all_end_lon), 1);

for j = 1:n_locations
    if mod(j, 1000) == 0
        fprintf('  Checked %d/%d polygons...\n', j, n_locations);
    end
    in = inpolygon(all_end_lon, all_end_lat, Xs(j,:), Ys(j,:));
    dest_idx(in) = j;
end

%% Build connectivity matrix
valid = dest_idx > 0 & all_source_idx > 0;
src = all_source_idx(valid);
dst = dest_idx(valid);

fprintf('\nValid connections: %d / %d particles (%.1f%%)\n', ...
    sum(valid), numel(valid), 100*sum(valid)/numel(valid));

ConnMatrix = zeros(n_locations, n_locations);
for i = 1:numel(src)
    ConnMatrix(src(i), dst(i)) = ConnMatrix(src(i), dst(i)) + 1;
end

active_sources = unique(src);
active_dests = unique(dst);

fprintf('Active sources: %d\n', numel(active_sources));
fprintf('Active destinations: %d\n', numel(active_dests));
fprintf('Self-recruitment: %d particles (%.1f%%)\n', sum(src == dst), ...
    100*sum(src == dst)/numel(src));

%% Network metrics
out_degree = sum(ConnMatrix, 2);
in_degree = sum(ConnMatrix, 1)';
total_connectivity = out_degree + in_degree;

fprintf('\nTop 5 most connected locations:\n');
[~, sorted_idx] = sort(total_connectivity, 'descend');
for i = 1:min(5, nnz(total_connectivity))
    idx = sorted_idx(i);
    fprintf('  Location %d (ID=%d): Out=%d, In=%d, Total=%d\n', ...
        idx, unique_IDs(idx), out_degree(idx), in_degree(idx), ...
        total_connectivity(idx));
end

%% Visualization 1: Connectivity matrix
figure('Position', [100 100 800 600]);
active_locs = unique([active_sources; active_dests]);
sub_conn = ConnMatrix(active_locs, active_locs);
imagesc(sub_conn);
colorbar;
xlabel('Destination Polygon');
ylabel('Source Polygon');
title(sprintf('Connectivity Matrix - Event %d', event_idx));
axis square;

%% Visualization 2: Spatial connectivity
figure('Position', [100 100 1000 800]);
hold on;

% Plot active polygons
for i = 1:numel(active_locs)
    idx = active_locs(i);
    plot(Xs(idx,:), Ys(idx,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
end

% Plot connectivity at centroids
centroids_x = mean(Xs(:,1:4), 2);
centroids_y = mean(Ys(:,1:4), 2);

conn_normalized = total_connectivity(active_locs);
conn_normalized = conn_normalized / max(conn_normalized) * 200 + 20;

scatter(centroids_x(active_locs), centroids_y(active_locs), ...
    conn_normalized, total_connectivity(active_locs), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

colorbar;
xlabel('Longitude');
ylabel('Latitude');
title(sprintf('Reef Connectivity - Event %d (size & color = total degree)', event_idx));
axis equal tight;
grid on;

fprintf('\n=== ANALYSIS COMPLETE ===\n');



















clear; clc;
fprintf('=== PASSING-THROUGH CONNECTIVITY TEST ===\n');

%% Setup paths
projectPath = matlab.project.rootProject().RootFolder;
dataPath = fullfile(projectPath, 'data');
tempPath = fullfile(projectPath, 'temp');

%% Load reef centroids and build polygons
centroids = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
unique_IDs = centroids(:,1);
Xs = [centroids(:,8) centroids(:,10) centroids(:,12) centroids(:,14) centroids(:,8)];
Ys = [centroids(:,9) centroids(:,11) centroids(:,13) centroids(:,15) centroids(:,9)];
n_locations = size(centroids,1);

fprintf('Loaded %d reef polygons.\n', n_locations);

%% Load release file
releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019_Q1.txt'));
release_reef_IDs = releasefile(:,1);

%% Pick one trajectory file for testing
trajlist = dir(fullfile(tempPath, 'traj*.nc'));
file_idx = randi(numel(trajlist));
file = fullfile(tempPath, trajlist(file_idx).name);
fprintf('Using file: %s\n', trajlist(file_idx).name);

%% Load ONE trajectory file's data
lon = ncread(file, 'lon');
lat = ncread(file, 'lat');
location = ncread(file, 'location') + 1;

% Convert longitude
lon(lon > 180) = lon(lon > 180) - 360;

% Map release locations to source indices
reef_IDs_for_particles = release_reef_IDs(location);
[~, source_idx] = ismember(reef_IDs_for_particles, unique_IDs);

n_particles = size(lon, 2);
n_timesteps = size(lon, 1);

fprintf('Processing %d particles with %d timesteps each.\n', n_particles, n_timesteps);

%% Initialize connectivity matrix
ConnMatrix = zeros(n_locations, n_locations);

%% Simple approach: Check each timestep against all polygons
fprintf('\nChecking particle positions at each timestep...\n');

for t = 1:n_timesteps
    if mod(t, 10) == 0
        fprintf('  Timestep %d/%d...\n', t, n_timesteps);
    end
    
    % Get all particle positions at this timestep
    lon_t = lon(t, :)';
    lat_t = lat(t, :)';
    
    % Remove NaN positions (particles that have exited)
    valid = ~isnan(lon_t) & ~isnan(lat_t);
    lon_t = lon_t(valid);
    lat_t = lat_t(valid);
    src_t = source_idx(valid);
    
    if isempty(lon_t)
        continue;
    end
    
    % Find which polygon each particle is in at this timestep
    dest_idx = zeros(numel(lon_t), 1);
    
    for j = 1:n_locations
        in = inpolygon(lon_t, lat_t, Xs(j,:), Ys(j,:));
        dest_idx(in) = j;
    end
    
    % Record connections (source -> destination at this timestep)
    valid_conn = dest_idx > 0 & src_t > 0;
    for i = find(valid_conn)'
        ConnMatrix(src_t(i), dest_idx(i)) = ConnMatrix(src_t(i), dest_idx(i)) + 1;
    end
end

%% Summary
fprintf('\n=== RESULTS ===\n');
fprintf('Total connections recorded: %d\n', nnz(ConnMatrix));

% Get active sources and destinations
[src_indices, dst_indices] = find(ConnMatrix);
active_sources = unique(src_indices);
active_dests = unique(dst_indices);

fprintf('Active sources: %d\n', numel(active_sources));
fprintf('Active destinations: %d\n', numel(active_dests));

% Self-recruitment (diagonal)
self_recruit = sum(diag(ConnMatrix));
fprintf('Self-recruitment counts: %d\n', self_recruit);

%% Network metrics
out_degree = sum(ConnMatrix, 2);
in_degree = sum(ConnMatrix, 1)';
total_connectivity = out_degree + in_degree;

%% Visualizations
% 1. Connectivity matrix heatmap
figure('Position', [100 100 800 600]);
active_locs = unique([active_sources; active_dests]);
sub_conn = ConnMatrix(active_locs, active_locs);
imagesc(log10(sub_conn + 1));  % Log scale for visibility
colorbar;
xlabel('Destination Polygon');
ylabel('Source Polygon');
title('Passing-Through Connectivity (log10 scale)');
axis square;

% 2. Spatial plot
figure('Position', [100 100 1000 800]);
hold on;

% Plot active polygons
for i = 1:numel(active_locs)
    idx = active_locs(i);
    plot(Xs(idx,:), Ys(idx,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
end

% Plot connectivity at centroids
centroids_x = mean(Xs(:,1:4), 2);
centroids_y = mean(Ys(:,1:4), 2);

conn_normalized = total_connectivity(active_locs);
conn_normalized = conn_normalized / max(conn_normalized) * 200 + 20;

scatter(centroids_x(active_locs), centroids_y(active_locs), ...
    conn_normalized, total_connectivity(active_locs), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

colorbar;
xlabel('Longitude');
ylabel('Latitude');
title('Passing-Through Connectivity (size & color = total degree)');
axis equal tight;
grid on;

fprintf('\n=== COMPLETE ===\n');
fprintf('Note: This counts every timestep a particle is inside each polygon.\n');
fprintf('Higher counts = more time spent in that source->dest connection.\n');



















clear; clc;
fprintf('=== OPTIMIZED PASSING-THROUGH CONNECTIVITY ===\n');

%% Setup paths
projectPath = matlab.project.rootProject().RootFolder;
dataPath = fullfile(projectPath, 'data');
tempPath = fullfile(projectPath, 'temp');

%% Load reef centroids and build polygons
centroids = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
unique_IDs = centroids(:,1);
Xs = [centroids(:,8) centroids(:,10) centroids(:,12) centroids(:,14) centroids(:,8)];
Ys = [centroids(:,9) centroids(:,11) centroids(:,13) centroids(:,15) centroids(:,9)];
n_locations = size(centroids,1);

fprintf('Loaded %d reef polygons.\n', n_locations);

%% Load release file
releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019_Q1.txt'));
release_reef_IDs = releasefile(:,1);

%% Identify all release events across trajectory files
fprintf('\nScanning trajectory files for release events...\n');
trajlist = dir(fullfile(tempPath, 'traj*.nc'));
fprintf('Found %d trajectory files.\n', numel(trajlist));

all_release_dates = [];
for i = 1:numel(trajlist)
    filename = fullfile(tempPath, trajlist(i).name);
    try
        releasedate = ncread(filename, 'releasedate');
        all_release_dates = [all_release_dates; releasedate(:)];
    catch
        warning('Could not read %s', trajlist(i).name);
    end
end

unique_release_dates = unique(all_release_dates);
n_events = length(unique_release_dates);

fprintf('Found %d unique release events.\n', n_events);

%% Pick one event to analyze
event_idx = randi(n_events);
target_date = unique_release_dates(event_idx);
fprintf('\nSelected Event %d (Julian date %.2f)\n', event_idx, target_date);

%% Collect all particles from this event across ALL files
fprintf('\nCollecting particles from all files for this event...\n');

all_lon = [];
all_lat = [];
all_source_idx = [];

for i = 1:numel(trajlist)
    if mod(i, 10) == 0
        fprintf('  Loading file %d/%d...\n', i, numel(trajlist));
    end
    
    filename = fullfile(tempPath, trajlist(i).name);
    
    try
        % Read data
        lon = ncread(filename, 'lon');
        lat = ncread(filename, 'lat');
        location = ncread(filename, 'location') + 1;
        releasedate = ncread(filename, 'releasedate');
        
        % Filter to target event
        event_mask = abs(releasedate - target_date) < 0.01;
        if ~any(event_mask)
            continue;
        end
        
        % Convert longitude
        lon(lon > 180) = lon(lon > 180) - 360;
        
        % Get particles for this event
        lon_event = lon(:, event_mask);
        lat_event = lat(:, event_mask);
        loc_event = location(event_mask);
        
        % Map to source indices
        reef_IDs = release_reef_IDs(loc_event);
        [~, src_idx] = ismember(reef_IDs, unique_IDs);
        
        % Flatten and accumulate
        lon_flat = lon_event(:);
        lat_flat = lat_event(:);
        n_timesteps_file = size(lon_event, 1);
        src_flat = repmat(src_idx', n_timesteps_file, 1);
        src_flat = src_flat(:);
        
        all_lon = [all_lon; lon_flat];
        all_lat = [all_lat; lat_flat];
        all_source_idx = [all_source_idx; src_flat];
        
    catch ME
        warning('Error processing %s: %s', trajlist(i).name, ME.message);
    end
end

fprintf('Collected %d total particle positions for this event.\n', numel(all_lon));

%% Remove NaN positions
valid = ~isnan(all_lon) & ~isnan(all_lat) & all_source_idx > 0;
lon_all = all_lon(valid);
lat_all = all_lat(valid);
source_all = all_source_idx(valid);

fprintf('Valid positions after filtering: %d\n', numel(lon_all));

%% OPTIMIZATION 1: Pre-compute bounding boxes for quick rejection
fprintf('\nPre-computing polygon bounding boxes...\n');
bbox_xmin = min(Xs(:,1:4), [], 2);
bbox_xmax = max(Xs(:,1:4), [], 2);
bbox_ymin = min(Ys(:,1:4), [], 2);
bbox_ymax = max(Ys(:,1:4), [], 2);

%% OPTIMIZATION 2: Data already prepared
fprintf('Data already flattened during collection.\n');
fprintf('Total valid positions to check: %d\n', numel(lon_all));

%% Check for parallel pool
pool = gcp('nocreate');
if isempty(pool)
    fprintf('Starting parallel pool...\n');
    pool = parpool('local');
else
    fprintf('Using existing parallel pool with %d workers.\n', pool.NumWorkers);
end

%% Check if INPOLY2 is available (with helper functions)
use_inpoly = false;
if exist('inpoly2', 'file') == 2
    try
        % Test if inpoly2 works (has all required helper functions)
        test_result = inpoly2([0 0], [0 0; 1 0; 1 1; 0 1]);
        use_inpoly = true;
        fprintf('INPOLY2 detected and working - will use fast algorithm.\n');
    catch
        fprintf('INPOLY2 found but missing helper functions (inpoly2_mat.m).\n');
        fprintf('Download complete package: https://github.com/dengwirda/inpoly\n');
        fprintf('Falling back to standard inpolygon.\n');
    end
else
    fprintf('INPOLY2 not found - using standard inpolygon.\n');
    fprintf('For faster performance, download: https://github.com/dengwirda/inpoly\n');
end

%% OPTIMIZATION 3: Parallel processing with bounding box pre-filtering
fprintf('\nProcessing with optimized parallel method...\n');
tic;

% Process in chunks to manage memory
chunk_size = 50000;
n_chunks = ceil(numel(lon_all) / chunk_size);

% Pre-allocate cell array for parallel results
chunk_connections = cell(n_chunks, 1);

% Copy variables for parfor (ensures proper broadcasting)
Xs_par = Xs;
Ys_par = Ys;
bbox_xmin_par = bbox_xmin;
bbox_xmax_par = bbox_xmax;
bbox_ymin_par = bbox_ymin;
bbox_ymax_par = bbox_ymax;
n_locations_par = n_locations;
use_inpoly_par = use_inpoly;

fprintf('Starting parallel processing of %d chunks...\n', n_chunks);

% PARALLEL LOOP over chunks
parfor chunk = 1:n_chunks
    % Get chunk indices
    idx_start = (chunk-1)*chunk_size + 1;
    idx_end = min(chunk*chunk_size, numel(lon_all));
    
    lon_chunk = lon_all(idx_start:idx_end);
    lat_chunk = lat_all(idx_start:idx_end);
    src_chunk = source_all(idx_start:idx_end);
    
    n_pts = numel(lon_chunk);
    dest_chunk = zeros(n_pts, 1);
    
    % Loop over polygons with bounding box pre-filtering
    for j = 1:n_locations_par
        % FAST: Bounding box check (vectorized)
        in_bbox = lon_chunk >= bbox_xmin_par(j) & lon_chunk <= bbox_xmax_par(j) & ...
                  lat_chunk >= bbox_ymin_par(j) & lat_chunk <= bbox_ymax_par(j);
        
        if ~any(in_bbox)
            continue;  % Skip this polygon entirely
        end
        
        % Only test points that passed bounding box
        candidates = find(in_bbox);
        
        if use_inpoly_par
            % INPOLY2: Fast algorithm with O((N+M)*log(N)) complexity
            in_poly = inpoly2([lon_chunk(candidates), lat_chunk(candidates)], ...
                             [Xs_par(j,1:4)', Ys_par(j,1:4)']);
        else
            % Standard MATLAB inpolygon
            in_poly = inpolygon(lon_chunk(candidates), lat_chunk(candidates), ...
                               Xs_par(j,:), Ys_par(j,:));
        end
        
        % Mark destinations
        dest_chunk(candidates(in_poly)) = j;
    end
    
    % Store connections as [source, destination] pairs
    valid_conn = dest_chunk > 0;
    chunk_connections{chunk} = [src_chunk(valid_conn), dest_chunk(valid_conn)];
end

fprintf('Parallel processing complete. Aggregating results...\n');
ConnMatrix = zeros(n_locations, n_locations);

for chunk = 1:n_chunks
    connections = chunk_connections{chunk};
    for i = 1:size(connections, 1)
        ConnMatrix(connections(i,1), connections(i,2)) = ...
            ConnMatrix(connections(i,1), connections(i,2)) + 1;
    end
end

elapsed = toc;
fprintf('Processing complete in %.1f seconds.\n', elapsed);

%% Normalize connectivity matrix (row normalization)
fprintf('\nNormalizing connectivity matrix...\n');

% Raw counts matrix
ConnMatrix_raw = ConnMatrix;

% Row normalization: each row sums to 1
% Interpretation: fraction of particle-time from source i spent in destination j
row_sums = sum(ConnMatrix, 2);
row_sums(row_sums == 0) = 1;  % Avoid division by zero for sources with no connections

ConnMatrix_normalized = ConnMatrix ./ row_sums;

fprintf('Normalization complete.\n');
fprintf('  Each row now sums to 1 (or 0 if no connections from that source).\n');

%% Summary
fprintf('\n=== RESULTS ===\n');
fprintf('RAW COUNTS:\n');
fprintf('  Total connections recorded: %d\n', nnz(ConnMatrix_raw));

[src_indices, dst_indices] = find(ConnMatrix_raw);
active_sources = unique(src_indices);
active_dests = unique(dst_indices);

fprintf('  Active sources: %d\n', numel(active_sources));
fprintf('  Active destinations: %d\n', numel(active_dests));
fprintf('  Self-recruitment counts: %d\n', sum(diag(ConnMatrix_raw)));

fprintf('\nNORMALIZED MATRIX:\n');
fprintf('  Non-zero entries: %d\n', nnz(ConnMatrix_normalized));
fprintf('  Self-recruitment probabilities range: [%.4f, %.4f]\n', ...
    min(diag(ConnMatrix_normalized(active_sources, active_sources))), ...
    max(diag(ConnMatrix_normalized(active_sources, active_sources))));

%% Network metrics (using raw counts)
out_degree = sum(ConnMatrix_raw, 2);
in_degree = sum(ConnMatrix_raw, 1)';
total_connectivity = out_degree + in_degree;

fprintf('\nTop 5 most connected locations:\n');
[~, sorted_idx] = sort(total_connectivity, 'descend');
for i = 1:min(5, nnz(total_connectivity))
    idx = sorted_idx(i);
    fprintf('  Location %d (ID=%d): Out=%d, In=%d, Total=%d\n', ...
        idx, unique_IDs(idx), out_degree(idx), in_degree(idx), ...
        total_connectivity(idx));
end

%% Visualizations
% 1. Connectivity matrix heatmap (NORMALIZED)
figure('Position', [100 100 800 600]);
active_locs = unique([active_sources; active_dests]);
sub_conn = ConnMatrix_normalized(active_locs, active_locs);
imagesc(sub_conn);
colorbar;
caxis([0 1]);  % Probabilities from 0 to 1
xlabel('Destination Polygon');
ylabel('Source Polygon');
title('Normalized Passing-Through Connectivity (row-normalized probabilities)');
axis square;

% 2. Raw counts heatmap for comparison
figure('Position', [920 100 800 600]);
sub_conn_raw = ConnMatrix_raw(active_locs, active_locs);
imagesc(log10(sub_conn_raw + 1));
colorbar;
xlabel('Destination Polygon');
ylabel('Source Polygon');
title('Raw Connectivity Counts (log10 scale)');
axis square;

% 2. Raw counts heatmap for comparison
figure('Position', [920 100 800 600]);
sub_conn_raw = ConnMatrix_raw(active_locs, active_locs);
imagesc(log10(sub_conn_raw + 1));
colorbar;
xlabel('Destination Polygon');
ylabel('Source Polygon');
title('Raw Connectivity Counts (log10 scale)');
axis square;

% 3. Spatial plot (using raw counts for node size)
figure('Position', [100 750 1000 800]);
hold on;

for i = 1:numel(active_locs)
    idx = active_locs(i);
    plot(Xs(idx,:), Ys(idx,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
end

centroids_x = mean(Xs(:,1:4), 2);
centroids_y = mean(Ys(:,1:4), 2);

conn_normalized = total_connectivity(active_locs);
conn_normalized = conn_normalized / max(conn_normalized) * 200 + 20;

scatter(centroids_x(active_locs), centroids_y(active_locs), ...
    conn_normalized, total_connectivity(active_locs), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

colorbar;
xlabel('Longitude');
ylabel('Latitude');
title('Reef Connectivity - Raw Counts (size & color = total degree)');
axis equal tight;
grid on;

% 4. Spatial plot (using normalized probabilities)
figure('Position', [920 750 1000 800]);
hold on;

for i = 1:numel(active_locs)
    idx = active_locs(i);
    plot(Xs(idx,:), Ys(idx,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
end

% Calculate normalized connectivity metrics
out_degree_norm = sum(ConnMatrix_normalized, 2);  % Sum of outgoing probabilities
in_degree_norm = sum(ConnMatrix_normalized, 1)';  % Sum of incoming probabilities
total_connectivity_norm = out_degree_norm + in_degree_norm;

conn_size_norm = total_connectivity_norm(active_locs);
conn_size_norm = conn_size_norm / max(conn_size_norm) * 200 + 20;

scatter(centroids_x(active_locs), centroids_y(active_locs), ...
    conn_size_norm, total_connectivity_norm(active_locs), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

colorbar;
xlabel('Longitude');
ylabel('Latitude');
title('Reef Connectivity - Normalized (size & color = total probability)');
axis equal tight;
grid on;

fprintf('\n=== COMPLETE ===\n');
fprintf('Optimizations applied:\n');
fprintf('  1. Bounding box pre-filtering\n');
fprintf('  2. Vectorized trajectory data\n');
fprintf('  3. Parallel processing (%d workers)\n', pool.NumWorkers);
if use_inpoly
    fprintf('  4. INPOLY2 fast algorithm (O((N+M)*log(N)))\n');
else
    fprintf('  4. Standard inpolygon (consider installing INPOLY for more speed)\n');
end
















%%% SEEMS TO BE WORKING WELL IN PARALLEL WITH LOTS OF DATA %%

clear; clc;
fprintf('=== OPTIMIZED PASSING-THROUGH CONNECTIVITY ===\n');

%% Setup paths
projectPath = matlab.project.rootProject().RootFolder;
dataPath = fullfile(projectPath, 'data');
tempPath = fullfile(projectPath, 'temp');

%% Load reef centroids and build polygons
centroids = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
unique_IDs = centroids(:,1);
Xs = [centroids(:,8) centroids(:,10) centroids(:,12) centroids(:,14) centroids(:,8)];
Ys = [centroids(:,9) centroids(:,11) centroids(:,13) centroids(:,15) centroids(:,9)];
n_locations = size(centroids,1);

fprintf('Loaded %d reef polygons.\n', n_locations);

%% Load release file
releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019_Q1.txt'));
release_reef_IDs = releasefile(:,1);

%% Identify all release events across trajectory files
fprintf('\nScanning trajectory files for release events...\n');
trajlist = dir(fullfile(tempPath, 'traj*.nc'));
fprintf('Found %d trajectory files.\n', numel(trajlist));

% Use cell array to avoid growing array in loop
release_date_cells = cell(numel(trajlist), 1);

for i = 1:numel(trajlist)
    filename = fullfile(tempPath, trajlist(i).name);
    try
        releasedate = ncread(filename, 'releasedate');
        release_date_cells{i} = releasedate(:);
    catch
        warning('Could not read %s', trajlist(i).name);
    end
end

% Concatenate all at once
all_release_dates = vertcat(release_date_cells{:});
unique_release_dates = unique(all_release_dates);
n_events = length(unique_release_dates);

fprintf('Found %d unique release events.\n', n_events);

%% Pick one event to analyze
event_idx = randi(n_events);
target_date = unique_release_dates(event_idx);
fprintf('\nSelected Event %d (Julian date %.2f)\n', event_idx, target_date);

%% Collect all particles from this event across ALL files
fprintf('\nCollecting particles from all files for this event...\n');

% Use cell arrays to avoid growing arrays in loop
lon_cells = cell(numel(trajlist), 1);
lat_cells = cell(numel(trajlist), 1);
src_cells = cell(numel(trajlist), 1);

for i = 1:numel(trajlist)
    if mod(i, 10) == 0
        fprintf('  Loading file %d/%d...\n', i, numel(trajlist));
    end
    
    filename = fullfile(tempPath, trajlist(i).name);
    
    try
        % Read data
        lon = ncread(filename, 'lon');
        lat = ncread(filename, 'lat');
        location = ncread(filename, 'location') + 1;
        releasedate = ncread(filename, 'releasedate');
        
        % Filter to target event
        event_mask = abs(releasedate - target_date) < 0.01;
        if ~any(event_mask)
            continue;
        end
        
        % Convert longitude
        lon(lon > 180) = lon(lon > 180) - 360;
        
        % Get particles for this event
        lon_event = lon(:, event_mask);
        lat_event = lat(:, event_mask);
        loc_event = location(event_mask);
        
        % Map to source indices
        reef_IDs = release_reef_IDs(loc_event);
        [~, src_idx] = ismember(reef_IDs, unique_IDs);
        
        % Flatten
        lon_flat = lon_event(:);
        lat_flat = lat_event(:);
        n_timesteps_file = size(lon_event, 1);
        src_flat = repmat(src_idx', n_timesteps_file, 1);
        src_flat = src_flat(:);
        
        % Store in cell array (NO GROWING!)
        lon_cells{i} = lon_flat;
        lat_cells{i} = lat_flat;
        src_cells{i} = src_flat;
        
    catch ME
        warning('Error processing %s: %s', trajlist(i).name, ME.message);
    end
end

% Concatenate all at once (MUCH faster than growing in loop)
fprintf('Concatenating data from all files...\n');
all_lon = vertcat(lon_cells{:});
all_lat = vertcat(lat_cells{:});
all_source_idx = vertcat(src_cells{:});

fprintf('Collected %d total particle positions for this event.\n', numel(all_lon));

%% Remove NaN positions
valid = ~isnan(all_lon) & ~isnan(all_lat) & all_source_idx > 0;
lon_all = all_lon(valid);
lat_all = all_lat(valid);
source_all = all_source_idx(valid);

fprintf('Valid positions after filtering: %d\n', numel(lon_all));

%% OPTIMIZATION 1: Pre-compute bounding boxes for quick rejection
fprintf('\nPre-computing polygon bounding boxes...\n');
bbox_xmin = min(Xs(:,1:4), [], 2);
bbox_xmax = max(Xs(:,1:4), [], 2);
bbox_ymin = min(Ys(:,1:4), [], 2);
bbox_ymax = max(Ys(:,1:4), [], 2);

%% OPTIMIZATION 2: Data already prepared
fprintf('Data already flattened during collection.\n');
fprintf('Total valid positions to check: %d\n', numel(lon_all));

%% Check for parallel pool
pool = gcp('nocreate');
if isempty(pool)
    fprintf('Starting parallel pool...\n');
    pool = parpool('local');
else
    fprintf('Using existing parallel pool with %d workers.\n', pool.NumWorkers);
end

%% Check if INPOLY2 is available (with helper functions)
use_inpoly = false;
if exist('inpoly2', 'file') == 2
    try
        % Test if inpoly2 works (has all required helper functions)
        test_result = inpoly2([0 0], [0 0; 1 0; 1 1; 0 1]);
        use_inpoly = true;
        fprintf('INPOLY2 detected and working - will use fast algorithm.\n');
    catch
        fprintf('INPOLY2 found but missing helper functions (inpoly2_mat.m).\n');
        fprintf('Download complete package: https://github.com/dengwirda/inpoly\n');
        fprintf('Falling back to standard inpolygon.\n');
    end
else
    fprintf('INPOLY2 not found - using standard inpolygon.\n');
    fprintf('For faster performance, download: https://github.com/dengwirda/inpoly\n');
end

%% OPTIMIZATION 3: Parallel processing with bounding box pre-filtering
fprintf('\nProcessing with optimized parallel method...\n');
tic;

% Process in chunks to manage memory
chunk_size = 50000;
n_chunks = ceil(numel(lon_all) / chunk_size);

% Pre-allocate cell array for parallel results
chunk_connections = cell(n_chunks, 1);

% Copy variables for parfor (ensures proper broadcasting)
Xs_par = Xs;
Ys_par = Ys;
bbox_xmin_par = bbox_xmin;
bbox_xmax_par = bbox_xmax;
bbox_ymin_par = bbox_ymin;
bbox_ymax_par = bbox_ymax;
n_locations_par = n_locations;
use_inpoly_par = use_inpoly;

fprintf('Starting parallel processing of %d chunks...\n', n_chunks);

% PARALLEL LOOP over chunks
parfor chunk = 1:n_chunks
    % Get chunk indices
    idx_start = (chunk-1)*chunk_size + 1;
    idx_end = min(chunk*chunk_size, numel(lon_all));
    
    % Slice data for this chunk
    lon_chunk = lon_all(idx_start:idx_end);
    lat_chunk = lat_all(idx_start:idx_end);
    src_chunk = source_all(idx_start:idx_end);
    
    n_pts = numel(lon_chunk);
    dest_chunk = zeros(n_pts, 1);
    
    % Loop over polygons with bounding box pre-filtering
    for j = 1:n_locations_par
        % FAST: Bounding box check (vectorized)
        in_bbox = lon_chunk >= bbox_xmin_par(j) & lon_chunk <= bbox_xmax_par(j) & ...
                  lat_chunk >= bbox_ymin_par(j) & lat_chunk <= bbox_ymax_par(j);
        
        if ~any(in_bbox)
            continue;  % Skip this polygon entirely
        end
        
        % Only test points that passed bounding box
        candidates = find(in_bbox);
        
        if use_inpoly_par
            % INPOLY2: Fast algorithm with O((N+M)*log(N)) complexity
            in_poly = inpoly2([lon_chunk(candidates), lat_chunk(candidates)], ...
                             [Xs_par(j,1:4)', Ys_par(j,1:4)']);
        else
            % Standard MATLAB inpolygon
            in_poly = inpolygon(lon_chunk(candidates), lat_chunk(candidates), ...
                               Xs_par(j,:), Ys_par(j,:));
        end
        
        % Mark destinations
        dest_chunk(candidates(in_poly)) = j;
    end
    
    % Store connections as [source, destination] pairs
    valid_conn = dest_chunk > 0;
    chunk_connections{chunk} = [src_chunk(valid_conn), dest_chunk(valid_conn)];
end

fprintf('Parallel processing complete. Aggregating results...\n');
ConnMatrix = zeros(n_locations, n_locations);

for chunk = 1:n_chunks
    connections = chunk_connections{chunk};
    for i = 1:size(connections, 1)
        ConnMatrix(connections(i,1), connections(i,2)) = ...
            ConnMatrix(connections(i,1), connections(i,2)) + 1;
    end
end

elapsed = toc;
fprintf('Processing complete in %.1f seconds.\n', elapsed);

%% Normalize connectivity matrix (row normalization)
fprintf('\nNormalizing connectivity matrix...\n');

% Raw counts matrix
ConnMatrix_raw = ConnMatrix;

% Row normalization: each row sums to 1
% Interpretation: fraction of particle-time from source i spent in destination j
row_sums = sum(ConnMatrix, 2);
row_sums(row_sums == 0) = 1;  % Avoid division by zero for sources with no connections

ConnMatrix_normalized = ConnMatrix ./ row_sums;

fprintf('Normalization complete.\n');
fprintf('  Each row now sums to 1 (or 0 if no connections from that source).\n');

%% Summary
fprintf('\n=== RESULTS ===\n');
fprintf('RAW COUNTS:\n');
fprintf('  Total connections recorded: %d\n', nnz(ConnMatrix_raw));

[src_indices, dst_indices] = find(ConnMatrix_raw);
active_sources = unique(src_indices);
active_dests = unique(dst_indices);

fprintf('  Active sources: %d\n', numel(active_sources));
fprintf('  Active destinations: %d\n', numel(active_dests));
fprintf('  Self-recruitment counts: %d\n', sum(diag(ConnMatrix_raw)));

fprintf('\nNORMALIZED MATRIX:\n');
fprintf('  Non-zero entries: %d\n', nnz(ConnMatrix_normalized));
fprintf('  Self-recruitment probabilities range: [%.4f, %.4f]\n', ...
    min(diag(ConnMatrix_normalized(active_sources, active_sources))), ...
    max(diag(ConnMatrix_normalized(active_sources, active_sources))));

%% Network metrics (using raw counts)
out_degree = sum(ConnMatrix_raw, 2);
in_degree = sum(ConnMatrix_raw, 1)';
total_connectivity = out_degree + in_degree;

fprintf('\nTop 5 most connected locations:\n');
[~, sorted_idx] = sort(total_connectivity, 'descend');
for i = 1:min(5, nnz(total_connectivity))
    idx = sorted_idx(i);
    fprintf('  Location %d (ID=%d): Out=%d, In=%d, Total=%d\n', ...
        idx, unique_IDs(idx), out_degree(idx), in_degree(idx), ...
        total_connectivity(idx));
end

%% Visualizations
% 1. Connectivity matrix heatmap (NORMALIZED)
figure('Position', [100 100 800 600]);
active_locs = unique([active_sources; active_dests]);
sub_conn = ConnMatrix_normalized(active_locs, active_locs);
imagesc(sub_conn);
colorbar;
caxis([0 1]);  % Probabilities from 0 to 1
xlabel('Destination Polygon');
ylabel('Source Polygon');
title('Normalized Passing-Through Connectivity (row-normalized probabilities)');
axis square;

% 2. Raw counts heatmap for comparison
figure('Position', [920 100 800 600]);
sub_conn_raw = ConnMatrix_raw(active_locs, active_locs);
imagesc(log10(sub_conn_raw + 1));
colorbar;
xlabel('Destination Polygon');
ylabel('Source Polygon');
title('Raw Connectivity Counts (log10 scale)');
axis square;

% 2. Raw counts heatmap for comparison
figure('Position', [920 100 800 600]);
sub_conn_raw = ConnMatrix_raw(active_locs, active_locs);
imagesc(log10(sub_conn_raw + 1));
colorbar;
xlabel('Destination Polygon');
ylabel('Source Polygon');
title('Raw Connectivity Counts (log10 scale)');
axis square;

% 3. Spatial plot (using raw counts for node size)
figure('Position', [100 750 1000 800]);
hold on;

for i = 1:numel(active_locs)
    idx = active_locs(i);
    plot(Xs(idx,:), Ys(idx,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
end

centroids_x = mean(Xs(:,1:4), 2);
centroids_y = mean(Ys(:,1:4), 2);

conn_normalized = total_connectivity(active_locs);
conn_normalized = conn_normalized / max(conn_normalized) * 200 + 20;

scatter(centroids_x(active_locs), centroids_y(active_locs), ...
    conn_normalized, total_connectivity(active_locs), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

colorbar;
xlabel('Longitude');
ylabel('Latitude');
title('Reef Connectivity - Raw Counts (size & color = total degree)');
axis equal tight;
grid on;

% 4. Spatial plot (using normalized probabilities)
figure('Position', [920 750 1000 800]);
hold on;

for i = 1:numel(active_locs)
    idx = active_locs(i);
    plot(Xs(idx,:), Ys(idx,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
end

% Calculate normalized connectivity metrics
out_degree_norm = sum(ConnMatrix_normalized, 2);  % Sum of outgoing probabilities
in_degree_norm = sum(ConnMatrix_normalized, 1)';  % Sum of incoming probabilities
total_connectivity_norm = out_degree_norm + in_degree_norm;

conn_size_norm = total_connectivity_norm(active_locs);
conn_size_norm = conn_size_norm / max(conn_size_norm) * 200 + 20;

scatter(centroids_x(active_locs), centroids_y(active_locs), ...
    conn_size_norm, total_connectivity_norm(active_locs), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

colorbar;
xlabel('Longitude');
ylabel('Latitude');
title('Reef Connectivity - Normalized (size & color = total probability)');
axis equal tight;
grid on;

fprintf('\n=== COMPLETE ===\n');
fprintf('Optimizations applied:\n');
fprintf('  1. Bounding box pre-filtering\n');
fprintf('  2. Vectorized trajectory data\n');
fprintf('  3. Parallel processing (%d workers)\n', pool.NumWorkers);
if use_inpoly
    fprintf('  4. INPOLY2 fast algorithm (O((N+M)*log(N)))\n');
else
    fprintf('  4. Standard inpolygon (consider installing INPOLY for more speed)\n');
end