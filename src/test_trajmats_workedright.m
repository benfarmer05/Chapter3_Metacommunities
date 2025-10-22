% % Load the event file
% load('D:\Dissertation\CMS_traj\output\CMS_traj\Q1_2019\event_2019_Jan01_120000.mat')

% Randomly select a source reef from this event
rand_idx = randi(length(event_data));
selected_source = event_data(rand_idx);

fprintf('\n=== SOURCE REEF VERIFICATION ===\n');
fprintf('Selected source reef index: %d\n', selected_source.source_reef_idx);
fprintf('Selected source reef ID: %d\n', selected_source.source_reef_id);

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

%% INDEXING DIAGNOSTIC
fprintf('\n=== INDEXING DIAGNOSTIC ===\n');
fprintf('selected_source.source_reef_idx = %d\n', selected_source.source_reef_idx);
fprintf('selected_source.source_reef_id = %d\n', selected_source.source_reef_id);
fprintf('unique_IDs(selected_source.source_reef_idx) = %d\n', ...
    unique_IDs(selected_source.source_reef_idx));

if unique_IDs(selected_source.source_reef_idx) == selected_source.source_reef_id
    fprintf('✓ Index matches reef ID correctly\n');
else
    warning('MISMATCH! The index does not correspond to the stored reef ID!');
end

%% VERIFICATION 1: Check that the reef ID exists in centroids data
centroid_match = find(unique_IDs == selected_source.source_reef_id);
if ~isempty(centroid_match)
    fprintf('✓ Reef ID %d found in centroids data at index %d\n', ...
        selected_source.source_reef_id, centroid_match);
else
    warning('✗ Reef ID %d NOT found in centroids data!', selected_source.source_reef_id);
end

%% VERIFICATION 2: Check that the reef ID exists in release file
release_match = find(release_reef_IDs == selected_source.source_reef_id);
if ~isempty(release_match)
    fprintf('✓ Reef ID %d found in release file (%d occurrences)\n', ...
        selected_source.source_reef_id, length(release_match));
else
    warning('✗ Reef ID %d NOT found in release file!', selected_source.source_reef_id);
end

%% Get the correct polygon using reef ID (not index)
% CRITICAL: Always look up polygon by reef ID, not by index
correct_idx = find(unique_IDs == selected_source.source_reef_id);
poly_x = Xs(correct_idx, :);
poly_y = Ys(correct_idx, :);

%% VERIFICATION 3: Check that release position is inside the correct polygon
release_lon = selected_source.lon(1,1);  % First timestep, first particle
release_lat = selected_source.lat(1,1);

is_inside = inpolygon(release_lon, release_lat, poly_x, poly_y);

if is_inside
    fprintf('✓ Release location (%.4f, %.4f) is inside polygon for reef ID %d\n', ...
        release_lon, release_lat, selected_source.source_reef_id);
else
    warning('✗ Release location (%.4f, %.4f) is NOT inside polygon for reef ID %d!', ...
        release_lon, release_lat, selected_source.source_reef_id);
end

%% RELEASE FILE VERIFICATION
fprintf('\n=== RELEASE FILE VERIFICATION ===\n');
release_indices = find(release_reef_IDs == selected_source.source_reef_id);
fprintf('Found %d release points for reef ID %d\n', length(release_indices), selected_source.source_reef_id);

if size(releasefile, 2) >= 3
    release_file_lons = releasefile(release_indices, 2);
    release_file_lats = releasefile(release_indices, 3);
    
    % Convert longitude if needed
    if min(release_file_lons) > 180
        release_file_lons = release_file_lons - 360;
    end
    
    fprintf('Release file coordinates for this reef:\n');
    fprintf('  Lon range: [%.4f, %.4f]\n', min(release_file_lons), max(release_file_lons));
    fprintf('  Lat range: [%.4f, %.4f]\n', min(release_file_lats), max(release_file_lats));
    
    % Check if these are inside the polygon
    in_poly = inpolygon(release_file_lons, release_file_lats, poly_x, poly_y);
    fprintf('  %d/%d release points are inside polygon\n', sum(in_poly), length(in_poly));
end

%% Plot trajectories with reef polygons
figure('Position', [100 100 1200 800]);
hold on;

% Plot all reef polygons (gray)
for i = 1:size(Xs, 1)
    plot(Xs(i,:), Ys(i,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
end

% Plot trajectories from selected source
plot(selected_source.lon, selected_source.lat, 'b-', 'LineWidth', 0.5, 'Color', [0 0.4470 0.7410]);

% Highlight the source reef polygon in red
plot(poly_x, poly_y, 'r-', 'LineWidth', 2.5);

% Mark release location with a star
plot(release_lon, release_lat, 'r*', 'MarkerSize', 15, 'LineWidth', 2);

% Plot release file points if available
if exist('release_file_lons', 'var')
    plot(release_file_lons, release_file_lats, 'go', 'MarkerSize', 8, 'LineWidth', 2);
    legend('All reefs', 'Particle trajectories', 'Source reef', 'Release location', ...
        'Release file points', 'Location', 'best');
else
    legend('All reefs', 'Particle trajectories', 'Source reef', 'Release location', ...
        'Location', 'best');
end

xlabel('Longitude');
ylabel('Latitude');
title(sprintf('Trajectories from Reef ID %d - %d particles', ...
    selected_source.source_reef_id, selected_source.n_particles));
axis equal tight;
grid on;

fprintf('\n=== PLOT COMPLETE ===\n');







%% Test the preprocessing indexing
fprintf('\n=== TESTING PREPROCESSING INDEXING ===\n');

% Load a random trajectory file
trajlist = dir(fullfile(tempPath, 'traj*.nc'));
test_file = fullfile(tempPath, trajlist(1).name);

% Read location as-is (no modification)
location_raw = ncread(test_file, 'location');
lon = ncread(test_file, 'lon');
lat = ncread(test_file, 'lat');

if min(lon(:)) > 180
    lon = lon - 360;
end

% Get first unique location
test_loc = location_raw(1);
fprintf('Raw location value: %d\n', test_loc);

% Test WITHOUT adding 1 (correct way)
reef_id_correct = release_reef_IDs(test_loc);
fprintf('Reef ID (no +1): %d\n', reef_id_correct);

% Test WITH adding 1 (preprocessing way - WRONG)
reef_id_wrong = release_reef_IDs(test_loc + 1);
fprintf('Reef ID (with +1): %d\n', reef_id_wrong);

% Check which polygon contains the actual trajectory start
start_lon = lon(1, 1);
start_lat = lat(1, 1);

% Check correct reef ID
idx_correct = find(unique_IDs == reef_id_correct);
in_correct = inpolygon(start_lon, start_lat, Xs(idx_correct,:), Ys(idx_correct,:));

% Check wrong reef ID
idx_wrong = find(unique_IDs == reef_id_wrong);
in_wrong = inpolygon(start_lon, start_lat, Xs(idx_wrong,:), Ys(idx_wrong,:));

fprintf('\nStart point (%.4f, %.4f):\n', start_lon, start_lat);
fprintf('  Inside polygon for reef %d (no +1): %d\n', reef_id_correct, in_correct);
fprintf('  Inside polygon for reef %d (with +1): %d\n', reef_id_wrong, in_wrong);






























%% Check time format in trajectory files
projectPath = matlab.project.rootProject().RootFolder;
dataPath = fullfile(projectPath, 'data');
tempPath = fullfile('D:\Dissertation\CMS_traj', 'Q1_2019');

trajlist = dir(fullfile(tempPath, 'traj*.nc'));
test_file = fullfile(tempPath, trajlist(1).name);

% Read time variable
time = ncread(test_file, 'time');
releasedate = ncread(test_file, 'releasedate');

fprintf('=== TIME FORMAT CHECK ===\n');
fprintf('Release date (Julian): %.6f\n', releasedate(1));
fprintf('Release date as datetime: %s\n', datestr(datetime(releasedate(1), 'ConvertFrom', 'juliandate')));
fprintf('\nTime variable:\n');
fprintf('  First value: %.6f\n', time(1,1));
fprintf('  Last value: %.6f\n', time(end,1));
fprintf('  Range: %.6f\n', time(end,1) - time(1,1));
fprintf('  Range in days: %.2f\n', (time(end,1) - time(1,1))/86400);

% Try different interpretations
fprintf('\nTrying different time interpretations:\n');
fprintf('1. If time is POSIX (seconds since 1970):\n');
fprintf('   Start: %s\n', datestr(datetime(time(1,1), 'ConvertFrom', 'posixtime')));
fprintf('   End: %s\n', datestr(datetime(time(end,1), 'ConvertFrom', 'posixtime')));

fprintf('\n2. If time is Julian date:\n');
fprintf('   Start: %s\n', datestr(datetime(time(1,1), 'ConvertFrom', 'juliandate')));
fprintf('   End: %s\n', datestr(datetime(time(end,1), 'ConvertFrom', 'juliandate')));

fprintf('\n3. If time is seconds since release (CORRECT FORMAT):\n');
fprintf('   Duration: %.2f days\n', (time(end,1) - time(1,1))/86400);
fprintf('   Absolute start: %s\n', datestr(datetime(releasedate(1), 'ConvertFrom', 'juliandate')));
fprintf('   Absolute end: %s\n', datestr(datetime(releasedate(1) + time(end,1)/86400, 'ConvertFrom', 'juliandate')));








% Check what's in your release file
releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019_Q1.txt'));

% If release file has more columns, check if there's a time column
fprintf('Release file columns: %d\n', size(releasefile, 2));
fprintf('First few rows:\n');
disp(releasefile(1:5, :));

% Also check the trajectory file's releasedate variable
fprintf('\nUnique release dates in this trajectory file:\n');
unique_releasedates = unique(releasedate);
for i = 1:length(unique_releasedates)
    fprintf('  %.6f = %s\n', unique_releasedates(i), ...
        datestr(datetime(unique_releasedates(i), 'ConvertFrom', 'juliandate')));
end







% Find nest files
nest_files = dir(fullfile(dataPath, 'nest_1_2019*.nc'));

if isempty(nest_files)
    fprintf('No nest*.nc files found in %s\n', dataPath);
else
    fprintf('Found %d nest files\n', length(nest_files));
    
    % Pick a random one
    rand_idx = randi(length(nest_files));
    nest_file = fullfile(dataPath, nest_files(rand_idx).name);
    
    fprintf('Selected file: %s\n', nest_files(rand_idx).name);
    
    % Get info about the file
    file_info = ncinfo(nest_file);
    
    fprintf('\n=== FILE VARIABLES ===\n');
    for i = 1:length(file_info.Variables)
        fprintf('  %s: %s\n', file_info.Variables(i).Name, mat2str(file_info.Variables(i).Size));
    end
    
    % Check if 'time' variable exists
    if any(strcmp({file_info.Variables.Name}, 'time'))
        time_hydro = ncread(nest_file, 'time');
        
        fprintf('\n=== HYDRODYNAMIC TIME VARIABLE ===\n');
        fprintf('Number of timesteps: %d\n', length(time_hydro));
        fprintf('First time value: %.6f\n', time_hydro(1));
        fprintf('Last time value: %.6f\n', time_hydro(end));
        fprintf('Time range: %.6f\n', time_hydro(end) - time_hydro(1));
        
        % Check for time attributes
        time_var_idx = find(strcmp({file_info.Variables.Name}, 'time'));
        if ~isempty(file_info.Variables(time_var_idx).Attributes)
            fprintf('\nTime variable attributes:\n');
            for j = 1:length(file_info.Variables(time_var_idx).Attributes)
                attr = file_info.Variables(time_var_idx).Attributes(j);
                fprintf('  %s: %s\n', attr.Name, mat2str(attr.Value));
            end
        end
        
        % Try interpretations
        fprintf('\n=== TIME INTERPRETATIONS ===\n');
        fprintf('If Julian date:\n');
        fprintf('  Start: %s\n', datestr(datetime(time_hydro(1), 'ConvertFrom', 'juliandate')));
        fprintf('  End: %s\n', datestr(datetime(time_hydro(end), 'ConvertFrom', 'juliandate')));
        
        fprintf('\nIf days since some epoch:\n');
        fprintf('  Duration: %.2f days\n', time_hydro(end) - time_hydro(1));
        
        fprintf('\nIf seconds since some epoch:\n');
        fprintf('  Duration: %.2f days\n', (time_hydro(end) - time_hydro(1))/86400);
    else
        fprintf('\nNo "time" variable found in this file.\n');
    end
end



time_hydro = ncread(nest_file, 'Time');
fprintf('Time value in this hydrodynamic file: %.6f hours\n', time_hydro);
fprintf('Which is: %.6f days since midnight Jan 1, 2019\n', time_hydro/24);
fprintf('Absolute datetime: %s\n', datestr(datetime(2019,1,1,0,0,0) + hours(time_hydro)));