%% Script to create a release file for the Connectivity Modeling System
%   simulating the dispersal and connectivity of SCTLD in the Virgin
%   Islands & Puerto Rico
%   
%   Modified to match Dobbelaere et al. 2020 methodology:
%   - Monthly releases (1st of each month)
%   - 150 particles per polygon per release
%   - 30-day tracking duration (set in CMS, not in this file)
%   - 30-day half-life for particle mass decay (set in CMS)
%
%   1 Oct 2025 - Modified for monthly releases

clear;clc

%% setup

% Get the project root directory
projectPath = matlab.project.rootProject().RootFolder;

% Define paths relative to the project root
dataPath = fullfile(projectPath, 'data');
outputPath = fullfile(projectPath, 'output');

%%

% Parameters for 2019 monthly releases
% Release on the 1st of each month throughout 2019
release_year = 2019;
release_months = 1:12;  % All 12 months
release_day = 1;        % 1st of each month
release_hour = 0;       % Midnight (00:00)

% Create datetime array for all monthly releases
release_datetimes = datetime(release_year, release_months, release_day, release_hour, 0, 0);
num_releases = length(release_datetimes);

fprintf('========== RELEASE SCHEDULE ==========\n');
fprintf('Release frequency: Monthly (1st of each month)\n');
fprintf('Year: %d\n', release_year);
fprintf('Total releases: %d\n', num_releases);
fprintf('First release: %s\n', datestr(release_datetimes(1)));
fprintf('Last release: %s\n', datestr(release_datetimes(end)));
fprintf('======================================\n\n');

%% MODIFIED: Load depth assignments from depth assignment script

fprintf('========== LOADING DEPTH ASSIGNMENTS ==========\n');

% Find the most recent depth assignment file
depthFiles = dir(fullfile(outputPath, 'depth_assignments_*.mat'));
if isempty(depthFiles)
    error('No depth assignment files found in output folder. Run the depth assignment script first!');
end

% Sort by date and get the most recent
[~, idx] = max([depthFiles.datenum]);
depthFile = fullfile(outputPath, depthFiles(idx).name);
fprintf('Loading depth assignments from: %s\n', depthFiles(idx).name);

% Load the depth assignments
load(depthFile, 'depth_assignments', 'IDs_release', 'deepest_depth');

% Extract components: [ID, Lon, Lat, Depth]
depth_IDs = depth_assignments(:,1);
depth_lons = depth_assignments(:,2);
depth_lats = depth_assignments(:,3);
assigned_depths = depth_assignments(:,4);

% CRITICAL: Filter out points with NaN depths (>60m or no ocean)
valid_depth_mask = ~isnan(assigned_depths);
fprintf('Total points in depth file: %d\n', length(assigned_depths));
fprintf('Points with valid depths: %d\n', sum(valid_depth_mask));
fprintf('Points excluded (NaN depths): %d\n', sum(~valid_depth_mask));

% Keep only valid points
IDs = depth_IDs(valid_depth_mask);
longitudes = depth_lons(valid_depth_mask) + 360;  % Add 360 for CMS format
latitudes = depth_lats(valid_depth_mask);
depths_assigned = assigned_depths(valid_depth_mask);

numpoints = length(IDs);  % Updated number of valid points (polygons)

fprintf('Using %d release points (polygons) with valid depths\n', numpoints);
fprintf('Depth range: %.2f to %.2f m\n', min(depths_assigned), max(depths_assigned));
fprintf('===============================================\n\n');

%% Create release file with 150 particles per polygon

fprintf('========== CREATING RELEASE DATA ==========\n');

% MODIFIED: 150 particles per polygon per release (matching Dobbelaere 2020 density)
num_particles = 150;

fprintf('Particles per polygon per release: %d\n', num_particles);
fprintf('Total particles per release: %d\n', numpoints * num_particles);
fprintf('Total particles over all releases: %d\n', numpoints * num_particles * num_releases);
fprintf('==========================================\n\n');

% Extract year, month, day, and hour from release datetimes
release_years = year(release_datetimes);
release_months_vals = month(release_datetimes);
release_days = day(release_datetimes);
release_hours = hour(release_datetimes);

% Create arrays properly aligned for each combination of point and time
% For each release datetime, repeat all points (ONE ROW PER RELEASE EVENT)
IDs_CMS = repmat(IDs, num_releases, 1);
longitudes_CMS = repmat(longitudes, num_releases, 1);
latitudes_CMS = repmat(latitudes, num_releases, 1);

% Use assigned depths for each release point
depths_CMS = repmat(depths_assigned, num_releases, 1);

% 150 particles per release event at each location
particles_CMS = repmat(num_particles, numpoints * num_releases, 1);

% Pre-allocate temporal arrays
total_rows = numpoints * num_releases;
days_CMS = zeros(total_rows, 1);
months_CMS = zeros(total_rows, 1);
years_CMS = zeros(total_rows, 1);
hours_CMS = zeros(total_rows, 1);

% Fill temporal arrays efficiently using vectorized operations
fprintf('Filling temporal arrays...\n');
for i = 1:num_releases
    start_idx = (i-1) * numpoints + 1;
    end_idx = i * numpoints;
    
    days_CMS(start_idx:end_idx) = release_days(i);
    months_CMS(start_idx:end_idx) = release_months_vals(i);
    years_CMS(start_idx:end_idx) = release_years(i);
    hours_CMS(start_idx:end_idx) = release_hours(i);
    
    if i == num_releases
        fprintf('Progress: %d/%d releases (100%%)\n', i, num_releases);
    end
end

% Convert hours to seconds for CMS format
time_CMS = hours_CMS * 3600;  % Convert hours to seconds

% Create matrix for sorting: [Year, Month, Day, Time, ID, Longitude, Latitude, Depth, Particles]
fprintf('Creating and sorting data matrix...\n');
sortMatrix = [years_CMS, months_CMS, days_CMS, time_CMS, IDs_CMS, longitudes_CMS, latitudes_CMS, depths_CMS, particles_CMS];

% Sort by Year, Month, Day, Time, then ID
fprintf('Sorting %d rows...\n', size(sortMatrix, 1));
sortMatrix = sortrows(sortMatrix, [1, 2, 3, 4, 5]);
fprintf('Sorting complete.\n');

% Extract sorted data
year_CMS = sortMatrix(:, 1);
months_CMS = sortMatrix(:, 2);
days_CMS = sortMatrix(:, 3);
time_CMS = sortMatrix(:, 4);
IDs_CMS = sortMatrix(:, 5);
longitudes_CMS = sortMatrix(:, 6);
latitudes_CMS = sortMatrix(:, 7);
depths_CMS = sortMatrix(:, 8);
particles_CMS = sortMatrix(:, 9);

% Open a file for writing
currentDateTime = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
currentDateTimeStr = string(currentDateTime);
fileName = sprintf('ReleaseFile_USVI_%d_monthly_150ppp_%s.txt', release_year, currentDateTimeStr);
fileID = fopen(fullfile(outputPath, fileName), 'w');

if fileID == -1
    error('Could not open file for writing: %s', fullfile(outputPath, fileName));
end

% Write data in chunks to handle large datasets
fprintf('Writing data to file: %s\n', fileName);
fprintf('Writing %d rows of data...\n', length(IDs_CMS));

chunk_size = 50000; % Write 50k rows at a time
total_rows_actual = length(IDs_CMS);
num_chunks = ceil(total_rows_actual / chunk_size);

tic;  % Start timer
for chunk = 1:num_chunks
    start_idx = (chunk-1) * chunk_size + 1;
    end_idx = min(chunk * chunk_size, total_rows_actual);
    
    % Create data matrix for this chunk
    chunk_data = [IDs_CMS(start_idx:end_idx), ...
                  longitudes_CMS(start_idx:end_idx), ...
                  latitudes_CMS(start_idx:end_idx), ...
                  depths_CMS(start_idx:end_idx), ...
                  particles_CMS(start_idx:end_idx), ...
                  year_CMS(start_idx:end_idx), ...
                  months_CMS(start_idx:end_idx), ...
                  days_CMS(start_idx:end_idx), ...
                  time_CMS(start_idx:end_idx)];
    
    % Check for any invalid values in this chunk
    if any(isnan(chunk_data(:))) || any(isinf(chunk_data(:)))
        warning('Invalid values detected in chunk %d', chunk);
    end
    
    % Write chunk to file
    fprintf(fileID, '%-4.0f %-15.9f %-15.9f %-2.0f %-3.0f %-6.0f %-2.0f %-3.0f %-6.0f\n', chunk_data');
    
    % Progress indicator
    if mod(chunk, 5) == 0 || chunk == num_chunks
        fprintf('Progress: %d/%d chunks (%.1f%%)\n', chunk, num_chunks, 100*chunk/num_chunks);
    end
end

elapsed_time = toc;  % End timer
fprintf('File writing completed in %.2f seconds.\n', elapsed_time);

fclose(fileID);

% Verify file was written correctly
file_info = dir(fullfile(outputPath, fileName));
fprintf('File size: %.2f MB\n', file_info.bytes / (1024*1024));

% Display confirmation message
fprintf('\n========== RELEASE FILE SUMMARY ==========\n');
fprintf('File: %s\n', fileName);
fprintf('Release schedule: Monthly (1st of each month)\n');
fprintf('Year: %d\n', release_year);
fprintf('Number of releases: %d\n', num_releases);
fprintf('Polygons per release: %d\n', numpoints);
fprintf('Particles per polygon: %d\n', num_particles);
fprintf('Total particles per release: %d\n', numpoints * num_particles);
fprintf('Total release events written: %d\n', total_rows_actual);
fprintf('==========================================\n\n');

fprintf('NOTE: Set CMS to track particles for 30 days with 30-day half-life\n');
fprintf('      to match Dobbelaere et al. 2020 methodology.\n\n');

%% Visualizations of release points

fprintf('========== CREATING VISUALIZATIONS ==========\n');

% Get unique release points (remove temporal replication)
unique_idx = 1:numpoints;
lon_unique = longitudes(unique_idx);
lat_unique = latitudes(unique_idx);
depth_unique = depths_assigned(unique_idx);

% 2D Map View
figure('Position', [100 100 1200 500]);

subplot(1,2,1);
scatter(lon_unique, lat_unique, 30, depth_unique, 'filled');
colormap(parula);
c = colorbar;
c.Label.String = 'Depth (m)';
xlabel('Longitude (°E)');
ylabel('Latitude (°N)');
title(sprintf('Release Point Locations (n=%d polygons)', numpoints));
grid on;
axis equal tight;

subplot(1,2,2);
histogram(depth_unique, 'BinWidth', 2, 'FaceColor', [0.3 0.6 0.9]);
xlabel('Depth (m)');
ylabel('Number of Release Points');
title('Distribution of Release Depths');
grid on;

% 3D Visualization (oceanographic convention: depth positive downward)
figure('Position', [100 100 900 700]);
scatter3(lon_unique, lat_unique, depth_unique, 20, depth_unique, 'filled');
colormap(parula);
c = colorbar;
c.Label.String = 'Depth below surface (m)';
xlabel('Longitude (°E)');
ylabel('Latitude (°N)');
zlabel('Depth below surface (m)');
title(sprintf('3D Release Point Space (n=%d polygons, 150 particles each)', numpoints));
grid on;
view(3);
set(gca, 'ZDir', 'reverse'); % Depth increases downward
zlim([0, max(depth_unique)+5]);

fprintf('Visualizations complete.\n');
fprintf('=============================================\n');

%% Validation Test - Verify depths in release file

fprintf('\n========== VALIDATION TEST ==========\n');

% Read release file and get unique points
releaseData = readmatrix(fullfile(outputPath, fileName));
[test_IDs, test_idx] = unique(releaseData(:,1), 'stable');
test_lons = releaseData(test_idx, 2);
test_lats = releaseData(test_idx, 3);
test_depths = releaseData(test_idx, 4);
test_particles = releaseData(test_idx, 5);

% Verify particle counts
fprintf('Particle count verification:\n');
fprintf('Expected particles per polygon: %d\n', num_particles);
fprintf('Unique particle counts in file: %s\n', mat2str(unique(test_particles)'));
if all(test_particles == num_particles)
    fprintf('✓ All polygons have correct particle count!\n');
else
    warning('✗ Inconsistent particle counts detected!');
end

% Sample 5 random points for validation
n = min(5, length(test_IDs));
sample_idx = randperm(length(test_IDs), n);

fprintf('\nSampled %d points:\n', n);
for i = 1:n
    fprintf('ID %d: Depth = %.2f m, Particles = %d\n', ...
        test_IDs(sample_idx(i)), test_depths(sample_idx(i)), test_particles(sample_idx(i)));
end

% Plot all points with samples highlighted
figure('Position', [100 100 1200 500]);
subplot(1,2,1);
scatter(test_lons, test_lats, 20, test_depths, 'filled', 'MarkerFaceAlpha', 0.3);
hold on;
scatter(test_lons(sample_idx), test_lats(sample_idx), 150, 'r', 'LineWidth', 2);
% Add ID labels next to red circles
for i = 1:length(sample_idx)
    text(test_lons(sample_idx(i)), test_lats(sample_idx(i)), ...
        sprintf(' ID %d', test_IDs(sample_idx(i))), ...
        'Color', 'red', 'FontWeight', 'bold', 'FontSize', 10);
end
colorbar; xlabel('Longitude'); ylabel('Latitude'); 
title('Sampled Points (red circles)'); axis equal tight; grid on;

subplot(1,2,2);
scatter3(test_lons, test_lats, test_depths, 20, test_depths, 'filled', 'MarkerFaceAlpha', 0.3);
hold on;
scatter3(test_lons(sample_idx), test_lats(sample_idx), test_depths(sample_idx), 150, 'r', 'LineWidth', 2);
% Add ID labels next to red circles in 3D plot
for i = 1:length(sample_idx)
    text(test_lons(sample_idx(i)), test_lats(sample_idx(i)), test_depths(sample_idx(i)), ...
        sprintf(' ID %d', test_IDs(sample_idx(i))), ...
        'Color', 'red', 'FontWeight', 'bold', 'FontSize', 10);
end
colorbar; xlabel('Longitude'); ylabel('Latitude'); zlabel('Depth (m)');
title('3D View'); view(3); set(gca, 'ZDir', 'reverse'); grid on;

fprintf('=====================================\n');

%% Summary Statistics

fprintf('\n========== FINAL SUMMARY ==========\n');
fprintf('Methodology: Dobbelaere et al. 2020\n');
fprintf('Release frequency: Monthly (1st of month)\n');
fprintf('Particles per polygon: %d\n', num_particles);
fprintf('Number of polygons: %d\n', numpoints);
fprintf('Particles per release: %d\n', numpoints * num_particles);
fprintf('Number of releases: %d\n', num_releases);
fprintf('Total particle releases: %d\n', numpoints * num_particles * num_releases);
fprintf('\nParticle density: %.1f particles/km²\n', num_particles / 0.4225);
fprintf('(Dobbelaere 2020 used 356.6 particles/km²)\n');
fprintf('\nREMINDER: Configure CMS for:\n');
fprintf('  - 30-day tracking duration\n');
fprintf('  - 30-day half-life (gamma parameter)\n');
fprintf('  - Write frequency: 15-20 minutes recommended\n');
fprintf('===================================\n');