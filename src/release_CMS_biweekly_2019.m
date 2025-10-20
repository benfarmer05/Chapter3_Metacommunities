%% Script to create a release file for the Connectivity Modeling System
%   simulating the dispersal and connectivity of SCTLD in the Virgin
%   Islands & Puerto Rico
%   
%   Modified to match Dobbelaere et al. 2020 methodology:
%   - Biweekly releases (every 14 days)
%   - 150 particles per polygon per release
%   - 30-day tracking duration (set in CMS, not in this file)
%   - 30-day half-life for particle mass decay (set in CMS)
%
%   Modified for 14-day release frequency
%   ADDED: ID validation against original release points file
%   ADDED: Quarterly file writing option

clear;clc

%% setup

% Get the project root directory
projectPath = matlab.project.rootProject().RootFolder;

% Define paths relative to the project root
dataPath = fullfile(projectPath, 'data');
outputPath = fullfile(projectPath, 'output');

%% TOGGLES: Control file creation behavior

% Toggle 1: Separate files per release event (original behavior)
CREATE_SEPARATE_FILES = false;  % Set to true for HPC batch processing (one file per release)
                                 % Set to false for combined files

% Toggle 2: Quarterly grouping (only applies when CREATE_SEPARATE_FILES = false)
CREATE_QUARTERLY_FILES = true;   % Set to true for 4 quarterly files (~6-7 releases each)
                                 % Set to false for single combined file (all releases)

%%

% Parameters for 2019 biweekly releases (every 14 days)
release_year = 2019;
start_date = datetime(2019, 1, 1, 0, 0, 0);
end_date = datetime(2019, 12, 25, 23, 59, 59); %December 25th set as final date, since USCROMS ends a bit early in 2019

% Create datetime array for all biweekly releases
release_datetimes = start_date:days(14):end_date;
num_releases = length(release_datetimes);

fprintf('========== RELEASE SCHEDULE ==========\n');
fprintf('Release frequency: Every 14 days\n');
fprintf('Year: %d\n', release_year);
fprintf('Total releases: %d\n', num_releases);
fprintf('First release: %s\n', datestr(release_datetimes(1)));
fprintf('Last release: %s\n', datestr(release_datetimes(end)));
fprintf('======================================\n\n');

%% Load original release points for validation

fprintf('========== LOADING ORIGINAL RELEASE POINTS ==========\n');
relpoints_file = fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv');

if ~exist(relpoints_file, 'file')
    error('Original release points file not found: %s', relpoints_file);
end

relpoints = readmatrix(relpoints_file);
fprintf('Loaded original release points file: centroids_vertices_FINALFORCMS.csv\n');
fprintf('Total points in original file: %d\n', size(relpoints, 1));

% Extract IDs from original file (assuming ID is in first column)
original_IDs = relpoints(:, 1);
fprintf('ID range in original file: %d to %d\n', min(original_IDs), max(original_IDs));
fprintf('=====================================================\n\n');

%% Load depth assignments from depth assignment script

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

%% VALIDATION: Check IDs match between depth data and original release points

fprintf('========== ID VALIDATION CHECK ==========\n');

% Sort both ID arrays for comparison
depth_IDs_sorted = sort(IDs);
original_IDs_sorted = sort(original_IDs);

% Check if all depth IDs exist in original file
depth_in_original = ismember(depth_IDs_sorted, original_IDs_sorted);
all_depth_match = all(depth_in_original);

% Check if all original IDs exist in depth file (accounting for NaN exclusions)
original_in_depth = ismember(original_IDs_sorted, depth_IDs_sorted);
num_original_matched = sum(original_in_depth);

fprintf('Validation Results:\n');
fprintf('-------------------\n');
if all_depth_match
    fprintf('✓ SUCCESS: All %d IDs in depth data exist in original release points file\n', numpoints);
else
    fprintf('✗ WARNING: %d IDs in depth data NOT found in original file\n', sum(~depth_in_original));
    missing_IDs = depth_IDs_sorted(~depth_in_original);
    fprintf('  Missing IDs: %s\n', mat2str(missing_IDs(1:min(10, length(missing_IDs)))));
end

fprintf('\nOriginal file coverage:\n');
fprintf('  Total IDs in original file: %d\n', length(original_IDs));
fprintf('  IDs with valid depths: %d (%.1f%%)\n', num_original_matched, 100*num_original_matched/length(original_IDs));
fprintf('  IDs excluded (NaN depths): %d (%.1f%%)\n', ...
    length(original_IDs) - num_original_matched, ...
    100*(length(original_IDs) - num_original_matched)/length(original_IDs));

% Identify which original IDs were excluded
if num_original_matched < length(original_IDs)
    excluded_IDs = original_IDs_sorted(~original_in_depth);
    fprintf('\nFirst 10 excluded IDs (depth too deep or no ocean): %s\n', ...
        mat2str(excluded_IDs(1:min(10, length(excluded_IDs)))));
end

fprintf('=========================================\n\n');

% Halt if validation fails
if ~all_depth_match
    error('ID validation failed! Depth data contains IDs not in original release points file.');
end

fprintf('✓ ID validation passed! Proceeding with release file creation...\n\n');

%% Create release file with 150 particles per polygon

fprintf('========== CREATING RELEASE DATA ==========\n');

% MODIFIED: 150 particles per polygon per release (matching Dobbelaere 2020 density)
num_particles = 65;

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

%% Write release file(s)

if CREATE_SEPARATE_FILES
    %% OPTION 1: Create separate file for each release event
    fprintf('Creating separate files for each release event...\n');
    
    for rel_idx = 1:num_releases
        % Extract data for this release
        start_idx = (rel_idx-1) * numpoints + 1;
        end_idx = rel_idx * numpoints;
        
        % Create filename with parseable date format: YYYY-MM-DD
        rel_date = release_datetimes(rel_idx);
        date_str = sprintf('%04d-%02d-%02d', ...
            year(rel_date), month(rel_date), day(rel_date));
        
        fileName = sprintf('ReleaseFile_USVI_%s.txt', date_str);
        fileID = fopen(fullfile(outputPath, fileName), 'w');
        
        if fileID == -1
            error('Could not open file for writing: %s', fullfile(outputPath, fileName));
        end
        
        % Write data for this release
        chunk_data = [IDs_CMS(start_idx:end_idx), ...
                      longitudes_CMS(start_idx:end_idx), ...
                      latitudes_CMS(start_idx:end_idx), ...
                      depths_CMS(start_idx:end_idx), ...
                      particles_CMS(start_idx:end_idx), ...
                      year_CMS(start_idx:end_idx), ...
                      months_CMS(start_idx:end_idx), ...
                      days_CMS(start_idx:end_idx), ...
                      time_CMS(start_idx:end_idx)];
        
        fprintf(fileID, '%-4.0f %-15.9f %-15.9f %-2.0f %-3.0f %-6.0f %-2.0f %-3.0f %-6.0f\n', chunk_data');
        fclose(fileID);
        
        if mod(rel_idx, 5) == 0 || rel_idx == num_releases
            fprintf('Progress: %d/%d files (%.1f%%)\n', rel_idx, num_releases, 100*rel_idx/num_releases);
        end
    end
    
    fprintf('\nCreated %d separate release files\n', num_releases);
    fprintf('File naming format: ReleaseFile_USVI_YYYY-MM-DD.txt\n');
    fprintf('Example: ReleaseFile_USVI_2019-01-01.txt\n');
    
elseif CREATE_QUARTERLY_FILES
    %% OPTION 2: Create quarterly release files
    fprintf('Creating quarterly release files...\n');
    
    % Define quarters - which release events belong to each quarter
    % Based on the month of each release datetime
    quarter_assignments = ceil(release_months_vals / 3);  % 1-3 months -> Q1, 4-6 -> Q2, etc.
    quarter_names = {'Q1', 'Q2', 'Q3', 'Q4'};
    
    fprintf('DEBUG: Number of unique release months: %d\n', length(unique(release_months_vals)));
    fprintf('DEBUG: Quarter assignments: %s\n', mat2str(unique(quarter_assignments)));
    fprintf('DEBUG: numpoints = %d\n', numpoints);
    
    for q = 1:4
        % Find which release events are in this quarter
        releases_in_quarter = find(quarter_assignments == q);
        
        if isempty(releases_in_quarter)
            fprintf('No releases found for Q%d\n', q);
            continue;
        end
        
        fprintf('DEBUG Q%d: Found %d releases: %s\n', q, length(releases_in_quarter), mat2str(releases_in_quarter'));
        
        % For each release in this quarter, we need all polygons (numpoints rows)
        quarter_indices = [];
        for i = 1:length(releases_in_quarter)
            rel_idx = releases_in_quarter(i);
            
            % Ensure rel_idx is a scalar
            if ~isscalar(rel_idx)
                error('rel_idx is not scalar! Value: %s', mat2str(rel_idx));
            end
            
            % Ensure numpoints is scalar
            if ~isscalar(numpoints)
                error('numpoints is not scalar! Value: %s', mat2str(numpoints));
            end
            
            start_idx = (rel_idx - 1) * numpoints + 1;
            end_idx = rel_idx * numpoints;
            
            % Validate indices are scalar real numbers
            if ~isscalar(start_idx) || ~isreal(start_idx)
                error('start_idx problem! rel_idx=%s, numpoints=%s, start_idx=%s', ...
                    mat2str(rel_idx), mat2str(numpoints), mat2str(start_idx));
            end
            if ~isscalar(end_idx) || ~isreal(end_idx)
                error('end_idx problem! rel_idx=%s, numpoints=%s, end_idx=%s', ...
                    mat2str(rel_idx), mat2str(numpoints), mat2str(end_idx));
            end
            
            % Create index range safely
            idx_range = (start_idx:end_idx)';
            quarter_indices = [quarter_indices; idx_range];
        end
        
        fprintf('DEBUG Q%d: Total indices to write: %d\n', q, length(quarter_indices));
        
        % Create filename
        fileName = sprintf('ReleaseFile_USVI_%d_%s.txt', release_year, quarter_names{q});
        fileID = fopen(fullfile(outputPath, fileName), 'w');
        
        if fileID == -1
            error('Could not open file for writing: %s', fileName);
        end
        
        % Write data for this quarter
        quarter_data = [IDs_CMS(quarter_indices), ...
                        longitudes_CMS(quarter_indices), ...
                        latitudes_CMS(quarter_indices), ...
                        depths_CMS(quarter_indices), ...
                        particles_CMS(quarter_indices), ...
                        year_CMS(quarter_indices), ...
                        months_CMS(quarter_indices), ...
                        days_CMS(quarter_indices), ...
                        time_CMS(quarter_indices)];
        
        fprintf(fileID, '%-4.0f %-15.9f %-15.9f %-2.0f %-3.0f %-6.0f %-2.0f %-3.0f %-6.0f\n', quarter_data');
        fclose(fileID);
        
        num_releases_in_q = length(releases_in_quarter);
        fprintf('Created %s with %d releases (%d total rows)\n', ...
            fileName, num_releases_in_q, length(quarter_indices));
    end
    
    fprintf('\nCreated 4 quarterly release files\n');
    fprintf('File naming format: ReleaseFile_USVI_YYYY_QX.txt\n');
    fprintf('Example: ReleaseFile_USVI_2019_Q1.txt\n');
    
else
    %% OPTION 3: Create single combined file (original behavior)
    currentDateTime = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
    currentDateTimeStr = string(currentDateTime);
    fileName = sprintf('ReleaseFile_USVI_%d_biweekly_%s.txt', release_year, currentDateTimeStr);
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
end

% Display confirmation message
fprintf('\n========== RELEASE FILE SUMMARY ==========\n');
if CREATE_SEPARATE_FILES
    fprintf('Mode: SEPARATE FILES (HPC batch mode)\n');
    fprintf('Files created: %d\n', num_releases);
    fprintf('Naming format: ReleaseFile_USVI_YYYY-MM-DD.txt\n');
elseif CREATE_QUARTERLY_FILES
    fprintf('Mode: QUARTERLY FILES\n');
    fprintf('Files created: 4 (Q1-Q4)\n');
    fprintf('Releases per file: ~6-7\n');
    fprintf('Naming format: ReleaseFile_USVI_YYYY_QX.txt\n');
else
    fprintf('Mode: SINGLE COMBINED FILE\n');
    fprintf('File: %s\n', fileName);
end
fprintf('Release schedule: Every 14 days\n');
fprintf('Year: %d\n', release_year);
fprintf('Number of releases: %d\n', num_releases);
fprintf('Polygons per release: %d\n', numpoints);
fprintf('Particles per polygon: %d\n', num_particles);
fprintf('Total particles per release: %d\n', numpoints * num_particles);
fprintf('Total release events written: %d\n', num_releases);
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

if CREATE_SEPARATE_FILES
    % Read first release file for validation
    first_date = release_datetimes(1);
    date_str = sprintf('%04d-%02d-%02d', ...
        year(first_date), month(first_date), day(first_date));
    fileName = sprintf('ReleaseFile_USVI_%s.txt', date_str);
    fprintf('Validating first file: %s\n', fileName);
elseif CREATE_QUARTERLY_FILES
    % Read Q1 file for validation
    fileName = sprintf('ReleaseFile_USVI_%d_Q1.txt', release_year);
    fprintf('Validating first quarterly file: %s\n', fileName);
else
    % Use the single combined file
    fprintf('Validating file: %s\n', fileName);
end

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
fprintf('Methodology: Dobbelaere et al. 2020 (modified)\n');
fprintf('Release frequency: Every 14 days\n');
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