%% Script to create a release file for the Connectivity Modeling System
%   simulating the dispersal and connectivity of SCTLD in the Virgin
%   Islands & Puerto Rico
%   25 Sep 2025

clear;clc

%% STOPPING POINT - 21 May 2024
%
% 3. use nearest neighbor / interpolation to select the closest reasonable
% depth to the release points
%
% 4. Dan 21 may 2024 Teams: My probably wrong prediction is that you're going to want daily releases, weekly matrices, and break up the analysis by month

%%
% .nc file format: %nest_nestnumber_yyyymmddhhmmss [.nc]
%   here, I'm working with 'nest_1_20190101000000' through
%                          'nest_1_20190101200000'
%   so,                    'nest_1_2019_January_1_00:00' through
%                          'nest_1_2019_January_1_20:00'

%release file format:
%   Polygon Longitude Latitude Depth Number Year Month Day Second
%   1       277.2     24.66    0     10     2016 1     1   0

%% setup

% Get the project root directory
projectPath = matlab.project.rootProject().RootFolder;

% Define paths relative to the project root
dataPath = fullfile(projectPath, 'data');

outputPath = fullfile(projectPath, 'output');

%%

%parameters for start and end dates
startDate = datetime(2019, 1, 1); % Example: Start date
endDate = datetime(2019, 12, 31);  % Example: End date

%read in release points from GIS output and ensure they are sorted by
% their unique ID
relpoints = readmatrix(fullfile(dataPath, 'points_650_none-on-land.csv'));
% relpoints = relpoints(:, 10:12); %for QGIS
% relpoints = sortrows(relpoints, 1); %for QGIS

IDs = relpoints(:,1);
longitudes = relpoints(:,2) + 360;
latitudes = relpoints(:,3);

numpoints = size(relpoints, 1);  % #/points

%generate daily release dates
release_dates = startDate:endDate;  % Creates daily dates from start to end
release_dates = release_dates';     % Convert to column vector

num_releases = length(release_dates);

% %define the release structure
% num_particles = 10; % #/particles released per lat/lon (line)
% 
% % Extract year, month, and day from release dates
% release_years = year(release_dates);
% release_months = month(release_dates);
% release_days = day(release_dates);
% 
% % Create arrays properly aligned for each combination of point and time
% % For each release day, repeat all points
% IDs_CMS = repmat(IDs, num_releases, 1);
% longitudes_CMS = repmat(longitudes, num_releases, 1);
% latitudes_CMS = repmat(latitudes, num_releases, 1);
% depths_CMS = ones(numpoints * num_releases, 1); %this will actually be reef depth (and thus affix directly to 'ID', so will need to be defined and applied 'repmat')
% particles_CMS = repmat(num_particles, numpoints * num_releases, 1); %this could also vary by ID / time / location, etc.
% 
% % OPTIMIZED: Pre-allocate temporal arrays
% total_rows = numpoints * num_releases;
% days_CMS = zeros(total_rows, 1);
% months_CMS = zeros(total_rows, 1);
% years_CMS = zeros(total_rows, 1);
% 
% % Fill temporal arrays efficiently using vectorized operations
% fprintf('Filling temporal arrays...\n');
% for i = 1:num_releases
%     if mod(i, 50) == 0 || i == num_releases
%         fprintf('Progress: %d/%d releases (%.1f%%)\n', i, num_releases, 100*i/num_releases);
%     end
% 
%     start_idx = (i-1) * numpoints + 1;
%     end_idx = i * numpoints;
% 
%     days_CMS(start_idx:end_idx) = release_days(i);
%     months_CMS(start_idx:end_idx) = release_months(i);
%     years_CMS(start_idx:end_idx) = release_years(i);
% end
% 
% time_CMS = zeros(numpoints * num_releases, 1); %can choose to release at different times of day if desired

%define the release structure
num_particles = 1; % #/particles released per release event (not per line)

% Extract year, month, and day from release dates
release_years = year(release_dates);
release_months = month(release_dates);
release_days = day(release_dates);

% Create arrays properly aligned for each combination of point and time
% For each release day, repeat all points (ONE ROW PER RELEASE EVENT)
IDs_CMS = repmat(IDs, num_releases, 1);
longitudes_CMS = repmat(longitudes, num_releases, 1);
latitudes_CMS = repmat(latitudes, num_releases, 1);
depths_CMS = ones(numpoints * num_releases, 1); %this will actually be reef depth
particles_CMS = repmat(num_particles, numpoints * num_releases, 1); % particles per release event

% OPTIMIZED: Pre-allocate temporal arrays
total_rows = numpoints * num_releases; % This should be 650 * 365 = 237,250
days_CMS = zeros(total_rows, 1);
months_CMS = zeros(total_rows, 1);
years_CMS = zeros(total_rows, 1);

% Fill temporal arrays efficiently using vectorized operations
fprintf('Filling temporal arrays...\n');
for i = 1:num_releases
    if mod(i, 50) == 0 || i == num_releases
        fprintf('Progress: %d/%d releases (%.1f%%)\n', i, num_releases, 100*i/num_releases);
    end
    
    start_idx = (i-1) * numpoints + 1;
    end_idx = i * numpoints;
    
    days_CMS(start_idx:end_idx) = release_days(i);
    months_CMS(start_idx:end_idx) = release_months(i);
    years_CMS(start_idx:end_idx) = release_years(i);
end

time_CMS = zeros(numpoints * num_releases, 1); %can choose to release at different times of day if desired

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

% % Open a file for writing
% currentDateTime = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
% currentDateTimeStr = string(currentDateTime);
% fileName = "ReleaseFile_USVI_2019_daily_" + currentDateTimeStr + ".txt";
% fileID = fopen(fullfile(outputPath, fileName), 'w');
% 
% % FAST VERSION: Write directly to file without cell array conversion
% fprintf('Writing data to file: %s\n', fileName);
% fprintf('Writing %d rows of data...\n', length(IDs_CMS));
% tic;  % Start timer
% fprintf(fileID, '%-4d %-15.9f %-15.9f %-2d %-3d %-6d %-2d %-3d %-2d\n', ...
%     [IDs_CMS, longitudes_CMS, latitudes_CMS, depths_CMS, particles_CMS, ...
%      year_CMS, months_CMS, days_CMS, time_CMS]');
% elapsed_time = toc;  % End timer
% fprintf('File writing completed in %.2f seconds.\n', elapsed_time);
% 
% fclose(fileID);
% 
% % Display confirmation message
% fprintf('Release file successfully created: %s\n', fileName);
% fprintf('Total number of release dates: %d (daily releases)\n', num_releases);

% Open a file for writing
currentDateTime = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
currentDateTimeStr = string(currentDateTime);
fileName = "ReleaseFile_USVI_2019_daily_" + currentDateTimeStr + ".txt";
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
    fprintf(fileID, '%-4.0f %-15.9f %-15.9f %-2.0f %-3.0f %-6.0f %-2.0f %-3.0f %-2.0f\n', chunk_data');
    
    % Progress indicator
    if mod(chunk, 10) == 0 || chunk == num_chunks
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
fprintf('Release file successfully created: %s\n', fileName);
fprintf('Total number of release dates: %d (daily releases)\n', num_releases);
fprintf('Total number of release points: %d\n', numpoints);
fprintf('Total release events written: %d\n', total_rows_actual);