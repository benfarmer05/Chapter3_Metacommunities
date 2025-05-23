%% Script to create a release file for the Connectivity Modeling System
%   simulating the dispersal and connectivity of SCTLD in the Virgin
%   Islands & Puerto Rico
%   22 May 2025

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
relpoints = relpoints(:, 10:12);
relpoints = sortrows(relpoints, 1);

IDs = relpoints(:,1);
longitudes = relpoints(:,2) + 360;
latitudes = relpoints(:,3);

numpoints = size(relpoints, 1);  % #/points

% Generate release dates for 1st and 15th of each month
release_dates = [];
current_date = startDate;

while current_date <= endDate
    % Add 1st of the month if within date range
    first_of_month = datetime(year(current_date), month(current_date), 1);
    if first_of_month >= startDate && first_of_month <= endDate
        release_dates = [release_dates; first_of_month];
    end
    
    % Add 15th of the month if within date range
    fifteenth_of_month = datetime(year(current_date), month(current_date), 15);
    if fifteenth_of_month >= startDate && fifteenth_of_month <= endDate
        release_dates = [release_dates; fifteenth_of_month];
    end
    
    % Move to next month
    if month(current_date) == 12
        current_date = datetime(year(current_date) + 1, 1, 1);
    else
        current_date = datetime(year(current_date), month(current_date) + 1, 1);
    end
end

% Sort release dates
release_dates = sort(release_dates);
num_releases = length(release_dates);

%define the release structure
num_particles = 1; % #/particles released per line

% Extract year, month, and day from release dates
release_years = year(release_dates);
release_months = month(release_dates);
release_days = day(release_dates);

% Create arrays properly aligned for each combination of point and time
% For each release day, repeat all points
IDs_CMS = repmat(IDs, num_releases, 1);
longitudes_CMS = repmat(longitudes, num_releases, 1);
latitudes_CMS = repmat(latitudes, num_releases, 1);
depths_CMS = ones(numpoints * num_releases, 1); %this will actually be reef depth (and thus affix directly to 'ID', so will need to be defined and applied 'repmat')
particles_CMS = repmat(num_particles, numpoints * num_releases, 1); %this could also vary by ID / time / location, etc.

% For temporal data, repeat each date for all points
days_CMS = [];
months_CMS = [];
years_CMS = [];
for i = 1:num_releases
    days_CMS = [days_CMS; repmat(release_days(i), numpoints, 1)];
    months_CMS = [months_CMS; repmat(release_months(i), numpoints, 1)];
    years_CMS = [years_CMS; repmat(release_years(i), numpoints, 1)];
end

time_CMS = zeros(numpoints * num_releases, 1); %can choose to release at different times of day if desired

% Create matrix for sorting: [Year, Month, Day, Time, ID, Longitude, Latitude, Depth, Particles]
sortMatrix = [years_CMS, months_CMS, days_CMS, time_CMS, IDs_CMS, longitudes_CMS, latitudes_CMS, depths_CMS, particles_CMS];

% Sort by Year, Month, Day, Time, then ID
sortMatrix = sortrows(sortMatrix, [1, 2, 3, 4, 5]);

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

% Create a cell array to store all data as strings with proper formatting
numRows = length(IDs_CMS);
outputData = cell(numRows, 9);

% Convert each column to properly formatted strings
for i = 1:numRows
    outputData{i, 1} = sprintf('%d', IDs_CMS(i));                   % ID - integer
    outputData{i, 2} = sprintf('%.9f', longitudes_CMS(i));          % Longitude - float with 9 decimal places
    outputData{i, 3} = sprintf('%.9f', latitudes_CMS(i));           % Latitude - float with 9 decimal places
    outputData{i, 4} = sprintf('%d', depths_CMS(i));                % Depth - integer
    outputData{i, 5} = sprintf('%d', particles_CMS(i));             % Particles - integer
    outputData{i, 6} = sprintf('%d', year_CMS(i));                  % Year - integer
    outputData{i, 7} = sprintf('%d', months_CMS(i));                % Month - integer
    outputData{i, 8} = sprintf('%d', days_CMS(i));                  % Day - integer
    outputData{i, 9} = sprintf('%d', time_CMS(i));                  % Time - integer
end

% Open a file for writing
currentDateTime = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
currentDateTimeStr = string(currentDateTime);
fileName = "ReleaseFile_USVI_2019_1st15th_" + currentDateTimeStr + ".txt";
fileID = fopen(fullfile(outputPath, fileName), 'w');

% Write formatted data to file with consistent spacing
for i = 1:numRows
    fprintf(fileID, '%-4s %-15s %-15s %-2s %-3s %-6s %-2s %-3s %-2s\n', ...
        outputData{i, 1}, ...         % ID (Polygon)
        outputData{i, 2}, ...         % Longitude
        outputData{i, 3}, ...         % Latitude
        outputData{i, 4}, ...         % Depth
        outputData{i, 5}, ...         % Particles (Number)
        outputData{i, 6}, ...         % Year
        outputData{i, 7}, ...         % Month
        outputData{i, 8}, ...         % Day
        outputData{i, 9});            % Second (Time)
end

% Close the file
fclose(fileID);

% Display confirmation message
fprintf('Release file successfully created: %s\n', fileName);
fprintf('Total number of release dates: %d (1st and 15th of each month)\n', num_releases);