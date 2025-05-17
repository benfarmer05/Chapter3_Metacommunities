%% Script to create a release file for the Connectivity Modeling System
%   simulating the dispersal and connectivity of SCTLD in the Virgin
%   Islands & Puerto Rico
%   17 May 2025

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
startDate = datetime(2019, 9, 1); % Example: Start date
endDate = datetime(2019, 9, 30);  % Example: End date
dateRange = startDate:endDate;

%read in release points from GIS output and ensure they are sorted by
% their unique ID
relpoints = readmatrix(fullfile(dataPath, 'points_650_none-on-land.csv'));
relpoints = relpoints(:, 10:12);
relpoints = sortrows(relpoints, 1);

IDs = relpoints(:,1);
longitudes = relpoints(:,2) + 360;
latitudes = relpoints(:,3);

numpoints = size(relpoints, 1);  % #/points

%define the release structure
repvalue = 1; % in days
num_releases = ceil(days(dateRange(end) - dateRange(1) + 1) / repvalue); % Calculate total # of releases based on repvalue

% Generate months and days for the specified date range
months = month(dateRange)';
days = day(dateRange)';

% Generate a continuous day count for the specified date range
day_of_year = (1:length(dateRange))';

%filter to select every repvalue-th day and month
filtered_day_index = mod(day_of_year - 1, repvalue) == 0;
filtered_days = days(filtered_day_index);
filtered_months = months(filtered_day_index);

%repeat the filtered days and months data to match the number of points
days_CMS = repmat(filtered_days, numpoints, 1);
months_CMS = repmat(filtered_months, numpoints, 1);

%generate IDs to match the length of days and months
IDs_CMS = repmat(IDs, num_releases, 1);
longitudes_CMS = repmat(longitudes, num_releases, 1);
latitudes_CMS = repmat(latitudes, num_releases, 1);
depths_CMS = ones(numpoints * num_releases, 1); %this will actually be reef depth (and thus affix directly to 'ID', so will need to be defined and applied 'repmat')
particles_CMS = repmat(10, numpoints * num_releases, 1); %this could also vary by ID / time / location, etc.
year_CMS = repmat(year(startDate), numpoints * num_releases, 1); %if more years are desired, should define and 'repmat' or similar above
time_CMS = zeros(numpoints * num_releases, 1); %can choose to release at different times of day if desired


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
fileName = "ReleaseFile_USVI_2019_" + currentDateTimeStr + ".txt";
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