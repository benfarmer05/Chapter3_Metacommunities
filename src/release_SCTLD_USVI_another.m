%% Script to create a release file for the Connectivity Modeling System
%   simulating the dispersal and connectivity of SCTLD in the Virgin
%   Islands & Puerto Rico
%   22 May 2024

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

% Convert numeric arrays to strings
data = {IDs_CMS, longitudes_CMS, latitudes_CMS, depths_CMS, particles_CMS, year_CMS, months_CMS, days_CMS, time_CMS};
data_str = cellfun(@string, data, 'UniformOutput', false);

% Function to pad strings with leading zeros to the maximum length in the array
pad_with_zeros = @(str_array) compose("%0" + max(strlength(str_array)) + "d", str2double(str_array));

% Apply zero-padding to each column and convert to string arrays
data_str_padded = cellfun(@(col) string(pad_with_zeros(col)), data_str, 'UniformOutput', false);

% Manually handle longitudes and latitudes to ensure they are not in scientific notation
longitudes_str = arrayfun(@(x) sprintf('%.16g', x), longitudes_CMS, 'UniformOutput', false);
latitudes_str = arrayfun(@(x) sprintf('%.16g', x), latitudes_CMS, 'UniformOutput', false);
% longitudes_str_padded = cellfun(@(col) string(pad_with_zeros(col)), longitudes_str, 'UniformOutput', false);
% latitudes_str_padded = cellfun(@(col) string(pad_with_zeros(col)), latitudes_str, 'UniformOutput', false);

% Determine the maximum string length for longitudes and latitudes
max_length_longitude = max(strlength(longitudes_str));
max_length_latitude = max(strlength(latitudes_str));

% Apply zero-padding to longitudes and latitudes at the end and convert to string arrays
longitudes_str_padded = cellfun(@(col) sprintf('%s%0.*d', col, max_length_longitude - numel(col), 0), longitudes_str, 'UniformOutput', false);
latitudes_str_padded = cellfun(@(col) sprintf('%s%0.*d', col, max_length_latitude - numel(col), 0), latitudes_str, 'UniformOutput', false);

% Replace the longitudes and latitudes in the padded data
data_str_padded{2} = string(longitudes_str_padded);
data_str_padded{3} = string(latitudes_str_padded);

% Combine padded columns into the release matrix
release = [data_str_padded{:}];

%write the matrix to the file with the constructed file name
currentDateTime = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
currentDateTimeStr = string(currentDateTime);
fileName = "ReleaseFile_USVI_2019_" + currentDateTimeStr + ".txt";
writematrix(release, fileName, 'delimiter', '\t');


% STOPPING POINT 16 may 2025 - adjusting the below to match the above, but
% with better pathing


% Create the full path for your file
fileName = fullfile(outputPath, "ReleaseFile_USVI_2019_" + currentDateTimeStr + ".txt");

% Write your data
writematrix(release, fileName, 'delimiter', '\t');