%% Script to create a release file for the CMS, simulating the dispersal and connectivity of SCTLD in the Virgin Islands & Puerto Rico
% 21 May 2024

clear;clc

%% test

% .nc file format: %nest_nestnumber_yyyymmddhhmmss [.nc]
%   here, I'm working with 'nest_1_20190101000000' through
%                          'nest_1_20190101200000'
%   so,                    'nest_1_2019_January_1_00:00' through
%                          'nest_1_2019_January_1_20:00'

% Release file format:
%   Polygon Longitude Latitude Depth Number Year Month Day Second
%   1       277.2     24.66    0     10     2016 1     1   0

scriptDrive = '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/Release';

cd(scriptDrive)

% Read in release points from GIS output and ensure they are sorted by
% their unique ID
relpoints = readmatrix('points_650_none-on-land.csv');
relpoints = relpoints(:, 10:12);
relpoints = sortrows(relpoints, 1);

IDs = relpoints(:,1);
longitudes = relpoints(:,2) + 360;
latitudes = relpoints(:,3);

numpoints = size(relpoints, 1);  % Number of points

% Construct bones of the release structure
repvalue = 7; %in days
num_releases = ceil(365 / repvalue); % Calculate total number of releases based on repvalue

% Generate days and months for a non-leap year (2019)
days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
num_days = sum(days_in_month); % Total number of days in a year

% Pre-initialize months and days arrays
months = zeros(num_days, 1);  % Pre-initialize months as a column vector of zeros
days = zeros(num_days, 1);    % Pre-initialize days as a column vector of zeros
index = 1; % Initialize index to keep track of the position in the days/months arrays

% Fill in the days and months arrays:
for month = 1:length(days_in_month)
    days(index:index + days_in_month(month) - 1) = 1:days_in_month(month); % Generate days for the month
    months(index:index + days_in_month(month) - 1) = month; % Set the month number
    index = index + days_in_month(month); % Update index for the next month
end

% Generate a continuous day count for the entire year
day_of_year = (1:num_days)';

% Filter to select every repvalue-th day
filtered_day_index = mod(day_of_year - 1, repvalue) == 0;
filtered_days = days(filtered_day_index);

% Corresponding months for the filtered days
filtered_months = months(filtered_day_index);



%%
longitudes_CMS = repmat(longitudes, num_releases, 1);
latitudes_CMS = repmat(latitudes, num_releases, 1);



% Repeat the filtered days and months data to match the number of points
days_CMS = repmat(filtered_days, numpoints, 1);
months_CMS = repmat(filtered_months, numpoints, 1);

%%

% Convert IDs to strings
IDs_str = string(IDs);

% Repeat the padded IDs to match the length of days and months
IDs_CMS = repmat(IDs_str, num_releases, 1);

% Convert numeric data to strings with consistent formatting
longitudes_str = string(longitudes_CMS);
latitudes_str = string(latitudes_CMS);

% Initialize release matrix
release = cell(numpoints * num_releases, 9);

% Fill in release matrix with pre-allocated columns
for i = 1:size(release, 1)
    release{i, 1} = IDs_CMS(i); % IDs
    release{i, 2} = longitudes_str(i); % Longitudes
    release{i, 3} = latitudes_str(i); % Latitudes
    release{i, 4} = "1"; % In meters. Releasing at surface
    release{i, 5} = "10"; % #/particles
    release{i, 6} = "2019"; % Year
    release{i, 7} = string(months_CMS(i)); % Months
    release{i, 8} = string(days_CMS(i)); % Days
    release{i, 9} = "0"; % In seconds. Releasing at midnight
end

% Write the matrix to the file with the constructed file name
currentDateTime = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
currentDateTimeStr = string(currentDateTime);
fileName = "ReleaseFile_tester_" + currentDateTimeStr + ".txt";

% Open the file for writing
fileID = fopen(fileName, 'w');

% Write data to the file
for i = 1:size(release, 1)
    % Write each column to the file
    for j = 1:size(release, 2)
        fprintf(fileID, '%s', release{i, j});
        
        % Add tab delimiter between columns (except for the last column)
        if j < size(release, 2)
            fprintf(fileID, '\t');
        end
    end
    % Add newline character at the end of each row
    fprintf(fileID, '\n');
end

% Close the file
fclose(fileID);