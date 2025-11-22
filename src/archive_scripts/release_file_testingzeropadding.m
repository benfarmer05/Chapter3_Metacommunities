
clear;clc

scriptDrive = '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/Release';

cd(scriptDrive)

% Read in release points from GIS output and ensure they are sorted by
% their unique ID
relpoints = readmatrix('points_650_none-on-land.csv');
relpoints = relpoints(:, 10:12);
relpoints = sortrows(relpoints, 1);

IDs = relpoints(:,1);

%%
% Convert IDs to strings
IDs_str = string(IDs);

% Determine the maximum length of IDs
max_length = max(strlength(IDs_str));

% Pad IDs with leading zeros 
IDs_str_padded = cellstr(IDs_str); % Convert to cell array of strings
for i = 1:numel(IDs_str_padded)
    IDs_str_padded{i} = char(sprintf('%0*d', max_length, str2double(IDs_str_padded{i})));
end
IDs_str_padded = string(IDs_str_padded); % Convert back to string array

IDs = IDs_str_padded;

%%
longitudes = relpoints(:,2) + 360;
latitudes = relpoints(:,3);


% Fill in release file matrices with pre-allocated columns
release(:,1) = IDs;
release(:,2) = longitudes;
release(:,3) = latitudes;

% Write the matrix to the file with the constructed file name
currentDateTime = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
currentDateTimeStr = string(currentDateTime);
fileName = "ReleaseFile_tester_" + currentDateTimeStr + ".txt";
writematrix(release, fileName, 'delimiter', '\t');