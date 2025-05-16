%% Script (requires output from 'grid_disassembly_VIRRS_PREDICT.m') that
%   produces USCROMS bottom-layer u- and v-velocities in CSV format for
%   the NSF-funded UVI VIRRS study
% 6 June 2024

clear; clc;

% Define directories
inputDirs = {'/Volumes/UVI_Hydro_2019_1-of-2', '/Volumes/UVI_Hydro_2019-2020/USCROMS_inputs'};
outputDir = '/Users/benja';
scriptDrive = '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/Release';
indicesFile = 'USCROMS_VIRRS_indices.mat';
errorLogFile = fullfile(outputDir, 'error_log.txt');
outputCSV_u = fullfile(outputDir, 'velocity_data_u_2019.csv');
outputCSV_v = fullfile(outputDir, 'velocity_data_v_2019.csv');

% Load the indices from the .mat file
cd(scriptDrive)
load(indicesFile, 'idx_rho_USCROMS_VIRRS', 'idx_u_USCROMS_VIRRS', 'idx_v_USCROMS_VIRRS');

% Initialize parameters
numLocations = numel(idx_rho_USCROMS_VIRRS);
numVariables = 2; % u, v (excluding w)
maxHours = 8760;

% Helper function to extract the hour from the file name
extractHourFromFileName = @(filename) str2double(filename(11:15));

% Initialize CSV files if they don't exist
if ~isfile(outputCSV_u) || ~isfile(outputCSV_v)
    % Initialize the headers and NaN-filled matrix
    header_u = ["Hour", strcat("Loc_", string(1:numLocations), "_u")];
    header_v = ["Hour", strcat("Loc_", string(1:numLocations), "_v")];

    data_u = array2table(NaN(maxHours, numLocations + 1), 'VariableNames', header_u);
    data_v = array2table(NaN(maxHours, numLocations + 1), 'VariableNames', header_v);

    % Initialize the hour column starting from 0
    data_u.Hour = (0:maxHours-1)';
    data_v.Hour = (0:maxHours-1)';

    % Write the initialized tables to the CSV files
    writetable(data_u, outputCSV_u);
    writetable(data_v, outputCSV_v);
end

% Load existing data to find the last written hour
data_u = readtable(outputCSV_u);
data_v = readtable(outputCSV_v);

% Determine the last written hour by finding the last non-NaN row
lastFilledHour_u = find(~isnan(data_u{:, 2}), 1, 'last');
lastFilledHour_v = find(~isnan(data_v{:, 2}), 1, 'last');

% Adjust the hour counter and file hour to start from 0
hourCounter = max([lastFilledHour_u, lastFilledHour_v, 0]);

% Loop through each input directory
for dirIdx = 1:numel(inputDirs)
    inputDir = inputDirs{dirIdx};
    cd(inputDir)

    % Get list of all ROMS files in the directory
    romsFiles = dir(fullfile(inputDir, 'croco_his.*.nc'));
    romsFiles = {romsFiles.name};

    % Sort files based on the hour extracted from file names
    [~, sortedIdx] = sort(cellfun(extractHourFromFileName, romsFiles));
    romsFiles = romsFiles(sortedIdx);

    % Loop through each ROMS file
    for i = 1:numel(romsFiles)
        try
            % Extract the hour from the file name starting from 0
           fileHour = extractHourFromFileName(romsFiles{i});

            % Skip if this hour has already been processed
            if fileHour < hourCounter
                continue;
            end

            % Display progress
            fprintf('Processing file %d of %d in directory %d: %s\n', i, numel(romsFiles), dirIdx, romsFiles{i});

            % Construct full file path
            filePath = fullfile(inputDir, romsFiles{i});

            % Read the dimensions of the u and v velocities
            uvel = ncread(filePath, 'u');
            vvel = ncread(filePath, 'v');

            % Get the number of time steps in the current file
            numTimeSteps = size(uvel, 4);

            % Loop through each time step
            for t = 1:numTimeSteps
                % Extract the bottom layer velocities for each time step
                uBottom = squeeze(uvel(:,:,1,t));
                vBottom = squeeze(vvel(:,:,1,t));

                % Prepare rows for the current time step
                row_u = NaN(1, numLocations);
                row_v = NaN(1, numLocations);

                % Extract the velocities at the specified indices
                for loc = 1:numLocations
                    uIdx = idx_u_USCROMS_VIRRS(loc);
                    vIdx = idx_v_USCROMS_VIRRS(loc);

                    if ~isnan(uIdx)
                        row_u(loc) = uBottom(uIdx);
                    end
                    if ~isnan(vIdx)
                        row_v(loc) = vBottom(vIdx);
                    end
                end

                % Write the current rows to the CSV files starting from 0
                data_u{hourCounter + 1, 2:end} = row_u;
                data_v{hourCounter + 1, 2:end} = row_v;

                % Increment the hour counter
                hourCounter = hourCounter + 1;

                % Check if the maximum hours have been reached
                if hourCounter >= maxHours
                    break;
                end
            end

            % Write the updated data to the CSV files after processing each file
            writetable(data_u, outputCSV_u);
            writetable(data_v, outputCSV_v);

            % Check if the maximum hours have been reached
            if hourCounter >= maxHours
                break;
            end

        catch ME
            % Log the error message to a text file
            fid = fopen(errorLogFile, 'a');
            if fid == -1
                error('Cannot open log file for writing: %s', errorLogFile);
            end
            fprintf(fid, 'Error processing file %s: %s\n', filePath, ME.message);
            fclose(fid);
            % Continue with the next file
            continue;
        end
    end

    % Check if the maximum hours have been reached
    if hourCounter >= maxHours
        break;
    end
end

% Display completion message
fprintf('Processing complete. Data saved to %s and %s\n', outputCSV_u, outputCSV_v);

%% testing that the output worked

% clear%;clc
% 
% % Define directories
% inputDirs = {'/Volumes/UVI_Hydro_2019_1-of-2', '/Volumes/UVI_Hydro_2019-2020/USCROMS_inputs'};
% outputDir = '/Users/benja';
% scriptDrive = '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/Release';
% indicesFile = 'USCROMS_VIRRS_indices.mat';
% 
% % Define paths to CSV files
% outputCSV_u = '/Users/benja/velocity_data_u_2019.csv';
% outputCSV_v = '/Users/benja/velocity_data_v_2019.csv';
% 
% % Load data from CSV files
% data_u = readtable(outputCSV_u);
% data_v = readtable(outputCSV_v);
% 
% % Find rows with no NaN values in the velocity data (excluding the first row which is the header)
% validRows_u = ~any(isnan(data_u{1:end, :}), 2);
% validRows_v = ~any(isnan(data_v{1:end, :}), 2);
% 
% % Get the valid rows (where both u and v data have no NaNs)
% validRows = validRows_u & validRows_v;
% 
% % Get the number of locations (excluding the first column which is the hour)
% numLocations = size(data_u, 2) - 1; % Subtract 1 for the hour column
% 
% % Randomly select a location
% randomLocation = randi(numLocations);
% 
% % Randomly select a row index from the valid rows
% validIndices = find(validRows);
% randomRowIndex = validIndices(randi(length(validIndices)));
% 
% % %test a specific hour & location:
% % testHour = 0; %8633 for the last one
% % randomRowIndex = testHour + 1; %account for zero-indexing of hours
% % randomLocation = 1;
% % %test a specific hour & location:
% 
% % Extract the hour value from the randomly selected row
% randomHour = data_u{randomRowIndex, 1}; % Add 1 to skip the header row
% 
% % Calculate the NetCDF file index based on the hour (divide by 3 and take the floor)
% netcdfIndex = floor(randomHour / 3) * 3;
% hoursequence = randomHour - netcdfIndex + 1;
% 
% % Read the velocity values from the CSV files
% actual_u = data_u{randomRowIndex, randomLocation + 1};
% actual_v = data_v{randomRowIndex, randomLocation + 1};
% 
% % Read the velocity values from the associated NetCDF file
% netcdfFileFound = false;
% for dirIndex = 1:length(inputDirs)
%     cd(inputDirs{dirIndex});
%     netcdfFilename = sprintf('croco_his.%05d.nc', netcdfIndex);
%     if exist(netcdfFilename, 'file') == 2
%         netcdfFileFound = true;
%         break; % Exit loop if file is found
%     end
% end
% 
% if ~netcdfFileFound
%     error('NetCDF file not found in any of the specified directories.');
% end
% 
% u_data = ncread(netcdfFilename, 'u');
% v_data = ncread(netcdfFilename, 'v');
% u_val = squeeze(u_data(:,:,1,hoursequence));
% v_val = squeeze(v_data(:,:,1,hoursequence));
% 
% % Load the indices from the .mat file
% cd(scriptDrive)
% load(indicesFile, 'idx_rho_USCROMS_VIRRS', 'idx_u_USCROMS_VIRRS', 'idx_v_USCROMS_VIRRS');
% 
% uIdx_rand = idx_u_USCROMS_VIRRS(randomLocation);
% vIdx_rand = idx_v_USCROMS_VIRRS(randomLocation);
% expected_u = u_val(uIdx_rand);
% expected_v = v_val(vIdx_rand);
% 
% % Determine the precision of actual_u and actual_v
% actual_u_precision = 10^(-floor(log10(abs(actual_u))) + 1);
% actual_v_precision = 10^(-floor(log10(abs(actual_v))) + 1);
% 
% % Compare the values
% fprintf('Random Hour: %d, Random Location: %d\n', randomHour, randomLocation);
% fprintf('Expected u: %f, Actual u: %f\n', expected_u, actual_u);
% fprintf('Expected v: %f, Actual v: %f\n', expected_v, actual_v);
% 
% % Check if the values match
% if abs(expected_u - actual_u) < eps && abs(expected_v - actual_v) < eps
%     disp('Values match!');
% else
%     disp('Values do not match!');
% end