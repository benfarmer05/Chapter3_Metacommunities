% Get the project root directory
projectPath = matlab.project.rootProject().RootFolder;

% Define paths relative to the project root
dataPath = fullfile(projectPath, 'data');
srcPath = fullfile(projectPath, 'src');
outputPath = fullfile(projectPath, 'output');
tempPath = fullfile(projectPath, 'temp');

% Get the original file (using your existing logic)
files = dir(fullfile(tempPath, 'HPC*'));
randomIndex = randi(length(files));
baseFileName_CMS_HPC = files(randomIndex).name;
originalFile = fullfile(tempPath, baseFileName_CMS_HPC);

% Create new filename
[~, name, ext] = fileparts(baseFileName_CMS_HPC);
newFileName = [name, '_renamed_dims', ext];
newFile = fullfile(tempPath, newFileName);

% Read all data from original file
fprintf('Reading original file: %s\n', originalFile);

% Read dimensions and variables
lon = ncread(originalFile, 'Longitude');
lat = ncread(originalFile, 'Latitude');
depth = ncread(originalFile, 'Depth');
time = ncread(originalFile, 'Time');

% Read all data variables
zu = ncread(originalFile, 'zu');
zv = ncread(originalFile, 'zv');
zw = ncread(originalFile, 'zw');
zt = ncread(originalFile, 'zt');
zs = ncread(originalFile, 'zs');

% Read global attributes
ncinfo_orig = ncinfo(originalFile);
global_attrs = ncinfo_orig.Attributes;

fprintf('Creating new file: %s\n', newFile);

% Delete file if it already exists
if exist(newFile, 'file')
    delete(newFile);
end

% Create new NetCDF file with renamed dimensions
nccreate(newFile, 'Longitude', 'Dimensions', {'Longitude', length(lon)}, 'Datatype', 'double');
nccreate(newFile, 'Latitude', 'Dimensions', {'Latitude', length(lat)}, 'Datatype', 'double');
nccreate(newFile, 'Depth', 'Dimensions', {'Depth', length(depth)}, 'Datatype', 'double');
nccreate(newFile, 'Time', 'Dimensions', {'Time', length(time)}, 'Datatype', 'double');

% Create data variables with new dimension names
nccreate(newFile, 'zu', 'Dimensions', {'Longitude', 'Latitude', 'Depth', 'Time'}, 'Datatype', 'single');
nccreate(newFile, 'zv', 'Dimensions', {'Longitude', 'Latitude', 'Depth', 'Time'}, 'Datatype', 'single');
nccreate(newFile, 'zw', 'Dimensions', {'Longitude', 'Latitude', 'Depth', 'Time'}, 'Datatype', 'single');
nccreate(newFile, 'zt', 'Dimensions', {'Longitude', 'Latitude', 'Depth', 'Time'}, 'Datatype', 'single');
nccreate(newFile, 'zs', 'Dimensions', {'Longitude', 'Latitude', 'Depth', 'Time'}, 'Datatype', 'single');

% Write coordinate variables
ncwrite(newFile, 'Longitude', lon);
ncwrite(newFile, 'Latitude', lat);
ncwrite(newFile, 'Depth', depth);
ncwrite(newFile, 'Time', time);

% Write data variables
ncwrite(newFile, 'zu', zu);
ncwrite(newFile, 'zv', zv);
ncwrite(newFile, 'zw', zw);
ncwrite(newFile, 'zt', zt);
ncwrite(newFile, 'zs', zs);

% Add global attributes
for i = 1:length(global_attrs)
    ncwriteatt(newFile, '/', global_attrs(i).Name, global_attrs(i).Value);
end

% Add variable attributes
% Get original variable info
vars_info = ncinfo_orig.Variables;

for i = 1:length(vars_info)
    var_name = vars_info(i).Name;
    var_attrs = vars_info(i).Attributes;
    
    for j = 1:length(var_attrs)
        ncwriteatt(newFile, var_name, var_attrs(j).Name, var_attrs(j).Value);
    end
end

fprintf('Successfully created new file with renamed dimensions!\n');
fprintf('Original file: %s\n', originalFile);
fprintf('New file: %s\n', newFile);

% Verify the new file
fprintf('\nVerifying new file structure:\n');
ncdisp(newFile);