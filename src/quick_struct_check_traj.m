%% Load one random trajectory file and explore its structure

% Get list of trajectory files
trajlist = dir(fullfile(tempPath,'traj*.nc'));

% Select one random file
random_idx = randi(length(trajlist));
random_file = fullfile(tempPath, trajlist(random_idx).name);

fprintf('Loading: %s\n', trajlist(random_idx).name);

% Get file info
file_info = ncinfo(random_file);

% Create structure with all variables
traj = struct();

for i = 1:length(file_info.Variables)
    var_name = file_info.Variables(i).Name;
    try
        traj.(var_name) = ncread(random_file, var_name);
        fprintf('Loaded: %s (size: %s)\n', var_name, mat2str(size(traj.(var_name))));
    catch
        fprintf('Could not load: %s\n', var_name);
    end
end

% Display the structure
disp(traj);