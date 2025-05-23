%% HOW I WANT TO TALK TO THE DATA
% Get the project root directory
projectPath = matlab.project.rootProject().RootFolder;
% Define paths relative to the project root
tempPath = fullfile(projectPath, 'temp');
outputPath = fullfile(projectPath, 'output'); % or wherever your shapefile is located

%% SECTION THAT NEEDS TO BE UPDATED FOR THE ABOVE
%%Convert netcdf data to structure
trajlist = dir(fullfile(tempPath,'traj*'));
% Initialize the structure
bigstruct_0 = struct();

% Determine the appropriate format based on number of files
if length(trajlist) >= 100
    format_str = '%03d'; % 3-digit format for 100+ files
else
    format_str = '%02d'; % 2-digit format for <100 files
end

% Loop over the range of NetCDF files
for i = 1:length(trajlist)
    % Construct the filename with full path using dynamic format
    filename = fullfile(tempPath, sprintf(['traj_file_' format_str '.nc'], i));
    
    % Check if the file exists
    if isfile(filename)
        % Read the variables from the NetCDF file
        bigstruct_0(i).time = ncread(filename, 'time');
        bigstruct_0(i).location = ncread(filename, 'location');
        bigstruct_0(i).lon = ncread(filename, 'lon');
        bigstruct_0(i).lat = ncread(filename, 'lat');
        bigstruct_0(i).depth = ncread(filename, 'depth');
        bigstruct_0(i).distance = ncread(filename, 'distance');
        bigstruct_0(i).exitcode = ncread(filename, 'exitcode');
        bigstruct_0(i).releasedate = ncread(filename, 'releasedate');
    else
        warning('File %s does not exist.', filename);
    end
end

% bigstruct = bigstruct_0
%Can create bigstruct immediatly if the number of release lines = number of
%traj files (not the case most of the time, so use the next loop)
%Depending on how many nodes/cores on the HPC, may need to adjust
%structure. This makes it so that each release line is a line in the
%structure.
bigstruct = struct();
% bigstruct = struct('time', {}, 'location', {}, 'lon', {},'lat', {},'depth', {}, 'distance', {}, 'exitcode', {}, 'releasedate', {});
row_no = 0;
for i = 1:length(bigstruct_0)
    locs_in_file = unique(bigstruct_0(i).location);
    for j = 1:length(locs_in_file)
        row_no = row_no + 1;
        r_index = bigstruct_0(i).location == locs_in_file(j);
        bigstruct(row_no).time = bigstruct_0(i).time;
        bigstruct(row_no).location = bigstruct_0(i).location(r_index);
        bigstruct(row_no).lon = bigstruct_0(i).lon(:,r_index);
        bigstruct(row_no).lat = bigstruct_0(i).lat(:,r_index);
        bigstruct(row_no).depth = bigstruct_0(i).depth(:,r_index);
        bigstruct(row_no).distance = bigstruct_0(i).distance(:,r_index);
        bigstruct(row_no).exitcode = bigstruct_0(i).exitcode(r_index);
        bigstruct(row_no).releasedate = bigstruct_0(i).releasedate(r_index);
    end
end
clear bigstruct_0
% Save to temp directory using the project-relative path
save(fullfile(tempPath, 'traj_all.mat'), 'bigstruct', '-v7.3') %just saves the structure

%% Quick figure to visualize trajectory files with particle skipping
% SET PARTICLE SKIP PARAMETER HERE
skip_particles = 100; % Plot every 100th particle (adjust as needed)

figure()
hold on;
% Read the shapefile and plot landmask first
landmask = shaperead(fullfile(outputPath, 'landmask_dissolved.shp'));
% Plot the landmask
mapshow(landmask);

% Calculate actual number of particles
total_release_locations = length(bigstruct);
particles_per_location = size(bigstruct(1).lon, 2); % Second dimension = number of particles
total_particles = total_release_locations * particles_per_location;

% Plot trajectory data with skipping - keep in same coordinate system as landmask (0-360)
fprintf('Total release locations: %d\n', total_release_locations);
fprintf('Particles per location: %d\n', particles_per_location);
fprintf('Total particles: %d\n', total_particles);
fprintf('Plotting every %dth release location (each contains %d particles)\n', skip_particles, particles_per_location);

locations_plotted = 0;
for i = 1:skip_particles:length(bigstruct) % Skip release locations based on skip_particles parameter
    plot([bigstruct(i).lon], [bigstruct(i).lat]);
    locations_plotted = locations_plotted + 1;
end
particles_plotted = locations_plotted * particles_per_location;
fprintf('Release locations plotted: %d\n', locations_plotted);
fprintf('Total particles plotted: %d\n', particles_plotted);

axis equal
axis([294, 296, 17.6, 19.2]) % Adjusted to match the coordinate system shown in your images
title(sprintf('Trajectory Data with Landmask (Every %dth release location)', skip_particles));
hold off;


%% Figure colored by release month
figure()
hold on;
% Plot the landmask first
mapshow(landmask);

% Create a colormap with 12 colors showing progression from "old" to "new"
colors = parula(12); % Parula goes from dark blue (old) to bright yellow (new)

% Extract release months from all trajectories
release_months = zeros(length(bigstruct), 1);
for i = 1:length(bigstruct)
    % Convert Julian day to month
    % Julian day 1 = January 1, 0000 (proleptic Gregorian calendar)
    % MATLAB's datenum uses January 1, 0000 as day 1
    julian_day = bigstruct(i).releasedate(1); % Take first value (should be same for all particles in this release)
    
    % Convert Julian day to MATLAB datenum and extract month
    % Note: You may need to adjust the Julian day conversion based on your specific Julian day reference
    % Common conversions: Julian Day Number (JDN) typically starts from January 1, 4713 BCE
    % If your Julian days are from a different epoch, adjust accordingly
    
    % Assuming standard Julian Day Number (adjust if needed)
    matlab_datenum = julian_day - 1721425.5; % Convert JDN to MATLAB datenum
    date_vec = datevec(matlab_datenum);
    release_months(i) = date_vec(2); % Extract month (1-12)
end

% Plot trajectories colored by release month with particle skipping
fprintf('Plotting trajectories colored by release month...\n');
locations_plotted = 0;
for i = 1:skip_particles:length(bigstruct)
    month = release_months(i);
    color = colors(month, :);
    plot([bigstruct(i).lon], [bigstruct(i).lat], 'Color', color, 'LineWidth', 0.5);
    locations_plotted = locations_plotted + 1;
end

% Create colorbar with January at top, December at bottom
colormap(colors);
cb = colorbar;
caxis([1 12]); % Set color axis limits
cb.Ticks = 1:12; % Set ticks at integer month values
cb.TickLabels = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
                 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
cb.Label.String = 'Release Month';

axis equal
axis([294, 296, 17.6, 19.2])
title(sprintf('Trajectory Data by Release Month (Every %dth release location)', skip_particles));
hold off;



%% Figure showing particle movement through depth
figure()
hold on;
% Plot the landmask first
mapshow(landmask);

% Create a bathymetry-inspired colormap (similar to cmocean depth)
% Deep blue for deep water, light blue/cyan for shallow water
n_depth_colors = 256;
% Create custom colormap: deep blue to light cyan
deep_color = [0.05, 0.1, 0.4];    % Deep blue for deep water  
shallow_color = [0.7, 0.9, 1.0];  % Light cyan for shallow water
depth_colormap = interp1([0 1], [deep_color; shallow_color], linspace(0, 1, n_depth_colors));

% % Find depth range across all trajectories for color scaling
% all_depths = [];
% for i = 1:length(bigstruct)
%     depths = bigstruct(i).depth(:);
%     depths = depths(~isnan(depths)); % Remove NaN values
%     all_depths = [all_depths; depths];
% end

% METHOD 1: Pre-allocate and use cell arrays (Recommended)
% First pass to count total non-NaN depths
total_depth_count = 0;
for i = 1:length(bigstruct)
    depths = bigstruct(i).depth(:);
    total_depth_count = total_depth_count + sum(~isnan(depths));
end

% Pre-allocate array
all_depths = zeros(total_depth_count, 1);
idx = 1;

% Second pass to fill the array
for i = 1:length(bigstruct)
    depths = bigstruct(i).depth(:);
    valid_depths = depths(~isnan(depths));
    
    if ~isempty(valid_depths)
        end_idx = idx + length(valid_depths) - 1;
        all_depths(idx:end_idx) = valid_depths;
        idx = end_idx + 1;
    end
end

depth_min = min(all_depths);
depth_max = max(all_depths);

fprintf('Depth range: %.2f to %.2f meters\n', depth_min, depth_max);
fprintf('Plotting trajectories colored by depth (Every %dth release location)...\n', skip_particles);

% % Plot trajectories with depth-based coloring
% locations_plotted = 0;
% for i = 1:skip_particles:length(bigstruct)
%     lons = bigstruct(i).lon;
%     lats = bigstruct(i).lat;
%     depths = bigstruct(i).depth;
% 
%     % Remove time steps where any coordinate is NaN
%     valid_idx = ~isnan(lons) & ~isnan(lats) & ~isnan(depths);
%     lons = lons(valid_idx);
%     lats = lats(valid_idx);
%     depths = depths(valid_idx);
% 
%     if length(lons) > 1
%         % Plot trajectory segments with colors based on depth
%         for j = 1:(length(lons)-1)
%             % Use depth at current point for segment color
%             depth_norm = (depths(j) - depth_min) / (depth_max - depth_min);
%             depth_norm = max(0, min(1, depth_norm)); % Clamp to [0,1]
% 
%             % Get color from colormap
%             color_idx = round(depth_norm * (n_depth_colors - 1)) + 1;
%             color_idx = max(1, min(n_depth_colors, color_idx));
%             segment_color = depth_colormap(color_idx, :);
% 
%             % Plot segment
%             plot([lons(j), lons(j+1)], [lats(j), lats(j+1)], ...
%                  'Color', segment_color, 'LineWidth', 1);
%         end
%     end
%     locations_plotted = locations_plotted + 1;
% end

% Pre-allocate arrays for all trajectory data
all_lons = [];
all_lats = [];
all_colors = [];

for i = 1:skip_particles:length(bigstruct)
    lons = bigstruct(i).lon;
    lats = bigstruct(i).lat;
    depths = bigstruct(i).depth;
    
    % Remove invalid points
    valid_idx = ~isnan(lons) & ~isnan(lats) & ~isnan(depths);
    lons = lons(valid_idx);
    lats = lats(valid_idx);
    depths = depths(valid_idx);
    
    if length(lons) > 1
        % Add NaN at the end to separate trajectories
        lons = [lons; NaN];
        lats = [lats; NaN];
        depths = [depths; NaN];
        
        % Accumulate data
        all_lons = [all_lons; lons];
        all_lats = [all_lats; lats];
        all_colors = [all_colors; depths];
    end
end

% Plot all trajectories at once using scatter plot with lines
% Method A: Use surface plot for colored lines
figure()
hold on;
mapshow(landmask);

% Create line plot with color data
surface([all_lons'; all_lons'], [all_lats'; all_lats'], ...
        zeros(2, length(all_lons)), [all_colors'; all_colors'], ...
        'EdgeColor', 'flat', 'FaceColor', 'none', 'LineWidth', 1);
% Set up colorbar for depth
colormap(depth_colormap);
caxis([depth_min, depth_max]);
cb = colorbar;
cb.Label.String = 'Depth (m)';

% Add some formatting for depth values in colorbar
if abs(depth_max - depth_min) > 1000
    cb.TickLabels = arrayfun(@(x) sprintf('%.0f', x), cb.Ticks, 'UniformOutput', false);
else
    cb.TickLabels = arrayfun(@(x) sprintf('%.1f', x), cb.Ticks, 'UniformOutput', false);
end

axis equal
axis([294, 296, 17.6, 19.2])
title(sprintf('Particle Trajectories Colored by Depth (Every %dth release location)', skip_particles));
xlabel('Longitude');
ylabel('Latitude');
hold off;