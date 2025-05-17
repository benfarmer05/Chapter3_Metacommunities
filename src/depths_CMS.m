%% Script to assign depths to lat/lon reef points, that are sensible for
%   the CMS. may be adopted to jitter points off of land (if there is room
%   in their containing polygon square) or otherwise that will become part
%   of a separate script. once no points are remaining on land, then it is
%   a matter of locating their deepest nearest u/v/rho-grid depth that is
%   not below the seafloor. the final step is deciding the best way to
%   settle on the most useful depth - the CMS likely interpolates the
%   "seafloor" using its available grid system, so I probably want to
%   use some kind of interpolation scheme. still testing that out
% 17 May 2025

clear;clc

%% setup

% Get the project root directory
projectPath = matlab.project.rootProject().RootFolder;

% Define paths relative to the project root
dataPath = fullfile(projectPath, 'data');
srcPath = fullfile(projectPath, 'src');
outputPath = fullfile(projectPath, 'output');
tempPath = fullfile(projectPath, 'temp');

%%

% % NOTE - may need to return to this
% CMSfilesDrive = '/Users/benja/testCMS/tests_for_depth-assignment'; %'/Volumes/UVI_Hydro_2019-2020/CMS_inputs_nozeta'; %

% NOTE - may need to adjust first value [currently 0.26] to reflect what
% was done with roms2cms_A_grid.m
zlevels = [0.26 1 2 4 6 8 10 12 14 16 18 20 22 24 26 29 32 35 38 42 46 50 60 70 85 100 120 140 160 190 220 250]'; %32 layers down to 250 m depth, mid-depth hi-res

% NOTE / STOPPING POINT - for the below, was previously just pulling a
% random hydro file I'd created for the CMS. just need to adopt this
% approach with proper pathing, and also double-check how I created my
% hydro files

%randomly select a CMS-ready file
files = dir(fullfile(tempPath, 'nest_1_*'));
randomIndex = randi(length(files)); % generate a random index between 1 and the number of files
baseFileName = files(randomIndex).name;
ncdisp(fullfile(tempPath, baseFileName))  % Pass the full path to ncdisp

%CMS coordinate systems
longitudes_CMS = ncread(fullfile(tempPath, baseFileName), 'Longitude');
latitudes_CMS = ncread(fullfile(tempPath, baseFileName), 'Latitude');

% Create a grid of all longitude and latitude pairs
[lat_grid, lon_grid] = meshgrid(latitudes_CMS, longitudes_CMS);

% Flatten the grid matrices
lon_grid_indices = lon_grid(:);
lat_grid_indices = lat_grid(:);

%CMS seafloor masks
fill_value = single(1.2676506e+30);
mask_u_seafloor = ncread(fullfile(tempPath, baseFileName), 'zu');
mask_u_seafloor(mask_u_seafloor==fill_value) = 0; %seafloor is 0's
mask_u_seafloor(mask_u_seafloor ~= 0) = 1; %the ocean is 1's
mask_v_seafloor = ncread(fullfile(tempPath, baseFileName), 'zv'); 
mask_v_seafloor(mask_v_seafloor==fill_value) = 0;
mask_v_seafloor(mask_v_seafloor ~= 0) = 1;
mask_rho_seafloor = ncread(fullfile(tempPath, baseFileName), 'zw'); 
mask_rho_seafloor(mask_rho_seafloor==fill_value) = 0; 
mask_rho_seafloor(mask_rho_seafloor ~= 0) = 1;

%extract release point coordinates from GIS output and ensure they are
% sorted by their unique ID
relpoints = readmatrix(fullfile(dataPath, 'points_650_none-on-land.csv'));
relpoints = relpoints(:, 10:12);
relpoints = sortrows(relpoints, 1);

IDs_release = relpoints(:,1);
longitudes_release = relpoints(:,2) + 360;
latitudes_release = relpoints(:,3);

% %% testing - difference in seafloor masks
% 
% depthslice = 15;
% 
% mask_u_seafloor = mask_u_seafloor(:,:,depthslice);
% figure();imagescn(lon_grid',lat_grid',mask_u_seafloor');axis equal;hold on;plot(relpoints(:,2) + 360,relpoints(:,3), '*r');xlim([360-65.5, 360-65]);ylim([18, 18.4])
% 
% mask_v_seafloor = mask_v_seafloor(:,:,depthslice);
% figure();imagescn(lon_grid',lat_grid',mask_v_seafloor');axis equal;hold on;plot(relpoints(:,2) + 360,relpoints(:,3), '*r');xlim([360-65.5, 360-65]);ylim([18, 18.4])
% 
% mask_rho_seafloor = mask_rho_seafloor(:,:,depthslice);
% figure();imagescn(lon_grid',lat_grid',mask_rho_seafloor');axis equal;hold on;plot(relpoints(:,2) + 360,relpoints(:,3), '*r');xlim([360-65.5, 360-65]);ylim([18, 18.4])
% 
% %% testing - difference in seafloor masks

% Locate CMS grid coordinates nearest to PREDICT points
idx_CMS_release = knnsearch([lon_grid_indices(:), lat_grid_indices(:)], [longitudes_release, latitudes_release]);

% Initialize the seamasks for each depth level
seamask_u_release = zeros(length(idx_CMS_release), 32);
seamask_v_release = zeros(length(idx_CMS_release), 32);
seamask_rho_release = zeros(length(idx_CMS_release), 32);

% Iterate over each depth level
for depth = 1:32
    % Extract the mask for the current depth level
    mask_u_depth = mask_u_seafloor(:,:,depth);
    mask_v_depth = mask_v_seafloor(:,:,depth);
    mask_rho_depth = mask_rho_seafloor(:,:,depth);
    
    % Apply the mask for the current depth level
    seamask_u_release(:, depth) = mask_u_depth(idx_CMS_release);
    seamask_v_release(:, depth) = mask_v_depth(idx_CMS_release);
    seamask_rho_release(:, depth) = mask_rho_depth(idx_CMS_release);
end

% Find the deepest ocean depth for each point
deepest_depth_u = zeros(length(idx_CMS_release), 1);
deepest_depth_v = zeros(length(idx_CMS_release), 1);
deepest_depth_rho = zeros(length(idx_CMS_release), 1);

for i = 1:length(idx_CMS_release) %test: i = 2495 to find a NaN 'u'
    deepest_u_idx = find(seamask_u_release(i, :) == 1, 1, 'last');
    deepest_v_idx = find(seamask_v_release(i, :) == 1, 1, 'last');
    deepest_rho_idx = find(seamask_rho_release(i, :) == 1, 1, 'last');
    
    if ~isempty(deepest_u_idx)
        deepest_depth_u(i) = zlevels(deepest_u_idx);
    else
        deepest_depth_u(i) = NaN;
    end
    
    if ~isempty(deepest_v_idx)
        deepest_depth_v(i) = zlevels(deepest_v_idx);
    else
        deepest_depth_v(i) = NaN;
    end
    
    if ~isempty(deepest_rho_idx)
        deepest_depth_rho(i) = zlevels(deepest_rho_idx);
    else
        deepest_depth_rho(i) = NaN;
    end
end

%find the shallowest nearest-neighbor depth across each grid
combined_depths = [deepest_depth_u, deepest_depth_v, deepest_depth_rho];
combined_depths(isnan(combined_depths)) = inf;
deepest_depth = min(combined_depths, [], 2);
deepest_depth(deepest_depth == inf) = NaN;

%% find deepest ocean depths for entire CMS grid, and then use it to 
% interpolate depths for release points

% Find the deepest ocean depth for each point
deepest_depth_u_CMS = zeros(size(lon_grid, 1), size(lon_grid, 2), 1);
deepest_depth_v_CMS = zeros(size(lon_grid, 1), size(lon_grid, 2), 1);
deepest_depth_rho_CMS = zeros(size(lon_grid, 1), size(lon_grid, 2), 1);

for i = 1:size(lon_grid, 1) %row 364
    for j = 1:size(lon_grid, 2) %column 543
        
        deepest_u_idx = find(squeeze(mask_u_seafloor(i,j,:)) == 1, 1, 'last');
        deepest_v_idx = find(squeeze(mask_v_seafloor(i,j,:)) == 1, 1, 'last');
        deepest_rho_idx = find(squeeze(mask_rho_seafloor(i,j,:)) == 1, 1, 'last');

        if ~isempty(deepest_u_idx)
            deepest_depth_u_CMS(i,j) = zlevels(deepest_u_idx);
        else
            deepest_depth_u_CMS(i,j) = NaN;
        end
        
        if ~isempty(deepest_v_idx)
            deepest_depth_v_CMS(i,j) = zlevels(deepest_v_idx);
        else
            deepest_depth_v_CMS(i,j) = NaN;
        end
        
        if ~isempty(deepest_rho_idx)
            deepest_depth_rho_CMS(i,j) = zlevels(deepest_rho_idx);
        else
            deepest_depth_rho_CMS(i,j) = NaN;
        end
    end
end

% Calculate the minimum of the three depth matrices at each grid point
deepest_depth_CMS = min(cat(3, deepest_depth_u_CMS, deepest_depth_v_CMS, deepest_depth_rho_CMS), [], 3);

% %% testing - difference in seafloor masks
% 
% figure();imagescn(lon_grid',lat_grid',deepest_depth_u_CMS');axis equal;hold on;plot(relpoints(:,2) + 360,relpoints(:,3), '*r');xlim([360-65.5, 360-65]);ylim([18, 18.4])
% figure();imagescn(lon_grid',lat_grid',deepest_depth_v_CMS');axis equal;hold on;plot(relpoints(:,2) + 360,relpoints(:,3), '*r');xlim([360-65.5, 360-65]);ylim([18, 18.4])
% figure();imagescn(lon_grid',lat_grid',deepest_depth_rho_CMS');axis equal;hold on;plot(relpoints(:,2) + 360,relpoints(:,3), '*r');xlim([360-65.5, 360-65]);ylim([18, 18.4])
% 
% %% testing - difference in seafloor masks

% Interpolate the depths at the release points using the deepest_depth_CMS matrix
release_depths = interp2(lat_grid, lon_grid, deepest_depth_CMS, latitudes_release, longitudes_release);

% Display the results
disp('Interpolated depths at release points:');
disp(release_depths);

% plot(release_depths)

% STOPPING POINT
%   11 June 2024
%    - still need to address NaNs in 'release_depths' (nnz(isnan(release_depths))) and the few cases
%    where it actually produces a point *deeper* rather than shallower than
%    the shotgun approach above (which produced 'deepest_depth')
%    - Then, address jittering the points within their release polygon (and
%    at that stage, check any issues with CRS misalignment b/t MATLAB &
%    QGIS). Any points that can't reach ocean without exiting their polygon
%    should be thrown out
%    - And voila, time to officially export these depths to the release
%    file and move on with life!

% Plotting the deepest depths and release points for visualization
figure();
imagescn(longitudes_CMS', latitudes_CMS', -deepest_depth_CMS');
axis equal;
hold on;
plot(longitudes_release, latitudes_release, '*r');
xlim([min(longitudes_CMS), max(longitudes_CMS)]);
ylim([min(latitudes_CMS), max(latitudes_CMS)]);
clim([-50 0]);
colorbar;
title('Deepest Depths and Release Points');


%%

% % Ensure that both release_depths and deepest_depth have the same length
% if length(release_depths) ~= length(deepest_depth)
%     error('release_depths and deepest_depth must have the same length');
% end
% 
% % Find indices where release_depths have NaNs but deepest_depth does not
% nan_in_release_not_deepest = find(isnan(release_depths) & ~isnan(deepest_depth));
% 
% % Find indices where deepest_depth has NaNs but release_depths does not
% nan_in_deepest_not_release = find(~isnan(release_depths) & isnan(deepest_depth));
% 
% % Display the indices
% disp('Indices where release_depths have NaNs but deepest_depth does not:');
% disp(nan_in_release_not_deepest);
% 
% disp('Indices where deepest_depth has NaNs but release_depths does not:');
% disp(nan_in_deepest_not_release);
% 
% % Optionally, display the values at those indices for verification
% disp('Values at those indices for release_depths having NaNs but deepest_depth does not:');
% disp(table(nan_in_release_not_deepest, deepest_depth(nan_in_release_not_deepest)));
% 
% disp('Values at those indices for deepest_depth having NaNs but release_depths does not:');
% disp(table(nan_in_deepest_not_release, release_depths(nan_in_deepest_not_release)));
% 
% % Insert values from deepest_depth into release_depths where nan_in_release_not_deepest is true
% release_depths(nan_in_release_not_deepest) = deepest_depth(nan_in_release_not_deepest);
% 
% % Verify the insertion
% disp('Updated release_depths:');
% disp(release_depths);
% 
% find(isnan(release_depths))
% find(isnan(deepest_depth))



%%

% % Ensure necessary variables are available
% if ~exist('release_depths', 'var') || ~exist('deepest_depth_CMS', 'var') || ~exist('latitudes_release', 'var') || ~exist('longitudes_release', 'var')
%     error('Required variables are missing. Ensure release_depths, deepest_depth_CMS, latitudes_release, and longitudes_release are defined.');
% end
% 
% % Plot the depth contours
% figure;
% imagescn(longitudes_CMS', latitudes_CMS', -deepest_depth_CMS');
% axis equal;
% hold on;
% 
% % Plot the release points with depths and add black borders for better visibility
% scatter(longitudes_release, latitudes_release, 50, -release_depths, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% colorbar;
% clim([-250 0]); % Adjust color limits according to the depth range
% xlabel('Longitude');
% ylabel('Latitude');
% title('Release Depths over Deepest Depth Contours');
% hold off;
% 
% % Optionally, enhance the plot
% % Add grid, labels, etc.
% grid on;
% set(gca, 'FontSize', 12);

% Ensure necessary variables are available
if ~exist('release_depths', 'var') || ~exist('deepest_depth_CMS', 'var') || ~exist('latitudes_release', 'var') || ~exist('longitudes_release', 'var')
    error('Required variables are missing. Ensure release_depths, deepest_depth_CMS, latitudes_release, and longitudes_release are defined.');
end

% Ensure that both release_depths and deepest_depth have the same length
if length(release_depths) ~= length(deepest_depth)
    error('release_depths and deepest_depth must have the same length');
end

% Find indices where release_depths have NaNs but deepest_depth does not
nan_in_release_not_deepest = find(isnan(release_depths) & ~isnan(deepest_depth));

% Insert values from deepest_depth into release_depths where nan_in_release_not_deepest is true
release_depths(nan_in_release_not_deepest) = deepest_depth(nan_in_release_not_deepest);

% Find remaining NaN values in release_depths
remaining_nans = isnan(release_depths);

% Plot the depth contours
figure;
imagescn(longitudes_CMS', latitudes_CMS', -deepest_depth_CMS');
axis equal;
hold on;

% Plot the release points with depths
scatter(longitudes_release, latitudes_release, 50, -release_depths, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% Highlight the points where values were inserted from deepest_depth
scatter(longitudes_release(nan_in_release_not_deepest), latitudes_release(nan_in_release_not_deepest), 100, -release_depths(nan_in_release_not_deepest), 'filled', 'MarkerEdgeColor', 'c', 'LineWidth', 1.5);

% Highlight the points with remaining NaNs with red circles
scatter(longitudes_release(remaining_nans), latitudes_release(remaining_nans), 50, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% Add a color bar and adjust color limits according to the depth range
colorbar;
clim([-250 0]); % Adjust color limits according to the depth range
xlabel('Longitude');
ylabel('Latitude');
title('Release Depths over Deepest Depth Contours');

% Optionally, enhance the plot
% Add grid, labels, etc.
grid on;
set(gca, 'FontSize', 12);

% Add a legend to explain the colors
legend({'Original Release Depths', 'Inserted from Deepest Depth', 'Remaining NaNs'}, 'Location', 'best');

hold off;



%%

% Ensure necessary variables are available
if ~exist('deepest_depth', 'var') || ~exist('deepest_depth_CMS', 'var') || ~exist('latitudes_release', 'var') || ~exist('longitudes_release', 'var')
    error('Required variables are missing. Ensure deepest_depth, deepest_depth_CMS, latitudes_release, and longitudes_release are defined.');
end

% Plot the depth contours
figure;
imagescn(longitudes_CMS', latitudes_CMS', -deepest_depth_CMS');
axis equal;
hold on;

% Plot the release points with deepest depths
scatter(longitudes_release, latitudes_release, 50, -deepest_depth, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% Enhance the plot
colorbar;
clim([-250 0]); % Adjust color limits according to the depth range
xlabel('Longitude');
ylabel('Latitude');
title('Deepest Depths over Deepest Depth Contours');
grid on;
set(gca, 'FontSize', 12);
hold off;

%%

% Ensure necessary variables are available
if ~exist('release_depths', 'var') || ~exist('deepest_depth', 'var') || ~exist('deepest_depth_CMS', 'var') || ~exist('latitudes_release', 'var') || ~exist('longitudes_release', 'var')
    error('Required variables are missing. Ensure release_depths, deepest_depth, deepest_depth_CMS, latitudes_release, and longitudes_release are defined.');
end

% Set the focus area limits
xlim_focus = [360-64.85, 360-64.65];
ylim_focus = [18.25, 18.45];
depthlimit = -60;
zlim_focus = [depthlimit 0];

% xlim_focus = [360-65.3, 360-64];
% ylim_focus = [17.2, 19.1];
% depthlimit = -60;
% zlim_focus = [depthlimit 0];

% % Create a grid of all longitude and latitude pairs
% [lon_grid, lat_grid] = meshgrid(longitudes_CMS, latitudes_CMS);

% Create a mask for the focus area
focus_mask = lon_grid >= xlim_focus(1) & lon_grid <= xlim_focus(2) & lat_grid >= ylim_focus(1) & lat_grid <= ylim_focus(2);

% Extract the relevant grid and depth data for the focus area
lon_grid_focus = lon_grid(focus_mask);
lat_grid_focus = lat_grid(focus_mask);
deepest_depth_CMS_focus = deepest_depth_CMS(focus_mask);

% Create 3D plot
figure;

% Plot the depth contours as a surface
scatter3(lon_grid_focus, lat_grid_focus, -deepest_depth_CMS_focus, 50, -deepest_depth_CMS_focus, 'filled');
clim([depthlimit 0]);colorbar
hold on;

% Plot the release depths as 3D scatter points
scatter3(longitudes_release, latitudes_release, -release_depths, 10, -release_depths, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5); 
% xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45])

% Plot the deepest depths as 3D scatter points
scatter3(longitudes_release, latitudes_release, -deepest_depth, 10, -deepest_depth, 'filled', 'MarkerEdgeColor', 'r', 'LineWidth', 1.5);

% Set the view and limits
view(3);
colorbar;
xlim(xlim_focus);
ylim(ylim_focus);
zlim(zlim_focus);
xlabel('Longitude');
ylabel('Latitude');
zlabel('Depth (negative values)');
title('3D Plot of Depth Contours and Release Points');
grid on;
set(gca, 'FontSize', 12);

hold off;