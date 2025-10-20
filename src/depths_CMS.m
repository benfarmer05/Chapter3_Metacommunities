%% Script to assign depths to lat/lon reef points, that are sensible for
%   the CMS. may be adopted to jitter points off of land (if there is room
%   in their containing polygon square) or otherwise that will become part
%   of a separate script. once no points are remaining on land, then it is
%   a matter of locating their deepest nearest u/v/rho-grid depth that is
%   not below the seafloor. the final step is deciding the best way to
%   settle on the most useful depth - the CMS likely interpolates the
%   "seafloor" using its available grid system, so I probably want to
%   use some kind of interpolation scheme. still testing that out
% 19 Oct 2025

%% Conservative depth assignment checking all three masks

clear; clc;

% Setup paths
projectPath = matlab.project.rootProject().RootFolder;
tempPath = fullfile(projectPath, 'temp');
dataPath = fullfile(projectPath, 'data');
outputPath = fullfile(projectPath, 'output');

% Load one hydro file
files = dir(fullfile(tempPath, 'drive*.nc'));
hydroFile = fullfile(tempPath, files(1).name);

% Read coordinates
lon_grid = ncread(hydroFile, 'Longitude');
lat_grid = ncread(hydroFile, 'Latitude');
[lat_mesh, lon_mesh] = meshgrid(lat_grid, lon_grid);

% Read ALL THREE masks to find the shallowest valid depth
fill_value = single(1.2676506e+30);
zu = ncread(hydroFile, 'zu');
zv = ncread(hydroFile, 'zv');
zw = ncread(hydroFile, 'zw');

mask_zu = (zu ~= fill_value);
mask_zv = (zv ~= fill_value);
mask_zw = (zw ~= fill_value);

% Define depth levels
zlevels = [0.26 1 2 4 6 8 10 12 14 16 18 20 22 24 26 29 32 35 38 42 46 50 60 70 85 100 120 140 160 190 220 250]';

% Load release points
% relpoints = readmatrix(fullfile(dataPath, 'points_650_none-on-land.csv'));
relpoints = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
lon_release = relpoints(:,2) + 360;
lat_release = relpoints(:,3);

% Find nearest grid point
idx_nearest = knnsearch([lon_mesh(:), lat_mesh(:)], [lon_release, lat_release]);
[row_idx, col_idx] = ind2sub(size(lon_mesh), idx_nearest);

% Find the absolute deepest depth for each mask at each release point
deepest_idx_zu = zeros(length(idx_nearest), 1);
deepest_idx_zv = zeros(length(idx_nearest), 1);
deepest_idx_zw = zeros(length(idx_nearest), 1);

for i = 1:length(idx_nearest)
    % zu mask
    mask_profile_zu = squeeze(mask_zu(row_idx(i), col_idx(i), :));
    idx_zu = find(mask_profile_zu == 1, 1, 'last');
    if ~isempty(idx_zu)
        deepest_idx_zu(i) = idx_zu;
    else
        deepest_idx_zu(i) = 0;
    end
    
    % zv mask
    mask_profile_zv = squeeze(mask_zv(row_idx(i), col_idx(i), :));
    idx_zv = find(mask_profile_zv == 1, 1, 'last');
    if ~isempty(idx_zv)
        deepest_idx_zv(i) = idx_zv;
    else
        deepest_idx_zv(i) = 0;
    end
    
    % zw mask
    mask_profile_zw = squeeze(mask_zw(row_idx(i), col_idx(i), :));
    idx_zw = find(mask_profile_zw == 1, 1, 'last');
    if ~isempty(idx_zw)
        deepest_idx_zw(i) = idx_zw;
    else
        deepest_idx_zw(i) = 0;
    end
end

% Find the SHALLOWEST of the three deepest indices (most conservative)
% Handle empty results (0) by treating them as inf so they don't affect min
deepest_idx_zu(deepest_idx_zu == 0) = inf;
deepest_idx_zv(deepest_idx_zv == 0) = inf;
deepest_idx_zw(deepest_idx_zw == 0) = inf;

shallowest_deepest_idx = min([deepest_idx_zu, deepest_idx_zv, deepest_idx_zw], [], 2);

% Apply safety margin and depth cutoff
deepest_depth = zeros(length(idx_nearest), 1);
deepest_depth_absolute = zeros(length(idx_nearest), 1);
safety_margin = 1; % Go 1 level shallower than absolute deepest

for i = 1:length(idx_nearest)
    if ~isinf(shallowest_deepest_idx(i))
        % Store absolute depth (no safety margin)
        deepest_depth_absolute(i) = zlevels(shallowest_deepest_idx(i));
        
        % Apply safety margin
        safe_idx = max(1, shallowest_deepest_idx(i) - safety_margin);
        depth_value = zlevels(safe_idx);
        
        % Set to NaN if ABSOLUTE depth is deeper than 60m
        if deepest_depth_absolute(i) > 60
            deepest_depth(i) = NaN;
        else
            deepest_depth(i) = depth_value;
        end
    else
        deepest_depth(i) = NaN;
        deepest_depth_absolute(i) = NaN;
    end
end

% Summary statistics
fprintf('\n========== DEPTH ASSIGNMENT SUMMARY ==========\n\n');
fprintf('Points with valid depths (<=60m): %d / %d\n', sum(~isnan(deepest_depth)), length(deepest_depth));
fprintf('Points excluded (>60m): %d\n', sum(isnan(deepest_depth)));
fprintf('Depth range: %.2f to %.2f m\n', min(deepest_depth), max(deepest_depth));

% Check how often each mask was the shallowest
valid_points = ~isinf(shallowest_deepest_idx);

zu_matches = (shallowest_deepest_idx == deepest_idx_zu) & valid_points;
zv_matches = (shallowest_deepest_idx == deepest_idx_zv) & valid_points;
zw_matches = (shallowest_deepest_idx == deepest_idx_zw) & valid_points;

% Count unique situations
all_three_same = zu_matches & zv_matches & zw_matches;
zu_zv_same = zu_matches & zv_matches & ~zw_matches;
zu_zw_same = zu_matches & zw_matches & ~zv_matches;
zv_zw_same = zv_matches & zw_matches & ~zu_matches;
only_zu = zu_matches & ~zv_matches & ~zw_matches;
only_zv = zv_matches & ~zu_matches & ~zw_matches;
only_zw = zw_matches & ~zu_matches & ~zv_matches;

fprintf('\nWhich mask(s) determined the final depth?\n');
fprintf('  All three identical:     %d points\n', sum(all_three_same));
fprintf('  Only zu was shallowest:  %d points\n', sum(only_zu));
fprintf('  Only zv was shallowest:  %d points\n', sum(only_zv));
fprintf('  Only zw was shallowest:  %d points\n', sum(only_zw));
fprintf('  zu & zv tied (not zw):   %d points\n', sum(zu_zv_same));
fprintf('  zu & zw tied (not zv):   %d points\n', sum(zu_zw_same));
fprintf('  zv & zw tied (not zu):   %d points\n', sum(zv_zw_same));

% Histogram
figure;
histogram(deepest_depth, 'BinWidth', 5);
xlabel('Depth (m)'); ylabel('Count');
title('Distribution of conservative release depths (<=60m only)');

% 3D visualization of release points (oceanographic convention: depth positive downward)
figure;

% Only plot valid (non-NaN) points
valid_idx = ~isnan(deepest_depth);

scatter3(lon_release(valid_idx), lat_release(valid_idx), deepest_depth(valid_idx), ...
    10, deepest_depth(valid_idx), 'filled');

colormap(parula);
c = colorbar;
c.Label.String = 'Depth below surface (m)';
xlabel('Longitude (°E)');
ylabel('Latitude (°N)');
zlabel('Depth below surface (m)');
title('3D Release Point Space (<=60m only, shallowest of zu/zv/zw)');
grid on;
view(3);
set(gca, 'ZDir', 'reverse'); % Depth increases downward
zlim([0, 65]);
hold off;

% Calculate the depth difference (only for valid points)
depth_difference = deepest_depth_absolute - deepest_depth;

% 3D visualization showing the gap between conservative and absolute depths
figure;
hold on;

% Only plot valid points
valid_for_comparison = ~isnan(deepest_depth) & ~isnan(deepest_depth_absolute);

% Draw lines first (without adding to legend)
first_line = true;
for i = 1:length(lon_release)
    if valid_for_comparison(i)
        h_line = plot3([lon_release(i), lon_release(i)], [lat_release(i), lat_release(i)], ...
              [deepest_depth(i), deepest_depth_absolute(i)], 'k-', 'LineWidth', 0.5);
        if ~first_line
            set(get(get(h_line,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
        first_line = false;
    end
end

% Plot conservative release depths (your assigned depths)
scatter3(lon_release(valid_for_comparison), lat_release(valid_for_comparison), ...
    deepest_depth(valid_for_comparison), 15, 'b', 'filled', 'DisplayName', 'Conservative depth');

% Plot absolute deepest depths
scatter3(lon_release(valid_for_comparison), lat_release(valid_for_comparison), ...
    deepest_depth_absolute(valid_for_comparison), 15, 'r', 'filled', 'DisplayName', 'Absolute deepest');

xlabel('Longitude (°E)');
ylabel('Latitude (°N)');
zlabel('Depth below surface (m)');
title(sprintf('Conservative vs Absolute Depths (mean difference: %.2f m)', mean(depth_difference, 'omitnan')));
legend('Location', 'best');
grid on;
view(3);
set(gca, 'ZDir', 'reverse');
zlim([0, 65]);
hold off;

fprintf('\nMean depth difference (safety margin): %.2f m\n', mean(depth_difference, 'omitnan'));
fprintf('Max depth difference: %.2f m\n', max(depth_difference(valid_for_comparison)));

fprintf('\n==============================================\n\n');

%% Verify consistency across all drive_*.nc files

fprintf('\n========== CONSISTENCY CHECK ACROSS FILES ==========\n\n');

% Get all drive_*.nc files
files = dir(fullfile(tempPath, 'drive*.nc'));

if length(files) < 2
    fprintf('Only %d file found. Need at least 2 files to compare.\n', length(files));
else
    fprintf('Found %d files to compare:\n', length(files));
    for f = 1:length(files)
        fprintf('  %d. %s\n', f, files(f).name);
    end
    fprintf('\n');
    
    % Store masks from all files
    all_masks_zu = cell(length(files), 1);
    all_masks_zv = cell(length(files), 1);
    all_masks_zw = cell(length(files), 1);
    
    % Read masks from each file
    for f = 1:length(files)
        hydroFile_test = fullfile(tempPath, files(f).name);
        
        zu_test = ncread(hydroFile_test, 'zu');
        zv_test = ncread(hydroFile_test, 'zv');
        zw_test = ncread(hydroFile_test, 'zw');
        
        all_masks_zu{f} = (zu_test ~= fill_value);
        all_masks_zv{f} = (zv_test ~= fill_value);
        all_masks_zw{f} = (zw_test ~= fill_value);
    end
    
    % Compare all files to the first file
    all_identical = true;
    
    for f = 2:length(files)
        zu_match = isequal(all_masks_zu{1}, all_masks_zu{f});
        zv_match = isequal(all_masks_zv{1}, all_masks_zv{f});
        zw_match = isequal(all_masks_zw{1}, all_masks_zw{f});
        
        fprintf('Comparing %s to %s:\n', files(1).name, files(f).name);
        fprintf('  zu masks identical: %s\n', mat2str(zu_match));
        fprintf('  zv masks identical: %s\n', mat2str(zv_match));
        fprintf('  zw masks identical: %s\n', mat2str(zw_match));
        fprintf('\n');
        
        if ~(zu_match && zv_match && zw_match)
            all_identical = false;
        end
    end
    
    if all_identical
        fprintf('✓ ALL FILES HAVE IDENTICAL MASKS - Safe to use any file!\n');
    else
        fprintf('⚠ WARNING: Files have different masks - investigate further!\n');
    end
end

fprintf('\n====================================================\n\n');

%% Export depth assignments for use in release file creation

fprintf('\n========== EXPORTING DEPTH ASSIGNMENTS ==========\n\n');

% Extract IDs from the release points file (assuming first column is ID)
IDs_release = relpoints(:,1);

% Create output structure with ID, Lon, Lat, and assigned depth
% Longitude is stored WITHOUT the +360 offset to match original data
depth_assignments = [IDs_release, relpoints(:,2), relpoints(:,3), deepest_depth];

% Create output filename with timestamp
currentDateTime = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
currentDateTimeStr = string(currentDateTime);
depthFileName = sprintf('depth_assignments_%s.csv', currentDateTimeStr);
depthFilePath = fullfile(outputPath, depthFileName);

% Write to CSV file with header
fprintf('Writing depth assignments to: %s\n', depthFileName);
fid = fopen(depthFilePath, 'w');
fprintf(fid, 'ID,Longitude,Latitude,Depth\n');
fclose(fid);
writematrix(depth_assignments, depthFilePath, 'WriteMode', 'append');

% Also save as .mat file for faster loading in MATLAB
matFileName = sprintf('depth_assignments_%s.mat', currentDateTimeStr);
matFilePath = fullfile(outputPath, matFileName);
save(matFilePath, 'depth_assignments', 'deepest_depth', 'deepest_depth_absolute', ...
     'lon_release', 'lat_release', 'IDs_release', 'zlevels', 'safety_margin');

fprintf('Depth assignments saved to:\n');
fprintf('  CSV: %s\n', depthFileName);
fprintf('  MAT: %s\n', matFileName);
fprintf('\nSummary of exported data:\n');
fprintf('  Total points: %d\n', length(deepest_depth));
fprintf('  Valid depths (<=60m): %d\n', sum(~isnan(deepest_depth)));
fprintf('  Excluded (>60m or no ocean): %d\n', sum(isnan(deepest_depth)));
fprintf('  Depth range: %.2f to %.2f m\n', min(deepest_depth), max(deepest_depth));

fprintf('\n====================================================\n\n');


% %below were earlier attempts with interpolation - decided to do a cleaner
% % "move one up" approach where the next highest depth level from deepest
% % common depth across zu, zv, and zw is chosen for CMS release.
% 
% clear;clc
% 
% %% setup
% 
% % Get the project root directory
% projectPath = matlab.project.rootProject().RootFolder;
% 
% % Define paths relative to the project root
% dataPath = fullfile(projectPath, 'data');
% srcPath = fullfile(projectPath, 'src');
% outputPath = fullfile(projectPath, 'output');
% tempPath = fullfile(projectPath, 'temp');
% 
% %%
% 
% % % NOTE - may need to return to this
% % CMSfilesDrive = '/Users/benja/testCMS/tests_for_depth-assignment'; %'/Volumes/UVI_Hydro_2019-2020/CMS_inputs_nozeta'; %
% 
% % NOTE - may need to adjust first value [currently 0.26] to reflect what
% % was done with roms2cms_A_grid.m
% zlevels = [0.26 1 2 4 6 8 10 12 14 16 18 20 22 24 26 29 32 35 38 42 46 50 60 70 85 100 120 140 160 190 220 250]'; %32 layers down to 250 m depth, mid-depth hi-res
% 
% % NOTE / STOPPING POINT - for the below, was previously just pulling a
% % random hydro file I'd created for the CMS. just need to adopt this
% % approach with proper pathing, and also double-check how I created my
% % hydro files
% 
% %randomly select a CMS-ready ROMS file from the HPC
% % files = dir(fullfile(tempPath, 'nest_1_2019*'));
% files = dir(fullfile(tempPath, 'HPC*'));
% randomIndex = randi(length(files)); % generate a random index between 1 and the number of files
% baseFileName_CMS_HPC = files(randomIndex).name;
% ncdisp(fullfile(tempPath, baseFileName_CMS_HPC))  % Pass the full path to ncdisp
% 
% %randomly select a CMS-ready ROMS file from my hard drive
% % files = dir(fullfile(tempPath, 'nest_1_2019*'));
% files = dir(fullfile(tempPath, 'drive*'));
% randomIndex = randi(length(files)); % generate a random index between 1 and the number of files
% baseFileName_CMS_drive = files(randomIndex).name;
% ncdisp(fullfile(tempPath, baseFileName_CMS_drive))  % Pass the full path to ncdisp
% 
% %randomly select a CMS-ready HYCOM file
% files = dir(fullfile(tempPath, 'nest_1_2023*'));
% randomIndex = randi(length(files)); % generate a random index between 1 and the number of files
% baseFileName_HYCOM = files(randomIndex).name;
% ncdisp(fullfile(tempPath, baseFileName_HYCOM))  % Pass the full path to ncdisp
% 
% %CMS coordinate systems
% longitudes_CMS = ncread(fullfile(tempPath, baseFileName_HYCOM), 'Longitude');
% latitudes_CMS = ncread(fullfile(tempPath, baseFileName_HYCOM), 'Latitude');
% 
% % Create a grid of all longitude and latitude pairs
% [lat_grid, lon_grid] = meshgrid(latitudes_CMS, longitudes_CMS);
% 
% % Flatten the grid matrices
% lon_grid_indices = lon_grid(:);
% lat_grid_indices = lat_grid(:);
% 
% %CMS seafloor masks
% fill_value = single(1.2676506e+30);
% mask_u_seafloor = ncread(fullfile(tempPath, baseFileName_CMS_HPC), 'zu');
% mask_u_seafloor(mask_u_seafloor==fill_value) = 0; %seafloor is 0's
% mask_u_seafloor(mask_u_seafloor ~= 0) = 1; %the ocean is 1's
% mask_v_seafloor = ncread(fullfile(tempPath, baseFileName_CMS_HPC), 'zv'); 
% mask_v_seafloor(mask_v_seafloor==fill_value) = 0;
% mask_v_seafloor(mask_v_seafloor ~= 0) = 1;
% mask_rho_seafloor = ncread(fullfile(tempPath, baseFileName_CMS_HPC), 'zw'); 
% mask_rho_seafloor(mask_rho_seafloor==fill_value) = 0; 
% mask_rho_seafloor(mask_rho_seafloor ~= 0) = 1;
% 
% %extract release point coordinates from GIS output and ensure they are
% % sorted by their unique ID
% relpoints = readmatrix(fullfile(dataPath, 'points_650_none-on-land.csv'));
% relpoints = relpoints(:, 10:12);
% relpoints = sortrows(relpoints, 1);
% 
% IDs_release = relpoints(:,1);
% longitudes_release = relpoints(:,2) + 360;
% latitudes_release = relpoints(:,3);
% 
% % %% testing - difference in seafloor masks
% % 
% % depthslice = 15;
% % 
% % mask_u_seafloor = mask_u_seafloor(:,:,depthslice);
% % figure();imagescn(lon_grid',lat_grid',mask_u_seafloor');axis equal;hold on;plot(relpoints(:,2) + 360,relpoints(:,3), '*r');xlim([360-65.5, 360-65]);ylim([18, 18.4])
% % 
% % mask_v_seafloor = mask_v_seafloor(:,:,depthslice);
% % figure();imagescn(lon_grid',lat_grid',mask_v_seafloor');axis equal;hold on;plot(relpoints(:,2) + 360,relpoints(:,3), '*r');xlim([360-65.5, 360-65]);ylim([18, 18.4])
% % 
% % mask_rho_seafloor = mask_rho_seafloor(:,:,depthslice);
% % figure();imagescn(lon_grid',lat_grid',mask_rho_seafloor');axis equal;hold on;plot(relpoints(:,2) + 360,relpoints(:,3), '*r');xlim([360-65.5, 360-65]);ylim([18, 18.4])
% % 
% % %% testing - difference in seafloor masks
% 
% % Locate CMS grid coordinates nearest to PREDICT points
% idx_CMS_release = knnsearch([lon_grid_indices(:), lat_grid_indices(:)], [longitudes_release, latitudes_release]);
% 
% % Initialize the seamasks for each depth level
% seamask_u_release = zeros(length(idx_CMS_release), 32);
% seamask_v_release = zeros(length(idx_CMS_release), 32);
% seamask_rho_release = zeros(length(idx_CMS_release), 32);
% 
% % Iterate over each depth level
% for depth = 1:32
%     % Extract the mask for the current depth level
%     mask_u_depth = mask_u_seafloor(:,:,depth);
%     mask_v_depth = mask_v_seafloor(:,:,depth);
%     mask_rho_depth = mask_rho_seafloor(:,:,depth);
% 
%     % Apply the mask for the current depth level
%     seamask_u_release(:, depth) = mask_u_depth(idx_CMS_release);
%     seamask_v_release(:, depth) = mask_v_depth(idx_CMS_release);
%     seamask_rho_release(:, depth) = mask_rho_depth(idx_CMS_release);
% end
% 
% % Find the deepest ocean depth for each point
% deepest_depth_u = zeros(length(idx_CMS_release), 1);
% deepest_depth_v = zeros(length(idx_CMS_release), 1);
% deepest_depth_rho = zeros(length(idx_CMS_release), 1);
% 
% for i = 1:length(idx_CMS_release) %test: i = 2495 to find a NaN 'u'
%     deepest_u_idx = find(seamask_u_release(i, :) == 1, 1, 'last');
%     deepest_v_idx = find(seamask_v_release(i, :) == 1, 1, 'last');
%     deepest_rho_idx = find(seamask_rho_release(i, :) == 1, 1, 'last');
% 
%     if ~isempty(deepest_u_idx)
%         deepest_depth_u(i) = zlevels(deepest_u_idx);
%     else
%         deepest_depth_u(i) = NaN;
%     end
% 
%     if ~isempty(deepest_v_idx)
%         deepest_depth_v(i) = zlevels(deepest_v_idx);
%     else
%         deepest_depth_v(i) = NaN;
%     end
% 
%     if ~isempty(deepest_rho_idx)
%         deepest_depth_rho(i) = zlevels(deepest_rho_idx);
%     else
%         deepest_depth_rho(i) = NaN;
%     end
% end
% 
% %find the shallowest nearest-neighbor depth across each grid
% combined_depths = [deepest_depth_u, deepest_depth_v, deepest_depth_rho];
% combined_depths(isnan(combined_depths)) = inf;
% deepest_depth = min(combined_depths, [], 2);
% deepest_depth(deepest_depth == inf) = NaN;
% 
% %% find deepest ocean depths for entire CMS grid, and then use it to 
% % interpolate depths for release points
% 
% % Find the deepest ocean depth for each point
% deepest_depth_u_CMS = zeros(size(lon_grid, 1), size(lon_grid, 2), 1);
% deepest_depth_v_CMS = zeros(size(lon_grid, 1), size(lon_grid, 2), 1);
% deepest_depth_rho_CMS = zeros(size(lon_grid, 1), size(lon_grid, 2), 1);
% 
% for i = 1:size(lon_grid, 1) %row 364
%     for j = 1:size(lon_grid, 2) %column 543
% 
%         deepest_u_idx = find(squeeze(mask_u_seafloor(i,j,:)) == 1, 1, 'last');
%         deepest_v_idx = find(squeeze(mask_v_seafloor(i,j,:)) == 1, 1, 'last');
%         deepest_rho_idx = find(squeeze(mask_rho_seafloor(i,j,:)) == 1, 1, 'last');
% 
%         if ~isempty(deepest_u_idx)
%             deepest_depth_u_CMS(i,j) = zlevels(deepest_u_idx);
%         else
%             deepest_depth_u_CMS(i,j) = NaN;
%         end
% 
%         if ~isempty(deepest_v_idx)
%             deepest_depth_v_CMS(i,j) = zlevels(deepest_v_idx);
%         else
%             deepest_depth_v_CMS(i,j) = NaN;
%         end
% 
%         if ~isempty(deepest_rho_idx)
%             deepest_depth_rho_CMS(i,j) = zlevels(deepest_rho_idx);
%         else
%             deepest_depth_rho_CMS(i,j) = NaN;
%         end
%     end
% end
% 
% % Calculate the minimum of the three depth matrices at each grid point
% deepest_depth_CMS = min(cat(3, deepest_depth_u_CMS, deepest_depth_v_CMS, deepest_depth_rho_CMS), [], 3);
% 
% % %% testing - difference in seafloor masks
% % 
% % figure();imagescn(lon_grid',lat_grid',deepest_depth_u_CMS');axis equal;hold on;plot(relpoints(:,2) + 360,relpoints(:,3), '*r');xlim([360-65.5, 360-65]);ylim([18, 18.4])
% % figure();imagescn(lon_grid',lat_grid',deepest_depth_v_CMS');axis equal;hold on;plot(relpoints(:,2) + 360,relpoints(:,3), '*r');xlim([360-65.5, 360-65]);ylim([18, 18.4])
% % figure();imagescn(lon_grid',lat_grid',deepest_depth_rho_CMS');axis equal;hold on;plot(relpoints(:,2) + 360,relpoints(:,3), '*r');xlim([360-65.5, 360-65]);ylim([18, 18.4])
% % 
% % %% testing - difference in seafloor masks
% 
% % Interpolate the depths at the release points using the deepest_depth_CMS matrix
% release_depths = interp2(lat_grid, lon_grid, deepest_depth_CMS, latitudes_release, longitudes_release);
% 
% % Display the results
% disp('Interpolated depths at release points:');
% disp(release_depths);
% 
% % plot(release_depths)
% 
% % STOPPING POINT
% %   11 June 2024
% %    - still need to address NaNs in 'release_depths' (nnz(isnan(release_depths))) and the few cases
% %    where it actually produces a point *deeper* rather than shallower than
% %    the shotgun approach above (which produced 'deepest_depth')
% %    - Then, address jittering the points within their release polygon (and
% %    at that stage, check any issues with CRS misalignment b/t MATLAB &
% %    QGIS). Any points that can't reach ocean without exiting their polygon
% %    should be thrown out
% %    - And voila, time to officially export these depths to the release
% %    file and move on with life!
% 
% % Plotting the deepest depths and release points for visualization
% figure();
% imagescn(longitudes_CMS', latitudes_CMS', -deepest_depth_CMS');
% axis equal;
% hold on;
% plot(longitudes_release, latitudes_release, '*r');
% xlim([min(longitudes_CMS), max(longitudes_CMS)]);
% ylim([min(latitudes_CMS), max(latitudes_CMS)]);
% clim([-50 0]);
% colorbar;
% title('Deepest Depths and Release Points');
% 
% 
% %%
% 
% % % Ensure that both release_depths and deepest_depth have the same length
% % if length(release_depths) ~= length(deepest_depth)
% %     error('release_depths and deepest_depth must have the same length');
% % end
% % 
% % % Find indices where release_depths have NaNs but deepest_depth does not
% % nan_in_release_not_deepest = find(isnan(release_depths) & ~isnan(deepest_depth));
% % 
% % % Find indices where deepest_depth has NaNs but release_depths does not
% % nan_in_deepest_not_release = find(~isnan(release_depths) & isnan(deepest_depth));
% % 
% % % Display the indices
% % disp('Indices where release_depths have NaNs but deepest_depth does not:');
% % disp(nan_in_release_not_deepest);
% % 
% % disp('Indices where deepest_depth has NaNs but release_depths does not:');
% % disp(nan_in_deepest_not_release);
% % 
% % % Optionally, display the values at those indices for verification
% % disp('Values at those indices for release_depths having NaNs but deepest_depth does not:');
% % disp(table(nan_in_release_not_deepest, deepest_depth(nan_in_release_not_deepest)));
% % 
% % disp('Values at those indices for deepest_depth having NaNs but release_depths does not:');
% % disp(table(nan_in_deepest_not_release, release_depths(nan_in_deepest_not_release)));
% % 
% % % Insert values from deepest_depth into release_depths where nan_in_release_not_deepest is true
% % release_depths(nan_in_release_not_deepest) = deepest_depth(nan_in_release_not_deepest);
% % 
% % % Verify the insertion
% % disp('Updated release_depths:');
% % disp(release_depths);
% % 
% % find(isnan(release_depths))
% % find(isnan(deepest_depth))
% 
% 
% 
% %%
% 
% % % Ensure necessary variables are available
% % if ~exist('release_depths', 'var') || ~exist('deepest_depth_CMS', 'var') || ~exist('latitudes_release', 'var') || ~exist('longitudes_release', 'var')
% %     error('Required variables are missing. Ensure release_depths, deepest_depth_CMS, latitudes_release, and longitudes_release are defined.');
% % end
% % 
% % % Plot the depth contours
% % figure;
% % imagescn(longitudes_CMS', latitudes_CMS', -deepest_depth_CMS');
% % axis equal;
% % hold on;
% % 
% % % Plot the release points with depths and add black borders for better visibility
% % scatter(longitudes_release, latitudes_release, 50, -release_depths, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% % colorbar;
% % clim([-250 0]); % Adjust color limits according to the depth range
% % xlabel('Longitude');
% % ylabel('Latitude');
% % title('Release Depths over Deepest Depth Contours');
% % hold off;
% % 
% % % Optionally, enhance the plot
% % % Add grid, labels, etc.
% % grid on;
% % set(gca, 'FontSize', 12);
% 
% % Ensure necessary variables are available
% if ~exist('release_depths', 'var') || ~exist('deepest_depth_CMS', 'var') || ~exist('latitudes_release', 'var') || ~exist('longitudes_release', 'var')
%     error('Required variables are missing. Ensure release_depths, deepest_depth_CMS, latitudes_release, and longitudes_release are defined.');
% end
% 
% % Ensure that both release_depths and deepest_depth have the same length
% if length(release_depths) ~= length(deepest_depth)
%     error('release_depths and deepest_depth must have the same length');
% end
% 
% % Find indices where release_depths have NaNs but deepest_depth does not
% nan_in_release_not_deepest = find(isnan(release_depths) & ~isnan(deepest_depth));
% 
% % Insert values from deepest_depth into release_depths where nan_in_release_not_deepest is true
% release_depths(nan_in_release_not_deepest) = deepest_depth(nan_in_release_not_deepest);
% 
% % Find remaining NaN values in release_depths
% remaining_nans = isnan(release_depths);
% 
% % Plot the depth contours
% figure;
% imagescn(longitudes_CMS', latitudes_CMS', -deepest_depth_CMS');
% axis equal;
% hold on;
% 
% % Plot the release points with depths
% scatter(longitudes_release, latitudes_release, 50, -release_depths, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% 
% % Highlight the points where values were inserted from deepest_depth
% scatter(longitudes_release(nan_in_release_not_deepest), latitudes_release(nan_in_release_not_deepest), 100, -release_depths(nan_in_release_not_deepest), 'filled', 'MarkerEdgeColor', 'c', 'LineWidth', 1.5);
% 
% % Highlight the points with remaining NaNs with red circles
% scatter(longitudes_release(remaining_nans), latitudes_release(remaining_nans), 50, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% 
% % Add a color bar and adjust color limits according to the depth range
% colorbar;
% clim([-250 0]); % Adjust color limits according to the depth range
% xlabel('Longitude');
% ylabel('Latitude');
% title('Release Depths over Deepest Depth Contours');
% 
% % Optionally, enhance the plot
% % Add grid, labels, etc.
% grid on;
% set(gca, 'FontSize', 12);
% 
% % Add a legend to explain the colors
% legend({'Original Release Depths', 'Inserted from Deepest Depth', 'Remaining NaNs'}, 'Location', 'best');
% 
% hold off;
% 
% 
% 
% %%
% 
% % Ensure necessary variables are available
% if ~exist('deepest_depth', 'var') || ~exist('deepest_depth_CMS', 'var') || ~exist('latitudes_release', 'var') || ~exist('longitudes_release', 'var')
%     error('Required variables are missing. Ensure deepest_depth, deepest_depth_CMS, latitudes_release, and longitudes_release are defined.');
% end
% 
% % Plot the depth contours
% figure;
% imagescn(longitudes_CMS', latitudes_CMS', -deepest_depth_CMS');
% axis equal;
% hold on;
% 
% % Plot the release points with deepest depths
% scatter(longitudes_release, latitudes_release, 50, -deepest_depth, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% 
% % Enhance the plot
% colorbar;
% clim([-250 0]); % Adjust color limits according to the depth range
% xlabel('Longitude');
% ylabel('Latitude');
% title('Deepest Depths over Deepest Depth Contours');
% grid on;
% set(gca, 'FontSize', 12);
% hold off;
% 
% %%
% 
% % Ensure necessary variables are available
% if ~exist('release_depths', 'var') || ~exist('deepest_depth', 'var') || ~exist('deepest_depth_CMS', 'var') || ~exist('latitudes_release', 'var') || ~exist('longitudes_release', 'var')
%     error('Required variables are missing. Ensure release_depths, deepest_depth, deepest_depth_CMS, latitudes_release, and longitudes_release are defined.');
% end
% 
% % Set the focus area limits
% xlim_focus = [360-64.85, 360-64.65];
% ylim_focus = [18.25, 18.45];
% depthlimit = -60;
% zlim_focus = [depthlimit 0];
% 
% % xlim_focus = [360-65.3, 360-64];
% % ylim_focus = [17.2, 19.1];
% % depthlimit = -60;
% % zlim_focus = [depthlimit 0];
% 
% % % Create a grid of all longitude and latitude pairs
% % [lon_grid, lat_grid] = meshgrid(longitudes_CMS, latitudes_CMS);
% 
% % Create a mask for the focus area
% focus_mask = lon_grid >= xlim_focus(1) & lon_grid <= xlim_focus(2) & lat_grid >= ylim_focus(1) & lat_grid <= ylim_focus(2);
% 
% % Extract the relevant grid and depth data for the focus area
% lon_grid_focus = lon_grid(focus_mask);
% lat_grid_focus = lat_grid(focus_mask);
% deepest_depth_CMS_focus = deepest_depth_CMS(focus_mask);
% 
% % Create 3D plot
% figure;
% 
% % Plot the depth contours as a surface
% scatter3(lon_grid_focus, lat_grid_focus, -deepest_depth_CMS_focus, 50, -deepest_depth_CMS_focus, 'filled');
% clim([depthlimit 0]);colorbar
% hold on;
% 
% % Plot the release depths as 3D scatter points
% scatter3(longitudes_release, latitudes_release, -release_depths, 10, -release_depths, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5); 
% % xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45])
% 
% % Plot the deepest depths as 3D scatter points
% scatter3(longitudes_release, latitudes_release, -deepest_depth, 10, -deepest_depth, 'filled', 'MarkerEdgeColor', 'r', 'LineWidth', 1.5);
% 
% % Set the view and limits
% view(3);
% colorbar;
% xlim(xlim_focus);
% ylim(ylim_focus);
% zlim(zlim_focus);
% xlabel('Longitude');
% ylabel('Latitude');
% zlabel('Depth (negative values)');
% title('3D Plot of Depth Contours and Release Points');
% grid on;
% set(gca, 'FontSize', 12);
% 
hold off;