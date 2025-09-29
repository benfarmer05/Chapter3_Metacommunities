clear;clc

%% Initialize paths
projectPath = matlab.project.rootProject().RootFolder;
tempPath = fullfile(projectPath, 'temp');
srcPath = fullfile(projectPath, 'src');
dataPath = fullfile(projectPath, 'data');

%% Load grid and setup
cd(srcPath) % For functions
addpath(dataPath) % For vigrid.nc

% Load ROMS grid
h = ncread(fullfile(dataPath, 'vigrid.nc'), 'h');
theta_s = ncread(fullfile(dataPath, 'vigrid.nc'), 'theta_s');
theta_b = ncread(fullfile(dataPath, 'vigrid.nc'), 'theta_b');
hc = ncread(fullfile(dataPath, 'vigrid.nc'), 'hc');

% Get ROMS grid coordinates
lon_rho = ncread(fullfile(dataPath, 'vigrid.nc'), 'lon_rho');
lat_rho = ncread(fullfile(dataPath, 'vigrid.nc'), 'lat_rho');
lon_u = ncread(fullfile(dataPath, 'vigrid.nc'), 'lon_u');
lat_u = ncread(fullfile(dataPath, 'vigrid.nc'), 'lat_u');

% Convert longitudes to positive values (if needed)
if any(lon_rho(:) < 0)
    lon_rho = 360 + lon_rho; % Convert from negative to positive longitude
    lon_u = 360 + lon_u;
end

% Find first CROCO file
crocolist = dir(fullfile(tempPath,'croco_his*.nc'));
if isempty(crocolist), error('No CROCO files found'); end
testFile = fullfile(tempPath, crocolist(1).name);

% Get test data
u_data = ncread(testFile, 'u');
layers = size(u_data,3);
zeta = zeros(size(h));
uvel = u_data(:,:,:,1); % First timestep only

% Calculate u-grid depths
u_zlev = set_depth(2,4,theta_s,theta_b,hc,layers,3,h,zeta,0);
u_zlev = u_zlev(:,:,layers:-1:1);

% Preprocess u-velocity
uvel = uvel(:,:,layers:-1:1);
uvel(uvel==0) = NaN;

% Define z-levels (identical to main script)
max_u_zlev = max(max(abs(squeeze(u_zlev(:,:,1)))));
zlevels = [max_u_zlev 1 2 4 6 8 10 12 14 16 18 20 22 24 26 29 32 35 38 42 46 50 60 70 85 100 120 140 160 190 220 250]';
numlevels = length(zlevels);

fprintf('Testing masking comparison on full domain\n');
fprintf('Using %d z-levels from %.2f to %.0f meters\n', numlevels, zlevels(1), zlevels(end));

%% Run comparison on full domain
[nrows, ncols, ~] = size(uvel);
u_rho_mask = NaN(nrows, ncols, numlevels);
u_grid_mask = NaN(nrows, ncols, numlevels);

fprintf('Processing %dx%d grid points\n', nrows, ncols);

for i = 1:nrows
    if mod(i, 50) == 0
        fprintf('  Processing row %d of %d\n', i, nrows);
    end
    
    for j = 1:ncols
        % Get depths at this grid point
        depth_profile = abs(squeeze(u_zlev(i,j,:)));
        
        % Interpolate u-velocity to z-levels
        u_interp = interp1(depth_profile, squeeze(uvel(i,j,:)), zlevels, 'linear');
        
        % Method 1: Mask by rho-grid bathymetry (current approach)
        u_rho = u_interp;
        for k = 1:numlevels
            if zlevels(k) > h(i,j)
                u_rho(k) = NaN;
            end
        end
        u_rho_mask(i,j,:) = u_rho;
        
        % Method 2: Mask by actual u-grid bathymetry (proposed approach)
        u_grid = u_interp;
        for k = 1:numlevels
            if zlevels(k) > depth_profile(end)
                u_grid(k) = NaN;
            end
        end
        u_grid_mask(i,j,:) = u_grid;
    end
end

%% Analyze results
rho_valid = ~isnan(u_rho_mask);
grid_valid = ~isnan(u_grid_mask);

total_rho = sum(rho_valid(:));
total_grid = sum(grid_valid(:));
net_preserved = total_grid - total_rho;
additional_data = sum(~rho_valid(:) & grid_valid(:));

fprintf('\n=== MASKING COMPARISON RESULTS ===\n');
fprintf('Rho-grid masking (current): %d valid points\n', total_rho);
fprintf('Grid-based masking (proposed): %d valid points\n', total_grid);
fprintf('Net data preserved: %d points (%.2f%% increase)\n', net_preserved, 100*net_preserved/total_rho);
fprintf('Additional valid data points: %d\n', additional_data);

%% Load PR bathymetry TIFF for comparison
try
    tiff_file = fullfile(dataPath, 'PR_wide_90m.tiff');
    fprintf('\nLoading PR bathymetry from: %s\n', tiff_file);
    
    % Read TIFF file
    [bathy_data, bathy_ref] = readgeoraster(tiff_file);
    
    % Get geographic coordinates from the raster
    [bathy_lat, bathy_lon] = geographicGrid(bathy_ref);
    
    % Convert bathymetry longitudes to match ROMS grid if needed
    if any(bathy_lon(:) < 0) && any(lon_u(:) > 180)
        bathy_lon = 360 + bathy_lon; % Convert to 0-360 if ROMS uses positive
    end
    
    fprintf('Bathymetry TIFF loaded: %dx%d points\n', size(bathy_data));
    fprintf('Bathy domain: %.3f°-%.3f° lon, %.3f°-%.3f° lat\n', ...
            min(bathy_lon(:)), max(bathy_lon(:)), min(bathy_lat(:)), max(bathy_lat(:)));
    
    % Create landmask from bathymetry (positive values = land)
    tiff_landmask = bathy_data > 0;
    
    % Interpolate TIFF landmask to ROMS u-grid
    u_grid_landmask = interp2(bathy_lon, bathy_lat, double(tiff_landmask), lon_u, lat_u, 'nearest', 0);
    u_grid_landmask = logical(u_grid_landmask);
    
    fprintf('TIFF landmask interpolated to ROMS u-grid\n');
    
catch ME
    fprintf('Error loading TIFF bathymetry: %s\n', ME.message);
    fprintf('Using model-derived landmask\n');
    u_grid_landmask = [];
end

%% Create focused plot at 50m depth with geographic context
target_depth = 50;
[~, depth_idx] = min(abs(zlevels - target_depth));
actual_depth = zlevels(depth_idx);

% Get masks at target depth
rho_mask_depth = ~isnan(u_rho_mask(:,:,depth_idx));
grid_mask_depth = ~isnan(u_grid_mask(:,:,depth_idx));
lost_data = grid_mask_depth & ~rho_mask_depth;

% Get landmask (TIFF if available, otherwise model-derived)
if ~isempty(u_grid_landmask)
    landmask_to_use = u_grid_landmask;
    landmask_source = 'TIFF';
else
    surface_mask = ~isnan(u_rho_mask(:,:,1));
    landmask_to_use = ~surface_mask;
    landmask_source = 'Model';
end

% Create figure with geographic coordinates
figure('Position', [200, 200, 800, 600]);

% Create composite map
composite_map = zeros(size(rho_mask_depth));
composite_map(landmask_to_use) = 0;       % Land = 0 (dark gray)
composite_map(~rho_mask_depth & ~landmask_to_use) = 1; % Deep seafloor = 1 (black)  
composite_map(rho_mask_depth) = 2;        % Valid data = 2 (blue)
composite_map(lost_data) = 3;             % Lost data = 3 (red)

% Plot with geographic coordinates
pcolor(lon_u, lat_u, composite_map);
colormap([0.3 0.3 0.3; 0 0 0; 0 0.5 1; 1 0 0]);
shading flat;
axis equal; axis tight;

% Add geographic labels
xlabel('Longitude (°E)');
ylabel('Latitude (°N)');
title(sprintf('Data Loss Analysis at %.1fm Depth - Lost Points: %d\n%s Landmask: Dark Gray=Land, Black=Deep Seafloor, Blue=Valid Data, Red=Lost Data', ...
              actual_depth, sum(lost_data(:)), landmask_source));

% Add colorbar
h = colorbar('Ticks', [0.375, 1.125, 1.875, 2.625], 'TickLabels', {'Land', 'Deep Seafloor', 'Valid Data', 'Lost Data'});
h.Label.String = 'Data Status';

% Print statistics
fprintf('\nFocused plot created showing data loss at %.1fm depth\n', actual_depth);
fprintf('Lost data points: %d (%.2f%% of total valid grid data)\n', ...
        sum(lost_data(:)), 100*sum(lost_data(:))/sum(grid_mask_depth(:)));
fprintf('Using %s landmask for reference\n', landmask_source);

if ~isempty(u_grid_landmask)
    % Compare landmasks
    model_land = ~(~isnan(u_rho_mask(:,:,1)));
    tiff_land_coverage = sum(u_grid_landmask(:));
    model_land_coverage = sum(model_land(:));
    fprintf('\nLandmask comparison:\n');
    fprintf('  TIFF landmask: %d land points\n', tiff_land_coverage);
    fprintf('  Model landmask: %d land points\n', model_land_coverage);
    fprintf('  Difference: %d points\n', tiff_land_coverage - model_land_coverage);
end

fprintf('\nComparison complete!\n');




%% Figure 2: TIFF bathymetry base with ROMS data overlay (properly cropped)
target_depth = 40;
[~, depth_idx] = min(abs(zlevels - target_depth));

% Load TIFF bathymetry
tiff_file = fullfile(dataPath, 'PR_wide_90m.tiff');
[bathy_data, bathy_ref] = readgeoraster(tiff_file);
[bathy_lat, bathy_lon] = geographicGrid(bathy_ref);
if any(bathy_lon(:) < 0), bathy_lon = 360 + bathy_lon; end

% Get ROMS domain bounds
roms_lon_bounds = [min(lon_u(:)) max(lon_u(:))];
roms_lat_bounds = [min(lat_u(:)) max(lat_u(:))];

% Crop TIFF data to ROMS domain
lon_mask = bathy_lon >= roms_lon_bounds(1) & bathy_lon <= roms_lon_bounds(2);
lat_mask = bathy_lat >= roms_lat_bounds(1) & bathy_lat <= roms_lat_bounds(2);
tiff_crop_mask = lon_mask & lat_mask;

bathy_data_crop = bathy_data;
bathy_data_crop(~tiff_crop_mask) = NaN;

% Create TIFF base layer
tiff_land = bathy_data_crop > 0;
tiff_shallow = (bathy_data_crop <= 0) & (bathy_data_crop >= -target_depth);
tiff_deep = bathy_data_crop < -target_depth;

composite_base = NaN(size(bathy_data_crop));
composite_base(tiff_deep) = 1;    % Deep water = 1 (dark blue)
composite_base(tiff_shallow) = 2; % Shallow seafloor = 2 (light blue)  
composite_base(tiff_land) = 3;    % Land = 3 (brown)

% Plot TIFF base layer
figure('Position', [300, 200, 900, 700]);
pcolor(bathy_lon, bathy_lat, composite_base);
colormap([0 0 0.5; 0.5 0.8 1; 0.8 0.6 0.4]); % dark blue, light blue, brown
shading flat; hold on;

% Create ROMS seafloor mask and overlay as semi-transparent layer
rho_mask_depth = ~isnan(u_rho_mask(:,:,depth_idx));
grid_mask_depth = ~isnan(u_grid_mask(:,:,depth_idx));
lost_data = grid_mask_depth & ~rho_mask_depth;
roms_seafloor_at_depth = ~rho_mask_depth;

% Overlay ROMS seafloor as semi-transparent black mask
roms_overlay = zeros(size(roms_seafloor_at_depth));
roms_overlay(roms_seafloor_at_depth) = 1;
[C1, h1] = contourf(lon_u, lat_u, roms_overlay, [0.5 1.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
h1.FaceColor = 'k';

% Overlay lost data as red mask
lost_overlay = zeros(size(lost_data));
lost_overlay(lost_data) = 1;
[C2, h2] = contourf(lon_u, lat_u, lost_overlay, [0.5 1.5], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
h2.FaceColor = 'r';

axis equal tight;
xlabel('Longitude (°E)'); ylabel('Latitude (°N)');
title(sprintf('TIFF Bathymetry + ROMS Seafloor at %.0fm\nBlack=ROMS Seafloor, Red=Lost Points (%d)', ...
              target_depth, sum(lost_data(:))));

h = colorbar('Ticks', [1.33, 2, 2.67], ...
             'TickLabels', {sprintf('Deep (>%dm)', target_depth), sprintf('Shallow (0-%dm)', target_depth), 'Land'});



%% Figure 3: TIFF vs ROMS bathymetry difference with 60m contour
figure('Position', [400, 300, 900, 700]);

% Define the plot domain (Vieques region)
plot_lon_range = [294.5 296.2];
plot_lat_range = [17.4 19.2];

% Extract TIFF data in plot domain
% Find indices for lon and lat separately (they're 1D in geographicGrid output)
tiff_lon_vec = bathy_lon(1,:);  % Longitude is along columns
tiff_lat_vec = bathy_lat(:,1);  % Latitude is along rows
tiff_lon_idx = tiff_lon_vec >= plot_lon_range(1) & tiff_lon_vec <= plot_lon_range(2);
tiff_lat_idx = tiff_lat_vec >= plot_lat_range(1) & tiff_lat_vec <= plot_lat_range(2);

% Subset TIFF data
tiff_lon_plot = bathy_lon(tiff_lat_idx, tiff_lon_idx);
tiff_lat_plot = bathy_lat(tiff_lat_idx, tiff_lon_idx);
tiff_bathy_plot = bathy_data(tiff_lat_idx, tiff_lon_idx);

% Extract ROMS data in plot domain
roms_bathy = ncread(fullfile(dataPath, 'vigrid.nc'), 'h');

% Find bounding box for ROMS data in plot domain
roms_in_domain = lon_rho >= plot_lon_range(1) & lon_rho <= plot_lon_range(2) & ...
                 lat_rho >= plot_lat_range(1) & lat_rho <= plot_lat_range(2);
[rows, cols] = find(roms_in_domain);
row_range = min(rows):max(rows);
col_range = min(cols):max(cols);

% Subset ROMS data
lon_rho_plot = lon_rho(row_range, col_range);
lat_rho_plot = lat_rho(row_range, col_range);
roms_depth_plot = roms_bathy(row_range, col_range);

% Interpolate TIFF to ROMS grid in plot domain (no extrapolation)
tiff_on_roms_plot = interp2(tiff_lon_plot, tiff_lat_plot, double(tiff_bathy_plot), ...
                            lon_rho_plot, lat_rho_plot, 'linear', NaN);

% Convert TIFF to positive depths and mask invalid regions
tiff_depth_plot = -tiff_on_roms_plot;
tiff_depth_plot(tiff_on_roms_plot > 0) = NaN;  % Mask land
tiff_depth_plot(isnan(tiff_on_roms_plot)) = NaN;  % Mask where TIFF has no data

% Calculate difference
bathy_diff_plot = tiff_depth_plot - roms_depth_plot;

% Clip depths below 250m AND mask where TIFF is invalid
bathy_diff_plot(tiff_depth_plot > 250 | roms_depth_plot > 250 | isnan(tiff_depth_plot)) = NaN;

% Plot difference
pcolor(lon_rho_plot, lat_rho_plot, bathy_diff_plot);
shading flat; hold on;


colormap(jet);

clim([-50 50]);

% Add 60m contour from TIFF in plot domain - white with thicker line
contour(tiff_lon_plot, tiff_lat_plot, tiff_bathy_plot, [-60 -60], 'w-', 'LineWidth', 2);

axis equal tight;
xlabel('Longitude (°E)'); ylabel('Latitude (°N)');
title('Bathymetry Difference: TIFF - ROMS with 60m Contour\nPositive = TIFF Deeper, Negative = ROMS Deeper, White Line = TIFF 60m, Black = >250m depth');

h = colorbar;
h.Label.String = 'Depth Difference (m)';

% Print statistics
valid_diff = bathy_diff_plot(~isnan(bathy_diff_plot));
fprintf('Bathymetry difference statistics (0-250m depth):\n');
fprintf('  Mean difference: %.1f m\n', mean(valid_diff));
fprintf('  Max difference: %.1f m\n', max(abs(valid_diff)));
fprintf('  Max TIFF depth in valid region: %.1f m\n', max(tiff_depth_plot(~isnan(tiff_depth_plot))));



%% Figure 4: ROMS Bathymetry
figure('Position', [500, 400, 900, 700]);

% Load ROMS bathymetry and coordinates
lon_rho = ncread(fullfile(dataPath, 'vigrid.nc'), 'lon_rho');
lat_rho = ncread(fullfile(dataPath, 'vigrid.nc'), 'lat_rho');
h = ncread(fullfile(dataPath, 'vigrid.nc'), 'h');

% Convert longitudes if needed
if any(lon_rho(:) < 0)
    lon_rho = 360 + lon_rho;
end

% Make depths negative and mask land (0 or NaN)
h_plot = -h;
h_plot(h == 0 | isnan(h)) = NaN;

pcolor(lon_rho, lat_rho, h_plot);
shading flat;
colormap(flip(parula));
clim([-250 0]);

% Make NaN values appear as grey
set(gca, 'Color', [0.7 0.7 0.7]);

axis equal tight;
xlabel('Longitude (°E)'); 
ylabel('Latitude (°N)');
title('ROMS Bathymetry (h)');
hbar = colorbar;
hbar.Label.String = 'Depth (m)';

fprintf('\nROMS bathymetry: Min=%.1fm, Max=%.1fm, Mean=%.1fm\n', ...
        min(h_plot,[],'all'), max(h_plot,[],'all'), mean(h_plot,'all'));