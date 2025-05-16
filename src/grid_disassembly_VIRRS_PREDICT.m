%% Script to assign break apart, plot, and analyze the different ROMS grids in order
%   to produce non-land coordinate indicees for 'VIRRS_velocities.m'
% 6 June 2024

clear;clc

inputDrive = '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023';
scriptDrive = '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/Release';
writeDrive = '/Users/benja/Downloads';

%extract USCROMS bathymetry and rho-coordinates
cd(inputDrive)
USCROMSgrid = 'vigrid.nc';
% ncdisp(USCROMSgrid)
bathymetry = ncread(USCROMSgrid, 'h');
longitudes_rho_USCROMS = ncread(USCROMSgrid, 'lon_rho') + 360;
latitudes_rho_USCROMS = ncread(USCROMSgrid, 'lat_rho');
longitudes_u_USCROMS = ncread(USCROMSgrid, 'lon_u') + 360;
latitudes_u_USCROMS = ncread(USCROMSgrid, 'lat_u');
longitudes_v_USCROMS = ncread(USCROMSgrid, 'lon_v') + 360;
latitudes_v_USCROMS = ncread(USCROMSgrid, 'lat_v');
mask = ncread(USCROMSgrid, 'mask_rho');

%make sure land is ignored
bathymetry(bathymetry==0) = NaN;
bathymetry(mask == 0) = NaN;

%construct landmasks for u- and v-grids
cd(scriptDrive)
USCROMSocean = 'croco_his.03969.nc';
wvel = ncread(USCROMSocean, 'w');
uvel = ncread(USCROMSocean, 'u');
uvel = uvel(:,:,:,1);
uvel = uvel(:,:,size(uvel, 3):-1:1); %flip so velocities are surface to bottom
uvel = uvel(:,:,1); %extract surface layer
mask_u = uvel;
mask_u(mask_u ~= 0) = 1; %the ocean is 1's, land is 0's
vvel = ncread(USCROMSocean, 'v');
vvel = vvel(:,:,:,1);
vvel = vvel(:,:,size(vvel, 3):-1:1); %flip so velocities are surface to bottom
vvel = vvel(:,:,1); %extract surface layer
mask_v = vvel;
mask_v(mask_v ~= 0) = 1;

% uvel(idx_u_USCROMS_VIRRS)
% vvel(idx_v_USCROMS_VIRRS)

%extract release point coordinates from GIS output and ensure they are
% sorted by their unique ID
relpoints = readmatrix('points_650_none-on-land.csv');
relpoints = relpoints(:, 10:12);
relpoints = sortrows(relpoints, 1);

IDs_PREDICT = relpoints(:,1);
longitudes_PREDICT = relpoints(:,2) + 360;
latitudes_PREDICT = relpoints(:,3);

% Load VIRRS points, then generate and add IDs to VIRRSpoints and its metadata, then sort
VIRRSpoints = readmatrix('RRS_siteMaster_abbreviated_18July2023.csv');
VIRRSpoints_metadata = readtable('RRS_siteMaster_abbreviated_18July2023.csv');
IDs_VIRRS = (0:size(VIRRSpoints, 1) - 1)'; % Assuming each row gets a unique ID
VIRRSpoints = [IDs_VIRRS, VIRRSpoints];
VIRRSpoints_metadata = addvars(VIRRSpoints_metadata, IDs_VIRRS, 'Before', 1, 'NewVariableNames', 'ID');
VIRRSpoints = sortrows(VIRRSpoints, 1);
VIRRSpoints_metadata = sortrows(VIRRSpoints_metadata, 'ID');
VIRRSpoints_metadata.lon = VIRRSpoints_metadata.lon+360;

longitudes_VIRRS = VIRRSpoints_metadata.lon;
latitudes_VIRRS = VIRRSpoints_metadata.lat;
depths_VIRRS = VIRRSpoints_metadata.depth_m;
maxdepth_VIRRS = max(depths_VIRRS);

% Locate USCROMS rho-grid points (and associated bathymetric values)
% nearest to PREDICT points, as well as u- and v-grid points
idx_rho_USCROMS_PREDICT = knnsearch([longitudes_rho_USCROMS(:), latitudes_rho_USCROMS(:)], [longitudes_PREDICT, latitudes_PREDICT]);
idx_u_USCROMS_PREDICT = knnsearch([longitudes_u_USCROMS(:), latitudes_u_USCROMS(:)], [longitudes_PREDICT, latitudes_PREDICT]);
idx_v_USCROMS_PREDICT = knnsearch([longitudes_v_USCROMS(:), latitudes_v_USCROMS(:)], [longitudes_PREDICT, latitudes_PREDICT]);

% Locate USCROMS rho-grid points (and associated bathymetric values) nearest to VIRRS points
idx_rho_USCROMS_VIRRS = knnsearch([longitudes_rho_USCROMS(:), latitudes_rho_USCROMS(:)], [longitudes_VIRRS, latitudes_VIRRS]);
idx_u_USCROMS_VIRRS = knnsearch([longitudes_u_USCROMS(:), latitudes_u_USCROMS(:)], [longitudes_VIRRS, latitudes_VIRRS]);
idx_v_USCROMS_VIRRS = knnsearch([longitudes_v_USCROMS(:), latitudes_v_USCROMS(:)], [longitudes_VIRRS, latitudes_VIRRS]);

% Identify points in each grid that are NaNs on their respective landmasks
land_rho_VIRRS = mask(idx_rho_USCROMS_VIRRS); %land is '0s', ocean is '1's
land_u_VIRRS = mask_u(idx_u_USCROMS_VIRRS);
land_v_VIRRS = mask_v(idx_v_USCROMS_VIRRS);

%% Construct the landmask polygons (boundaries)

%rho
% Binary landmask contours
mask_plot = mask;
mask_plot(mask_plot == 0) = NaN;
mask_plot(~isnan(mask_plot)) = 0;
mask_plot(isnan(mask_plot)) = 1;

% Lon/lat vectors
lat = latitudes_rho_USCROMS(1, :)';
lon = longitudes_rho_USCROMS(:, 1);

% Make the landmask
trace = bwboundaries(mask_plot); 
for i = 1:length(trace)
    trace{i}(:, 1) = lon(trace{i}(:, 1));
    trace{i}(:, 2) = lat(trace{i}(:, 2));
end

% Create land polygons
land = struct();
for k = 1:length(trace)
    land(k).Geometry = 'Polygon';
    Xs = trace{k}(:, 1); % Might have to subtract 360 at the end
    Xs(end + 1) = NaN;
    land(k).X = Xs';
    Ys = trace{k}(:, 2);
    Ys(end + 1) = NaN;
    land(k).Y = Ys';
end

%u
% Binary landmask contours
mask_plot = mask_u;
mask_plot(mask_plot == 0) = NaN;
mask_plot(~isnan(mask_plot)) = 0;
mask_plot(isnan(mask_plot)) = 1;

% Lon/lat vectors
lat = latitudes_u_USCROMS(1, :)';
lon = longitudes_u_USCROMS(:, 1);

% Make the landmask
trace = bwboundaries(mask_plot); 
for i = 1:length(trace)
    trace{i}(:, 1) = lon(trace{i}(:, 1));
    trace{i}(:, 2) = lat(trace{i}(:, 2));
end

% Create land polygons
land_u = struct();
for k = 1:length(trace)
    land_u(k).Geometry = 'Polygon';
    Xs = trace{k}(:, 1); % Might have to subtract 360 at the end
    Xs(end + 1) = NaN;
    land_u(k).X = Xs';
    Ys = trace{k}(:, 2);
    Ys(end + 1) = NaN;
    land_u(k).Y = Ys';
end

%v
% Binary landmask contours
mask_plot = mask_v;
mask_plot(mask_plot == 0) = NaN;
mask_plot(~isnan(mask_plot)) = 0;
mask_plot(isnan(mask_plot)) = 1;

% Lon/lat vectors
lat = latitudes_v_USCROMS(1, :)';
lon = longitudes_v_USCROMS(:, 1);

% Make the landmask
trace = bwboundaries(mask_plot); 
for i = 1:length(trace)
    trace{i}(:, 1) = lon(trace{i}(:, 1));
    trace{i}(:, 2) = lat(trace{i}(:, 2));
end

% Create land polygons
land_v = struct();
for k = 1:length(trace)
    land_v(k).Geometry = 'Polygon';
    Xs = trace{k}(:, 1); % Might have to subtract 360 at the end
    Xs(end + 1) = NaN;
    land_v(k).X = Xs';
    Ys = trace{k}(:, 2);
    Ys(end + 1) = NaN;
    land_v(k).Y = Ys';
end

%% Identify coordinates in each grid system that are the closest to VIRRS points and also in the ocean

% Define the Haversine formula for distance calculation
haversine = @(lat1, lon1, lat2, lon2) 2 * 6371 * asin(sqrt(sin((deg2rad(lat2) - deg2rad(lat1)) / 2).^2 + ...
    cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin((deg2rad(lon2) - deg2rad(lon1)) / 2).^2));

for i = 1:numel(land_rho_VIRRS)

    %rho-grid
    % Check if the grid-specific nearest-neighbor USCROMS coordinate is on
    % land ('0')
    if land_rho_VIRRS(i) == 0

        % Initialize the minimum distance and index
        min_distance = Inf;
        min_idx = -1;

        % Iterate through USCROMS coordinates to find the nearest one that is not on land
        for j = 1:numel(longitudes_rho_USCROMS)

            % Calculate the distance between the VIRRS point and the USCROMS point using the Haversine formula
            d = haversine(latitudes_VIRRS(i), longitudes_VIRRS(i), latitudes_rho_USCROMS(j), longitudes_rho_USCROMS(j));

            % Check if the USCROMS point is not on land and closer than the current minimum distance
            if mask(j) && d < min_distance
                min_distance = d;
                min_idx = j;
            end
        end

        % Assign the nearest USCROMS index for the VIRRS point
        idx_rho_USCROMS_VIRRS(i) = min_idx;
    end

    %u-grid
    if land_u_VIRRS(i) == 0

        % Initialize the minimum distance and index
        min_distance = Inf;
        min_idx = -1;

        % Iterate through USCROMS coordinates to find the nearest one that is not on land
        for j = 1:numel(mask_u)

            % Calculate the distance between the VIRRS point and the USCROMS point using the Haversine formula
            d = haversine(latitudes_VIRRS(i), longitudes_VIRRS(i), latitudes_u_USCROMS(j), longitudes_u_USCROMS(j));

            % Check if the USCROMS point is not on land and closer than the current minimum distance
            if mask_u(j) && d < min_distance
                min_distance = d;
                min_idx = j;
            end
        end

        % Assign the nearest USCROMS index for the VIRRS point
        idx_u_USCROMS_VIRRS(i) = min_idx;
    end

    %v-grid
    if land_v_VIRRS(i) == 0

        % Initialize the minimum distance and index
        min_distance = Inf;
        min_idx = -1;

        % Iterate through USCROMS coordinates to find the nearest one that is not on land
        for j = 1:numel(mask_v)

            % Calculate the distance between the VIRRS point and the USCROMS point using the Haversine formula
            d = haversine(latitudes_VIRRS(i), longitudes_VIRRS(i), latitudes_v_USCROMS(j), longitudes_v_USCROMS(j));

            % Check if the USCROMS point is not on land and closer than the current minimum distance
            if mask_v(j) && d < min_distance
                min_distance = d;
                min_idx = j;
            end
        end

        % Assign the nearest USCROMS index for the VIRRS point
        idx_v_USCROMS_VIRRS(i) = min_idx;
    end
end

% Save the idx_rho_USCROMS_VIRRS, idx_u_USCROMS_VIRRS, and idx_v_USCROMS_VIRRS vectors to a .mat file
save('USCROMS_VIRRS_indices.mat', 'idx_rho_USCROMS_VIRRS', 'idx_u_USCROMS_VIRRS', 'idx_v_USCROMS_VIRRS');

%% Compare the depths of nearest-neighbor ROMS points with real depths of VIRRS LTR sites

% Ensure the indices match correctly (if not already)
depths_VIRRS = depths_VIRRS(IDs_VIRRS + 1);

bathymetry_rho_VIRRS = bathymetry(idx_rho_USCROMS_VIRRS);

% Create a new table comparing bathymetry_VIRRS with the depth values from VIRRSpoints
comparison_table = table(IDs_VIRRS, depths_VIRRS, bathymetry_rho_VIRRS, ...
    'VariableNames', {'ID', 'Real_Depth_VIRRS', 'USCROMS_Bathymetry_VIRRS'});

% Display the first few rows of the comparison table
disp(comparison_table(1:10, :));

% Save the table to a CSV file if needed
cd(writeDrive)
writetable(comparison_table, 'VIRRS_Bathymetry_Comparison.csv');
cd(scriptDrive)

%% plot the nearest neighbor-adjusted VIRRS points

% Visualize data
figure;
hold on;

% Plot the land mask separately without adding it to the legend directly
for k = 1:length(land)
    plot(land(k).X, land(k).Y, 'k'); % Black lines for the land mask
end

% Plot VIRRS points
h1 = plot(longitudes_VIRRS, latitudes_VIRRS, '*g'); % Green stars for VIRRS points
h2 = plot(longitudes_rho_USCROMS(idx_rho_USCROMS_VIRRS), latitudes_rho_USCROMS(idx_rho_USCROMS_VIRRS), 'om'); % Magenta circles for nearest USCROMS coordinates for VIRRS

% Indicate VIRRS points on land in ROMS
onland_idx = logical(1 - land_rho_VIRRS); % '1' now indicates land
% nan_bathymetry_idx = isnan(bathymetry_rho_VIRRS);

% h3 = plot(longitudes_VIRRS(onland_idx), latitudes_VIRRS(onland_idx), 'xr', 'LineWidth', 2, 'MarkerSize', 10); % Red crosses for NaN bathymetry VIRRS points

% plot(longitudes_VIRRS, latitudes_VIRRS, '*g');hold on
% plot(longitudes_rho_USCROMS(idx_rho_USCROMS_VIRRS), latitudes_rho_USCROMS(idx_rho_USCROMS_VIRRS), 'om');hold on
% Draw lines between VIRRS points and nearest USCROMS coordinates
for i = 1:numel(longitudes_VIRRS)
    % Check if the nearest USCROMS coordinate is in the ocean
    % if onland_idx(i)
        % Plot a line between the VIRRS point and its nearest USCROMS point
        plot([longitudes_VIRRS(i), longitudes_rho_USCROMS(idx_rho_USCROMS_VIRRS(i))], ...
             [latitudes_VIRRS(i), latitudes_rho_USCROMS(idx_rho_USCROMS_VIRRS(i))], 'm--');
    % end
end

% Add the legend using the handles
legend([h1, h2], ...
    {'VIRRS Points', 'Nearest USCROMS for VIRRS'}, ...
    'Location', 'best');

title('ROMS Coordinates and VIRRS Points with NaN Bathymetry Indication');
xlabel('Longitude');
ylabel('Latitude');
grid on;
hold off;

% Set aspect ratio to be equal
daspect([1, 1, 1]);

%% plot the nearest neighbor-adjusted VIRRS points, with the USCROMS grid as well

figure;
hold on;

% Plot the land mask separately without adding it to the legend directly
for k = 1:length(land)
    plot(land(k).X, land(k).Y, 'k'); % Black lines for the land mask
end

% Plot VIRRS points
h1 = plot(longitudes_VIRRS, latitudes_VIRRS, '*g'); % Green stars for VIRRS points

% Plot nearest USCROMS coordinates for VIRRS
h2 = plot(longitudes_rho_USCROMS(idx_rho_USCROMS_VIRRS), latitudes_rho_USCROMS(idx_rho_USCROMS_VIRRS), 'om'); % Magenta circles for nearest USCROMS coordinates for VIRRS

% Plot every single USCROMS coordinate
h3 = plot(longitudes_rho_USCROMS(:), latitudes_rho_USCROMS(:), 'or', 'MarkerSize', 2); % Red circles for every single USCROMS coordinate

% Draw lines between VIRRS points and nearest USCROMS coordinates
for i = 1:numel(longitudes_VIRRS)
    plot([longitudes_VIRRS(i), longitudes_rho_USCROMS(idx_rho_USCROMS_VIRRS(i))], ...
     [latitudes_VIRRS(i), latitudes_rho_USCROMS(idx_rho_USCROMS_VIRRS(i))], 'm--');
end

% Add the legend using the handles
legend([h1, h2, h3], ...
    {'VIRRS Points', 'Nearest USCROMS for VIRRS', 'USCROMS Coordinates'}, ...
    'Location', 'best');

% Set x and y limits for St John
xlim([360-64.85, 360-64.65]); % Longitude range for St John
ylim([18.25, 18.45]); % Latitude range for St John

title('ROMS Coordinates, VIRRS Points, and Every USCROMS Coordinate');
xlabel('Longitude');
ylabel('Latitude');
grid on;
hold off;

% Set aspect ratio to be equal
daspect([1, 1, 1]);

%% plot the u,v, and rho-grids separately, with rho-landmask

% Plot the land mask separately without adding it to the legend directly
figure;
hold on;
for k = 1:length(land)
    plot(land(k).X, land(k).Y, 'k'); % Black lines for the land mask
end

% Plot every single USCROMS coordinate
% h_rho = plot(longitudes_rho_USCROMS(:), latitudes_rho_USCROMS(:), 'or', 'MarkerSize', 3);xlim([360-65.1, 360-64.8]);ylim([18.2, 18.5]);hold on
% h_u = plot(longitudes_u_USCROMS(:), latitudes_u_USCROMS(:), 'ob', 'MarkerSize', 3);xlim([360-65.1, 360-64.8]);ylim([18.2, 18.5])
% h_v = plot(longitudes_v_USCROMS(:), latitudes_v_USCROMS(:), 'og', 'MarkerSize', 3);xlim([360-65.1, 360-64.8]);ylim([18.2, 18.5])

h_rho = plot(longitudes_rho_USCROMS(:), latitudes_rho_USCROMS(:), 'or', 'MarkerSize', 3);%xlim([360-65.1, 360-64.8]);ylim([18.2, 18.5]);hold on
h_u = plot(longitudes_u_USCROMS(:), latitudes_u_USCROMS(:), 'ob', 'MarkerSize', 3);%xlim([360-65.1, 360-64.8]);ylim([18.2, 18.5])
h_v = plot(longitudes_v_USCROMS(:), latitudes_v_USCROMS(:), 'og', 'MarkerSize', 3);%xlim([360-65.1, 360-64.8]);ylim([18.2, 18.5])

% Set aspect ratio to be equal
daspect([1, 1, 1]);

% Add legend
legend([h_rho, h_u, h_v], {'Rho-coordinates', 'U-coordinates', 'V-coordinates'}, 'Location', 'best');

%% plot the u,v, and rho-grids separately, with all landmasks


% Plot the land mask separately without adding it to the legend directly
figure;
hold on;

% Plot land mask for rho-grid
for k = 1:length(land)
    plot(land(k).X, land(k).Y, 'r'); % Red lines for the land mask for rho-grid
end

% Plot land mask for u-grid
for k = 1:length(land_u)
    plot(land_u(k).X, land_u(k).Y, 'b'); % Blue lines for the land mask for u-grid
end

% Plot land mask for v-grid
for k = 1:length(land_v)
    plot(land_v(k).X, land_v(k).Y, 'g'); % Green lines for the land mask for v-grid
end

% Plot every single USCROMS coordinate
h_rho = plot(longitudes_rho_USCROMS(:), latitudes_rho_USCROMS(:), 'or', 'MarkerSize', 3); % Red circles for every single rho-coordinate
h_u = plot(longitudes_u_USCROMS(:), latitudes_u_USCROMS(:), 'ob', 'MarkerSize', 3); % Blue circles for every single u-coordinate
h_v = plot(longitudes_v_USCROMS(:), latitudes_v_USCROMS(:), 'og', 'MarkerSize', 3); % Green circles for every single v-coordinate

% Set aspect ratio to be equal
daspect([1, 1, 1]);

% Add legend
legend([h_rho, h_u, h_v], {'Rho-coordinates', 'U-coordinates', 'V-coordinates'}, 'Location', 'best');

% Set x and y limits for St John
xlim([360-64.85, 360-64.65]); % Longitude range for St John
ylim([18.25, 18.45]); % Latitude range for St John

% Add title and labels
title('USCROMS Coordinates for St John, USVI');
xlabel('Longitude');
ylabel('Latitude');

%% Plot the nearest neighbor-adjusted VIRRS points on rho-grid, u-grid, and v-grid

% Visualize data
figure;
hold on;

% Plot the land mask separately without adding it to the legend directly
for k = 1:length(land_v)
    plot(land_v(k).X, land_v(k).Y, 'k'); % Black lines for the land mask
end

% Plot VIRRS points
h1 = plot(longitudes_VIRRS, latitudes_VIRRS, '*g'); % Green stars for VIRRS points

% Plot nearest USCROMS coordinates for VIRRS on rho-grid
h2 = plot(longitudes_rho_USCROMS(idx_rho_USCROMS_VIRRS), latitudes_rho_USCROMS(idx_rho_USCROMS_VIRRS), 'om'); % Magenta circles for rho-grid

% Plot nearest USCROMS coordinates for VIRRS on u-grid
h3 = plot(longitudes_u_USCROMS(idx_u_USCROMS_VIRRS), latitudes_u_USCROMS(idx_u_USCROMS_VIRRS), 'sc'); % Cyan squares for u-grid

% Plot nearest USCROMS coordinates for VIRRS on v-grid
h4 = plot(longitudes_v_USCROMS(idx_v_USCROMS_VIRRS), latitudes_v_USCROMS(idx_v_USCROMS_VIRRS), '^b'); % Blue triangles for v-grid

% Draw lines between VIRRS points and nearest USCROMS coordinates on rho-grid
for i = 1:numel(longitudes_VIRRS)
    plot([longitudes_VIRRS(i), longitudes_rho_USCROMS(idx_rho_USCROMS_VIRRS(i))], ...
         [latitudes_VIRRS(i), latitudes_rho_USCROMS(idx_rho_USCROMS_VIRRS(i))], 'm--');
end

% Draw lines between VIRRS points and nearest USCROMS coordinates on u-grid
for i = 1:numel(longitudes_VIRRS)
    plot([longitudes_VIRRS(i), longitudes_u_USCROMS(idx_u_USCROMS_VIRRS(i))], ...
         [latitudes_VIRRS(i), latitudes_u_USCROMS(idx_u_USCROMS_VIRRS(i))], 'c--');
end

% Draw lines between VIRRS points and nearest USCROMS coordinates on v-grid
for i = 1:numel(longitudes_VIRRS)
    plot([longitudes_VIRRS(i), longitudes_v_USCROMS(idx_v_USCROMS_VIRRS(i))], ...
         [latitudes_VIRRS(i), latitudes_v_USCROMS(idx_v_USCROMS_VIRRS(i))], 'b--');
end

% Add legend
legend([h1, h2, h3, h4], ...
    {'VIRRS Points', 'Nearest USCROMS for VIRRS (rho-grid)', 'Nearest USCROMS for VIRRS (u-grid)', 'Nearest USCROMS for VIRRS (v-grid)'}, ...
    'Location', 'best');

title('ROMS Coordinates and VIRRS Points with NaN Bathymetry Indication');
xlabel('Longitude');
ylabel('Latitude');
grid on;
hold off;

% Set aspect ratio to be equal
daspect([1, 1, 1]);



%%
figure;
hold on;

% Plot the land mask for rho-grid
for k = 1:length(land)
    plot(land(k).X, land(k).Y, 'r', 'LineWidth', 1); % Red lines for the land mask for rho-grid
end

% Plot the land mask for u-grid
for k = 1:length(land_u)
    plot(land_u(k).X, land_u(k).Y, 'b', 'LineWidth', 1); % Blue lines for the land mask for u-grid
end

% Plot the land mask for v-grid
for k = 1:length(land_v)
    plot(land_v(k).X, land_v(k).Y, 'g', 'LineWidth', 1); % Green lines for the land mask for v-grid
end

% Plot every single USCROMS coordinate
h_rho = plot(longitudes_rho_USCROMS(:), latitudes_rho_USCROMS(:), 'or', 'MarkerSize', 4, 'LineWidth', 1); % Red circles for every single rho-coordinate
h_u = plot(longitudes_u_USCROMS(:), latitudes_u_USCROMS(:), 'ob', 'MarkerSize', 4, 'LineWidth', 1); % Blue circles for every single u-coordinate
h_v = plot(longitudes_v_USCROMS(:), latitudes_v_USCROMS(:), 'og', 'MarkerSize', 4, 'LineWidth', 1); % Green circles for every single v-coordinate

% Plot VIRRS points
h1 = plot(longitudes_VIRRS, latitudes_VIRRS, '*k', 'MarkerSize', 6, 'LineWidth', 1); % Black stars for VIRRS points

% Plot nearest USCROMS coordinates for VIRRS for each grid
h2_rho = plot(longitudes_rho_USCROMS(idx_rho_USCROMS_VIRRS), latitudes_rho_USCROMS(idx_rho_USCROMS_VIRRS), 'om', 'MarkerSize', 5, 'LineWidth', 1); % Magenta circles for nearest rho-coordinates for VIRRS
h2_u = plot(longitudes_u_USCROMS(idx_u_USCROMS_VIRRS), latitudes_u_USCROMS(idx_u_USCROMS_VIRRS), 'oc', 'MarkerSize', 5, 'LineWidth', 1); % Cyan circles for nearest u-coordinates for VIRRS
h2_v = plot(longitudes_v_USCROMS(idx_v_USCROMS_VIRRS), latitudes_v_USCROMS(idx_v_USCROMS_VIRRS), 'oy', 'MarkerSize', 5, 'LineWidth', 1); % Yellow circles for nearest v-coordinates for VIRRS

% Draw lines between VIRRS points and nearest USCROMS coordinates for each grid
for i = 1:numel(longitudes_VIRRS)
    plot([longitudes_VIRRS(i), longitudes_rho_USCROMS(idx_rho_USCROMS_VIRRS(i))], ...
        [latitudes_VIRRS(i), latitudes_rho_USCROMS(idx_rho_USCROMS_VIRRS(i))], 'r--', 'LineWidth', 1); % Lines for rho grid
    plot([longitudes_VIRRS(i), longitudes_u_USCROMS(idx_u_USCROMS_VIRRS(i))], ...
        [latitudes_VIRRS(i), latitudes_u_USCROMS(idx_u_USCROMS_VIRRS(i))], 'b--', 'LineWidth', 1); % Lines for u grid
    plot([longitudes_VIRRS(i), longitudes_v_USCROMS(idx_v_USCROMS_VIRRS(i))], ...
        [latitudes_VIRRS(i), latitudes_v_USCROMS(idx_v_USCROMS_VIRRS(i))], 'g--', 'LineWidth', 1); % Lines for v grid
end

% Add the legend using the handles
legend([h_rho, h_u, h_v, h1, h2_rho, h2_u, h2_v], ...
    {'Rho-coordinates', 'U-coordinates', 'V-coordinates', 'VIRRS Points', 'Nearest Rho for VIRRS', 'Nearest U for VIRRS', 'Nearest V for VIRRS'}, ...
    'Location', 'best');

% Set x and y limits for St John
xlim([360-64.85, 360-64.65]); % Longitude range for St John
ylim([18.25, 18.45]); % Latitude range for St John

title('Nearest neighbor-adjusted VIRRS points');
xlabel('Longitude');
ylabel('Latitude');
grid on;

% Set aspect ratio to be equal
daspect([1, 1, 1]);

hold off;


%%

% Determine the corners of the grid based on the input data
min_lon = min(longitudes_rho_USCROMS(:));
max_lon = max(longitudes_rho_USCROMS(:));
min_lat = min(latitudes_rho_USCROMS(:));
max_lat = max(latitudes_rho_USCROMS(:));

% Define a zoom factor (adjust this value to control how zoomed-in the corners are)
zoom_factor = 0.05; 

% Calculate the x and y limits for each corner
corner_limits = [
    min_lon, min_lon + zoom_factor, min_lat, min_lat + zoom_factor; % Bottom-left corner
    min_lon, min_lon + zoom_factor, max_lat - zoom_factor, max_lat; % Top-left corner
    max_lon - zoom_factor, max_lon, min_lat, min_lat + zoom_factor; % Bottom-right corner
    max_lon - zoom_factor, max_lon, max_lat - zoom_factor, max_lat; % Top-right corner
];

% Create figure
figure;

% Titles and annotations for each corner
corner_titles = {'Bottom-left corner', 'Top-left corner', 'Bottom-right corner', 'Top-right corner'};
annotation_positions = [0.1, 0.1; 0.1, 0.9; 0.9, 0.1; 0.9, 0.9];

% Main plot showing the full grid with boundaries
subplot(2, 2, [1 2 3 4])
hold on;
for k = 1:length(land)
    plot(land(k).X, land(k).Y, 'k'); % Black lines for the land mask
end
% Plot the boundary of the main grid
plot([min_lon, max_lon, max_lon, min_lon, min_lon], [min_lat, min_lat, max_lat, max_lat, min_lat], 'k--');


% Highlight each corner with a rectangle
for i = 1:4
    rectangle('Position', [corner_limits(i, 1), corner_limits(i, 3), ...
        corner_limits(i, 2) - corner_limits(i, 1), corner_limits(i, 4) - corner_limits(i, 3)], ...
        'EdgeColor', 'r', 'LineStyle', '--');
end
title('Main Grid with Zoomed-in Areas');
daspect([1, 1, 1]);

%
%
%

% Create subplots for each corner with zoomed-in view
figure();
subplot_positions = [3, 1, 4, 2]; % Match corners: BL, TL, BR, TR in subplot positions
for i = 1:4
    subplot(2, 2, subplot_positions(i))
    hold on;
    for k = 1:length(land)
        plot(land(k).X, land(k).Y, 'k'); % Black lines for the land mask
    end
    % Plot points with respective colors
    p1 = plot(longitudes_rho_USCROMS(:), latitudes_rho_USCROMS(:), 'or', 'MarkerSize', 3);
    p2 = plot(longitudes_u_USCROMS(:), latitudes_u_USCROMS(:), 'ob', 'MarkerSize', 3);
    p3 = plot(longitudes_v_USCROMS(:), latitudes_v_USCROMS(:), 'og', 'MarkerSize', 3);
    xlim(corner_limits(i, 1:2));
    ylim(corner_limits(i, 3:4));
    daspect([1, 1, 1]);
    title(corner_titles{i});
    % Add background color to distinguish each corner
    set(gca, 'Color', [0.9, 0.9, 0.9]);
    % Store legend labels for each subplot
    if i == 1 % Store legend labels only once
        legends = {p1, p2, p3};
    end
end

% Combine legends into a single legend for the entire figure
legend_labels = {'Rho points', 'U points', 'V points'};
legend([legends{:}], legend_labels, 'Location', 'best');
