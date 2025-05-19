%% Script to create a landmask from USCROMS hydrodynamic files
%   19 May 2025

clear;clc

%% setup

% Get the project root directory
projectPath = matlab.project.rootProject().RootFolder;

% Define paths relative to the project root
dataPath = fullfile(projectPath, 'data');
srcPath = fullfile(projectPath, 'src');
outputPath = fullfile(projectPath, 'output');

%%

%extract USCROMS bathymetry and rho-coordinates
USCROMSgrid = 'vigrid.nc';
% ncdisp(fullfile(dataPath, USCROMSgrid))  % Pass the full path to ncdisp
bathymetry = ncread(fullfile(dataPath, USCROMSgrid), 'h');
longitudes_rho_USCROMS = ncread(fullfile(dataPath, USCROMSgrid), 'lon_rho') + 360;
latitudes_rho_USCROMS = ncread(fullfile(dataPath, USCROMSgrid), 'lat_rho');
longitudes_u_USCROMS = ncread(fullfile(dataPath, USCROMSgrid), 'lon_u') + 360;
latitudes_u_USCROMS = ncread(fullfile(dataPath, USCROMSgrid), 'lat_u');
longitudes_v_USCROMS = ncread(fullfile(dataPath, USCROMSgrid), 'lon_v') + 360;
latitudes_v_USCROMS = ncread(fullfile(dataPath, USCROMSgrid), 'lat_v');
mask = ncread(fullfile(dataPath, USCROMSgrid), 'mask_rho');

%make sure land is ignored
bathymetry(bathymetry==0) = NaN;
bathymetry(mask == 0) = NaN;

%construct landmasks for u- and v-grids
USCROMSocean = 'croco_his.03969.nc';
wvel = ncread(fullfile(dataPath, USCROMSocean), 'w');
uvel = ncread(fullfile(dataPath, USCROMSocean), 'u');
uvel = uvel(:,:,:,1);
uvel = uvel(:,:,size(uvel, 3):-1:1); %flip so velocities are surface to bottom
uvel = uvel(:,:,1); %extract surface layer
mask_u = uvel;
mask_u(mask_u ~= 0) = 1; %the ocean is 1's, land is 0's
vvel = ncread(fullfile(dataPath, USCROMSocean), 'v');
vvel = vvel(:,:,:,1);
vvel = vvel(:,:,size(vvel, 3):-1:1); %flip so velocities are surface to bottom
vvel = vvel(:,:,1); %extract surface layer
mask_v = vvel;
mask_v(mask_v ~= 0) = 1;

%extract release point coordinates from GIS output and ensure they are
% sorted by their unique ID
points_650 = 'points_650_none-on-land.csv';
relpoints = readmatrix(fullfile(dataPath, points_650));
relpoints = relpoints(:, 10:12);
relpoints = sortrows(relpoints, 1);

IDs_PREDICT = relpoints(:,1);
longitudes_PREDICT = relpoints(:,2) + 360;
latitudes_PREDICT = relpoints(:,3);

% Locate USCROMS rho-grid points (and associated bathymetric values)
% nearest to PREDICT points, as well as u- and v-grid points
idx_rho_USCROMS_PREDICT = knnsearch([longitudes_rho_USCROMS(:), latitudes_rho_USCROMS(:)], [longitudes_PREDICT, latitudes_PREDICT]);
idx_u_USCROMS_PREDICT = knnsearch([longitudes_u_USCROMS(:), latitudes_u_USCROMS(:)], [longitudes_PREDICT, latitudes_PREDICT]);
idx_v_USCROMS_PREDICT = knnsearch([longitudes_v_USCROMS(:), latitudes_v_USCROMS(:)], [longitudes_PREDICT, latitudes_PREDICT]);

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

% % Plot every single USCROMS coordinate
% h_rho = plot(longitudes_rho_USCROMS(:), latitudes_rho_USCROMS(:), 'or', 'MarkerSize', 3); % Red circles for every single rho-coordinate
% h_u = plot(longitudes_u_USCROMS(:), latitudes_u_USCROMS(:), 'ob', 'MarkerSize', 3); % Blue circles for every single u-coordinate
% h_v = plot(longitudes_v_USCROMS(:), latitudes_v_USCROMS(:), 'og', 'MarkerSize', 3); % Green circles for every single v-coordinate

% Set aspect ratio to be equal
daspect([1, 1, 1]);

% % Add legend
% legend([h_rho, h_u, h_v], {'Rho-coordinates', 'U-coordinates', 'V-coordinates'}, 'Location', 'best');

% % Set x and y limits for St John
% xlim([360-64.85, 360-64.65]); % Longitude range for St John
% ylim([18.25, 18.45]); % Latitude range for St John

% Add title and labels
% title('USCROMS Coordinates for St John, USVI');
title('USCROMS Landmask');
xlabel('Longitude');
ylabel('Latitude');


%% write to shapefile

% % Write the shapefile
% shapewrite(land, fullfile(outputPath, 'landmask_justrho.shp'));
% 
% % Write the .prj file
% fid = fopen(fullfile(outputPath, 'landmask_justrho.prj'), 'w');
% fwrite(fid, prjText);
% fclose(fid);
% 
% 
% shapewrite(land, fullfile(outputPath, 'landmask_justrho.shp'));
% 
% % Read the shapefile
% landmask_justrho = shaperead(fullfile(outputPath, 'landmask_justrho.shp'));
% 
% % Plot it
% figure;
% mapshow(landmask_justrho);
% title('Land Mask');

% Define .prj file contents for EPSG:4269 (NAD83)
prjText = ['GEOGCS["NAD83",', ...
    'DATUM["North_American_Datum_1983",', ...
    'SPHEROID["GRS 1980",6378137,298.257222101]],', ...
    'PRIMEM["Greenwich",0],', ...
    'UNIT["degree",0.0174532925199433],', ...
    'AXIS["Latitude",NORTH],', ...
    'AXIS["Longitude",EAST]]'];

% Combine all into one array of structs
land_all = [land(:); land_u(:); land_v(:)];

% Write to a single shapefile
shapewrite(land_all, fullfile(outputPath, 'landmask_merged.shp'));

% Write the .prj file
fid = fopen(fullfile(outputPath, 'landmask_merged.prj'), 'w');
fwrite(fid, prjText);
fclose(fid);

% Read the shapefile
landmask = shaperead(fullfile(outputPath, 'landmask_merged.shp'));

% Plot it
figure;
mapshow(landmask);
title('Land Mask');

%% dissolve shapefile

% Step 1: Read merged shapefile
S = shaperead(fullfile(outputPath, 'landmask_merged.shp'));

% Step 2: Convert to polyshapes
P = polyshape.empty(length(S), 0);
for k = 1:length(S)
    P(k) = polyshape(S(k).X, S(k).Y);
end

% Step 3: Union all shapes (dissolve)
P_union = P(1);
for k = 2:length(P)
    P_union = union(P_union, P(k));
end

% Step 4: Convert back to struct
[x, y] = boundary(P_union);
S_dissolved = struct( ...
    'Geometry', 'Polygon', ...
    'X', x(:)', ...
    'Y', y(:)', ...
    'ID', 1); % Optional field

% Step 5: Write to new shapefile
shapewrite(S_dissolved, fullfile(outputPath, 'landmask_dissolved.shp'));

% Write the .prj file
fid = fopen(fullfile(outputPath, 'landmask_dissolved.prj'), 'w');
fwrite(fid, prjText);
fclose(fid);

% Read the shapefile
landmask_dissolved = shaperead(fullfile(outputPath, 'landmask_dissolved.shp'));

% Plot it
figure;
mapshow(landmask_dissolved);
title('Land Mask Dissolved');
