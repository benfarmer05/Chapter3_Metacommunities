%% Script to create a landmask from USCROMS hydrodynamic files
%   21 May 2025

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

% %extract USCROMS coordinates
% USCROMSgrid = 'vigrid.nc';
% % ncdisp(fullfile(dataPath, USCROMSgrid))  % Pass the full path to ncdisp
% longitudes_rho_USCROMS = ncread(fullfile(dataPath, USCROMSgrid), 'lon_rho') + 360;
% latitudes_rho_USCROMS = ncread(fullfile(dataPath, USCROMSgrid), 'lat_rho');
% longitudes_u_USCROMS = ncread(fullfile(dataPath, USCROMSgrid), 'lon_u') + 360;
% latitudes_u_USCROMS = ncread(fullfile(dataPath, USCROMSgrid), 'lat_u');
% longitudes_v_USCROMS = ncread(fullfile(dataPath, USCROMSgrid), 'lon_v') + 360;
% latitudes_v_USCROMS = ncread(fullfile(dataPath, USCROMSgrid), 'lat_v');
% mask = ncread(fullfile(dataPath, USCROMSgrid), 'mask_rho');
% 
% %construct landmasks for USCROMS u- and v-grids
% USCROMSocean = 'croco_his.04608.nc';
% % ncdisp(fullfile(dataPath, USCROMSocean))  % Pass the full path to ncdisp
% uvel_USCROMS = ncread(fullfile(dataPath, USCROMSocean), 'u');
% uvel_USCROMS = uvel_USCROMS(:,:,:,1);
% uvel_USCROMS = uvel_USCROMS(:,:,size(uvel_USCROMS, 3):-1:1); %flip so velocities are surface to bottom
% uvel_USCROMS = uvel_USCROMS(:,:,1); %extract surface layer
% mask_u_USCROMS = uvel_USCROMS;
% mask_u_USCROMS(mask_u_USCROMS ~= 0) = 1; %the ocean is 1's, land is 0's
% vvel_USCROMS = ncread(fullfile(dataPath, USCROMSocean), 'v');
% vvel_USCROMS = vvel_USCROMS(:,:,:,1);
% vvel_USCROMS = vvel_USCROMS(:,:,size(vvel_USCROMS, 3):-1:1); %flip so velocities are surface to bottom
% vvel_USCROMS = vvel_USCROMS(:,:,1); %extract surface layer
% mask_v_USCROMS = vvel_USCROMS;
% mask_v_USCROMS(mask_v_USCROMS ~= 0) = 1;

%extract random sigma-to-z converted CMS file and its coordinates
files = dir(fullfile(tempPath, 'nest_1_*'));
randomIndex = randi(length(files)); % generate a random index between 1 and the number of files
CMSname = files(randomIndex).name;
% ncdisp(fullfile(tempPath, CMSname))  % Pass the full path to ncdisp
templon_CMS = ncread(fullfile(tempPath, CMSname), 'Longitude');
templat_CMS = ncread(fullfile(tempPath, CMSname), 'Latitude');
% [longitudes_CMS, latitudes_CMS] = meshgrid(templat_CMS', templon_CMS); %convert to meshgrid format of USCROMS
[latitudes_CMS, longitudes_CMS] = meshgrid(templat_CMS', templon_CMS); %convert to meshgrid format of USCROMS

%construct landmask for CMS unified Arakawa grid (A-grid)
% NOTE - 1st layer here is the surface; 32nd layer is the seafloor
wvel_CMS = ncread(fullfile(tempPath, CMSname), 'zw');
uvel_CMS = ncread(fullfile(tempPath, CMSname), 'zu');
vvel_CMS = ncread(fullfile(tempPath, CMSname), 'zv');
wvel_CMS = wvel_CMS(:,:,1); %extract surface layer
uvel_CMS = uvel_CMS(:,:,1);
vvel_CMS = vvel_CMS(:,:,1);
mask_w_CMS = wvel_CMS;
mask_u_CMS = uvel_CMS;
mask_v_CMS = vvel_CMS;
mask_w_CMS(abs(mask_w_CMS) > 100) = 0; %the ocean is 1's, land is 0's
mask_u_CMS(abs(mask_u_CMS) > 100) = 0; %the ocean is 1's, land is 0's
mask_v_CMS(abs(mask_v_CMS) > 100) = 0; %the ocean is 1's, land is 0's
mask_w_CMS(abs(mask_w_CMS) ~= 0) = 1; %the ocean is 1's, land is 0's
mask_u_CMS(abs(mask_u_CMS) ~= 0) = 1; %the ocean is 1's, land is 0's
mask_v_CMS(abs(mask_v_CMS) ~= 0) = 1; %the ocean is 1's, land is 0's

%% THE BELOW IS FOR CMS

%% Construct the landmask polygons (boundaries)

% %plot masks
% figure('Position', [100, 100, 1000, 800]);
% subplot(2, 1, 1);
% pcolor(latitudes_CMS, longitudes_CMS, mask_u_CMS);
% shading flat;
% colorbar;
% set(gca, 'YDir', 'normal');  % 'normal' = bottom to top, 'reverse' = top to bottom
% title('Mask Data (pcolor)');
% xlabel('Longitude');
% ylabel('Latitude');
% axis tight;
% colormap jet;
% 
% figure('Position', [100, 100, 1000, 800]);
% subplot(2, 1, 1);
% pcolor(latitudes_CMS, longitudes_CMS, mask_v_CMS);
% shading flat;
% colorbar;
% set(gca, 'YDir', 'normal');  % 'normal' = bottom to top, 'reverse' = top to bottom
% title('Mask Data (pcolor)');
% xlabel('Longitude');
% ylabel('Latitude');
% axis tight;
% colormap jet;
% 
% figure('Position', [100, 100, 1000, 800]);
% subplot(2, 1, 1);
% pcolor(latitudes_CMS, longitudes_CMS, mask_w_CMS);
% shading flat;
% colorbar;
% set(gca, 'YDir', 'normal');  % 'normal' = bottom to top, 'reverse' = top to bottom
% title('Mask Data (pcolor)');
% xlabel('Longitude');
% ylabel('Latitude');
% axis tight;
% colormap jet;

%rho (w-velocity)
% Binary landmask contours
mask_plot = mask_w_CMS;
mask_plot(mask_plot == 0) = NaN;
mask_plot(~isnan(mask_plot)) = 0;
mask_plot(isnan(mask_plot)) = 1;

% Lon/lat vectors
lat = latitudes_CMS(1, :)';
lon = longitudes_CMS(:, 1);

% Make the landmask
trace = bwboundaries(mask_plot); 
for i = 1:length(trace)
    trace{i}(:, 1) = lon(trace{i}(:, 1));
    trace{i}(:, 2) = lat(trace{i}(:, 2));
end

% Create land polygons
land_CMS = struct();
for k = 1:length(trace)
    land_CMS(k).Geometry = 'Polygon';
    Xs = trace{k}(:, 1); % Might have to subtract 360 at the end
    Xs(end + 1) = NaN;
    land_CMS(k).X = Xs';
    Ys = trace{k}(:, 2);
    Ys(end + 1) = NaN;
    land_CMS(k).Y = Ys';
end

%u
% Binary landmask contours
mask_plot = mask_u_CMS;
mask_plot(mask_plot == 0) = NaN;
mask_plot(~isnan(mask_plot)) = 0;
mask_plot(isnan(mask_plot)) = 1;

% Lon/lat vectors
lat = latitudes_CMS(1, :)';
lon = longitudes_CMS(:, 1);

% Make the landmask
trace = bwboundaries(mask_plot); 
for i = 1:length(trace)
    trace{i}(:, 1) = lon(trace{i}(:, 1));
    trace{i}(:, 2) = lat(trace{i}(:, 2));
end

% Create land polygons
land_u_CMS = struct();
for k = 1:length(trace)
    land_u_CMS(k).Geometry = 'Polygon';
    Xs = trace{k}(:, 1); % Might have to subtract 360 at the end
    Xs(end + 1) = NaN;
    land_u_CMS(k).X = Xs';
    Ys = trace{k}(:, 2);
    Ys(end + 1) = NaN;
    land_u_CMS(k).Y = Ys';
end

%v
% Binary landmask contours
mask_plot = mask_v_CMS;
mask_plot(mask_plot == 0) = NaN;
mask_plot(~isnan(mask_plot)) = 0;
mask_plot(isnan(mask_plot)) = 1;

% Lon/lat vectors
lat = latitudes_CMS(1, :)';
lon = longitudes_CMS(:, 1);

% Make the landmask
trace = bwboundaries(mask_plot); 
for i = 1:length(trace)
    trace{i}(:, 1) = lon(trace{i}(:, 1));
    trace{i}(:, 2) = lat(trace{i}(:, 2));
end

% Create land polygons
land_v_CMS = struct();
for k = 1:length(trace)
    land_v_CMS(k).Geometry = 'Polygon';
    Xs = trace{k}(:, 1); % Might have to subtract 360 at the end
    Xs(end + 1) = NaN;
    land_v_CMS(k).X = Xs';
    Ys = trace{k}(:, 2);
    Ys(end + 1) = NaN;
    land_v_CMS(k).Y = Ys';
end

%% plot the u,v, and rho-grids separately, with all landmasks

% Plot the land mask separately without adding it to the legend directly
figure;
hold on;

% Plot land mask for rho-grid
for k = 1:length(land_CMS)
    plot(land_CMS(k).X, land_CMS(k).Y, 'r'); % Red lines for the land mask for rho-grid
end

% Plot land mask for u-grid
for k = 1:length(land_u_CMS)
    plot(land_u_CMS(k).X, land_u_CMS(k).Y, 'b'); % Blue lines for the land mask for u-grid
end

% Plot land mask for v-grid
for k = 1:length(land_v_CMS)
    plot(land_v_CMS(k).X, land_v_CMS(k).Y, 'g'); % Green lines for the land mask for v-grid
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
land_all = [land_CMS(:); land_u_CMS(:); land_v_CMS(:)];

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

%% write shapefile for entire hydrodynamic extent

% Assume longitudes_CMS and latitudes_CMS are 2D matrices
lon = longitudes_CMS;
lat = latitudes_CMS;

% Get domain bounds
minLon = min(lon(:));
maxLon = max(lon(:));
minLat = min(lat(:));
maxLat = max(lat(:));

% Create rectangular polygon coordinates
X = [minLon, maxLon, maxLon, minLon, minLon];
Y = [minLat, minLat, maxLat, maxLat, minLat];

% Create shapefile structure
domain_shape = struct( ...
    'Geometry', 'Polygon', ...
    'X', X, ...
    'Y', Y, ...
    'ID', 1); % Optional

% Define .prj file contents for EPSG:4269 (NAD83)
prjText = ['GEOGCS["NAD83",', ...
    'DATUM["North_American_Datum_1983",', ...
    'SPHEROID["GRS 1980",6378137,298.257222101]],', ...
    'PRIMEM["Greenwich",0],', ...
    'UNIT["degree",0.0174532925199433],', ...
    'AXIS["Latitude",NORTH],', ...
    'AXIS["Longitude",EAST]]'];

% Write shapefile
shapewrite(domain_shape, fullfile(outputPath, 'hydro_domain_extent.shp'));

% Write .prj file
fid = fopen(fullfile(outputPath, 'hydro_domain_extent.prj'), 'w');
fwrite(fid, prjText);
fclose(fid);

% Read the shapefile
hydro_domain_extent = shaperead(fullfile(outputPath, 'hydro_domain_extent.shp'));

% Plot it
figure;
mapshow(hydro_domain_extent);
title('Hydro Domain Extent');


%% THE BELOW IS FOR USCROMS - not needed if just require CMS landmask
% 
% %% Construct the landmask polygons (boundaries) for USCROMS
% 
% %rho
% % Binary landmask contours
% mask_plot = mask;
% mask_plot(mask_plot == 0) = NaN;
% mask_plot(~isnan(mask_plot)) = 0;
% mask_plot(isnan(mask_plot)) = 1;
% 
% % Lon/lat vectors
% lat = latitudes_rho_USCROMS(1, :)';
% lon = longitudes_rho_USCROMS(:, 1);
% 
% % Make the landmask
% trace = bwboundaries(mask_plot); 
% for i = 1:length(trace)
%     trace{i}(:, 1) = lon(trace{i}(:, 1));
%     trace{i}(:, 2) = lat(trace{i}(:, 2));
% end
% 
% % Create land polygons
% land_USCROMS = struct();
% for k = 1:length(trace)
%     land_USCROMS(k).Geometry = 'Polygon';
%     Xs = trace{k}(:, 1); % Might have to subtract 360 at the end
%     Xs(end + 1) = NaN;
%     land_USCROMS(k).X = Xs';
%     Ys = trace{k}(:, 2);
%     Ys(end + 1) = NaN;
%     land_USCROMS(k).Y = Ys';
% end
% 
% %u
% % Binary landmask contours
% mask_plot = mask_u_USCROMS;
% mask_plot(mask_plot == 0) = NaN;
% mask_plot(~isnan(mask_plot)) = 0;
% mask_plot(isnan(mask_plot)) = 1;
% 
% % Lon/lat vectors
% lat = latitudes_u_USCROMS(1, :)';
% lon = longitudes_u_USCROMS(:, 1);
% 
% % Make the landmask
% trace = bwboundaries(mask_plot); 
% for i = 1:length(trace)
%     trace{i}(:, 1) = lon(trace{i}(:, 1));
%     trace{i}(:, 2) = lat(trace{i}(:, 2));
% end
% 
% % Create land polygons
% land_u_USCROMS = struct();
% for k = 1:length(trace)
%     land_u_USCROMS(k).Geometry = 'Polygon';
%     Xs = trace{k}(:, 1); % Might have to subtract 360 at the end
%     Xs(end + 1) = NaN;
%     land_u_USCROMS(k).X = Xs';
%     Ys = trace{k}(:, 2);
%     Ys(end + 1) = NaN;
%     land_u_USCROMS(k).Y = Ys';
% end
% 
% %v
% % Binary landmask contours
% mask_plot = mask_v_USCROMS;
% mask_plot(mask_plot == 0) = NaN;
% mask_plot(~isnan(mask_plot)) = 0;
% mask_plot(isnan(mask_plot)) = 1;
% 
% % Lon/lat vectors
% lat = latitudes_v_USCROMS(1, :)';
% lon = longitudes_v_USCROMS(:, 1);
% 
% % Make the landmask
% trace = bwboundaries(mask_plot); 
% for i = 1:length(trace)
%     trace{i}(:, 1) = lon(trace{i}(:, 1));
%     trace{i}(:, 2) = lat(trace{i}(:, 2));
% end
% 
% % Create land polygons
% land_v_USCROMS = struct();
% for k = 1:length(trace)
%     land_v_USCROMS(k).Geometry = 'Polygon';
%     Xs = trace{k}(:, 1); % Might have to subtract 360 at the end
%     Xs(end + 1) = NaN;
%     land_v_USCROMS(k).X = Xs';
%     Ys = trace{k}(:, 2);
%     Ys(end + 1) = NaN;
%     land_v_USCROMS(k).Y = Ys';
% end
% 
% %% plot the u,v, and rho-grids separately, with all landmasks
% 
% % Plot the land mask separately without adding it to the legend directly
% figure;
% hold on;
% 
% % Plot land mask for rho-grid
% for k = 1:length(land_USCROMS)
%     plot(land_USCROMS(k).X, land_USCROMS(k).Y, 'r'); % Red lines for the land mask for rho-grid
% end
% 
% % Plot land mask for u-grid
% for k = 1:length(land_u_USCROMS)
%     plot(land_u_USCROMS(k).X, land_u_USCROMS(k).Y, 'b'); % Blue lines for the land mask for u-grid
% end
% 
% % Plot land mask for v-grid
% for k = 1:length(land_v_USCROMS)
%     plot(land_v_USCROMS(k).X, land_v_USCROMS(k).Y, 'g'); % Green lines for the land mask for v-grid
% end
% 
% % % Plot every single USCROMS coordinate
% % h_rho = plot(longitudes_rho_USCROMS(:), latitudes_rho_USCROMS(:), 'or', 'MarkerSize', 3); % Red circles for every single rho-coordinate
% % h_u = plot(longitudes_u_USCROMS(:), latitudes_u_USCROMS(:), 'ob', 'MarkerSize', 3); % Blue circles for every single u-coordinate
% % h_v = plot(longitudes_v_USCROMS(:), latitudes_v_USCROMS(:), 'og', 'MarkerSize', 3); % Green circles for every single v-coordinate
% 
% % Set aspect ratio to be equal
% daspect([1, 1, 1]);
% 
% % % Add legend
% % legend([h_rho, h_u, h_v], {'Rho-coordinates', 'U-coordinates', 'V-coordinates'}, 'Location', 'best');
% 
% % % Set x and y limits for St John
% % xlim([360-64.85, 360-64.65]); % Longitude range for St John
% % ylim([18.25, 18.45]); % Latitude range for St John
% 
% % Add title and labels
% % title('USCROMS Coordinates for St John, USVI');
% title('USCROMS Landmask');
% xlabel('Longitude');
% ylabel('Latitude');
% 
% %% write to shapefile
% 
% % % Write the shapefile
% % shapewrite(land, fullfile(outputPath, 'landmask_justrho.shp'));
% % 
% % % Write the .prj file
% % fid = fopen(fullfile(outputPath, 'landmask_justrho.prj'), 'w');
% % fwrite(fid, prjText);
% % fclose(fid);
% % 
% % 
% % shapewrite(land, fullfile(outputPath, 'landmask_justrho.shp'));
% % 
% % % Read the shapefile
% % landmask_justrho = shaperead(fullfile(outputPath, 'landmask_justrho.shp'));
% % 
% % % Plot it
% % figure;
% % mapshow(landmask_justrho);
% % title('Land Mask');
% 
% % Define .prj file contents for EPSG:4269 (NAD83)
% prjText = ['GEOGCS["NAD83",', ...
%     'DATUM["North_American_Datum_1983",', ...
%     'SPHEROID["GRS 1980",6378137,298.257222101]],', ...
%     'PRIMEM["Greenwich",0],', ...
%     'UNIT["degree",0.0174532925199433],', ...
%     'AXIS["Latitude",NORTH],', ...
%     'AXIS["Longitude",EAST]]'];
% 
% % Combine all into one array of structs
% land_all = [land(:); land_u(:); land_v(:)];
% 
% % Write to a single shapefile
% shapewrite(land_all, fullfile(outputPath, 'landmask_merged.shp'));
% 
% % Write the .prj file
% fid = fopen(fullfile(outputPath, 'landmask_merged.prj'), 'w');
% fwrite(fid, prjText);
% fclose(fid);
% 
% % Read the shapefile
% landmask = shaperead(fullfile(outputPath, 'landmask_merged.shp'));
% 
% % Plot it
% figure;
% mapshow(landmask);
% title('Land Mask');
% 
% %% dissolve shapefile
% 
% % Step 1: Read merged shapefile
% S = shaperead(fullfile(outputPath, 'landmask_merged.shp'));
% 
% % Step 2: Convert to polyshapes
% P = polyshape.empty(length(S), 0);
% for k = 1:length(S)
%     P(k) = polyshape(S(k).X, S(k).Y);
% end
% 
% % Step 3: Union all shapes (dissolve)
% P_union = P(1);
% for k = 2:length(P)
%     P_union = union(P_union, P(k));
% end
% 
% % Step 4: Convert back to struct
% [x, y] = boundary(P_union);
% S_dissolved = struct( ...
%     'Geometry', 'Polygon', ...
%     'X', x(:)', ...
%     'Y', y(:)', ...
%     'ID', 1); % Optional field
% 
% % Step 5: Write to new shapefile
% shapewrite(S_dissolved, fullfile(outputPath, 'landmask_dissolved.shp'));
% 
% % Write the .prj file
% fid = fopen(fullfile(outputPath, 'landmask_dissolved.prj'), 'w');
% fwrite(fid, prjText);
% fclose(fid);
% 
% % Read the shapefile
% landmask_dissolved = shaperead(fullfile(outputPath, 'landmask_dissolved.shp'));
% 
% % Plot it
% figure;
% mapshow(landmask_dissolved);
% title('Land Mask Dissolved');
