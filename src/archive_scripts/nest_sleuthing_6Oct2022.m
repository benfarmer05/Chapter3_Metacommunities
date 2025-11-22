%% notes

% 10 Nov 2021 (from Sonaljit Mukherjee, UVI):
% The model takes a few days to stabilize the flow. So, for your analysis,
% start from January 20th, 00:00:00. That would be croco_his.00480. The
% number on each file represents the number of hours from January 1,
% 00:00:00.

%% some notes 17 feb 2022

% 1.) Need to get Longitude, Latitude, Depth, and Time dimensions for a new
%format of .nc file. This is what the CMS can use. (update: got lon and
%lat moved from grid to hydro file. this is a good start)

% 2.) To do that, figure out what is missing in current .nc files (update:
% doesn't look like many or any NAs in Sonaljit's ROMS output)

% 3.) Then what - fill those missing data in. Interpolate? (update: maybe
% N/A

% 4.) Next steps: Time, but the velocity, temp, and salinity components
% are all there and may be fine or require minor editing hopefully.
% the real hard part - after time is sorted out, will need to load grid
% info (coordinate & depth arrays) into each of the new nydro files. this
% may require a for loop or could hopefully be done with a lovely MATLAB
% command.

%% notes for Dan meeting Friday 18 Feb 2022
% 1.) Depth: what stretching function do we want to use? see here:
%   https://www.myroms.org/wiki/Vertical_S-coordinate
%   https://www.myroms.org/forum/viewtopic.php?t=1170
% 2.) Time: We have three time steps per file...how many is this per day?
%   - And do we really want to day-average anyways? It might be nice to see
%   how coastal forcing from subdaily events like tide, storm surge, and
%   hurricanes play into connectivity. Dorian is an example.
% 3.) Format: How to rename all the files in the format we want for CMS
% 4.) Temperature: I need to know the best format for analyzing heat
%   content in 2019

%% notes from Dan meeting Friday 18 Feb 2022
% 1.) Time: we have hourly output data (each file is 3 hours apart and has
% 3 hours in it
% 2.) Loop: 3 parts?
% 3.) CMS user guide: do we have to interpolate horizontally ?


%% deprecated 
% %% Load in 2019 UVI .nc files
% %restart file
% rst = 'croco_rst.nc';ncdisp(rst)
% 
% %forcing 2019 - idealized?
% forcing = 'forcing-2019.nc';ncdisp(forcing)

%% Grid file sleuthing

% cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/VI_ROMS
% parentgrid = 'parent2km100.nc';ncdisp(parentgrid)

cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023'
% parentgrid = 'croco_his.00000.nc';
parentgrid = 'vigrid.nc';

%lon/lat and depth matrices
lons = ncread(parentgrid, 'lon_rho');
lats = ncread(parentgrid, 'lat_rho');
mask = ncread(parentgrid, 'mask_rho');

%depthbins represents fractions applied to the max depth scalar from ocean bottom to surface
depthbins = ncread(parentgrid, 's_rho'); %need to convert this to actual depths, to depth of interest
theta_s = ncread(parentgrid,'theta_s'); %S-coordinate surface control parameter (theta)
theta_b = ncread(parentgrid,'theta_b'); %S-coordinate bottom control parameter (b)
Tcline = ncread(parentgrid,'Tcline'); %S-coordinate surface/bottom layer width
hc = ncread(parentgrid,'hc'); %S-coordinate parameter, critical depth
% z = clc(depths,hc,theta_s,theta_b,length(depthbins));

%% Combine variables from grid and hydro file (child)

% childgrid = 'sttstj300m.nc';
% childhydro = 'nest_his.02604.nc';

childgrid = parentgrid;
childhydro = 'croco_his.00000.nc';

%lon/lat and depth matrices
lons = ncread(childgrid, 'lon_rho');
lats = ncread(childgrid, 'lat_rho');
mask = ncread(childgrid, 'mask_rho');

%bathymetry
depths = ncread(childgrid, 'h');
depthcut1 = depths;
depthcut1(depthcut1==0) = NaN;
% depthcut(depthcut > 40) = NaN;
depthcut1(mask == 0) = NaN;
% figure();surf(lons', lats', -depthcut');colorbar;%caxis([-160 0])
% title 'Bathymetry in the Virgin Islands';hold on
% cb = colorbar;ylabel(cb,'Depth (m)')
imagescn(lons', lats', -depthcut1');set(gca, 'ydir', 'normal');hold on;colorbar;clim([-50 -1])


%% figure out the grid resolution

% Specify the path to your ROMS .nc file
ncFilePath = childhydro;

% Read the lon_rho and lat_rho variables using ncread
lonRho = ncread(ncFilePath, 'lon_rho');
latRho = ncread(ncFilePath, 'lat_rho');

% Get the dimensions of the grid
etaRhoDim = size(lonRho, 1);
xiRhoDim = size(lonRho, 2);

% Define the UTM zone parameters for NAD83 / UTM zone 20N
utmZone = 20;
utmHemisphere = 'N';

% Define the projection parameters
projParams = projcrs('Authority', 'utm', 'Authority', utmZone, 'Authority', 'WGS84', 'Authority', 'WGS84', 'Authority', utmHemisphere);

% Convert the latitude and longitude values to projected coordinates (meters)
[xRho, yRho] = projfwd(projParams, latRho, lonRho);

%% release points
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/Release';     
[relpnts] = readmatrix('points.xlsx');
plot(relpnts(:,2),relpnts(:,3), '*r')


%%
cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023

%lon/lat and depth vectors
lat = lats(1,:)';
lon = lons(:,1);

%x-velocity of current?
u = ncread(childhydro,'u');
u(u==0) = NaN;

figure(2);subplot(1,2,1)
imagescn(lon', lat', u(:,:,1,1)');set(gca, 'ydir', 'normal');colorbar;axis square
pcolor(lons',lats',mask')
subplot(1,2,2)
imagesc(u(:,:,1,1)');set(gca,'ydir','normal');colorbar;axis square %demonstrates variables were extracted correctly

% v = ncread(childhydro,'v');
% w = ncread(childhydro,'w');
% temp = ncread(childhydro, 'temp'); %sea water potential temp (C)
% 
% 
% temp(temp==0) = NaN;
% subplot(1,2,1);imagesc(temp(:,:,1,1)'); set(gca, 'ydir', 'normal'); axis square; colorbar
% caxis([0 5])
% 
depths_temp = depths;
depths_temp(depths_temp > 100) = 100;
surf(lons', lats', -depths_temp')

% %% Make new .nc file / edit one
% % nccreate('nesttest.nc','Longitude'); 
% % nccreate('nesttest.nc','Latitude')
% % ncdisp('nesttest.nc')
% % 
% % % Create a netCDF file.
% % ncid = netcdf.create('foo.nc','NC_NOCLOBBER')
% % 
% % % Define a dimension.
% % lat_dimID = netcdf.defDim(ncid,'latitude',360);
% % 
% % % Define an unlimited dimension.
% % long_dimID = netcdf.defDim(ncid,'longitude',...
% % 		netcdf.getConstant('NC_UNLIMITED'));
% % 
% % %Create a new two-dimensional variable named peaks in a classic (NetCDF 3)
% % %format file named myncclassic.nc. Use the 'Dimensions' name-value pair
% % %argument to specify the names and lengths of the two dimensions. Use the
% % %'Format' name-value pair argument to specify the file format.
% % nccreate('myncclassic.nc','peaks',...
% %           'Dimensions',{'r',200,'c',200},...
% %           'Format','classic');
% %       
% % myncclassic = 'myncclassic.nc';ncdisp(myncclassic);
% % 
% % %Write data to the variable.
% % ncwrite(myncclassic,'peaks',peaks(200));
% 
% %append variables from grid file into hydro
% % this will permanently add variables to the original .nc file in the
% % MATLAB directory. be very careful with this kind of command. Only need
% % to/can run once
% % nccreate(childhydro,'lats','Dimensions',{'lat coordinate',182},'Format','64bit')
% % nccreate(childhydro,'lons','Dimensions',{'lon coordinate',314},'Format','64bit')
% % nccreate(childhydro,'depths','Dimensions',{'lons',314,'lats',182},'Format','64bit')
% 
% %write in coordinate/depth data from the grid to the hydro file
% ncwrite(childhydro, 'depths',depths)
% ncwrite(childhydro,'lats',lat) %see lines just above if this doesn't run in your env
% ncwrite(childhydro,'lons',lon)
% 
% newlats = ncread(childhydro, 'lats');
% newlons = ncread(childhydro, 'lons');
% newlat = newlats(1,:)';
% newlon = newlons(:,1);
% 
% %plot with new variables!
% figure(3)%;subplot(1,2,1)
% imagesc(newlon, newlat', u(:,:,1,1)');set(gca, 'ydir', 'normal');colorbar;axis square %plotting 1st sigma, 1st time step
% 
% depthcut1 = depths;
% depthcut1(depthcut1==0) = NaN;
% depthcut1(depthcut1 > 100) = 100;
% figure(1)%;subplot(1,2,1)
% surf(lons', lats', -depthcut1');colorbar

%% New grid from Sonaljit (15 Mar 2023)
clear;clc
cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/VI_ROMS
newgrid = 'bathy_2.nc'; %ncdisp(newgrid)
oldgrid = 'vigrid.nc';

%lon/lat and depth matrices
lons = ncread(newgrid, 'lon_rho');
lats = ncread(newgrid, 'lat_rho');
lon = lons(:,1);
lat = lats(1,:)';
mask = ncread(newgrid, 'mask_rho');
% depths1 = ncread(newgrid, 'hgauss1');
% depths2 = ncread(newgrid, 'hgauss2');
% depths3 = ncread(newgrid, 'hgauss3');
% depthscroix = ncread(newgrid, 'hcroix');
depthsnew = ncread(newgrid, 'hgauss');

%lon/lat and depth matrices for old grid (31 Aug 2022)
lonsold = ncread(oldgrid, 'lon_rho');
latsold = ncread(oldgrid, 'lat_rho');
lonold = latsold(:,1);
latold = lats(1,:)';
maskold = ncread(oldgrid, 'mask_rho');
depthsold = ncread(oldgrid, 'h');

%bathymetry
% depthcut1 = depths1;
% depthcut1(depthcut1==0) = NaN;
% % depthcut(depthcut > 40) = NaN;
% depthcut1(mask == 0) = NaN;
% % figure();surf(lons', lats', -depthcut');colorbar;%caxis([-160 0])
% % title 'Bathymetry in the Virgin Islands';hold on
% % cb = colorbar;ylabel(cb,'Depth (m)')
% % imagescn(lons', lats', -depthcut1');set(gca, 'ydir', 'normal');colorbar;clim([-50 -1]);axis equal
% % hold on
% 
% depthcut2 = depths2;
% depthcut2(depthcut2==0) = NaN;
% depthcut2(mask == 0) = NaN;
% % imagescn(lons', lats', -depthcut2');set(gca, 'ydir', 'normal');colorbar;clim([-50 -1]);axis equal
% % hold on
% 
% depthcut3 = depths3;
% depthcut3(depthcut3==0) = NaN;
% depthcut3(mask == 0) = NaN;
% % imagescn(lons', lats', -depthcut3');set(gca, 'ydir', 'normal');colorbar;clim([-50 -1]);axis equal
% % hold on
% 
% depthcutcroix = depthscroix;
% depthcutcroix(depthcutcroix==0) = NaN;
% depthcutcroix(mask == 0) = NaN;

depthcutnew = depthsnew;
depthcutnew(depthcutnew==0) = NaN;
depthcutnew(mask == 0) = NaN;

%bathymetry for old grid
depthcutold = depthsold;
depthcutold(depthcutold==0) = NaN;
depthcutold(maskold == 0) = NaN;

% figure();imagescn(lons', lats', -depthcut1'); title('Haidvogel: 0.22');set(gca, 'ydir', 'normal');colorbar;clim([-70 -1]);axis equal
% figure();imagescn(lons', lats', -depthcut2'); title('Haidvogel: 0.2?');set(gca, 'ydir', 'normal');colorbar;clim([-50 -1]);axis equal
% figure();imagescn(lons', lats', -depthcut3'); title('Haidvogel: 0.?');set(gca, 'ydir', 'normal');colorbar;clim([-50 -1]);axis equal
% figure();imagescn(lons', lats', -depthcutcroix'); title('Haidvogel: 0.?');set(gca, 'ydir', 'normal');colorbar;clim([-50 -1]);axis equal
figure();imagescn(lons', lats', -depthcutnew'); title('Haidvogel: 0.?');set(gca, 'ydir', 'normal');colorbar;clim([-50 -1]);axis equal

figure();imagescn(lonsold', latsold', -depthcutold');title('Old grid (sort of fixed but missing our more high-res bathy for STT/STJ)');set(gca, 'ydir', 'normal');colorbar;clim([-50 -1]);axis equal

%% release points
% clear;clc

cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/VI_ROMS/Release';     
[relpnts] = readmatrix('points_with_depth.xlsx');

relpnts = relpnts(relpnts(:,4)<0, :);
% relpnts = relpnts(relpnts(:,4)>-100, :);
relpnts = relpnts(relpnts(:,4)>-1.960081909179680e+03, :);

%bathymetry plot (using bathy from Q)
rel_lons = relpnts(:,2);
rel_lats = relpnts(:,3);
rel_depths = relpnts(:,4);
% bathdepthssort = sort(bathdepths);

% figure();imagescn(lons', lats', -depthcut1'); title('Updated grid (Haidvogel: 0.22)');set(gca, 'ydir', 'normal');colorbar;clim([-70 -1]);axis equal;hold on;plot(relpnts(:,2),relpnts(:,3), '*r');
% figure();imagescn(lons', lats', -depthcut2'); title('Updated grid (Haidvogel: 0.2?)');set(gca, 'ydir', 'normal');colorbar;clim([-50 -1]);axis equal; hold on;plot(relpnts(:,2),relpnts(:,3), '*r')
% figure();imagescn(lons', lats', -depthcut3'); title('Updated grid (Haidvogel: 0.?)');set(gca, 'ydir', 'normal');colorbar;clim([-50 -1]);axis equal; hold on;plot(relpnts(:,2),relpnts(:,3), '*r')
% figure();imagescn(lons', lats', -depthcutcroix'); title('Updated grid (Haidvogel: 0.?)');set(gca, 'ydir', 'normal');colorbar;clim([-50 -1]);axis equal; hold on;plot(relpnts(:,2),relpnts(:,3), '*r')
figure();imagescn(lons', lats', -depthcutnew'); title('Updated grid (Haidvogel: 0.?)');set(gca, 'ydir', 'normal');colorbar;clim([-50 -1]);axis equal; hold on;plot(relpnts(:,2),relpnts(:,3), '*r')

figure();imagescn(lonsold', latsold', -depthcutold');title('Old grid (sort of fixed but missing our more high-res bathy for STT/STJ)');set(gca, 'ydir', 'normal');colorbar;clim([-50 -1]);axis equal;hold on;plot(relpnts(:,2),relpnts(:,3), '*r')
    
%% Compare smoothed vs realistic release point depths
clc

%interpolating smoothed ROMS model bathymetry ('depths') at release point
%coordinates ('rel_lats' and 'rel_lons'
ROMSdepths_at_relpnts = griddata(lat,lon,depths,rel_lats,rel_lons); 
ROMSdepths_at_relpnts = ROMSdepths_at_relpnts*-1;
% relpnts(:,5) = ROMSdepths_at_relpnts; %add interpolated ROMS smoothed depths to tabular data

excessdepths_proportion = sum(ROMSdepths_at_relpnts<-50)/sum(ROMSdepths_at_relpnts>-50)*100;

figure
scatter(relpnts(:,4), ROMSdepths_at_relpnts); xlabel('Realistic depth'); ylabel('Smoothed depth')
h = lsline;
set(h(1),'color','r', 'LineWidth', 2)

figure
scatter(relpnts(:,4), ROMSdepths_at_relpnts); xlabel('Realistic depth'); ylabel('Smoothed depth'); ylim([-50 0]);title('Smoothed depths restricted to 0-50 m')
h = lsline;
set(h(1),'color','r', 'LineWidth', 2)

figure
scatter(-relpnts(:,4), -ROMSdepths_at_relpnts); xlabel('Realistic depth'); ylabel('Smoothed depth');title('Smoothed depths restricted to 20-50 m')
h = refline(1,0);
axis equal
axis([0 50 20 100])
set(h(1),'color','r', 'LineWidth', 2)

%% 3D map of real vs modeled depths

% Northern USVI Deep
figure
meshc(lons', lats', -depthcut1');colorbar;caxis([-350 -1]);xlim([-65.2 -64.1]);ylim([18.1 18.8]);zlim([-350 -1]);cmocean('-deep');hold on
plot3(rel_lons, rel_lats, rel_depths, '.'); hold on
plot3(rel_lons, rel_lats, ROMSdepths_at_relpnts, '*')
grid on
view(35,25)
title('Northern USVI Deep')

% Northern USVI to 50 m
figure
meshc(lons', lats', -depthcut1');colorbar;caxis([-50 -1]);xlim([-65.2 -64.1]);ylim([18.1 18.8]);zlim([-50 -1]);cmocean('-deep');hold on
plot3(rel_lons, rel_lats, rel_depths, '.'); hold on
plot3(rel_lons, rel_lats, ROMSdepths_at_relpnts, '*')
grid on
view(35,25)
title('Northern USVI to 50 m')

% St Croix Deep
figure
meshc(lons', lats', -depthcut1');colorbar;caxis([-350 -1]);xlim([-65 -64.4]);ylim([17.4 17.9]);zlim([-350 -1]);cmocean('-deep');hold on
plot3(rel_lons, rel_lats, rel_depths, '.'); hold on
plot3(rel_lons, rel_lats, ROMSdepths_at_relpnts, '*')
grid on
view(35,25)
title('St Croix Deep')

% St Croix to 50 m
figure
meshc(lons', lats', -depthcut1');colorbar;caxis([-50 -1]);xlim([-65 -64.4]);ylim([17.4 17.9]);zlim([-50 -1]);cmocean('-deep');hold on
plot3(rel_lons, rel_lats, rel_depths, '.'); hold on
plot3(rel_lons, rel_lats, ROMSdepths_at_relpnts, '*')
grid on
view(35,25)
title('St Croix to 50 m')

%% Creating merged bathymetry map from VI Shapes for Sonaljit
cd /Users/benja/Documents/Carib_Habitat/Bathy_merged_BF/

nu_bath = 'final_merge_6Oct2022_zeros-for-land.nc'; %came from a Python script
nu2_bath = 'final_merge_Sonaljit.nc'; %came straight from QGIS Translate function (possibly very similar)

ncdisp(nu_bath)
ncdisp(nu2_bath)
x_range = ncread(nu2_bath, 'x_range');
y_range = ncread(nu2_bath, 'y_range');
z_range = ncread(nu2_bath, 'z_range');
spacing = ncread(nu2_bath, 'spacing');
dimension = ncread(nu2_bath, 'dimension');

% I HAVE NO IDEA WHAT'S GOING ON LOL

% [Z,R] = readgeoraster('final_merge_6Oct2022_zeros-for-land.tif');
% [lat,lon] = geographicGrid(R);

%https://www.mathworks.com/matlabcentral/answers/38255-display-usgs-dem-using-geotiffread-and-mapshow
[A,R] = readgeoraster('final_merge_6Oct2022_zeros-for-land.tif');
mapshow(A,R,'DisplayType','surface')
A=double(A);
A(A <= -3.40282346638528e+38) = NaN;
mapshow(A,R)

%this one might actually be getting somewhere, but Puerto Rico is all
%messed up :(
[geotiff1, R] = geotiffread('final_merge_6Oct2022_zeros-for-land.tif');
geotiff1(geotiff1 <= -3.40282346638528e+38) = NaN;
imagesc(geotiff1);




%lon/lat and depth matrices
x = ncread(nu_bath, 'x');
y = ncread(nu_bath, 'y');
Band1 = ncread(nu_bath, 'Band1');
Band1(Band1 <= -3.40282346638528e+38) = NaN;
land = Band1(Band1 == 0);
figure();imagescn(Band1');set(gca, 'ydir', 'normal');colorbar;clim([-80 0])



%Chad Greene: https://www.mathworks.com/matlabcentral/answers/348056-how-to-use-qgis-geotiff-in-matlab
[lon,lat] = meshgrid(-180:0.25:180,90:-0.25,-90);
Z = imread('final_merge_6Oct2022_NAs-for-land.tif');
zi = interp2(lon,lat,Z,loni,lati);

%Working with TIFFs in QGIS - https://www.youtube.com/watch?v=VlUg1lQJhms





%% Other stuff

% xLon = linspace(min(rel_lons), max(rel_lons), 1E+3);
% yLat = linspace(min(rel_lats), max(rel_lats), 1E+3);
% [X,Y] = meshgrid(xLon, yLat);
% zDep = griddata(rel_lons, rel_lats, ROMSdepths_at_relpnts, X, Y);
% figure
% mesh(X, Y, zDep)
% grid on
% view(35,25)
% title('Mesh Plot')
% 
% figure
% mesh(X, Y, zDep)
% hold on
% contour3(X, Y, zDep, 20, 'k', 'LineWidth',2)
% hold off
% grid on
% view(35,25)
% title('Mesh Plot With Contours')





