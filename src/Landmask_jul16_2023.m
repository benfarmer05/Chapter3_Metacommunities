%% Make landmasks of nests
clear;clc
cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities

% some info about the real spatial extent of this landmask that took way too
%   long to figure out manually, using https://twcc.fr/en/#donate and
%   converting the max extents of this landmask (found in QGIS) from GPS
%   degrees (WGS84) --> NAD 1983 HARN Florida East, since I couldn't really
%   find one on their website for the VI/Puerto Rico.
%
% - Converted extents on the x-axis are 1865369.18 to 2026765.85 m
%       This comes out to 161,396.67 m width
% - Converted extents on the y-axis are -735859.42 to -475801.937 m
%       This comes out to 260,057.483 m height
% - After dividing the width by 536 (number of etaRho grid cells) and the
%       height by 716 (number of xiRho grid cells)...
% - Horizontal resolution is ~301.1131902985 m
% - Vertical resolution is ~363.7 m
%
% So, all in all, we could say the "resolution" is ~332.160982719 m. To
% have a habitat map that roughly doubles this, we'd be at ~664.321965438
% m. The conservative thing may be to round up to 700 m; slightly "riskier"
% decision might be to go with 650 m. Has to be a multiple of 50 m to match
% up well with the NCRMP grid. I personally think that 650 m is more than
% fair.

%% Smol
% childgrid = 'sttstj300m.nc';
childgrid = 'vigrid.nc';

%lon/lat and mask matrices
lons = ncread(childgrid, 'lon_rho');
lats = ncread(childgrid, 'lat_rho');
smolmask = ncread(childgrid, 'mask_rho');

%lon/lat vectors
lat = lats(1,:)';
lon = lons(:,1);

%% Domain boundary
trace = bwboundaries(smolmask); %trace region boundaries in binary image
% bwboundaries(BW) traces the exterior boundaries of objects, as well as boundaries of holes inside these objects, in the binary image BW. 
% bwboundaries also descends into the outermost objects (parents) and traces their children (objects completely enclosed by the parents). 
% Returns B, a cell array of boundary pixel locations.

for i = 1:length(trace)
    trace{i}(:,1) = lon(trace{i}(:,1));
    trace{i}(:,2) = lat(trace{i}(:,2));
end
 
for k = 1:length(trace)
    boundary = trace{k};
    plot(boundary(:,1),boundary(:,2),'r'); hold on
end

domain = struct();

for k = 1:length(trace)
    domain(k).Geometry = 'Polygon';
    Xs = trace{k}(:,1); %might have to subtract 360 at the end
    Xs(end + 1) = NaN;
    domain(k).X = Xs';
    Ys = trace{k}(:,2);
    Ys(end + 1) = NaN;
    domain(k).Y = Ys';
end

% cd '/Users/benja/MATLAB-DRIVE/USVI_release';
% load 'ReleaseFile.txt'
mapshow(domain)
% hold on;
% plot(ReleaseFile(:,2),ReleaseFile(:,3),'*')

shapewrite(domain,'smol_domain')

%% Landmask

%swap 1s and 0s so land becomes 1
smolmask(smolmask == 0) = NaN;
smolmask(~isnan(smolmask)) = 0;
smolmask(isnan(smolmask)) = 1;

%make the landmask
trace = bwboundaries(smolmask); 
for i = 1:length(trace)
    trace{i}(:,1) = lon(trace{i}(:,1));
    trace{i}(:,2) = lat(trace{i}(:,2));
end
 
for k = 1:length(trace)
    boundary = trace{k};
    plot(boundary(:,1),boundary(:,2),'r'); hold on
end

land = struct();
for k = 1:length(trace)
    land(k).Geometry = 'Polygon';
    Xs = trace{k}(:,1); %might have to subtract 360 at the end
    Xs(end + 1) = NaN;
    land(k).X = Xs';
    Ys = trace{k}(:,2);
    Ys(end + 1) = NaN;
    land(k).Y = Ys';
end

mapshow(land)
shapewrite(land,'smol_land')

% %% Big
% parentgrid = 'parent2km100.nc';
% 
% %lon/lat and mask matrices
% lons = ncread(parentgrid, 'lon_rho');
% lats = ncread(parentgrid, 'lat_rho');
% bigmask = ncread(parentgrid, 'mask_rho');
% 
% %lon/lat vectors
% lat = lats(1,:)';
% lon = lons(:,1);
% 
% %% Domain boundary
% trace = bwboundaries(bigmask); %trace region boundaries in binary image
% % bwboundaries(BW) traces the exterior boundaries of objects, as well as boundaries of holes inside these objects, in the binary image BW. 
% % bwboundaries also descends into the outermost objects (parents) and traces their children (objects completely enclosed by the parents). 
% % Returns B, a cell array of boundary pixel locations.
% 
% for i = 1:length(trace)
%     trace{i}(:,1) = lon(trace{i}(:,1));
%     trace{i}(:,2) = lat(trace{i}(:,2));
% end
% 
% for k = 1:length(trace)
%     boundary = trace{k};
%     plot(boundary(:,1),boundary(:,2),'r'); hold on
% end
% 
% domain = struct();
% 
% for k = 1:length(trace)
%     domain(k).Geometry = 'Polygon';
%     Xs = trace{k}(:,1); %might have to subtract 360 at the end
%     Xs(end + 1) = NaN;
%     domain(k).X = Xs';
%     Ys = trace{k}(:,2);
%     Ys(end + 1) = NaN;
%     domain(k).Y = Ys';
% end
% 
% mapshow(domain)
% 
% shapewrite(domain,'big_domain')
% 
% %% Landmask
% 
% %swap 1s and 0s so land becomes 1
% bigmask(bigmask == 0) = NaN;
% bigmask(~isnan(bigmask)) = 0;
% bigmask(isnan(bigmask)) = 1;
% 
% %make the landmask
% trace = bwboundaries(bigmask); 
% for i = 1:length(trace)
%     trace{i}(:,1) = lon(trace{i}(:,1));
%     trace{i}(:,2) = lat(trace{i}(:,2));
% end
% 
% for k = 1:length(trace)
%     boundary = trace{k};
%     plot(boundary(:,1),boundary(:,2),'r'); hold on
% end
% 
% land = struct();
% for k = 1:length(trace)
%     land(k).Geometry = 'Polygon';
%     Xs = trace{k}(:,1); %might have to subtract 360 at the end
%     Xs(end + 1) = NaN;
%     land(k).X = Xs';
%     Ys = trace{k}(:,2);
%     Ys(end + 1) = NaN;
%     land(k).Y = Ys';
% end
% 
% mapshow(land)
% shapewrite(land,'big_land')
% 