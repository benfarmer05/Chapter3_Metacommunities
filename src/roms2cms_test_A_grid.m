%%
% 1.) make sure the w-velocity interpolation isn't going haywire. CHECK
% 2.) fix the "valid range" issues. CHECK
% 3.) create a release file. CHECK
% 4.) test this bad boy out on the CMS. CHECK
% 5.) maybe increase depth range down to 300 m? NO
% 6.) does netcdf even work anymore?? figure out by trying a job, etc. CHECK

clear;clc

scriptDrive = '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023';
writeDrive = '/Users/benja/testCMS'; %'/Volumes/UVI_Hydro_2019-2020/CMS_inputs_nozeta'; %
readDrive = '/Users/benja/testCMS'; %'/Volumes/UVI_Hydro_2019_1-of-2';
readDrive2 = '/Users/benja/testCMS'; %'/Volumes/UVI_Hydro_2019-2020/USCROMS_inputs';

%% select an output file and check if it matches correctly with the corresponding input file

%randomly select an output file
cd(writeDrive)
myDir = pwd;
files = dir(fullfile(myDir,'nest_1_*'));
randomIndex = randi(length(files)); %generate a random index between 1 and the number of files
baseFileName = files(randomIndex).name;
fullFileName = fullfile(myDir, baseFileName);
[~, randomfile, ext] = fileparts(fullFileName); %extract the file name from the full path
randomfile = [randomfile ext];

ncdisp(randomfile)

%depth level information (guessing the first level; see 'roms2cms_A_grid'
%for exact first level reading)
zlevels = [0.26 1 2 4 6 8 10 12 14 16 18 20 22 24 26 29 32 35 38 42 46 50 60 70 85 100 120 140 160 190 220 250]'; %32 layers down to 250 m depth, mid-depth hi-res

% Extract year, month, day, and hour from the file name
year = str2double(randomfile(8:11));
month = str2double(randomfile(12:13));
day = str2double(randomfile(14:15));
hour = str2double(randomfile(16:17));

% Define the number of days in each month for a non-leap year
daysInMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

% Calculate the number of days elapsed up to the given month
daysElapsed = sum(daysInMonth(1:month-1)) + (day - 1);

% Calculate the total number of hours elapsed
hoursElapsed_net = daysElapsed * 24 + hour;

% Adjust hoursElapsed for the file naming system (each file contains 3 hours)
hoursElapsed_format = floor(hoursElapsed_net / 3) * 3;
hoursequence = hoursElapsed_net - hoursElapsed_format + 1;

% Format the hour count to a zero-padded string of length 5
hourCountStr = sprintf('%05d', hoursElapsed_format);

% Construct the new file name
newFileName = ['croco_his.', hourCountStr, '.nc'];

% Display the new file name
disp(newFileName);

% Assuming the directory is specified in 'directoryPath'
directoryPath = readDrive; % Update this to your directory path
directoryPath2 = readDrive2; %since not all the files fit onto on drive

% Construct the full file path
crocofile = fullfile(directoryPath, newFileName);
crocofile2 = fullfile(directoryPath2, newFileName);

% % Display the full file path
% disp(crocofile);
% disp(crocofile2);

% Check if either of the files exists in the directory
crocoreadfile = '';
if isfile(crocofile)
    disp(['File found: ', crocofile]);
    crocoreadfile = crocofile;
elseif isfile(crocofile2)
    disp(['File found: ', crocofile2]);
    crocoreadfile = crocofile2;
else
    disp('Neither file found.');
end

%% test how the interpolated variables in the output files look

varu = ncread(randomfile, 'zu');
varv = ncread(randomfile, 'zv');
varw = ncread(randomfile, 'zw');
vart = ncread(randomfile, 'zt');
vars = ncread(randomfile, 'zs');
timer = ncread(randomfile, 'Time');

lon = ncread(randomfile, 'Longitude');
lat = ncread(randomfile, 'Latitude');

uvel = ncread(crocoreadfile, 'u');
vvel = ncread(crocoreadfile, 'v');
wvel = ncread(crocoreadfile, 'w');
salt = ncread(crocoreadfile,'salt');
temp = ncread(crocoreadfile,'temp');
uvel = uvel(:,:,:,hoursequence);
vvel = vvel(:,:,:,hoursequence);
wvel = wvel(:,:,:,hoursequence);
temp = temp(:,:,:,hoursequence);
salt = salt(:,:,:,hoursequence);

layers = size(salt,3);

uvel = uvel(:,:,layers:-1:1);
vvel = vvel(:,:,layers:-1:1);
wvel = wvel(:,:,layers:-1:1);
temp = temp(:,:,layers:-1:1);
salt = salt(:,:,layers:-1:1);

uvel(uvel==0) = NaN;
vvel(vvel==0) = NaN;
wvel(wvel==0) = NaN;
temp(temp==0) = NaN;
salt(salt==0) = NaN;

% ROMS_grid = 'sttstj300m.nc';
% ROMS_grid = 'parent2km100.nc';
ROMS_grid = '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed/Testing/vigrid.nc';
ncdisp(ROMS_grid)
lon_u_ROMS = ncread(ROMS_grid, 'lon_u') + 360;
lat_u_ROMS = ncread(ROMS_grid, 'lat_u');
lon_v_ROMS = ncread(ROMS_grid, 'lon_v') + 360;
lat_v_ROMS = ncread(ROMS_grid, 'lat_v');
lon_rho_ROMS = ncread(ROMS_grid, 'lon_rho') + 360; %w-vel, temp, salinity, zeta, and rho use rho-grid
lat_rho_ROMS = ncread(ROMS_grid, 'lat_rho'); %w-vel, temp, salinity, zeta, and rho use rho-grid

%% release points

% cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/Release';     
% % [relpnts] = readmatrix('points.xlsx');
% [relpnts] = readmatrix('points_650_none-on-land.csv');
% plot(relpnts(:,2)+360,relpnts(:,3), '*r')

%%

fill_value = single(1.2676506e+30);

%replace fill values with NaNs for plotting
varu(varu==fill_value) = NaN;
varv(varv==fill_value) = NaN; %varv(varv==1.2676506e+30) = NaN;
varw(varw==fill_value) = NaN;
vart(vart==fill_value) = NaN;
vars(vars==fill_value) = NaN;

clc
max(max(max(varu)))
min(min(min(varu)))
max(max(max(uvel)))
min(min(min(uvel)))
maxvaru = max(max(max(varu))); minvaru = min(min(min(varu)));
plotu = max([abs(maxvaru) abs(minvaru)]);

clc
max(max(max(varv)))
min(min(min(varv)))
max(max(max(vvel)))
min(min(min(vvel)))
maxvarv = max(max(max(varv))); minvarv = min(min(min(varv)));
plotv = max([abs(maxvarv) abs(minvarv)]);

clc
max(max(max(varw)))
min(min(min(varw)))
max(max(max(wvel)))
min(min(min(wvel)))
maxvarw = max(max(varw(:,:,1))); minvarw = min(min(varw(:,:,1)));
plotw = max([abs(maxvarw) abs(minvarw)]);
% plotw = plotw/5; %the actual value here much higher than the usual. maybe some issues with USCROMS' w-velocities

clc
max(max(max(vart)))
min(min(min(vart)))
max(max(max(temp)))
min(min(min(temp)))
maxvart = max(max(vart(:,:,1))); minvart = min(min(vart(:,:,1)));

clc
max(max(max(vars)))
min(min(min(vars)))
max(max(max(salt)))
min(min(min(salt)))
maxvars = max(max(vars(:,:,1))); minvars = min(min(vars(:,:,1)));

%%
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023'

lons = ncread(ROMS_grid, 'lon_rho') + 360;
lats = ncread(ROMS_grid, 'lat_rho');
mask_USCROMS = ncread(ROMS_grid, 'mask_rho');

%USCROMS landmasks
mask_u_land_ROMS = uvel(:,:,1);
mask_u_land_ROMS(isnan(mask_u_land_ROMS)) = 0; %land is 0's
mask_u_land_ROMS(mask_u_land_ROMS ~= 0) = 1; %the ocean is 1's

mask_v_land_ROMS = vvel(:,:,1);
mask_v_land_ROMS(isnan(mask_v_land_ROMS)) = 0;
mask_v_land_ROMS(mask_v_land_ROMS ~= 0) = 1;

mask_rho_land_ROMS = wvel(:,:,1);
mask_rho_land_ROMS(isnan(mask_rho_land_ROMS)) = 0;
mask_rho_land_ROMS(mask_rho_land_ROMS ~= 0) = 1;

%USCROMS seafloor masks
mask_u_seafloor_ROMS = uvel(:,:,72);
mask_u_seafloor_ROMS(isnan(mask_u_seafloor_ROMS)) = 0; %seafloor is 0's
mask_u_seafloor_ROMS(mask_u_seafloor_ROMS ~= 0) = 1; %the ocean is 1's

mask_v_seafloor_ROMS = vvel(:,:,72);
mask_v_seafloor_ROMS(isnan(mask_v_seafloor_ROMS)) = 0;
mask_v_seafloor_ROMS(mask_v_seafloor_ROMS ~= 0) = 1;

mask_rho_seafloor_ROMS = wvel(:,:,72); 
mask_rho_seafloor_ROMS(isnan(mask_rho_seafloor_ROMS)) = 0;
mask_rho_seafloor_ROMS(mask_rho_seafloor_ROMS ~= 0) = 1;

%CMS landmasks
mask_u_land_CMS = varu(:,:,1); 
mask_u_land_CMS(isnan(mask_u_land_CMS)) = 0; %land is 0's
mask_u_land_CMS(mask_u_land_CMS ~= 0) = 1; %the ocean is 1's

mask_v_land_CMS = varv(:,:,1); 
mask_v_land_CMS(isnan(mask_v_land_CMS)) = 0; 
mask_v_land_CMS(mask_v_land_CMS ~= 0) = 1; 

mask_rho_land_CMS = varw(:,:,1);
mask_rho_land_CMS(isnan(mask_rho_land_CMS)) = 0; 
mask_rho_land_CMS(mask_rho_land_CMS ~= 0) = 1;

%CMS seafloor masks
desired_depthlevel = 17;
mask_u_seafloor_CMS = varu(:,:,desired_depthlevel); 
mask_u_seafloor_CMS(isnan(mask_u_seafloor_CMS)) = 0; %seafloor is 0's
mask_u_seafloor_CMS(mask_u_seafloor_CMS ~= 0) = 1; %the ocean is 1's

mask_v_seafloor_CMS = varv(:,:,desired_depthlevel); 
mask_v_seafloor_CMS(isnan(mask_v_seafloor_CMS)) = 0;
mask_v_seafloor_CMS(mask_v_seafloor_CMS ~= 0) = 1; 

mask_rho_seafloor_CMS = varw(:,:,desired_depthlevel); 
mask_rho_seafloor_CMS(isnan(mask_rho_seafloor_CMS)) = 0; 
mask_rho_seafloor_CMS(mask_rho_seafloor_CMS ~= 0) = 1; 

depths = ncread(ROMS_grid, 'h');
depthcut1 = depths;
depthcut1(depthcut1==0) = NaN;
depthcut1(mask_rho_land_CMS == 0) = NaN;
% figure();surf(lons', lats', -depthcut1');colorbar;clim([-30 0])

% figure();imagescn(lons', lats', -depthcut1');set(gca, 'ydir', 'normal');hold on;colorbar;axis equal; clim([-15 -1]); hold on
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/Release';     
[relpnts] = readmatrix('points_650_none-on-land.csv');
% plot(relpnts(:,11) + 360,relpnts(:,12), '*r')

%%

%seafloor mask comparisons
% figure();imagescn(lon',lat',mask_rho_seafloor_CMS');cmocean thermal;axis equal;xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Rho seafloor mask');
% figure();imagescn(lon',lat',mask_u_seafloor_CMS');cmocean thermal;axis equal;xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('U seafloor mask');
% figure();imagescn(lon',lat',mask_v_seafloor_CMS');cmocean thermal;axis equal;xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('V seafloor mask');

% %map mismatch between rho-grid and u-grid
% nan_diff_mask = mask_u_seafloor_CMS==0 & mask_rho_seafloor_CMS==1;
% figure();
% imagescn(lon', lat', nan_diff_mask');
% cmocean('thermal');
% axis equal;
% % xlim([360-65.5, 360-65]);ylim([18, 18.4])
% % xlim([360-64.85, 360-64.65]);
% % ylim([18.25, 18.45]);
% hold on;
% % plot(relpnts(:,11) + 360, relpnts(:,12), '*r');
% title('NaNs in U seafloor mask but not in Rho seafloor mask');
% 
% %map mismatch between rho-grid and v-grid
% nan_diff_mask = mask_v_seafloor_CMS==0 & mask_rho_seafloor_CMS==1;
% figure();
% imagescn(lon', lat', nan_diff_mask');
% cmocean('thermal');
% axis equal;
% % xlim([360-65.5, 360-65]);ylim([18, 18.4])
% % xlim([360-64.85, 360-64.65]);
% % ylim([18.25, 18.45]);
% hold on;
% % plot(relpnts(:,11) + 360, relpnts(:,12), '*r');
% title('NaNs in V seafloor mask but not in Rho seafloor mask');

%land mask comparisons
% figure();imagescn(lon',lat',mask_rho_land_CMS');cmocean thermal;axis equal;xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Rho land mask');
% figure();imagescn(lon',lat',mask_u_land_CMS');cmocean thermal;axis equal;xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('U land mask');
% figure();imagescn(lon',lat',mask_v_land_CMS');cmocean thermal;axis equal;xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('V land mask');

% %seafloor velocity comparisons
% figure();imagescn(lon',lat',varw(:,:,desired_depthlevel)');axis equal;clim([-plotw plotw]);xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('z-level Interpolated & rho-Regridded at-depth w-velocity');
% figure();imagescn(lon',lat',varu(:,:,desired_depthlevel)');axis equal;clim([-plotu plotu]);xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('z-level Interpolated & rho-Regridded surface u-velocity');
% figure();imagescn(lon',lat',varv(:,:,desired_depthlevel)');axis equal;clim([-plotv plotv]);xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('z-level Interpolated & rho-Regridded surface v-velocity');

%u
%interpolated
% figure(1);imagescn(lon',lat',varu(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-plotu plotu]);xlim([360-65.3, 360-64.3]);ylim([17.6, 18.7]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('z-level Interpolated & rho-Regridded surface u-velocity');
% figure(1);imagescn(lon',lat',varu(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-plotu plotu]);
% figure(2);imagescn(lon',lat',varu(:,:,32)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-plotu plotu]);
%original
% figure(2);imagescn(lon_u_ROMS',lat_u_ROMS',uvel(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-plotu plotu]);xlim([360-65.3, 360-64.3]);ylim([17.6, 18.7]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Original sigma-level surface u-velocity');
% figure();imagescn(lon_u_ROMS',lat_u_ROMS',uvel(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-plotu plotu]);
% figure();imagescn(lon_u_ROMS',lat_u_ROMS',uvel(:,:,72)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-plotu plotu]);

% %STJ
% figure();imagescn(lon_rho_ROMS',lat_rho_ROMS',mask_rho_land_ROMS');axis equal;xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Rho-grid mask');
% figure();imagescn(lon_u_ROMS',lat_u_ROMS',uvel(:,:,1)');axis equal;clim([-plotu plotu]);xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Original sigma-level surface u-velocity');
% figure();imagescn(lon',lat',varu(:,:,1)');axis equal;clim([-plotu plotu]);xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('z-level Interpolated & rho-Regridded surface u-velocity');

% %STX
% figure();imagescn(lon_rho_ROMS',lat_rho_ROMS',mask');axis equal;xlim([360-64.95, 360-64.5]);ylim([17.6, 17.9]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Rho-grid mask');
% figure();imagescn(lon_u_ROMS',lat_u_ROMS',uvel(:,:,1)');axis equal;clim([-plotu plotu]);xlim([360-64.95, 360-64.5]);ylim([17.6, 17.9]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Original sigma-level surface u-velocity');
% figure();imagescn(lon',lat',varu(:,:,1)');axis equal;clim([-plotu plotu]);xlim([360-64.95, 360-64.5]);ylim([17.6, 17.9]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('z-level Interpolated & rho-Regridded surface u-velocity');

% %TEST
% [row, column] = find (varu <= -10)
% [row, column] = find (uvel >= .1)
% %TEST

%v
%interpolated
% figure();imagescn(lon',lat',varv(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-plotv plotv]);xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('z-level Interpolated & rho-Regridded surface v-velocity');
% figure();imagescn(lon',lat',varv(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-plotv plotv]);
% figure();imagescn(lon',lat',varv(:,:,32)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-plotv plotv]);
%original
% figure();imagescn(lon_v_ROMS',lat_v_ROMS',vvel(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-plotv plotv]);xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Original sigma-level surface v-velocity');
% figure();imagescn(lon_v_ROMS',lat_v_ROMS',vvel(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-plotv plotv]);
% figure();imagescn(lon_v_ROMS',lat_v_ROMS',vvel(:,:,72)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-plotv plotv]);

% %all
% figure();imagescn(lon_rho_ROMS',lat_rho_ROMS',mask');axis equal;xlim([360-65.3, 360-64]);ylim([17.2, 19.1]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Rho-grid mask');
% figure();imagescn(lon_v_ROMS',lat_v_ROMS',vvel(:,:,1)');axis equal;clim([-plotv plotv]);xlim([360-65.3, 360-64]);ylim([17.2, 19.1]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Original sigma-level surface v-velocity');
% figure();imagescn(lon',lat',varv(:,:,1)');axis equal;clim([-plotv plotv]);xlim([360-65.3, 360-64]);ylim([17.2, 19.1]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('z-level Interpolated & rho-Regridded surface v-velocity');

% %STJ
% figure();imagescn(lon_rho_ROMS',lat_rho_ROMS',mask');axis equal;xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Rho-grid mask');
% figure();imagescn(lon_v_ROMS',lat_v_ROMS',vvel(:,:,1)');axis equal;clim([-plotv plotv]);xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Original sigma-level surface v-velocity');
% figure();imagescn(lon',lat',varv(:,:,1)');axis equal;clim([-plotv plotv]);xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('z-level Interpolated & rho-Regridded surface v-velocity');

% %STX
% figure();imagescn(lon_rho_ROMS',lat_rho_ROMS',mask');axis equal;xlim([360-64.95, 360-64.5]);ylim([17.6, 17.9]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Rho-grid mask');
% figure();imagescn(lon_v_ROMS',lat_v_ROMS',vvel(:,:,1)');axis equal;clim([-plotv plotv]);xlim([360-64.95, 360-64.5]);ylim([17.6, 17.9]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Original sigma-level surface v-velocity');
% figure();imagescn(lon',lat',varv(:,:,1)');axis equal;clim([-plotv plotv]);xlim([360-64.95, 360-64.5]);ylim([17.6, 17.9]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('z-level Interpolated & rho-Regridded surface v-velocity');

%w
% %interpolated
% figure(1);imagescn(lon',lat',varw(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;%clim([-plotw plotw]);
% figure(1);imagescn(lon',lat',varw(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;axis equal;%clim([-plotw plotw]);
% figure(1);imagescn(lon',lat',varw(:,:,32)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;axis equal;%clim([-plotw plotw]);
% %original
% figure(2);imagescn(lon_rho_ROMS',lat_rho_ROMS',wvel(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;%clim([-plotw plotw]);
% figure(4);imagescn(lon_rho_ROMS',lat_rho_ROMS',wvel(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;axis equal;%clim([-plotw plotw]);
% figure(5);imagescn(lon_rho_ROMS',lat_rho_ROMS',wvel(:,:,72)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;axis equal;%clim([-plotw plotw]);

% %all
% figure();imagescn(lon_rho_ROMS',lat_rho_ROMS',mask');axis equal;xlim([360-65.3, 360-64]);ylim([17.2, 19.1]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Rho-grid mask');
% figure();imagescn(lon_rho_ROMS',lat_rho_ROMS',wvel(:,:,depthlevel)');colorbar;axis equal;clim([-plotw plotw]);xlim([360-65.3, 360-64]);ylim([17.2, 19.1]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Original sigma-level surface w-velocity');
% figure();imagescn(lon',lat',varw(:,:,depthlevel)');colorbar;axis equal;clim([-plotw plotw]);xlim([360-65.3, 360-64]);ylim([17.2, 19.1]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('z-level Interpolated & rho-Regridded surface w-velocity');

% % STJ
% %desired_depthlevel
% figure();imagescn(lon_rho_ROMS',lat_rho_ROMS',mask_rho_land_ROMS');axis equal;xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Rho-grid mask');
% figure();imagescn(lon_rho_ROMS',lat_rho_ROMS',wvel(:,:,4)');colorbar;axis equal;cmocean speed;clim([-plotw plotw]);xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Original sigma-level surface w-velocity');
% figure();imagescn(lon',lat',varw(:,:,1)');colorbar;axis equal;cmocean speed;clim([-plotw plotw]);xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('z-level Interpolated & rho-Regridded surface w-velocity');

%t
% %interpolated
% figure(1);imagescn(lon',lat',vart(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean thermal;axis equal;clim([26 27.2]);xlim([360-65.3, 360-64.3]);ylim([17.6, 18.7]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('z-level Interpolated & rho-Regridded surface temperature');
% figure(1);imagescn(lon',lat',vart(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean thermal;axis equal;%clim([minvart maxvart]);
% figure(1);imagescn(lon',lat',vart(:,:,32)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean thermal;hold on;axis equal;%clim([minvart maxvart]);
% 
% %original
% figure(2);imagescn(lon_rho_ROMS',lat_rho_ROMS',temp(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean thermal;axis equal;clim([26 27.2]);xlim([360-65.3, 360-64.3]);ylim([17.6, 18.7]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('z-level Interpolated & rho-Regridded surface temperature');
% figure(2);imagescn(lon_rho_ROMS',lat_rho_ROMS',temp(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean thermal;axis equal;%clim([minvart maxvart]);
% figure(2);imagescn(lon_rho_ROMS',lat_rho_ROMS',temp(:,:,72)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean thermal;axis equal;%clim([minvart maxvart]);

%STJ
% figure();imagescn(lon_rho_ROMS',lat_rho_ROMS',mask_rho_land_ROMS');cmocean thermal;axis equal;xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Rho-grid mask');
% figure();imagescn(lon_rho_ROMS',lat_rho_ROMS',temp(:,:,1)');cmocean thermal;axis equal;clim([26 27.2]);xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Original sigma-level surface temperature');
% figure();imagescn(lon',lat',vart(:,:,1)');cmocean thermal;axis equal;clim([26 27.2]);xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('z-level Interpolated & rho-Regridded surface temperature');
% figure();imagescn(lon',lat',mask_rho_land_CMS');cmocean thermal;axis equal;xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('z-level Interpolated & rho-Regridded surface temperature');
% figure();imagescn(lon_rho_ROMS',lat_rho_ROMS',mask_rho_seafloor_ROMS');cmocean thermal;axis equal;xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Rho-grid mask');

%s
% %interpolated
% figure(1);imagescn(lon',lat',vars(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean haline;axis equal;  %caxis([-.6 .6]);
% figure(1);imagescn(lon',lat',vars(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean haline;axis equal; %caxis([-.015 .015]);
% figure(1);imagescn(lon',lat',vars(:,:,72)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean haline;axis equal; %caxis([-.015 .015]);
% %original
% figure(2);imagescn(lon_rho_ROMS',lat_rho_ROMS',salt(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean haline;axis equal;  %caxis([-.6 .6]);
% figure(2);imagescn(lon_rho_ROMS',lat_rho_ROMS',salt(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean haline;axis equal; %caxis([-.015 .015]);
% figure(2);imagescn(lon_rho_ROMS',lat_rho_ROMS',salt(:,:,72)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean haline;axis equal; %caxis([-.015 .015]);

% %all
% figure();imagescn(lon_rho_ROMS',lat_rho_ROMS',mask');cmocean haline;axis equal;xlim([360-65.3, 360-64]);ylim([17.2, 19.1]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Rho-grid mask');
% figure();imagescn(lon_rho_ROMS',lat_rho_ROMS',salt(:,:,72)');cmocean haline;colorbar;axis equal;clim([minvars maxvars]);xlim([360-65.3, 360-64]);ylim([17.2, 19.1]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Original sigma-level surface salinity');
% figure();imagescn(lon',lat',vars(:,:,32)');cmocean haline;colorbar;axis equal;clim([minvars maxvars]);xlim([360-65.3, 360-64]);ylim([17.2, 19.1]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('z-level Interpolated & rho-Regridded surface salinity');

% %STJ
% figure();imagescn(lon_rho_ROMS',lat_rho_ROMS',mask');cmocean haline;axis equal;xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Rho-grid mask');
% figure();imagescn(lon_rho_ROMS',lat_rho_ROMS',salt(:,:,1)');cmocean haline;colorbar;axis equal;clim([minvars maxvars]);xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('Original sigma-level surface salinity');
% figure();imagescn(lon',lat',vars(:,:,1)');cmocean haline;colorbar;axis equal;clim([minvars maxvars]);xlim([360-64.85, 360-64.65]);ylim([18.25, 18.45]);hold on;plot(relpnts(:,11) + 360,relpnts(:,12), '*r');title('z-level Interpolated & rho-Regridded surface salinity');

%% Testing different variable naming conventions for very small CMS run. C-grid was running into issues (24 August 2023)

% clc
% cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed/Testing
% tester = 'smallchungus.nc';
% ncdisp(tester)
% 
% % Specify the input and output file names
% inputFile = tester; % Replace with the name of your input NetCDF4 file
% outputFile = 'bigchungus.nc'; % Replace with the desired name of your output 64-bit NetCDF file
% 
% % Open the input NetCDF4 file
% ncid = netcdf.open(inputFile, 'NOWRITE');
% 
% % Get information about the input file
% [ndims, nvars, ngatts, unlimdimid] = netcdf.inq(ncid);
% 
% % Create the output NetCDF file
% ncid_out = netcdf.create(outputFile, 'CLOBBER');
% 
% % Define dimensions in the output file (copy from the input file)
% for dimid = 0:(ndims - 1)
%     [name, len] = netcdf.inqDim(ncid, dimid);
%     netcdf.defDim(ncid_out, name, len);
% end
% 
% % Define variables in the output file (copy from the input file)
% for varid = 0:(nvars - 1)
%     [varname, xtype, dimids, natts] = netcdf.inqVar(ncid, varid);
%     netcdf.defVar(ncid_out, varname, xtype, dimids);
% 
%     % Copy variable attributes
%     for attnum = 0:(natts - 1)
%         attname = netcdf.inqAttName(ncid, varid, attnum);
%         attvalue = netcdf.getAtt(ncid, varid, attname);
%         netcdf.putAtt(ncid_out, varid, attname, attvalue);
%     end
% end
% 
% % End definitions in the output file
% netcdf.endDef(ncid_out);
% 
% % Copy variable data from input to output
% for varid = 0:(nvars - 1)
%     varname = netcdf.inqVar(ncid, varid);
%     data = netcdf.getVar(ncid, varid);
%     netcdf.putVar(ncid_out, varid, data);
% end
% 
% % Close both input and output files
% netcdf.close(ncid);
% netcdf.close(ncid_out);
% 
% bigchungus = 'bigchungus.nc';
% ncdisp(bigchungus)
% 
% %FUNCTIONING MAR HYDRO
% MARtester = 'nest_1_20140301000000.nc';
% ncdisp(MARtester)
% 
% GOMtester = 'nest_1_20100101000000.nc';
% ncdisp(GOMtester)
% 
% testerdepth = ncread(tester, 'Depth');
% MARdepth = ncread(MARtester, 'Depth');
% 
% testerw = ncread(tester, 'w');
% testerw32 = testerw(:,:,32);
% testerw1 = testerw(:,:,1);
% 
% MARw = ncread(MARtester, 'zw');
% MARw21 = MARw(:,:,21);
% MARw1 = MARw(:,:,1);
% 
% 
% netcdf.renameVar(tester,'lons','longitude')
% 
% ncwrite(filename,varname,vardata)
% 
% ncwrite(tester,'lons',lon_u_convert(:,1))
% 
% 
% 
% 
% % Open netCDF file.
% % ncid = netcdf.open('my_example.nc','NC_WRITE')
% 
% % Put file in define mode.
% netcdf.redef(ncid)
% 
% % Get name of first variable
% [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,0);
% 
% % varname
% 
% 
% avagadros_number
% 
% % Rename the variable, using a capital letter to start the name.
% netcdf.renameVar(ncid,0,'Avagadros_number')

%% testing traj files

% %25 aug 2023: there was an issue with the data type or fill value mismatch
% cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed/traj_file_troubleshoot
% trajtest = 'traj_file_01.nc';
% ncdisp(trajtest)

%% extracting the grid dimensions (grid extent)

% cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023'
% clear;clc
% 
% % ROMS_grid = 'parent2km100.nc';
% ROMS_grid = 'vigrid.nc';
% ncdisp(ROMS_grid)
% lon_u = ncread(ROMS_grid, 'lon_u');
% lat_u = ncread(ROMS_grid, 'lat_u');
% lon_v = ncread(ROMS_grid, 'lon_v');
% lat_v = ncread(ROMS_grid, 'lat_v');
% lon_rho = ncread(ROMS_grid, 'lon_rho'); %w-vel, temp, salinity, zeta, and rho use rho-grid
% lat_rho = ncread(ROMS_grid, 'lat_rho'); %w-vel, temp, salinity, zeta, and rho use rho-grid
% 
% min_lon_u = min(min(lon_u));
% min_lon_v = min(min(lon_v));
% min_lon_rho = min(min(lon_rho));
% xstart = min([min_lon_u min_lon_v min_lon_rho]); % +360
% 
% max_lon_u = max(max(lon_u));
% max_lon_v = max(max(lon_v));
% max_lon_rho = max(max(lon_rho));
% xend = max([max_lon_u max_lon_v max_lon_rho]); % +360
% 
% min_lat_u = min(min(lat_u));
% min_lat_v = min(min(lat_v));
% min_lat_rho = min(min(lat_rho));
% ystart = min([min_lat_u min_lat_v min_lat_rho]);
% 
% max_lat_u = max(max(lat_u));
% max_lat_v = max(max(lat_v));
% max_lat_rho = max(max(lat_rho));
% yend = max([max_lat_u max_lat_v max_lat_rho]);


% % % some stuff for plotting Sonaljit's old parent USCROMS grid (2021) if it
% % % ends up being useful
% % cd /Users/benja/Downloads
% % ROMS_grid = 'parent2km100.nc';
% % 
% % lon_u_ROMS = ncread(ROMS_grid, 'lon_u');
% % lat_u_ROMS = ncread(ROMS_grid, 'lat_u');
% % lon_v_ROMS = ncread(ROMS_grid, 'lon_v');
% % lat_v_ROMS = ncread(ROMS_grid, 'lat_v');
% % lon_rho_ROMS = ncread(ROMS_grid, 'lon_rho'); %w-vel, temp, salinity, zeta, and rho use rho-grid
% % lat_rho_ROMS = ncread(ROMS_grid, 'lat_rho'); %w-vel, temp, salinity, zeta, and rho use rho-grid
% % 
% % nest = 'croco_his.04863.nc';
% % ncdisp(nest)
% % 
% % uvel = ncread(nest, 'u');
% % vvel = ncread(nest, 'v');
% % wvel = ncread(nest, 'w');
% % salt = ncread(nest,'salt');
% % temp = ncread(nest,'temp');
% % uvel = uvel(:,:,:,3);
% % vvel = vvel(:,:,:,3);
% % wvel = wvel(:,:,:,3);
% % temp = temp(:,:,:,3);
% % salt = salt(:,:,:,3);
% % 
% % layers = size(salt,3);
% % 
% % uvel = uvel(:,:,layers:-1:1);
% % vvel = vvel(:,:,layers:-1:1);
% % wvel = wvel(:,:,layers:-1:1);
% % temp = temp(:,:,layers:-1:1);
% % salt = salt(:,:,layers:-1:1);
% % 
% % uvel(uvel==0) = NaN;
% % vvel(vvel==0) = NaN;
% % wvel(wvel==0) = NaN;
% % temp(temp==0) = NaN;
% % salt(salt==0) = NaN;
% % 
% % 
% % figure(2);imagescn(lon_u_ROMS',lat_u_ROMS',uvel(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;%clim([-0.5 0.5]);
