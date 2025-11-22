%% Script to convert netCDF file output file from ROMS to a netCDF file that CMS can use
% 22 May 2024

tic
clear;clc

%% extract and define global variables
scriptDrive = '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023';
writeDrive = '/Users/benja/testCMS';%'/Volumes/UVI_Hydro_2019-2020/hydro_inputs/error'; %
readDrive = '/Users/benja/testCMS'; %'/Volumes/UVI_Hydro_2019-2020/hydro_inputs/error'; %

cd(scriptDrive)

ROMS_grid = 'vigrid.nc';
ncdisp(ROMS_grid)

theta_s = ncread(ROMS_grid, 'theta_s');
theta_b = ncread(ROMS_grid, 'theta_b');
hc = ncread(ROMS_grid, 'hc');
h = ncread(ROMS_grid, 'h');
lon_u = ncread(ROMS_grid, 'lon_u');
lat_u = ncread(ROMS_grid, 'lat_u');
lon_v = ncread(ROMS_grid, 'lon_v');
lat_v = ncread(ROMS_grid, 'lat_v');
lon_rho = ncread(ROMS_grid, 'lon_rho'); %w-vel, temp, salinity, zeta, and rho use rho-grid
lat_rho = ncread(ROMS_grid, 'lat_rho'); %w-vel, temp, salinity, zeta, and rho use rho-grid

%make longitudes positive & specify ranges for attributes
lon_u_convert = 360-abs(lon_u);
lon_u_range = [min(lon_u_convert(:,1)) max(lon_u_convert(:,1))];
lon_v_rho_convert = 360-abs(lon_v); %v & rho - they are on same grid for longitudes
lon_v_rho_range = [min(lon_v_rho_convert(:,1)) max(lon_v_rho_convert(:,1))];
lat_u_rho_range = [min(lat_u(1,:)) max(lat_u(1,:))]; %u & rho - they are on same grid for latitudes
lat_v_range = [min(lat_v(1,:)) max(lat_v(1,:))];

%chop off longitudes and latitudes
lon_v_rho_convert = lon_v_rho_convert(1:end-1, :);
lat_rho = lat_rho(1:end-1, 1:end-1);

%% calculate z-level depths from sigma-layers

%extract number of sigma layers from u-velocity in the first history file
% (this is arbitrary; the number of depth layers is identical across
% u,v,w,t,s in each file in the simulation)
cd(readDrive)
myDir = pwd;
files = dir(fullfile(myDir,'croco_his*'));
baseFileName = files(1).name; %name of the file itself
fullFileName = fullfile(myDir, baseFileName); %path to the file name
placeholder = ncread(fullFileName, 'u');
layers = size(placeholder,3);

%ensure that free surface zeta is ignored, as it can cause issues with
% 'set_depth'.
%
%  NOTE - reached out to Sonaljit 15 May 2024 to confirm whether this is
%  an impactful choice for the model
zeta = zeros(size(h));

%convert from sigma-layer to z-level depths at each lat/lon on each grid
% system
%
% NOTE - the rho-grid is here assumed as the ideal grid for w-velocity,
% temperature, and salinity. However, ROMS does have a w-grid and
% 'set_depth' documentation seems to relate it with w-velocity. That said,
% USCROMS does not appear to have been constructed in that way, given that
% w-velocity has a depth dimension of 's_rho' just like u,v,t,s. Rather
% than a depth dimension of 's_w', which appears to only be used by 'AKv',
% or vertical viscosity coefficient, and 'AKt', or temperature vertical
% diffusion coefficient. Reached out to Sonaljit 17 May 2024 to confirm
cd(scriptDrive)
u_zlev = set_depth(2,4,theta_s,theta_b,hc,layers,3,h,zeta,0); %'3' indicates u-point grid
v_zlev = set_depth(2,4,theta_s,theta_b,hc,layers,4,h,zeta,0); %'4' indicates v-point grid
rho_zlev = set_depth(2,4,theta_s,theta_b,hc,layers,1,h,zeta,0); %'1' indicates rho-point (density) grid

%flip the depth dimension so that it is ordered surface to bottom
u_zlev = u_zlev(:,:,layers:-1:1);
v_zlev = v_zlev(:,:,layers:-1:1);
rho_zlev = rho_zlev(:,:,layers:-1:1);

%calculate the maximum depth of the uppermost 'set-depth'-produced
% z-level. use this as the first depth to interpolate to (the first z-level
% has to be equal to or greater than the first 'set_depth'-produced level,
% or an NA will be produced in the resulting value (e.g., u-velocity) at
% that depth and location
max_u_zlev = max(max(abs(squeeze(u_zlev(:,:,1)))));
max_v_zlev = max(max(abs(squeeze(v_zlev(:,:,1)))));
max_rho_zlev = max(max(abs(squeeze(rho_zlev(:,:,1)))));
max_zlev = max([max_u_zlev, max_v_zlev, max_rho_zlev]);

%define depth levels in meters
zlevels = [max_zlev 1 2 4 6 8 10 12 14 16 18 20 22 24 26 29 32 35 38 42 46 50 60 70 85 100 120 140 160 190 220 250]'; %32 layers down to 250 m depth, mid-depth hi-res
% zlevels = [0 1 2 3 4 5 6 7 8 9 10 13 16 19 22 26 30 35 40 45 50 60 70 80 90 110 130 150 170 190 220 250]'; %32 layers down to 250 m depth, shallow hi-res
% zlevels = [0 .5 1 2 3 4 5 6 7 8 9 10 15 20 25 30 35 40 45 50 55 60 70 80 90 100 120 145 170 195 220 250]'; %32 layers down to 250 m depth, shallow lo-res
% zlevels = [0 .5 1 2 3 4 5 6 7 8 9 10 15 20 25 30 35 40 45 50 55 60 70 80 90 100 110 120 130 140 150 160]'; %32 layers down to 160 m depth

% zlevels = [0 .5 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 28 30 32 34 36 38 40 42 44 46 48 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120 125 130 135 140 145 150 155 160 165 170 180 190 200 210 220 230 240 250]'; %72 layers down to 250 m depth
% zlevels = [0 .5 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 28 30 32 34 36 38 40 42 44 46 48 50 53 56 59 62 65 68 71 74 77 80 83 86 89 92 95 98 100 102 104 107 110 113 117 121 125 130 135 140 145 150 155 160]'; %72 layers down to 160 m depth
depth_range = [min(zlevels) max(zlevels)];

% plot(zlevels, '-rd')

% test = rescale(zlevels160, 0, 250);
% plot(test, '-rd');hold on
% plot(zlevels250, '-gd')
% plot(zlevels160, '-bd')

% Some relevant literature to inform z-levels:
%   - https://www.int-res.com/articles/feature/a071p193.pdf
%   - https://link.springer.com/content/pdf/10.1007/BF00258221.pdf
%   - https://www.frontiersin.org/articles/10.3389/frans.2022.868515/full
%   - https://www.sciencedirect.com/science/article/pii/S0304380011003115
%   - https://link.springer.com/article/10.1007/s00338-010-0593-6
%
%   - Cherubin 2011 saw USVI neutrally buoyant particles downwell then
%   upwell in rapid succession, but did not discuss much past 200 m of
%   depth.
%
%   - Huck 2022 demonstrated stark differences in neutrally buoyant plastic
%   residence times between the upper 200 m of ocean and anything below
%   that. Perhaps a good sign that once particles exit the upper 200 m,
%   they do not often come back.
%
%   - Overall, I think it is fair to cut off our model at 250 m. This is
%   above the 200 m distinction made above, yet allows for corals at the
%   very deepest limits to disperse and receive particles. Agaricia has
%   been discovered as deep as 119 m in the Bahamas (Reed 1985), and is
%   one of our susceptible species. 150 m is the theoretical mesophotic
%   limit, and the euphotic zone is limited to 200 m. The average depth of
%   the mixed layer is also 200 m. Going a bit deeper than this provides
%   some wiggle room for particles to stay in the model that otherwise
%   wouldn't.

%% Read in all USCROMS files in a directory and rewrite them in another directory as CMS-ready netCDF files
cd(writeDrive)
error_log = 'error_log.txt';
warning_log = 'warning_log.text';

cd(readDrive)
myDir = pwd;
files = dir(fullfile(myDir,'croco_his*'));

%ensure the 'files' structure is sorted correctly, from hour 0 on January
% 1st 2019, to hour 8631 (happens to be the 359th day of the simulation run,
% when Dr. Sonaljit Mukherjee completed the simulation)
fileNumbers = cellfun(@(x) sscanf(x, 'croco_his.%05d.nc'), {files.name}); %extract numerical values from filenames
[~, sortedIndices] = sort(fileNumbers); %sort files based on numerical values
files = files(sortedIndices);

t = 0; %start at day 0 - if you are testing a simulation NOT starting on January 1st, 2019 at midnight, edit this
% % dummy stuff to test contents of loop:
% myDir = pwd;
% i=1;j=1;
% % files = dir(fullfile(myDir,'croco_his.03741.nc')); %create an index for every netcdf file in directory
% % j = 3;
% % dummy stuff to test contents of loop
for i = 1:length(files) %loop through each file
    baseFileName = files(i).name; %name of the file itself
    fullFileName = fullfile(myDir, baseFileName); %path to the file name
    fprintf(1, 'Now reading %s\n', fullFileName);
    
    for j = 1:3 %loop through each hour in each file. there are three hours per file
        
        % This tracks how far along the script is and monitors
         % output. if any files are missing in theenvironment, it will not
         % properly reflect the within-file hour sequence
        t = t+1;
        fprintf(1, 'Hour is %.1f \n', t);
                
        %create new .nc file names for CMS, for u,v,w
        cd(readDrive) %ensures still in the read drive when the loop resets
        secondsaccum = ncread(fullFileName, 'time'); %ROMS accumulates number of seconds since simulation initialization for its 'time' variable
        hoursaccum = secondsaccum/60/60; %each time-step is an hour interval
        startdate = datetime(2019,1,1,0,0,0); %beginning of ROMS output
        hoursequence = j;
        currentdate = startdate + hours(hoursaccum(hoursequence)); %calculate current date by adding #/hours passed since midnight January 1st
        [y,m,d] = ymd(currentdate);
        hr = hour(currentdate);
        cms_fname = strcat('nest_1_',num2str(y),sprintf('%02d',m),sprintf('%02d',d),sprintf('%02d',hr),'0000.nc'); %nest_nestnumber_yyyymmddhhmmss [.nc]

        %update output log
        cd(writeDrive)
        fileName = 'output_log.txt';
        fileID = fopen(fileName, 'a');
        fprintf(fileID, '%s: Time at end of job (whether from error or end of data environment).\n', char(datetime('now')));
        fprintf(fileID, 'Original .nc-file path: %s\n', fullFileName);
        fprintf(fileID, 'CMS-ready file name: %s\n\n', cms_fname);
        fclose(fileID);
 
        %create new .nc files and dummy variables in correct dimensions
        % fill_value = -30; % NOTE - could try 1.267650600228229e+30 (from MAR Dan run, March 2023 on HPC. CMS one: 1.2676506e+30
        % fill_value = 1.2676506e+30; % NOTE - testing!
        nccreate(cms_fname,'Longitude','Dimensions',{'Lon',length(lon_v_rho_convert(:,1))},'Format','netcdf4','Datatype','double')
        nccreate(cms_fname,'Latitude','Dimensions',{'Lat',length(lat_rho(1,:))},'Format','netcdf4','Datatype','double')
        nccreate(cms_fname,'Depth','Dimensions',{'Depth',length(zlevels)},'Format','netcdf4','Datatype','double')
        nccreate(cms_fname,'Time','Dimensions',{'Time',1},'Format','netcdf4','Datatype','double'); %include only 1 time slice per .nc output file
        nccreate(cms_fname,'zu','Dimensions',{'Lon','Lat','Depth','Time'},'Format','netcdf4','Datatype','single') %,'FillValue',fill_value
        nccreate(cms_fname,'zv','Dimensions',{'Lon','Lat','Depth','Time'},'Format','netcdf4','Datatype','single')
        nccreate(cms_fname,'zw','Dimensions',{'Lon','Lat','Depth','Time'},'Format','netcdf4','Datatype','single')
        nccreate(cms_fname,'zt','Dimensions',{'Lon','Lat','Depth','Time'},'Format','netcdf4','Datatype','single')
        nccreate(cms_fname,'zs','Dimensions',{'Lon','Lat','Depth','Time'},'Format','netcdf4','Datatype','single')
        
        %assign theta_s, theta_b and hc as global attributes
        ncwriteatt(cms_fname, '/', 'theta_s', theta_s)
        ncwriteatt(cms_fname, '/', 'theta_b', theta_b)
        ncwriteatt(cms_fname, '/', 'hc', hc)
        ncwriteatt(cms_fname, '/', 'source', fullFileName)

        %write variables
        ncwrite(cms_fname,'Longitude',lon_v_rho_convert(:,1))
        ncwrite(cms_fname,'Latitude',lat_rho(1,:))
        ncwrite(cms_fname,'Depth',zlevels)

        %extract local variables, then break into current time-slice
        cd(readDrive)
        uvel = ncread(fullFileName, 'u');
        vvel = ncread(fullFileName, 'v');
        wvel = ncread(fullFileName, 'w');
        salt = ncread(fullFileName,'salt');
        temp = ncread(fullFileName,'temp');
        ssh = ncread(fullFileName, 'zeta');
        uvel = uvel(:,:,:,j);
        vvel = vvel(:,:,:,j);
        wvel = wvel(:,:,:,j);
        temp = temp(:,:,:,j);
        salt = salt(:,:,:,j);
        ssh = ssh(:,:,j);
        
        %flip the depth dimension so that it is ordered surface to bottom,
        %with increasing number indicating deeper depth
        uvel = uvel(:,:,layers:-1:1);
        vvel = vvel(:,:,layers:-1:1);
        wvel = wvel(:,:,layers:-1:1);
        temp = temp(:,:,layers:-1:1);
        salt = salt(:,:,layers:-1:1);
        
        %replace zeros with NaNs for interpolation. this is identical to
        %converting the landmask to NaNs
        % NOTE - this is making "land" at actual land, seafloor, AND
        % anywhere that interpolation failed (if there was a failure -
        % which can be hard to detect). Pay careful attention to the
        % actions of 'convert_sigma_z'
        uvel(uvel==0) = NaN;
        vvel(vvel==0) = NaN;
        wvel(wvel==0) = NaN; 
        temp(temp==0) = NaN;
        salt(salt==0) = NaN;
       
        %interpolate u,v,w,t,s at z-levels
        cd(scriptDrive)
        try %error handler - loop will not break if interpolation fails for a file
            [u_new,v_new,w_new,t_new,s_new] = convert_sigma_z(uvel,vvel,wvel,temp,salt,h,u_zlev,v_zlev,rho_zlev,zlevels,scriptDrive); %layers,theta_s,theta_b,hc,h,ssh,warning_log,cms_fname,fullFileName
        catch ME
            %update the error log with the current errant file
            cd(writeDrive)
            fileID = fopen(error_log, 'a');
            fprintf(fileID, '%s: Time of error (due to an issue with interpolation).\n', char(datetime('now')));
            fprintf(fileID, 'Error type: %s\n', ME.message); %the error spit out directly by MATLAB
            fprintf(fileID, 'Original .nc-file path: %s\n', fullFileName);
            fprintf(fileID, 'CMS-ready file name: %s\n\n', cms_fname);
            fclose(fileID);

            continue %jump to the next iteration of j-loop
        end

        %convert to single format rather than double, lowers precision but
        % halves file size. CMS is okay with singles. may want to return to
        % this
        u_new = single(u_new);
        v_new = single(v_new);
        w_new = single(w_new);
        t_new = single(t_new);
        s_new = single(s_new);

        %calculate ranges for variables
        % NOTE: this needs to happen before replacing NaNs with CMS
        % fill value! Otherwise the max on the range is not useful
        % information (it'll be 1.2676506e+30)
        uvel_range = [min(min(min(u_new))) max(max(max(u_new)))];
        vvel_range = [min(min(min(v_new))) max(max(max(v_new)))];
        wvel_range = [min(min(min(w_new))) max(max(max(w_new)))];
        temp_range = [min(min(min(t_new))) max(max(max(t_new)))];
        salt_range = [min(min(min(s_new))) max(max(max(s_new)))];
        
        %replace NaNs with fill value for CMS
        fill_value = single(1.2676506e+30);
        % fill_value = 0; % NOTE - testing!
        u_new(isnan(u_new)) = fill_value;
        v_new(isnan(v_new)) = fill_value;
        w_new(isnan(w_new)) = fill_value;
        t_new(isnan(t_new)) = fill_value;
        s_new(isnan(s_new)) = fill_value;
        
        %write sigma-to-z converted variables to files
        cd(writeDrive)
        ncwrite(cms_fname, 'zu', u_new)
        ncwrite(cms_fname, 'zv', v_new)
        ncwrite(cms_fname, 'zw', w_new)
        ncwrite(cms_fname, 'zt', t_new)
        ncwrite(cms_fname, 'zs', s_new)
        
        %assign local attributes
        ncwriteatt(cms_fname, 'Longitude', 'units', 'degrees_east')
        ncwriteatt(cms_fname, 'Longitude', 'valid_range', lon_u_range)
        ncwriteatt(cms_fname, 'Latitude', 'units', 'degrees_north')
        ncwriteatt(cms_fname, 'Latitude', 'valid_range', lat_u_rho_range)        
        ncwriteatt(cms_fname, 'Depth', 'units', 'meter')
        ncwriteatt(cms_fname, 'Depth', 'positive', 'down')
        ncwriteatt(cms_fname, 'Depth', 'valid_range', depth_range)
        ncwriteatt(cms_fname, 'Time', 'units', 'hours since midnight January 1 2019')
        ncwriteatt(cms_fname, 'zu', 'units', 'm/s')
        ncwriteatt(cms_fname, 'zu', 'valid_range', uvel_range) %% NOTE: not sure if this is necessary / if I'm doing it correctly      
        ncwriteatt(cms_fname, 'zu', 'long_name', 'eastward_sea_water_velocity')
        ncwriteatt(cms_fname, 'zv', 'units', 'm/s')
        ncwriteatt(cms_fname, 'zv', 'valid_range', vvel_range) %% NOTE: not sure if this is necessary / if I'm doing it correctly      
        ncwriteatt(cms_fname, 'zv', 'long_name', 'northward_sea_water_velocity')
        ncwriteatt(cms_fname, 'zw', 'units', 'm/s')
        ncwriteatt(cms_fname, 'zw', 'valid_range', wvel_range) %% NOTE: not sure if this is necessary / if I'm doing it correctly      
        ncwriteatt(cms_fname, 'zw', 'long_name', 'downward sea velocity')
        ncwriteatt(cms_fname, 'zt', 'units', 'Celcius')
        ncwriteatt(cms_fname, 'zt', 'valid_range', temp_range) %% NOTE: not sure if this is necessary / if I'm doing it correctly      
        ncwriteatt(cms_fname, 'zt', 'long_name', 'temperature')
        ncwriteatt(cms_fname, 'zs', 'units', 'PSU')
        ncwriteatt(cms_fname, 'zs', 'valid_range', salt_range) %% NOTE: not sure if this is necessary / if I'm doing it correctly      
        ncwriteatt(cms_fname, 'zs', 'long_name', 'salinity')
    end
end

toc