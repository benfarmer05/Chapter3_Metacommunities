%% Script to convert netCDF file output file from HYCOM .nc-files to a netCDF file that CMS can use

%% extract and define global variables
clear;clc

cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed/Testing'

examplenest2 = 'nest_1_20021129080000.nc';
ncdisp(examplenest2)

theta_s = ncread(ROMS_grid, 'theta_s'); %this kind of stuff is standard for ROMS but may not be for HYCOM. they are called global attributes
theta_b = ncread(ROMS_grid, 'theta_b');
hc = ncread(ROMS_grid, 'hc');
h = ncread(ROMS_grid, 'h'); %bathymetry (2-D matrix of depth at every lat/lon)
lon_u = ncread(ROMS_grid, 'lon_u');
lat_u = ncread(ROMS_grid, 'lat_u');
lon_v = ncread(ROMS_grid, 'lon_v');
lat_v = ncread(ROMS_grid, 'lat_v');
lon_rho = ncread(ROMS_grid, 'lon_rho'); %w-vel, temp, salinity, zeta, and rho use rho-grid for ROMS. HYCOM is A-grid Arakawa so you won't have to deal with that issue
lat_rho = ncread(ROMS_grid, 'lat_rho'); %w-vel, temp, salinity, zeta, and rho use rho-grid

%make longitudes positive & specify ranges for attributes [HYCOM may be
%different]
lon_u_convert = 360-abs(lon_u);
lon_u_range = [min(lon_u_convert(:,1)) max(lon_u_convert(:,1))];
lon_v_rho_convert = 360-abs(lon_v); %v & rho - they are on same grid for longitudes
lon_v_rho_range = [min(lon_v_rho_convert(:,1)) max(lon_v_rho_convert(:,1))];
lat_u_rho_range = [min(lat_u(1,:)) max(lat_u(1,:))]; %u & rho - they are on same grid for latitudes
lat_v_range = [min(lat_v(1,:)) max(lat_v(1,:))];

%% Loop through files

myDir = pwd; %navigate manually or with 'cd' code to the directory of interest
files = dir(fullfile(myDir,'croco_his*')); %create an index for every netcdf file in directory. tailor this to the HYCOM naming conventions

for i = 1:length(files) %loop through each file
    baseFileName = files(i).name; %name of the file itself
    fullFileName = fullfile(myDir, baseFileName); %path to the file name
    fprintf(1, 'Now reading %s\n', fullFileName);
    
    for j = 1:3 %loop through each hour in each file [note - tailor to HYCOM. probably don't need this loop at all
        t = t+1;
        fprintf(1, 'Hour is %.1f \n', t);
        
        %create new .nc file names for CMS, for u,v,w [note - also do this
        %for whatever else you're interested in w/ HYCOM]
        secondsaccum = ncread(fullFileName, 'time'); %ROMS accumulates number of seconds since simulation initialization for its 'time' variable. HYCOM may be different
        hoursaccum = secondsaccum/60/60; %each time-step is an hour interval. again, probably unnecessary for HYCOM
        startdate = datetime(2019,1,1,0,0,0); %beginning of ROMS output [may have to tailor this by year for HYCOM]
        hoursequence = j;
        currentdate = startdate + hours(hoursaccum(hoursequence)); %calculate current date by adding #/hours passed since midnight January 1st
        [y,m,d] = ymd(currentdate);
        hr = hour(currentdate);
        cms_fname = strcat('nest_1_',num2str(y),sprintf('%02d',m),sprintf('%02d',d),sprintf('%02d',hr),'0000.nc'); %nest_nestnumber_yyyymmddhhmmss [.nc]
 
        %create new .nc files and dummy variables in correct dimensions
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
        
        %write variables [NOTE - once you do this, you cannot edit the
        %.nc-files you have created. to run this loop again, you would have
        %to delete all your output and start from scratch. netCDF is
        %finicky like that]
        ncwrite(cms_fname,'Longitude',lon_v_rho_convert(:,1))
        ncwrite(cms_fname,'Latitude',lat_rho(1,:))
        ncwrite(cms_fname,'Depth',zlevels)

        %extract local variables, then break into current time-slice
        uvel = ncread(fullFileName, 'u');
        vvel = ncread(fullFileName, 'v');
        wvel = ncread(fullFileName, 'w');
        salt = ncread(fullFileName,'salt');
        temp = ncread(fullFileName,'temp');
        ssh = ncread(fullFileName, 'zeta');
        uvel = uvel(:,:,:,j); %this may only be 3 dimensions for HYCOM, not 4 (the fourth for ROMS here was time)
        vvel = vvel(:,:,:,j);
        wvel = wvel(:,:,:,j);
        temp = temp(:,:,:,j);
        salt = salt(:,:,:,j);
        ssh = ssh(:,:,j);
        
        %flip the depth dimension so that it is ordered surface to bottom,
        %with increasing number indicating deeper depth [this may be
        %unnecessary for HYCOM)
        layers = size(salt,3); %salinity variable used arbitrarily, just to find the number of depth layers
        uvel = uvel(:,:,layers:-1:1);
        vvel = vvel(:,:,layers:-1:1);
        wvel = wvel(:,:,layers:-1:1);
        temp = temp(:,:,layers:-1:1);
        salt = salt(:,:,layers:-1:1);
        
        %replace zeros with NaNs for interpolation. this is identical to
        %converting the landmask to NaNs (consider whether this is
        %necessary for CMS given your application)
        uvel(uvel==0) = NaN;
        vvel(vvel==0) = NaN;
        wvel(wvel==0) = NaN;
        temp(temp==0) = NaN;
        salt(salt==0) = NaN;
        
        %interpolate u,v,w at z-levels
        % [NOTE - for HYCOM, you're not going to have to do this since it's
        % already on Arakawa A-grid]
        [u_new,v_new,w_new,t_new,s_new] = convert_sigma_z(uvel,vvel,wvel,temp,salt,theta_s,theta_b,hc,h,ssh,zlevels,layers);

        %convert to single format rather than double, lowers precision but
        %halves file size. CMS is okay with singles. may want to return to
        %this
        u_new = single(u_new); %again, this is just a ROMS thing. after interpolation I had to transform the variables back to singles to conserve space
        v_new = single(v_new); %    I strongly recommend investigating whether the HYCOM output is in 'single' or 'double' format.
        w_new = single(w_new); %    'double' format will double the storage requirements of all of your input files to CMS, when single precision...
        t_new = single(t_new); %    ...is likely all you need for the 1-km resolution with HYCOM
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
        
        %replace NaNs with fill value for CMS. [NOTE - fill value is
        %extremely finicky with CMS. I recommend testing a few different
        %options and making sure what definitely works without any "under
        %the hood" issues
        fill_value = single(1.2676506e+30);
        % fill_value = 0; % NOTE - testing!
        u_new(isnan(u_new)) = fill_value;
        v_new(isnan(v_new)) = fill_value;
        w_new(isnan(w_new)) = fill_value;
        t_new(isnan(t_new)) = fill_value;
        s_new(isnan(s_new)) = fill_value;
        
        %write sigma-to-z converted variables to files (in your case,
        %simply renamed files, not converted)
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