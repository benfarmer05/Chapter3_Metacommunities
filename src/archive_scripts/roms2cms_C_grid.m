%% Script to convert netCDF file output file from ROMS to a netCDF file that CMS can use
% 28 June 2023

%% IMPORTANT NOTE
% So, this doesn't work. The reason is that the CMS does not have any
% references to lon_nameS or lat_nameS. It only has C-grid name references
% for u,v,w, and temp. This is true across at least mod_getdata.f90,
% getnestinfo.f90, and getfilename.f90.


tic
%% Large grid
%% extract and define global variables
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023'
clear;clc

% ROMS_grid = 'parent2km100.nc';
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

%depth levels in meters
% zlevels = [0 .5 1 2 3 4 5 6 7 8 9 10 15 20 25 30 35 40 45 50 75 100 125 150 175 200 250 300 350 400 450 500]';
zlevels = [0 1 2 4 6 8 10 12 14 16 18 20 22 24 26 29 32 35 38 42 46 50 60 70 85 100 120 140 160 190 220 250]'; %32 layers down to 250 m depth, mid-depth hi-res
depth_range = [min(zlevels) max(zlevels)];
% NOTE - maybe try 300 m as a middle-ground. The reasoning: if a 0 m-depth
% coral could, in theory, get its larvae to a 150-m coral, then a 150-m
% coral could send its larvae another 150-m down and back up. Also 300 m is
% a nice clean number to simulate to.

%% Loop through files

% cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/VI_ROMS/CMS-formatting_Spring2022_to_Summer2023/2022_hydro'
% ncdisp("croco_his.04320.nc")
% cd '/ddnA/project/holstein/bfarme8/USVI/Nests/STT-330m'

% cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro
% cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed
cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed/Testing

myDir = pwd;
files = dir(fullfile(myDir,'croco_his*')); %create an index for every netcdf file in directory

t = 0;
% % dummy stuff to test contents of loop:
% cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed/Testing
% myDir = pwd;
% i=1;j=1;
% files = dir(fullfile(myDir,'croco_his.01224.nc')); %create an index for every netcdf file in directory
for i = 1:length(files) %loop through each file
    baseFileName = files(i).name; %name of the file itself
    fullFileName = fullfile(myDir, baseFileName); %path to the file name
    fprintf(1, 'Now reading %s\n', fullFileName);
    
    for j = 1:3 %loop through each hour in each file. there are three hours per file
        t = t+1;
        fprintf(1, 'Hour is %.1f \n', t);
        
        %create new .nc file names for CMS, for u,v,w
        secondsaccum = ncread(fullFileName, 'time'); %ROMS accumulates number of seconds since simulation initialization for its 'time' variable
        hoursaccum = secondsaccum/60/60; %each time-step is an hour interval
        startdate = datetime(2019,1,1,0,0,0); %beginning of ROMS output
        hoursequence = j;
        currentdate = startdate + hours(hoursaccum(hoursequence)); %calculate current date by adding #/hours passed since midnight January 1st
        [y,m,d] = ymd(currentdate);
        hr = hour(currentdate);
        CMS_fname = strcat('nest_1_',num2str(y),sprintf('%02d',m),sprintf('%02d',d),sprintf('%02d',hr),'0000'); %nest_nestnumber_yyyymmddhhmmss [.nc]
        cms_uname = strcat(CMS_fname,'u.nc');
        cms_vname = strcat(CMS_fname,'v.nc');
        cms_wname = strcat(CMS_fname,'w.nc');
        cms_tname = strcat(CMS_fname,'t.nc');
        cms_sname = strcat(CMS_fname,'s.nc');

        %create new .nc files and dummy variables in correct dimensions
        % fill_value = -30; % NOTE - could try 1.267650600228229e+30 (from MAR Dan run, March 2023 on HPC. CMS one: 1.2676506e+30
        nccreate(cms_uname,'lonu','Dimensions',{'Lon',length(lon_u(:,1))},'Format','netcdf4','Datatype','double')
        nccreate(cms_uname,'latu','Dimensions',{'Lat',length(lat_u(1,:))},'Format','netcdf4','Datatype','double')
        nccreate(cms_uname,'Depth','Dimensions',{'Depth',length(zlevels)},'Format','netcdf4','Datatype','double')
        nccreate(cms_uname,'Time','Dimensions',{'Time',1},'Format','netcdf4','Datatype','double'); %include only 1 time slice per .nc output file
        nccreate(cms_uname,'u','Dimensions',{'Lon','Lat','Depth','Time'},'Format','netcdf4','Datatype','single')

        nccreate(cms_vname,'lonv','Dimensions',{'Lon',length(lon_v(:,1))},'Format','netcdf4','Datatype','double')
        nccreate(cms_vname,'latv','Dimensions',{'Lat',length(lat_v(1,:))},'Format','netcdf4','Datatype','double')
        nccreate(cms_vname,'Depth','Dimensions',{'Depth',length(zlevels)},'Format','netcdf4','Datatype','double')
        nccreate(cms_vname,'Time','Dimensions',{'Time',1},'Format','netcdf4','Datatype','double')
        nccreate(cms_vname,'v','Dimensions',{'Lon','Lat','Depth','Time'},'Format','netcdf4','Datatype','single')

        nccreate(cms_wname,'lonw','Dimensions',{'Lon',length(lon_rho(:,1))},'Format','netcdf4','Datatype','double')
        nccreate(cms_wname,'latw','Dimensions',{'Lat',length(lat_rho(1,:))},'Format','netcdf4','Datatype','double')
        nccreate(cms_wname,'Depth','Dimensions',{'Depth',length(zlevels)},'Format','netcdf4','Datatype','double')
        nccreate(cms_wname,'Time','Dimensions',{'Time',1},'Format','netcdf4','Datatype','double')
        nccreate(cms_wname,'w','Dimensions',{'Lon','Lat','Depth','Time'},'Format','netcdf4','Datatype','single')
        
        nccreate(cms_tname,'lont','Dimensions',{'Lon',length(lon_rho(:,1))},'Format','netcdf4','Datatype','double')
        nccreate(cms_tname,'latt','Dimensions',{'Lat',length(lat_rho(1,:))},'Format','netcdf4','Datatype','double')
        nccreate(cms_tname,'Depth','Dimensions',{'Depth',length(zlevels)},'Format','netcdf4','Datatype','double')
        nccreate(cms_tname,'Time','Dimensions',{'Time',1},'Format','netcdf4','Datatype','double')
        nccreate(cms_tname,'t','Dimensions',{'Lon','Lat','Depth','Time'},'Format','netcdf4','Datatype','single')

        nccreate(cms_sname,'lons','Dimensions',{'Lon',length(lon_rho(:,1))},'Format','netcdf4','Datatype','double')
        nccreate(cms_sname,'lats','Dimensions',{'Lat',length(lat_rho(1,:))},'Format','netcdf4','Datatype','double')
        nccreate(cms_sname,'Depth','Dimensions',{'Depth',length(zlevels)},'Format','netcdf4','Datatype','double')
        nccreate(cms_sname,'Time','Dimensions',{'Time',1},'Format','netcdf4','Datatype','double')
        nccreate(cms_sname,'s','Dimensions',{'Lon','Lat','Depth','Time'},'Format','netcdf4','Datatype','single')

        %assign theta_s, theta_b and hc as global attributes
        ncwriteatt(cms_uname, '/', 'theta_s', theta_s)
        ncwriteatt(cms_uname, '/', 'theta_b', theta_b)
        ncwriteatt(cms_uname, '/', 'hc', hc)
        ncwriteatt(cms_uname, '/', 'source', fullFileName)
        ncwriteatt(cms_vname, '/', 'theta_s', theta_s)
        ncwriteatt(cms_vname, '/', 'theta_b', theta_b)
        ncwriteatt(cms_vname, '/', 'hc', hc)
        ncwriteatt(cms_vname, '/', 'source', fullFileName)
        ncwriteatt(cms_wname, '/', 'theta_s', theta_s)
        ncwriteatt(cms_wname, '/', 'theta_b', theta_b)
        ncwriteatt(cms_wname, '/', 'hc', hc)
        ncwriteatt(cms_wname, '/', 'source', fullFileName)
        ncwriteatt(cms_tname, '/', 'theta_s', theta_s)
        ncwriteatt(cms_tname, '/', 'theta_b', theta_b)
        ncwriteatt(cms_tname, '/', 'hc', hc)
        ncwriteatt(cms_tname, '/', 'source', fullFileName)
        ncwriteatt(cms_sname, '/', 'theta_s', theta_s)
        ncwriteatt(cms_sname, '/', 'theta_b', theta_b)
        ncwriteatt(cms_sname, '/', 'hc', hc)
        ncwriteatt(cms_sname, '/', 'source', fullFileName)
        
        %write variables
        ncwrite(cms_uname,'lonu',lon_u_convert(:,1))
        ncwrite(cms_uname,'latu',lat_u(1,:))
        ncwrite(cms_uname,'Depth',zlevels)
        ncwrite(cms_vname,'lonv',lon_v_rho_convert(:,1))
        ncwrite(cms_vname,'latv',lat_v(1,:))
        ncwrite(cms_vname,'Depth',zlevels)
        ncwrite(cms_wname,'lonw',lon_v_rho_convert(:,1))
        ncwrite(cms_wname,'latw',lat_rho(1,:))
        ncwrite(cms_wname,'Depth',zlevels)
        ncwrite(cms_tname,'lont',lon_v_rho_convert(:,1))
        ncwrite(cms_tname,'latt',lat_rho(1,:))
        ncwrite(cms_tname,'Depth',zlevels)
        ncwrite(cms_sname,'lons',lon_v_rho_convert(:,1))
        ncwrite(cms_sname,'lats',lat_rho(1,:))
        ncwrite(cms_sname,'Depth',zlevels)
        
        %extract local variables, then break into current time-slice
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
        layers = size(salt,3);
        uvel = uvel(:,:,layers:-1:1);
        vvel = vvel(:,:,layers:-1:1);
        wvel = wvel(:,:,layers:-1:1);
        temp = temp(:,:,layers:-1:1);
        salt = salt(:,:,layers:-1:1);
        
        %replace zeros with NaNs for interpolation. this is identical to
        %converting the landmask to NaNs
        uvel(uvel==0) = NaN;
        vvel(vvel==0) = NaN;
        wvel(wvel==0) = NaN; %NOTE - this is making "land" at actual land, seafloor, AND anywhere that interpolation failed
        temp(temp==0) = NaN;
        salt(salt==0) = NaN;
        
        cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023
        %interpolate u,v,w at z-levels
        tic
        [u_new,v_new,w_new,t_new,s_new] = convert_sigma_z(uvel,vvel,wvel,temp,salt,theta_s,theta_b,hc,h,ssh,zlevels,layers);
        toc

        % %
        % figure(1);imagescn(lon_v',lat_v',v_new(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean curl;clim([-1.3 1.3]);axis equal
        % figure(1);imagescn(lon_u',lat_u',u_new(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;clim([-.6 .6]);axis equal
        % figure(1);imagescn(lon_u',lat_u',u_new(:,:,32)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;clim([-.6 .6]);axis equal
        % %

        % cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro
        % cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed
        cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed/Testing

        %convert to single format rather than double, lowers precision but
        %halves file size. CMS is okay with singles. may want to return to
        %this
        u_new = single(u_new);
        v_new = single(v_new);
        w_new = single(w_new);
        t_new = single(t_new);
        s_new = single(s_new);

        %calculate ranges for variables
        % NOTE: I think this needs to happen before replacing NaNs with CMS
        % fill value! Otherwise the max on the range is not useful
        % information (it'll be 1.2676506e+30)
        uvel_range = [min(min(min(u_new))) max(max(max(u_new)))];
        vvel_range = [min(min(min(v_new))) max(max(max(v_new)))];
        wvel_range = [min(min(min(w_new))) max(max(max(w_new)))];
        temp_range = [min(min(min(t_new))) max(max(max(t_new)))];
        salt_range = [min(min(min(s_new))) max(max(max(s_new)))];
        
        %replace NaNs with fill value for CMS
        fill_value = single(1.2676506e+30);
        u_new(isnan(u_new)) = fill_value;
        v_new(isnan(v_new)) = fill_value;
        w_new(isnan(w_new)) = fill_value;
        t_new(isnan(t_new)) = fill_value;
        s_new(isnan(s_new)) = fill_value;

        %write sigma-to-z converted variables to files
        ncwrite(cms_uname, 'u', u_new)
        ncwrite(cms_vname, 'v', v_new)
        ncwrite(cms_wname, 'w', w_new)
        ncwrite(cms_tname, 't', t_new)
        ncwrite(cms_sname, 's', s_new)

        % STOPPING POINT - [SKIP TO BOTTOM UPDATE FOR THIS STOPPING POINT] - I am now thoroughly confused on fill values and
        % what the CMS needs. ncwrite is automatically converting all fill
        % values to NaNs, which I thought might have been what broke the
        % CMS. But, functioning MAR hydro files are working just fine with
        % only 0's and NaNs - no fill values, even though ncdisp shows in
        % the variable attributes that a fill value was specified. Who
        % knows

        %replace NaNs with fill value for CMS? this must happen AFTER
        %'ncwrite', because netcdf automatically converts that fill value
        %to NaNs which *may* break the CMS. but can't overwrite an existing
        %ncfile, by then it's read only
        %
        % update - I may have found a silly workaround where I set the
        % fillvalue attribute as a fake value (-30) in the netcdf file, so
        % then netcdf doesn't automatically yeet all the fill values that I
        % actually want for the CMS. not sure if this will break the CMS
        % too...alternatively, might be able to get away with simply not
        % including an attribute indicating a fill value in the first
        % place! try that next
        %
        % 8 SEP UPDATE: now, in a bit of a rut. spline version, either with
        % surface extrap on OR off, produces silly output with 15 m/s
        % current speeds and stuff. akima simply takes forever (2-3X as
        % long as linear or spline). linear IS working right now - so need
        % to figure out if it will work well with surface extrap on or not.
        % also, gonna test whether it makes sense to change the way the
        % interpolation is being done in the first place. i.e., only
        % allowing depths 160 m or shallower to be accounted for in the
        % interp1 function!
        %
        % update: no. I think how it's being done is fine, any effect on
        % speed of interpolation is definitely neglible and I don't see the
        % quality of the output changing. main thing is avoiding spline and
        % now trying surface extrap again
        
        %assign local attributes
        %u
        ncwriteatt(cms_uname, 'lonu', 'units', 'degrees_east')
        ncwriteatt(cms_uname, 'lonu', 'valid_range', lon_u_range)
        ncwriteatt(cms_uname, 'latu', 'units', 'degrees_north')
        ncwriteatt(cms_uname, 'latu', 'valid_range', lat_u_rho_range)        
        ncwriteatt(cms_uname, 'Depth', 'units', 'meter')
        ncwriteatt(cms_uname, 'Depth', 'positive', 'down')
        ncwriteatt(cms_uname, 'Depth', 'valid_range', depth_range)
        ncwriteatt(cms_uname, 'Time', 'units', 'hours since midnight January 1 2019')
        ncwriteatt(cms_uname, 'u', 'units', 'm/s')
        ncwriteatt(cms_uname, 'u', 'valid_range', uvel_range) %% NOTE: not sure if this is necessary / if I'm doing it correctly      
        ncwriteatt(cms_uname, 'u', 'long_name', 'eastward_sea_water_velocity')
        %v
        ncwriteatt(cms_vname, 'lonv', 'units', 'degrees_east')
        ncwriteatt(cms_vname, 'lonv', 'valid_range', lon_v_rho_range)
        ncwriteatt(cms_vname, 'latv', 'units', 'degrees_north')
        ncwriteatt(cms_vname, 'latv', 'valid_range', lat_v_range)        
        ncwriteatt(cms_vname, 'Depth', 'units', 'meter')
        ncwriteatt(cms_vname, 'Depth', 'positive', 'down')
        ncwriteatt(cms_vname, 'Depth', 'valid_range', depth_range)
        ncwriteatt(cms_vname, 'Time', 'units', 'hours since midnight January 1 2019')
        ncwriteatt(cms_vname, 'v', 'units', 'm/s')
        ncwriteatt(cms_vname, 'v', 'valid_range', vvel_range) %% NOTE: not sure if this is necessary / if I'm doing it correctly      
        ncwriteatt(cms_vname, 'v', 'long_name', 'northward_sea_water_velocity')
        %w
        ncwriteatt(cms_wname, 'lonw', 'units', 'degrees_east')
        ncwriteatt(cms_wname, 'lonw', 'valid_range', lon_v_rho_range)
        ncwriteatt(cms_wname, 'latw', 'units', 'degrees_north')
        ncwriteatt(cms_wname, 'latw', 'valid_range', lat_u_rho_range)        
        ncwriteatt(cms_wname, 'Depth', 'units', 'meter')
        ncwriteatt(cms_wname, 'Depth', 'positive', 'down')
        ncwriteatt(cms_wname, 'Depth', 'valid_range', depth_range)
        ncwriteatt(cms_wname, 'Time', 'units', 'hours since midnight January 1 2019')
        ncwriteatt(cms_wname, 'w', 'units', 'm/s')
        ncwriteatt(cms_wname, 'w', 'valid_range', wvel_range) %% NOTE: not sure if this is necessary / if I'm doing it correctly      
        ncwriteatt(cms_wname, 'w', 'long_name', 'downward sea velocity')
        %t
        ncwriteatt(cms_tname, 'lont', 'units', 'degrees_east')
        ncwriteatt(cms_tname, 'lont', 'valid_range', lon_v_rho_range)
        ncwriteatt(cms_tname, 'latt', 'units', 'degrees_north')
        ncwriteatt(cms_tname, 'latt', 'valid_range', lat_u_rho_range)        
        ncwriteatt(cms_tname, 'Depth', 'units', 'meter')
        ncwriteatt(cms_tname, 'Depth', 'positive', 'down')
        ncwriteatt(cms_tname, 'Depth', 'valid_range', depth_range)
        ncwriteatt(cms_tname, 'Time', 'units', 'hours since midnight January 1 2019')
        ncwriteatt(cms_tname, 't', 'units', 'Celcius')
        ncwriteatt(cms_tname, 't', 'valid_range', temp_range) %% NOTE: not sure if this is necessary / if I'm doing it correctly      
        ncwriteatt(cms_tname, 't', 'long_name', 'temperature')
        %s
        ncwriteatt(cms_sname, 'lons', 'units', 'degrees_east')
        ncwriteatt(cms_sname, 'lons', 'valid_range', lon_v_rho_range)
        ncwriteatt(cms_sname, 'lats', 'units', 'degrees_north')
        ncwriteatt(cms_sname, 'lats', 'valid_range', lat_u_rho_range)       
        ncwriteatt(cms_sname, 'Depth', 'units', 'meter')
        ncwriteatt(cms_sname, 'Depth', 'positive', 'down')
        ncwriteatt(cms_sname, 'Depth', 'valid_range', depth_range)
        ncwriteatt(cms_sname, 'Time', 'units', 'hours since midnight January 1 2019')
        ncwriteatt(cms_sname, 's', 'units', 'PSU')
        ncwriteatt(cms_sname, 's', 'valid_range', salt_range) %% NOTE: not sure if this is necessary / if I'm doing it correctly      
        ncwriteatt(cms_sname, 's', 'long_name', 'salinity')
    end
end

toc
% %% Smol grid (nest 2)
% %% extract and define global variables
% cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/VI_ROMS/CMS-formatting_Spring2022'
% clear;clc
% 
% ROMS_grid = 'sttstj300m.nc';
% theta_s = ncread(ROMS_grid, 'theta_s');
% theta_b = ncread(ROMS_grid, 'theta_b');
% hc = ncread(ROMS_grid, 'hc');
% h = ncread(ROMS_grid, 'h');
% lon_u = ncread(ROMS_grid, 'lon_u');
% lat_u = ncread(ROMS_grid, 'lat_u');
% lon_v = ncread(ROMS_grid, 'lon_v');
% lat_v = ncread(ROMS_grid, 'lat_v');
% lon_rho = ncread(ROMS_grid, 'lon_rho'); %w-vel, temp, salinity, zeta, and rho use rho-grid
% lat_rho = ncread(ROMS_grid, 'lat_rho'); %w-vel, temp, salinity, zeta, and rho use rho-grid
% 
% %make longitudes positive & specify ranges for attributes
% lon_u_convert = 360-abs(lon_u);
% lon_u_range = [min(lon_u_convert(:,1)) max(lon_u_convert(:,1))];
% lon_v_rho_convert = 360-abs(lon_v); %v & rho - they are on same grid for longitudes
% lon_v_rho_range = [min(lon_v_rho_convert(:,1)) max(lon_v_rho_convert(:,1))];
% lat_u_rho_range = [min(lat_u(1,:)) max(lat_u(1,:))]; %u & rho - they are on same grid for latitudes
% lat_v_range = [min(lat_v(1,:)) max(lat_v(1,:))];
% 
% %depth levels in meters
% % zlevels = [0 .5 1 2 3 4 5 6 7 8 9 10 15 20 25 30 35 40 45 50 75 100 125 150 175 200 250 300 350 400 450 500]';
% zlevels = [0 .5 1 2 3 4 5 6 7 8 9 10 15 20 25 30 35 40 45 50 55 60 70 80 90 100 110 120 130 140 150 160]';
% depth_range = [min(zlevels) max(zlevels)];
% 
% %% Loop through files
% cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/VI_ROMS/CMS-formatting_Spring2022'
% % cd '/ddnA/project/holstein/bfarme8/USVI/Nests/STT-330m'
% myDir = pwd;
% files = dir(fullfile(myDir,'nest_his*')); %create an index for every nc file in directory
% 
% % t = 0;
% % %dummy stuff to test contents of loop:
% % i=1;j=1;
% % files = dir(fullfile(myDir,'nest_his.04365.nc')); %create an index for every netcdf file in directory
% % for i = 1:length(files) %loop through each file
%     baseFileName = files(i).name; %name of the file itself
%     fullFileName = fullfile(myDir, baseFileName); %path to the file name
%     fprintf(1, 'Now reading %s\n', fullFileName);
% 
%     for j = 1:3 %loop through each hour in each file. there are three hours per file
% %         t = t+1;
% %         fprintf(1, 'Hour is %.1f \n', t);
% 
%         %create new .nc file names for CMS, for u,v,w
%         secondsaccum = ncread(fullFileName, 'time'); %ROMS accumulates number of seconds since simulation initialization for its 'time' variable
%         hoursaccum = secondsaccum/60/60; %each time-step is an hour interval
%         startdate = datetime(2019,1,1,0,0,0); %beginning of ROMS output
%         hoursequence = j;
%         currentdate = startdate + hours(hoursaccum(hoursequence)); %calculate current date by adding #/hours passed since midnight January 1st
%         [y,m,d] = ymd(currentdate);
%         hr = hour(currentdate);
%         CMS_fname = strcat('nest_2_',num2str(y),sprintf('%02d',m),sprintf('%02d',d),sprintf('%02d',hr),'0000'); %nest_nestnumber_yyyymmddhhmmss [.nc]
%         cms_uname = strcat(CMS_fname,'u.nc');
%         cms_vname = strcat(CMS_fname,'v.nc');
%         cms_wname = strcat(CMS_fname,'w.nc');
%         cms_tname = strcat(CMS_fname,'t.nc');
%         cms_sname = strcat(CMS_fname,'s.nc');
% 
%         %create new .nc files and dummy variables in correct dimensions
%         fill_value = 1.2676506e+30;
%         nccreate(cms_uname,'lonu','Dimensions',{'Lon',length(lon_u(:,1))},'Format','netcdf4','Datatype','double')
%         nccreate(cms_uname,'latu','Dimensions',{'Lat',length(lat_u(1,:))},'Format','netcdf4','Datatype','double')
%         nccreate(cms_uname,'Depth','Dimensions',{'Depth',length(zlevels)},'Format','netcdf4','Datatype','double')
%         nccreate(cms_uname,'Time','Dimensions',{'Time',1},'Format','netcdf4','Datatype','double'); %include only 1 time slice per .nc output file
%         nccreate(cms_uname,'u','Dimensions',{'Lon','Lat','Depth','Time'},'Format','netcdf4','Datatype','double','FillValue',fill_value)
% 
%         nccreate(cms_vname,'lonv','Dimensions',{'Lon',length(lon_v(:,1))},'Format','netcdf4','Datatype','double')
%         nccreate(cms_vname,'latv','Dimensions',{'Lat',length(lat_v(1,:))},'Format','netcdf4','Datatype','double')
%         nccreate(cms_vname,'Depth','Dimensions',{'Depth',length(zlevels)},'Format','netcdf4','Datatype','double')
%         nccreate(cms_vname,'Time','Dimensions',{'Time',1},'Format','netcdf4','Datatype','double')
%         nccreate(cms_vname,'v','Dimensions',{'Lon','Lat','Depth','Time'},'Format','netcdf4','Datatype','double','FillValue',fill_value)
% 
%         nccreate(cms_wname,'lonw','Dimensions',{'Lon',length(lon_rho(:,1))},'Format','netcdf4','Datatype','double')
%         nccreate(cms_wname,'latw','Dimensions',{'Lat',length(lat_rho(1,:))},'Format','netcdf4','Datatype','double')
%         nccreate(cms_wname,'Depth','Dimensions',{'Depth',length(zlevels)},'Format','netcdf4','Datatype','double')
%         nccreate(cms_wname,'Time','Dimensions',{'Time',1},'Format','netcdf4','Datatype','double')
%         nccreate(cms_wname,'w','Dimensions',{'Lon','Lat','Depth','Time'},'Format','netcdf4','Datatype','double','FillValue',fill_value)
% 
%         nccreate(cms_tname,'lont','Dimensions',{'Lon',length(lon_rho(:,1))},'Format','netcdf4','Datatype','double')
%         nccreate(cms_tname,'latt','Dimensions',{'Lat',length(lat_rho(1,:))},'Format','netcdf4','Datatype','double')
%         nccreate(cms_tname,'Depth','Dimensions',{'Depth',length(zlevels)},'Format','netcdf4','Datatype','double')
%         nccreate(cms_tname,'Time','Dimensions',{'Time',1},'Format','netcdf4','Datatype','double')
%         nccreate(cms_tname,'t','Dimensions',{'Lon','Lat','Depth','Time'},'Format','netcdf4','Datatype','double','FillValue',fill_value)
% 
%         nccreate(cms_sname,'lons','Dimensions',{'Lon',length(lon_rho(:,1))},'Format','netcdf4','Datatype','double')
%         nccreate(cms_sname,'lats','Dimensions',{'Lat',length(lat_rho(1,:))},'Format','netcdf4','Datatype','double')
%         nccreate(cms_sname,'Depth','Dimensions',{'Depth',length(zlevels)},'Format','netcdf4','Datatype','double')
%         nccreate(cms_sname,'Time','Dimensions',{'Time',1},'Format','netcdf4','Datatype','double')
%         nccreate(cms_sname,'s','Dimensions',{'Lon','Lat','Depth','Time'},'Format','netcdf4','Datatype','double','FillValue',fill_value)
% 
%         %assign theta_s, theta_b and hc as global attributes
%         ncwriteatt(cms_uname, '/', 'theta_s', theta_s)
%         ncwriteatt(cms_uname, '/', 'theta_b', theta_b)
%         ncwriteatt(cms_uname, '/', 'hc', hc)
%         ncwriteatt(cms_uname, '/', 'source', fullFileName)
%         ncwriteatt(cms_vname, '/', 'theta_s', theta_s)
%         ncwriteatt(cms_vname, '/', 'theta_b', theta_b)
%         ncwriteatt(cms_vname, '/', 'hc', hc)
%         ncwriteatt(cms_vname, '/', 'source', fullFileName)
%         ncwriteatt(cms_wname, '/', 'theta_s', theta_s)
%         ncwriteatt(cms_wname, '/', 'theta_b', theta_b)
%         ncwriteatt(cms_wname, '/', 'hc', hc)
%         ncwriteatt(cms_wname, '/', 'source', fullFileName)
%         ncwriteatt(cms_tname, '/', 'theta_s', theta_s)
%         ncwriteatt(cms_tname, '/', 'theta_b', theta_b)
%         ncwriteatt(cms_tname, '/', 'hc', hc)
%         ncwriteatt(cms_tname, '/', 'source', fullFileName)
%         ncwriteatt(cms_sname, '/', 'theta_s', theta_s)
%         ncwriteatt(cms_sname, '/', 'theta_b', theta_b)
%         ncwriteatt(cms_sname, '/', 'hc', hc)
%         ncwriteatt(cms_sname, '/', 'source', fullFileName)
% 
%         %write variables
%         ncwrite(cms_uname,'lonu',lon_u_convert(:,1))
%         ncwrite(cms_uname,'latu',lat_u(1,:))
%         ncwrite(cms_uname,'Depth',zlevels)
%         ncwrite(cms_vname,'lonv',lon_v_rho_convert(:,1))
%         ncwrite(cms_vname,'latv',lat_v(1,:))
%         ncwrite(cms_vname,'Depth',zlevels)
%         ncwrite(cms_wname,'lonw',lon_v_rho_convert(:,1))
%         ncwrite(cms_wname,'latw',lat_rho(1,:))
%         ncwrite(cms_wname,'Depth',zlevels)
%         ncwrite(cms_tname,'lont',lon_v_rho_convert(:,1))
%         ncwrite(cms_tname,'latt',lat_rho(1,:))
%         ncwrite(cms_tname,'Depth',zlevels)
%         ncwrite(cms_sname,'lons',lon_v_rho_convert(:,1))
%         ncwrite(cms_sname,'lats',lat_rho(1,:))
%         ncwrite(cms_sname,'Depth',zlevels)
% 
%         %extract local variables, then break into current time-slice
%         uvel = ncread(fullFileName, 'u');
%         vvel = ncread(fullFileName, 'v');
%         wvel = ncread(fullFileName, 'w');
%         salt = ncread(fullFileName,'salt');
%         temp = ncread(fullFileName,'temp');
%         ssh = ncread(fullFileName, 'zeta');
%         uvel = uvel(:,:,:,j);
%         vvel = vvel(:,:,:,j);
%         wvel = wvel(:,:,:,j);
%         temp = temp(:,:,:,j);
%         salt = salt(:,:,:,j);
%         ssh = ssh(:,:,j);
% 
%         %flip the depth dimension so that it is ordered surface to bottom,
%         %with increasing number indicating deeper depth
%         uvel = uvel(:,:,32:-1:1);
%         vvel = vvel(:,:,32:-1:1);
%         wvel = wvel(:,:,32:-1:1);
%         temp = temp(:,:,32:-1:1);
%         salt = salt(:,:,32:-1:1);
% 
%         %replace zeros with NaNs for interpolation. this is identical to
%         %converting the landmask to NaNs
%         uvel(uvel==0) = NaN;
%         vvel(vvel==0) = NaN;
%         wvel(wvel==0) = NaN;
%         temp(temp==0) = NaN;
%         salt(salt==0) = NaN;
% 
%         %interpolate u,v,w at z-levels
%         [u_new,v_new,w_new,t_new,s_new] = convert_sigma_z(uvel,vvel,wvel,temp,salt,theta_s,theta_b,hc,h,ssh,zlevels);
% 
%         %calculate ranges for variables
%         % NOTE: I think this needs to happen before replacing NaNs with CMS
%         % fill value! Otherwise the max on the range is not useful
%         % information (it'll be 1.2676506e+30)
%         uvel_range = [min(min(min(u_new))) max(max(max(u_new)))];
%         vvel_range = [min(min(min(v_new))) max(max(max(v_new)))];
%         wvel_range = [min(min(min(w_new))) max(max(max(w_new)))];
%         temp_range = [min(min(min(t_new))) max(max(max(t_new)))];
%         salt_range = [min(min(min(s_new))) max(max(max(s_new)))];
% 
%         %replace NaNs with fill value for CMS
%         u_new(isnan(u_new)) = fill_value;
%         v_new(isnan(v_new)) = fill_value;
%         w_new(isnan(w_new)) = fill_value;
%         t_new(isnan(t_new)) = fill_value;
%         s_new(isnan(s_new)) = fill_value;
% 
%         %write sigma-to-z converted variables to files
%         ncwrite(cms_uname, 'u', u_new)
%         ncwrite(cms_vname, 'v', v_new)
%         ncwrite(cms_wname, 'w', w_new)
%         ncwrite(cms_tname, 't', t_new)
%         ncwrite(cms_sname, 's', s_new)
% 
%         %assign local attributes
%         %u
%         ncwriteatt(cms_uname, 'lonu', 'units', 'degrees_east')
%         ncwriteatt(cms_uname, 'lonu', 'valid_range', lon_u_range)
%         ncwriteatt(cms_uname, 'latu', 'units', 'degrees_north')
%         ncwriteatt(cms_uname, 'latu', 'valid_range', lat_u_rho_range)        
%         ncwriteatt(cms_uname, 'Depth', 'units', 'meter')
%         ncwriteatt(cms_uname, 'Depth', 'positive', 'down')
%         ncwriteatt(cms_uname, 'Depth', 'valid_range', depth_range)
%         ncwriteatt(cms_uname, 'Time', 'units', 'hours since midnight January 1 2019')
%         ncwriteatt(cms_uname, 'u', 'units', 'm/s')
%         ncwriteatt(cms_uname, 'u', 'valid_range', uvel_range) %% NOTE: not sure if this is necessary / if I'm doing it correctly      
%         ncwriteatt(cms_uname, 'u', 'long_name', 'eastward_sea_water_velocity')
%         %v
%         ncwriteatt(cms_vname, 'lonv', 'units', 'degrees_east')
%         ncwriteatt(cms_vname, 'lonv', 'valid_range', lon_v_rho_range)
%         ncwriteatt(cms_vname, 'latv', 'units', 'degrees_north')
%         ncwriteatt(cms_vname, 'latv', 'valid_range', lat_v_range)        
%         ncwriteatt(cms_vname, 'Depth', 'units', 'meter')
%         ncwriteatt(cms_vname, 'Depth', 'positive', 'down')
%         ncwriteatt(cms_vname, 'Depth', 'valid_range', depth_range)
%         ncwriteatt(cms_vname, 'Time', 'units', 'hours since midnight January 1 2019')
%         ncwriteatt(cms_vname, 'v', 'units', 'm/s')
%         ncwriteatt(cms_vname, 'v', 'valid_range', vvel_range) %% NOTE: not sure if this is necessary / if I'm doing it correctly      
%         ncwriteatt(cms_vname, 'v', 'long_name', 'northward_sea_water_velocity')
%         %w
%         ncwriteatt(cms_wname, 'lonw', 'units', 'degrees_east')
%         ncwriteatt(cms_wname, 'lonw', 'valid_range', lon_v_rho_range)
%         ncwriteatt(cms_wname, 'latw', 'units', 'degrees_north')
%         ncwriteatt(cms_wname, 'latw', 'valid_range', lat_u_rho_range)        
%         ncwriteatt(cms_wname, 'Depth', 'units', 'meter')
%         ncwriteatt(cms_wname, 'Depth', 'positive', 'down')
%         ncwriteatt(cms_wname, 'Depth', 'valid_range', depth_range)
%         ncwriteatt(cms_wname, 'Time', 'units', 'hours since midnight January 1 2019')
%         ncwriteatt(cms_wname, 'w', 'units', 'm/s')
%         ncwriteatt(cms_wname, 'w', 'valid_range', wvel_range) %% NOTE: not sure if this is necessary / if I'm doing it correctly      
%         ncwriteatt(cms_wname, 'w', 'long_name', 'downward sea velocity')
%         %t
%         ncwriteatt(cms_tname, 'lont', 'units', 'degrees_east')
%         ncwriteatt(cms_tname, 'lont', 'valid_range', lon_v_rho_range)
%         ncwriteatt(cms_tname, 'latt', 'units', 'degrees_north')
%         ncwriteatt(cms_tname, 'latt', 'valid_range', lat_u_rho_range)        
%         ncwriteatt(cms_tname, 'Depth', 'units', 'meter')
%         ncwriteatt(cms_tname, 'Depth', 'positive', 'down')
%         ncwriteatt(cms_tname, 'Depth', 'valid_range', depth_range)
%         ncwriteatt(cms_tname, 'Time', 'units', 'hours since midnight January 1 2019')
%         ncwriteatt(cms_tname, 't', 'units', 'Celcius')
%         ncwriteatt(cms_tname, 't', 'valid_range', temp_range) %% NOTE: not sure if this is necessary / if I'm doing it correctly      
%         ncwriteatt(cms_tname, 't', 'long_name', 'temperature')
%         %s
%         ncwriteatt(cms_sname, 'lons', 'units', 'degrees_east')
%         ncwriteatt(cms_sname, 'lons', 'valid_range', lon_v_rho_range)
%         ncwriteatt(cms_sname, 'lats', 'units', 'degrees_north')
%         ncwriteatt(cms_sname, 'lats', 'valid_range', lat_u_rho_range)        
%         ncwriteatt(cms_sname, 'Depth', 'units', 'meter')
%         ncwriteatt(cms_sname, 'Depth', 'positive', 'down')
%         ncwriteatt(cms_sname, 'Depth', 'valid_range', depth_range)
%         ncwriteatt(cms_sname, 'Time', 'units', 'hours since midnight January 1 2019')
%         ncwriteatt(cms_sname, 's', 'units', 'PSU')
%         ncwriteatt(cms_sname, 's', 'valid_range', salt_range) %% NOTE: not sure if this is necessary / if I'm doing it correctly      
%         ncwriteatt(cms_sname, 's', 'long_name', 'salinity')
%     end
% end