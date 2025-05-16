

%% May be veru useful!: https://romsmatlab.tiddlyspot.com/
% John Wilkins made this with the help of Gordon Zhang from WHOI

% 1.) make sure the w-velocity interpolation isn't going haywire. CHECK
% 2.) fix the "valid range" issues. CHECK
% 3.) create a release file. SHIT
% 4.) test this bad boy out on the CMS
% 5.) maybe increase depth range down to 300 m?
% 6.) does netcdf even work anymore?? figure out by trying a job, etc.

clear;clc

%% tests
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed/Testing'
% examplenest = 'nest_1_20190630180000w.nc';
% examplenest = 'nest_2_20190701010000t.nc';
% exampleMAR = 'example_MAR_nest.nc';
% clc;ncdisp(examplenest);ncdisp(exampleMAR)
% MARtime = ncread(exampleMAR, 'Time');MARtime
% MARdepth = ncread(exampleMAR, 'Depth');MARdepth
% MARnestU = ncread(exampleMAR, 'zu');

% the thang
% examplenestu = 'nest_2_20190701210000u.nc';
% examplenestv = 'nest_2_20190701210000v.nc';
% examplenestw = 'nest_2_20190701210000w.nc';
% examplenestt = 'nest_2_20190701210000t.nc';
% examplenests = 'nest_2_20190701210000s.nc';

examplenestu = 'nest_1_20190221000000u.nc';
examplenestv = 'nest_1_20190221000000v.nc';
examplenestw = 'nest_1_20190221000000w.nc';
examplenestt = 'nest_1_20190221000000t.nc';
examplenests = 'nest_1_20190221000000s.nc';
ncdisp(examplenestt)

varu = ncread(examplenestu, 'u');
varv = ncread(examplenestv, 'v');
varw = ncread(examplenestw, 'w');
vart = ncread(examplenestt, 't');
vars = ncread(examplenests, 's');

lon_u = ncread(examplenestu, 'lonu');
lat_u = ncread(examplenestu, 'latu');
lon_v = ncread(examplenestv, 'lonv');
lat_v = ncread(examplenestv, 'latv');
lon_rho = ncread(examplenestw, 'lonw'); %w-vel, temp, salinity, zeta, and rho use rho-grid
lat_rho = ncread(examplenestw, 'latw'); %w-vel, temp, salinity, zeta, and rho use rho-grid

%%
% ROMS_examplenest = 'croco_his.04365.nc'; %nest 1; large one; 2km
% ROMS_examplenest = 'nest_his.04365.nc'; %nest 2; small one; 330m
ROMS_examplenest = 'croco_his.01224.nc';

uvel = ncread(ROMS_examplenest, 'u');
vvel = ncread(ROMS_examplenest, 'v');
wvel = ncread(ROMS_examplenest, 'w');
salt = ncread(ROMS_examplenest,'salt');
temp = ncread(ROMS_examplenest,'temp');
uvel = uvel(:,:,:,3);
vvel = vvel(:,:,:,3);
wvel = wvel(:,:,:,3);
temp = temp(:,:,:,3);
salt = salt(:,:,:,3);
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
ROMS_grid = 'vigrid.nc';
lon_u_ROMS = ncread(ROMS_grid, 'lon_u');
lat_u_ROMS = ncread(ROMS_grid, 'lat_u');
lon_v_ROMS = ncread(ROMS_grid, 'lon_v');
lat_v_ROMS = ncread(ROMS_grid, 'lat_v');
lon_rho_ROMS = ncread(ROMS_grid, 'lon_rho'); %w-vel, temp, salinity, zeta, and rho use rho-grid
lat_rho_ROMS = ncread(ROMS_grid, 'lat_rho'); %w-vel, temp, salinity, zeta, and rho use rho-grid

%% release points
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/Release';     
% [relpnts] = readmatrix('points.xlsx');
[relpnts] = readmatrix('points_650_none-on-land.csv');
plot(relpnts(:,2)+360,relpnts(:,3), '*r')

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

clc
max(max(max(varv)))
min(min(min(varv)))
max(max(max(vvel)))
min(min(min(vvel)))

clc
max(max(max(varw)))
min(min(min(varw)))
max(max(max(wvel)))
min(min(min(wvel)))

clc
max(max(max(vart)))
min(min(min(vart)))
max(max(max(temp)))
min(min(min(temp)))

clc
max(max(max(vars)))
min(min(min(vars)))
max(max(max(salt)))
min(min(min(salt)))


%%
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023'

%lon/lat and depth matrices
lons = ncread(ROMS_grid, 'lon_rho');
lats = ncread(ROMS_grid, 'lat_rho');
mask = ncread(ROMS_grid, 'mask_rho');

depths = ncread(ROMS_grid, 'h');
depthcut1 = depths;
depthcut1(depthcut1==0) = NaN;
% depthcut(depthcut > 40) = NaN;
depthcut1(mask == 0) = NaN;
figure();surf(lons', lats', -depthcut1');colorbar;%caxis([-160 0])
% title 'Bathymetry in the Virgin Islands';hold on
% cb = colorbar;ylabel(cb,'Depth (m)')
figure();imagescn(lons', lats', -depthcut1');set(gca, 'ydir', 'normal');hold on;colorbar;axis equal; clim([-15 -1]); hold on
% figure();imagescn(lons', lats', -depthcut1');set(gca, 'ydir', 'normal');hold on;colorbar;axis equal; clim([-3000 -1])
% figure();imagescn(lons', lats', -depthcut1');set(gca, 'ydir', 'normal');hold on;colorbar;axis equal; clim([-7300 -1]) %shows the tiny lil sliver of Puerto Rico trench!
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/Release';     
% [relpnts] = readmatrix('points.xlsx');
[relpnts] = readmatrix('points_650_none-on-land.csv');
plot(relpnts(:,8),relpnts(:,9), '*r')

%%
%u
%interpolated
figure(1);imagescn(lon_u',lat_u',varu(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-1.3 1.3]);title interpolated
figure(1);imagescn(lon_u',lat_u',varu(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;axis equal;clim([-.6 .6]);
figure(1);imagescn(lon_u',lat_u',varu(:,:,72)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;axis equal;clim([-.6 .6]);
%original
figure(2);imagescn(lon_u_ROMS',lat_u_ROMS',uvel(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-1.3 1.3]);title ROMS
figure(2);imagescn(lon_u_ROMS',lat_u_ROMS',uvel(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;axis equal;clim([-.6 .6]);
figure(2);imagescn(lon_u_ROMS',lat_u_ROMS',uvel(:,:,72)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;axis equal;clim([-.6 .6]);

% %TEST
% [row, column] = find (varu <= -10)
% [row, column] = find (uvel >= .1)
% %TEST

%v
%interpolated
figure(1);imagescn(lon_v',lat_v',varv(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal; clim([-.6 .6]);
figure(1);imagescn(lon_v',lat_v',varv(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal; clim([-.6 .6]);
figure(1);imagescn(lon_v',lat_v',varv(:,:,32)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal; clim([-.6 .6]);
%original
figure(2);imagescn(lon_v_ROMS',lat_v_ROMS',vvel(:,:,5)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-.6 .6]);
figure(2);imagescn(lon_v_ROMS',lat_v_ROMS',vvel(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal; clim([-.6 .6]);
figure(2);imagescn(lon_v_ROMS',lat_v_ROMS',vvel(:,:,32)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal; clim([-.6 .6]);

%w
%interpolated
figure(1);imagescn(lon_rho',lat_rho',varw(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean curl;axis equal; clim([-.07 .07]);
figure(1);imagescn(lon_rho',lat_rho',varw(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal; clim([-.05 .05]);
figure(1);imagescn(lon_rho',lat_rho',varw(:,:,32)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal; clim([-.2 .2]);
%original
figure(2);imagescn(lon_rho_ROMS',lat_rho_ROMS',wvel(:,:,2)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean curl;axis equal; clim([-.07 .07]);
figure(2);imagescn(lon_rho_ROMS',lat_rho_ROMS',wvel(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal; clim([-.05 .05]);
figure(2);imagescn(lon_rho_ROMS',lat_rho_ROMS',wvel(:,:,32)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal; clim([-.05 .05]);

%t
%interpolated
figure(1);imagescn(lon_rho',lat_rho',vart(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean thermal;%axis equal;  caxis([25 30]);
figure(1);imagescn(lon_rho',lat_rho',vart(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean thermal;axis equal; clim([0 30]);
figure(1);imagescn(lon_rho',lat_rho',vart(:,:,32)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean thermal;axis equal; clim([0 10]);
%original
figure(2);imagescn(lon_rho_ROMS',lat_rho_ROMS',temp(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean thermal;%axis equal;  caxis([25 30]);
figure(2);imagescn(lon_rho_ROMS',lat_rho_ROMS',temp(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean thermal;axis equal; clim([0 30]);
figure(2);imagescn(lon_rho_ROMS',lat_rho_ROMS',temp(:,:,32)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean thermal;axis equal; clim([0 10]);

%s
%interpolated
figure(1);imagescn(lon_rho',lat_rho',vars(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean haline;axis equal;clim([34 38])
figure(1);imagescn(lon_rho',lat_rho',vars(:,:,10)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean halinea;axis equal;clim([34 38])
figure(1);imagescn(lon_rho',lat_rho',vars(:,:,32)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean haline;axis equal;clim([34 38])
%original
figure(2);imagescn(lon_rho_ROMS',lat_rho_ROMS',salt(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean haline;axis equal;clim([34 38])
figure(2);imagescn(lon_rho_ROMS',lat_rho_ROMS',salt(:,:,10)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean haline;axis equal;clim([34 38])
figure(2);imagescn(lon_rho_ROMS',lat_rho_ROMS',salt(:,:,32)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean haline;axis equal;clim([34 38])



%% Testing different variable naming conventions for very small CMS run. C-grid was running into issues (24 August 2023)

clc
cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed/Testing
tester = 'smallchungus.nc';
ncdisp(tester)

% Specify the input and output file names
inputFile = tester; % Replace with the name of your input NetCDF4 file
outputFile = 'bigchungus.nc'; % Replace with the desired name of your output 64-bit NetCDF file

% Open the input NetCDF4 file
ncid = netcdf.open(inputFile, 'NOWRITE');

% Get information about the input file
[ndims, nvars, ngatts, unlimdimid] = netcdf.inq(ncid);

% Create the output NetCDF file
ncid_out = netcdf.create(outputFile, 'CLOBBER');

% Define dimensions in the output file (copy from the input file)
for dimid = 0:(ndims - 1)
    [name, len] = netcdf.inqDim(ncid, dimid);
    netcdf.defDim(ncid_out, name, len);
end

% Define variables in the output file (copy from the input file)
for varid = 0:(nvars - 1)
    [varname, xtype, dimids, natts] = netcdf.inqVar(ncid, varid);
    netcdf.defVar(ncid_out, varname, xtype, dimids);
    
    % Copy variable attributes
    for attnum = 0:(natts - 1)
        attname = netcdf.inqAttName(ncid, varid, attnum);
        attvalue = netcdf.getAtt(ncid, varid, attname);
        netcdf.putAtt(ncid_out, varid, attname, attvalue);
    end
end

% End definitions in the output file
netcdf.endDef(ncid_out);

% Copy variable data from input to output
for varid = 0:(nvars - 1)
    varname = netcdf.inqVar(ncid, varid);
    data = netcdf.getVar(ncid, varid);
    netcdf.putVar(ncid_out, varid, data);
end

% Close both input and output files
netcdf.close(ncid);
netcdf.close(ncid_out);

bigchungus = 'bigchungus.nc';
ncdisp(bigchungus)

%FUNCTIONING MAR HYDRO
MARtester = 'nest_1_20140301000000.nc';
ncdisp(MARtester)

fill_value = 1.267650600228229e+30; % NOTE - could try 1.267650600228229e+30 (from MAR Dan run, March 2023 on HPC. CMS one: 1.2676506e+30

varu = ncread(MARtester, 'zu');
yeet = varu(:,:,1);
max(max(yeet))
min(min(yeet))

balls = varu(isnan(varu));
balls2 = varu(varu==fill_value);




GOMtester = 'nest_1_20100101000000.nc';
ncdisp(GOMtester)

testerdepth = ncread(tester, 'Depth');
MARdepth = ncread(MARtester, 'Depth');

testerw = ncread(tester, 'w');
testerw32 = testerw(:,:,32);
testerw1 = testerw(:,:,1);

MARw = ncread(MARtester, 'zw');
MARw21 = MARw(:,:,21);
MARw1 = MARw(:,:,1);


netcdf.renameVar(tester,'lons','longitude')

ncwrite(filename,varname,vardata)

ncwrite(tester,'lons',lon_u_convert(:,1))




% Open netCDF file.
ncid = netcdf.open('my_example.nc','NC_WRITE');

% Put file in define mode.
netcdf.redef(ncid)

% Get name of first variable
[varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,0);

avagadros_number

% Rename the variable, using a capital letter to start the name.
netcdf.renameVar(ncid,0,'Avagadros_number')


%% testing traj files

%25 aug 2023: there was an issue with the data type or fill value mismatch
cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed/traj_file_troubleshoot
trajtest = 'traj_file_01.nc';
ncdisp(trajtest)

%% extracting the grid dimensions (grid extent)
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023'
clear;clc

% ROMS_grid = 'parent2km100.nc';
ROMS_grid = 'vigrid.nc';
ncdisp(ROMS_grid)
lon_u = ncread(ROMS_grid, 'lon_u');
lat_u = ncread(ROMS_grid, 'lat_u');
lon_v = ncread(ROMS_grid, 'lon_v');
lat_v = ncread(ROMS_grid, 'lat_v');
lon_rho = ncread(ROMS_grid, 'lon_rho'); %w-vel, temp, salinity, zeta, and rho use rho-grid
lat_rho = ncread(ROMS_grid, 'lat_rho'); %w-vel, temp, salinity, zeta, and rho use rho-grid

min_lon_u = min(min(lon_u));
min_lon_v = min(min(lon_v));
min_lon_rho = min(min(lon_rho));
xstart = min([min_lon_u min_lon_v min_lon_rho]); % +360

max_lon_u = max(max(lon_u));
max_lon_v = max(max(lon_v));
max_lon_rho = max(max(lon_rho));
xend = max([max_lon_u max_lon_v max_lon_rho]); % +360

min_lat_u = min(min(lat_u));
min_lat_v = min(min(lat_v));
min_lat_rho = min(min(lat_rho));
ystart = min([min_lat_u min_lat_v min_lat_rho]);

max_lat_u = max(max(lat_u));
max_lat_v = max(max(lat_v));
max_lat_rho = max(max(lat_rho));
yend = max([max_lat_u max_lat_v max_lat_rho]);


