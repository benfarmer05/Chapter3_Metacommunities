
clear;clc %clears all environment variables and command window

%% tests
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed/Testing'

% examplenest2 = 'nest_1_20190221000000.nc';
examplenest2 = 'nest_1_20021129080000.nc';
ncdisp(examplenest2)

varu = ncread(examplenest2, 'u');
varv = ncread(examplenest2, 'v');
varw = ncread(examplenest2, 'w_velocity');
vart = ncread(examplenest2, 'water_temp');
vars = ncread(examplenest2, 'salinity');
time = ncread(examplenest2, 'MT');

lon = ncread(examplenest2, 'lon');
lat = ncread(examplenest2, 'lat');

%%

fill_value = single(1.2676506e+30);

%replace fill values with NaNs for plotting, if relevant
varu(varu==fill_value) = NaN;
varv(varv==fill_value) = NaN; %varv(varv==1.2676506e+30) = NaN;
varw(varw==fill_value) = NaN;
vart(vart==fill_value) = NaN;
vars(vars==fill_value) = NaN;

% vars(isnan(vars))

clc
max(max(max(varu)))
min(min(min(varu)))

clc
max(max(max(varv)))
min(min(min(varv)))

clc
max(max(max(varw)))
min(min(min(varw)))

clc
max(max(max(vart)))
min(min(min(vart)))

clc
max(max(max(vars)))
min(min(min(vars)))

%%
HYCOM_2D = '010_archv.2002_001_00_2d.nc';
ncdisp(HYCOM_2D)

depths = ncread(examplenest2, 'depth');

%%
%u
figure(1);imagescn(lon',lat',varu(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-2.2 2.2]);hold on;plot(266.17022, 27.88414, 'r.', 'markers', 8)
figure(2);imagescn(lon',lat',varu(:,:,10)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-2.2 2.2]);hold on;plot(266.17022, 27.88414, 'r.', 'markers', 8); plot(266.40233, 27.91967, 'b.', 'markers', 8)
figure(3);imagescn(lon',lat',varu(:,:,17)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-2.2 2.2]);hold on;plot(266.17022, 27.88414, 'r.', 'markers', 8); plot(266.40233, 27.91967, 'b.', 'markers', 8)
figure(4);imagescn(lon',lat',varu(:,:,23)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-2.2 2.2]);hold on;plot(266.17022, 27.88414, 'r.', 'markers', 8)

%v
figure(1);imagescn(lon',lat',varv(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;axis equal; clim([-.6 .6]);
figure(1);imagescn(lon',lat',varv(:,:,10)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;axis equal; clim([-.6 .6]);
figure(1);imagescn(lon',lat',varv(:,:,23)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;axis equal; clim([-.6 .6]);

%w
figure(1);imagescn(lon',lat',varw(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean balance;axis equal;clim([-.00009 .00009]);
figure(1);imagescn(lon',lat',varw(:,:,10)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;axis equal;clim([-7e-04 4e-04]);
figure(1);imagescn(lon',lat',varw(:,:,23)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;axis equal; clim([-.2 .2]);

%t
figure(1);imagescn(lon',lat',vart(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean thermal;axis equal;%clim([24 27]);
figure(1);imagescn(lon',lat',vart(:,:,10)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean thermal;axis equal; %caxis([-.015 .015]);
figure(1);imagescn(lon',lat',vart(:,:,23)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean thermal;hold on;axis equal;%clim([0 10])

%s
figure(1);imagescn(lon',lat',vars(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean haline;axis equal;  %caxis([-.6 .6]);
figure(1);imagescn(lon',lat',vars(:,:,3)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean haline;axis equal; %caxis([-.015 .015]);
figure(1);imagescn(lon',lat',vars(:,:,23)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean haline;axis equal; %caxis([-.015 .015]);

% %% Testing different variable naming conventions for very small CMS run. C-grid was running into issues (24 August 2023)
% 
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
% 
% 
% %% testing traj files
% 
% %25 aug 2023: there was an issue with the data type or fill value mismatch
% cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed/traj_file_troubleshoot
% trajtest = 'traj_file_01.nc';
% ncdisp(trajtest)
% 
% %% extracting the grid dimensions (grid extent)
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
% 
% 
