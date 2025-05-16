%% Script to convert netCDF file output file from ROMS to a netCDF file that CMS can use
% NOTE: The ROMS output from Jim Hench and Walter Torres is grid C. When
% using grid C with the CMS, you need to provide a ROMS file for each
% variable. For example: nest_1_yyymmddhhmmssu.nc;
% nest_1_yyymmddhhmmssv.nc, nest_1_yyymmddhhmmssw.nc, etc.
%
% For each variable, you need to define different variable names for lon,
% lat, etc., set in nest_x.nml: lon_nameU, lon_nameY, etc.
%
% CMS can also handle multiple snapshots in a file, apparently. I don't
% think this will be necessary. We can just use two files and set them a
% month apart... 
%
% Some information about CMS nest files:
% Global attrivutes:
% theta_s; theta_b; hc
%
% Variables:
% Longitude
% Latitude
% Depth (this we may have to redo using convert_sigma_z)
% u
% v
% w
% h - how do we do these? Are we sure we need it? Pretty sure NaN in the
% velocity fields does the trick here... 
% zeta - how do we do these? Do we need zeta? Don't think so.

clear

% Specify files name of ROMS output file and name of new file to be created in CMS format
% ROMS_fname = 'ocean_his_stokes.nc';
ROMS_grid = 'sttstj300m.nc';
ROMS_hydro = 'nest_his.02604.nc';
CMS_fname = 'nest_1_20010101000000'; % BF - update as needed
cms_uname = strcat(CMS_fname,'u.nc');
cms_vname = strcat(CMS_fname,'v.nc');
cms_wname = strcat(CMS_fname,'w.nc');

depth = [.5 .75 1 1.5 2.5 3.5 5 6.5 8.5 10.5 13 16 20 25 30 40 50 60 70 80 90 100 ...
    120 140 160 180 200 220 240 260 280 300]; % BF - see convert_sigma_z.m; want to specify actual depths being used
 
% Read in u,v,w and stokes velocities
u_euler      = ncread(ROMS_hydro, 'u'); % This is usually just 'u'
v_euler      = ncread (ROMS_hydro, 'v'); % This is usuall just 'v'
w_euler      = ncread (ROMS_hydro, 'w');
u_stokes     = ncread (ROMS_hydro, 'u_stokes'); % BF - this likely needs to change
v_stokes     = ncread (ROMS_hydro, 'v_stokes');
w_stokes     = ncread (ROMS_hydro, 'w_stokes');

% Read in x and y UTM
x_w_UTM        = ncread (ROMS_grid, 'x_rho');
y_w_UTM        = ncread (ROMS_grid, 'y_rho');
x_u_UTM        = ncread (ROMS_grid, 'x_u');
y_u_UTM        = ncread (ROMS_grid, 'y_u');
x_v_UTM        = ncread (ROMS_grid, 'x_v');
y_v_UTM        = ncread (ROMS_grid, 'y_v');

% Add euler and stokes velocities
% BF - we will do something different
u = u_euler + u_stokes; % These don't line up. u_eular and u_stokes may be on different grids...
v = v_euler + v_stokes; % These don't line up. v_eular and v_stokes may be on different grids...
w = w_euler + w_stokes; % This does work. Makes me think the stokes vel for u and v are on the rho grid...

% Convert meters to decimal degrees USING COARSE ESTIMATE OF 111KM/degree
x_w = x_w_UTM / 111000;
y_w = y_w_UTM / 111000;
x_u = x_u_UTM / 111000;
y_u = y_u_UTM / 111000;
x_v = x_v_UTM / 111000;
y_v = y_v_UTM / 111000;

% Create new netcdfs and empty variables with the correct dims
nccreate(cms_uname,'lonu',...
    'Dimensions',{'Lon',length(x_u(:,1))},'Datatype','double');
nccreate(cms_uname,'latu',...
    'Dimensions',{'Lat',length(y_u(1,:))},'Datatype','double');
nccreate(cms_uname,'Depth',...
    'Dimensions',{'Depth',length(depth)},'Datatype','double');
nccreate(cms_uname,'Time',...
    'Dimensions',{'Time',1},'Datatype','double');
nccreate(cms_uname,'u',...
    'Dimensions',{'Lon','Lat','Depth','Time'},'Datatype','double');

nccreate(cms_vname,'lonv',...
    'Dimensions',{'Lon',length(x_v(:,1))},'Datatype','double');
nccreate(cms_vname,'latv',...
    'Dimensions',{'Lat',length(y_v(1,:))},'Datatype','double');
nccreate(cms_vname,'Depth',...
    'Dimensions',{'Depth',length(depth)},'Datatype','double');
nccreate(cms_vname,'Time',...
    'Dimensions',{'Time',1},'Datatype','double');
nccreate(cms_vname,'v',...
    'Dimensions',{'Lon','Lat','Depth','Time'},'Datatype','double');

nccreate(cms_wname,'lonw',...
    'Dimensions',{'Lon',length(x_w(:,1))},'Datatype','double');
nccreate(cms_wname,'latw',...
    'Dimensions',{'Lat',length(y_w(1,:))},'Datatype','double');
nccreate(cms_wname,'Depth',...
    'Dimensions',{'Depth',length(depth)},'Datatype','double');
nccreate(cms_wname,'Time',...
    'Dimensions',{'Time',1},'Datatype','double');
nccreate(cms_wname,'w',...
    'Dimensions',{'Lon','Lat','Depth','Time'},'Datatype','double');

% Need to move the variables theta_s, theta_b and hc to GLOBAL ATTRIBUTES
% for CMS

theta_s = ncread(ROMS_grid, 'theta_s');
theta_b = ncread(ROMS_grid, 'theta_b');
hc = ncread(ROMS_grid, 'hc');
h = ncread(ROMS_grid, 'h');
zeta = ncread(ROMS_hydro, 'zeta'); % BF - had to grab this from hydro

ncwriteatt(cms_uname, '/', 'theta_s', theta_s);
ncwriteatt(cms_vname, '/', 'theta_s', theta_s);
ncwriteatt(cms_wname, '/', 'theta_s', theta_s);
ncwriteatt(cms_uname, '/', 'theta_b', theta_b);
ncwriteatt(cms_vname, '/', 'theta_b', theta_b);
ncwriteatt(cms_wname, '/', 'theta_b', theta_b);
ncwriteatt(cms_uname, '/', 'hc', hc);
ncwriteatt(cms_vname, '/', 'hc', hc);
ncwriteatt(cms_wname, '/', 'hc', hc);

ncwrite(cms_uname,'lonu',x_u(:,1));
ncwrite(cms_uname,'latu',y_u(1,:));
ncwrite(cms_uname,'Depth',depth); % Must update

ncwrite(cms_vname,'lonv',x_v(:,1));
ncwrite(cms_vname,'latv',y_v(1,:));
ncwrite(cms_vname,'Depth',depth); % Must update

ncwrite(cms_wname,'lonw',x_w(:,1));
ncwrite(cms_wname,'latw',y_w(1,:));
ncwrite(cms_wname,'Depth',depth); % Must update

% Convert depth
u_new = convert_sigma_z_sep(u(:,:,:,4), theta_s, theta_b, hc, h, zeta, depth, 1:10, 3); % BF - I think rho-points (sigma-levels) need to be in here somewhere

ncwrite(cms_uname, 'u', u_new);
ncwrite(cms_vname, 'v', v_new);
ncwrite(cms_wname, 'w', w_new); % This has one more depth level... WILL NEED TO UPDATE THIS



% 
% % 2D depth averaged velocities
% ubar = ubar_euler + ubar_stokes;
% vbar = vbar_euler + vbar_stokes;
% 
% % __
% 
% % lon_u/lat_u, lon_v/lat_v, etc. matricies are different dimensions and locations
% % due to computational stencil. for now, will trim u so that it matches
% % dimensions of v. need to sort this out properly at some point ...
% lon_u_vec = lon_u (1:end, 1);
% lon_v_vec = lon_v (1:387, 1);
% lat_u_vec = lat_u (1,1:321);
% lat_v_vec = lat_v (1,1:end);
% 
% zu = ones  (387, 321, 1, 1);
% zv = ones  (387, 321, 1, 1);
% zw = ones  (387, 321, 1, 1);
% zu (:, :, 1, 1) = ubar (1:end, 1:321);
% zv (:, :, 1, 1) = vbar (1:387, 1:end);
% zw = zeros (387, 321, 1, 1);   % assume vertical velocity is zero for now ...
% 
% % __
% 
% % Specify ranges for attributes in CMS netCDF file
% 
% lon_u_vec_range = minmax ( lon_u_vec');
% lat_u_vec_range = minmax ( lat_u_vec);
% 
% % For now, use 2D depth-averaged velocities, so just one depth level
% depth_levels = [1]; % Dan changed from 1
% depth_level_range = [0 10]; % Dan changed from [0 10]
% 
% % Specify velocity ranges; vertical velocities hard-wired as zero for now
% u_vel_range = [min(min(zu)) max(max(zu))];
% v_vel_range = [min(min(zv)) max(max(zv))];
% w_vel_range = [ 0 0];
% 
% % Specify current time as days since 1900-12-31 
% % t0 = datenum (1900, 12, 31, 23, 59, 59); 
% % tc = datenum (2016, 10, 10, 10, 10, 10);
% time_0 = juliandate(datetime(2010,1,1)) - juliandate(datetime(1900,12,31));
% 
% % Specify grid dimensions; hard-wired for Moorea north shore grid
% nrow  = 387;
% ncol  = 321;
% ndep  =   1;
% ntime =   1;
% 
% % __

% In ROMS masked values are set equal to NaN; switch to mask value used in HYCOM example file 

fill_value = single (1.267650600228229e+30);

for i=1:387
    for j=1:321
        k = isnan ( zu (i,j) );
        if k > 0
           zu (i,j) = fill_value;
           zv (i,j) = fill_value;
           zw (i,j) = fill_value;
        end
    end
end

%pcolor (zv');   shading flat; colorbar; caxis ([-0.2 0.2]);

% __

% Open netCDF file to write in CMS format 
%   note: HYCOM example is 'Format: 64bit', but this causes some problems with the 'netcdf.defVarFill' function.
%         Will just use 'netcdf4' format for now and see if this is compatible with CMS.
%ncid = netcdf.create(['roms_cms_01.nc'], '64BIT_OFFSET');
ncid = netcdf.create(CMS_fname, 'netcdf4');

% Specify global attributes for CMS file 
varid0 = netcdf.getConstant('GLOBAL');
netcdf.putAtt (ncid, varid0, 'source', ROMS_fname);
 
% Specify dimensions for variables to be packed into CMS file
dimid_row = netcdf.defDim (ncid, 'Longitude', nrow);
dimid_col = netcdf.defDim (ncid, 'Latitude',  ncol);
dimid_dep = netcdf.defDim (ncid, 'Depth',     ndep);
dimid_tim = netcdf.defDim (ncid, 'Time',     ntime);

% Specify variables to pack into CMS file (declaring as single precision for now, but may need to use double ...)
varid1 = netcdf.defVar (ncid, 'Longitude', 'NC_FLOAT', [dimid_row]);
varid2 = netcdf.defVar (ncid, 'Latitude',  'NC_FLOAT', [dimid_col]);
varid3 = netcdf.defVar (ncid, 'Depth',     'NC_FLOAT', [dimid_dep]);
varid4 = netcdf.defVar (ncid, 'Time',      'NC_FLOAT', [dimid_tim]);
varid5 = netcdf.defVar (ncid, 'zu',        'NC_FLOAT', [dimid_row dimid_col dimid_dep dimid_tim]);
varid6 = netcdf.defVar (ncid, 'zv',        'NC_FLOAT', [dimid_row dimid_col dimid_dep dimid_tim]);
varid7 = netcdf.defVar (ncid, 'zw',        'NC_FLOAT', [dimid_row dimid_col dimid_dep dimid_tim]);

% Assign attributes to variables
%   Longitude
netcdf.putAtt (ncid, varid1, 'units',      'degrees_east');
netcdf.putAtt (ncid, varid1, 'valid_range', lon_u_vec_range);

%   Latitude
netcdf.putAtt (ncid, varid2, 'units',      'degrees_north');
netcdf.putAtt (ncid, varid2, 'valid_range', lat_u_vec_range);

%   Depth levels 
netcdf.putAtt (ncid, varid3, 'units',      'meter');
netcdf.putAtt (ncid, varid3, 'valid_range', depth_level_range);
netcdf.putAtt (ncid, varid3, 'positive',      'down');

%   Time
netcdf.putAtt (ncid, varid4, 'units',      'days since 1900-12-31');
netcdf.putAtt (ncid, varid4, 'calender',   'standard');

%   U-velocities 
netcdf.putAtt     (ncid, varid5, 'long_name',  'eastward_sea_water_velocity');
netcdf.putAtt     (ncid, varid5, 'units',      'm/s');
netcdf.defVarFill (ncid, varid5,  false,        fill_value);
netcdf.putAtt     (ncid, varid5, 'valid_range', u_vel_range);

%   V-velocities 
netcdf.putAtt     (ncid, varid6, 'long_name',  'northward_sea_water_velocity');
netcdf.putAtt     (ncid, varid6, 'units',      'm/s');
netcdf.defVarFill (ncid, varid6,  false,        fill_value);
netcdf.putAtt     (ncid, varid6, 'valid_range', v_vel_range);

%   W-velocities 
netcdf.putAtt     (ncid, varid7, 'long_name',  'downward_sea_water_velocity');
netcdf.putAtt     (ncid, varid7, 'units',      'm/s');
netcdf.defVarFill (ncid, varid7,  false,        fill_value);
netcdf.putAtt     (ncid, varid7, 'valid_range', w_vel_range);

% Done setting up CMS file
netcdf.endDef (ncid);

% __

% Write data to CMS file
netcdf.putVar (ncid, varid1, lon_u_vec);
netcdf.putVar (ncid, varid2, lat_u_vec);
netcdf.putVar (ncid, varid3, depth_levels);
netcdf.putVar (ncid, varid4, time_0);
netcdf.putVar (ncid, varid5, zu);
netcdf.putVar (ncid, varid6, zv);
netcdf.putVar (ncid, varid7, zw);

% Close up the shop
netcdf.close (ncid);

% __
% __

% Check new CMS file contents
ncdisp (CMS_fname);

cms_u  = ncread (CMS_fname, 'zu');
cms_v  = ncread (CMS_fname, 'zv');

subplot (211)
pcolor (zv');
shading flat;
colorbar;
caxis ([-0.2 0.2]);
axis equal
axis tight6
title ('ROMS v-vel')

subplot (212)
pcolor (cms_v');
shading flat;
colorbar;
caxis ([-0.2 0.2]);
axis equal
axis tight
title ('CMS v-vel')

