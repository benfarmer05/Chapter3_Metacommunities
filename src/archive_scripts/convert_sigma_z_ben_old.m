function [u_new, v_new, w_new] = convert_sigma_z(nest_name, Z) % Inputs will change
%% Convert ROMS daily nc files from sigma-levels to z-levels
% Requires set_depth.m and stretching.m - I found them online. NOTE: This
% script outputs NaN (benthos mask) in u, v, and w as the FILL VALUE in the
% nc file.

% UPDATES FOR BEN FARMER AND THE UPDATED USVI HYDRODYNAMICS 2/18/2022
cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/VI_Hydro-Connectivity-Modeling
nest_name = 'nest_his.02604.nc';

% Set z levels
zlevels = Z;

% Read in global attributes
% THESE ARE IN THE ORIGINAL NEST AND/OR GRID FILES
% YOU MAY NEED TO CHANGE SOME OF THESE VARIABLE NAMES

% We will be running this in a loop, and may have to bring these variables
% in explicitly 

%ncreadatt retrieves global attributes from each hydro file
theta_s = ncreadatt(nest_name, '/', 'theta_s'); %S-coordinate surface control parameter (theta)
theta_b = ncreadatt(nest_name, '/', 'theta_b'); %S-coordinate bottom control parameter (b)
hc = ncreadatt(nest_name, '/', 'hc'); %S-coordinate parameter, critical depth

% Read in variables
h = ncread(nest_name, 'h'); %BF - I think this comes from roms2cms.m
h = ncread(nest_name, 'depths'); %BF - this is a temporary solution
zeta = ncread(nest_name, 'zeta');
D = ncread(nest_name, 'Depth'); %BF - will I need to generate this? This should be the depth levels
% u_old = ncread(nest_name, 'zu');
% v_old = ncread(nest_name, 'zv');
% w_old = ncread(nest_name, 'zw');
u_old = ncread(nest_name, 'u'); % BF - are these the same as zu zv zw?
v_old = ncread(nest_name, 'v');
w_old = ncread(nest_name, 'w');

% Read in fill values
% zuFill = ncreadatt(nest_name,'zu','_FillValue'); % BF - no fill values present
% zvFill = ncreadatt(nest_name,'zv','_FillValue');
% zwFill = ncreadatt(nest_name,'zw','_FillValue');

% u_old(u_old == zuFill) = NaN;
% v_old(v_old == zvFill) = NaN;
% w_old(w_old == zwFill) = NaN;

if size(D,1) ~= size(zlevels,2) % BF - I believe 'D' here is sigma levels, and 'zlevels' is depth levels (m) at which other u,v,w are interpolated
    disp('This script currently only works when the number of zlevels = the number of sigma levels; Exiting...')
    return
end

% NEED TO UPDATE THIS - WE WILL BE INTERPOLATING AT C-GRID: U, V, & W
% The i-grid input of set_depth() allows for indicating what type of point
% you're interpolating at. So we will need to run this loop at least 3
% times for u, v and w. The result will be something like new_z_u, new_z_v,
% new_z_w....
new_z = set_depth(2,4,theta_s, theta_b, hc, size(D,1), 1,h,zeta,1); % Puts z's at sigma levels -- CHANGE NUMBER ('N' and 'igrid' may need to change)
for i = 1:size(h,1)
    for j = 1:size(h,2)
        new_z(i,j,:) = sort(-new_z(i,j,:)); % Reverse order
    end
end

new_z = double(new_z);

% BF edits. 22 Feb 2022
new_u = set_depth(2,4,theta_s,theta_b,hc,32,3,h,zeta,1); %changing 'igrid' here to match relevant velocity point ('u' in this case)
for i = 1:size(new_u,1) %this was originally indexing off of 'h', not 'new_u' - might cause problems. did this because lons (x, or xi) were 313 instead of 314
    for j = 1:size(new_u,2)
        new_u(i,j,:) = sort(-new_u(i,j,:));
    end
end % BF - error: Index in position 1 exceeds array bounds (must not exceed 313).

new_v = set_depth(2,4,theta_s,theta_b,hc,32,4,h,zeta,1);
for i = 1:size(new_v,1)
    for j = 1:size(new_v,2)
        new_v(i,j,:) = sort(-new_v(i,j,:));
    end
end

new_w = set_depth(2,4,theta_s,theta_b,hc,32,5,h,zeta,1); %Error in set_depth (line 259) ... z(:,:,k)=zetar+(zetar+hr).*z0;
for i = 1:size(new_w,1)
    for j = 1:size(new_w,2)
        new_w(i,j,:) = sort(-new_w(i,j,:));
    end
end



% Interpolate u, v, and w at new z-positions
%U
Z = 32; % BF - added this, levels should be something different than this though
u_new = zeros(size(h,1), size(h,2), size(zlevels,2));
u_new(u_new == 0) = NaN;


% REPLACE new_z with appropriate new_z_'s
for i = 1:size(h,1)
    for j = 1:size(h,2)
        %u_new(i,j,:) = interp1(test_z(i,j,:),squeeze(u_z(i,j,:)),zlevels); % Not sure what Z is - I think these depths are all equal
        vec = interp1(new_z(i,j,:), squeeze(u_old(i,j,:)),zlevels,'linear','extrap'); % solves the 0 level problem, but now no bottom.
        for k = 1:length(zlevels)
            if zlevels(k) > h(i,j)
                vec(k) = NaN;
            end
        end
        u_new(i,j,:) = vec;
    end
end 

%V
v_new = zeros(size(h,1), size(h,2), size(zlevels,2));
v_new(v_new == 0) = NaN;

for i = 1:size(h,1)
    for j = 1:size(h,2)
        %u_new(i,j,:) = interp1(test_z(i,j,:),squeeze(u_z(i,j,:)),zlevels); % Not sure what Z is - I think these depths are all equal
        vec = interp1(new_z(i,j,:), squeeze(v_old(i,j,:)),zlevels,'linear','extrap');
        for k = 1:length(zlevels)
            if zlevels(k) > h(i,j)
                vec(k) = NaN;
            end
        end
        v_new(i,j,:) = vec;
    end
end 

%W
w_new = zeros(size(h,1), size(h,2), size(zlevels,2));
w_new(w_new == 0) = NaN;

for i = 1:size(h,1)
    for j = 1:size(h,2)
        %u_new(i,j,:) = interp1(test_z(i,j,:),squeeze(u_z(i,j,:)),zlevels); % Not sure what Z is - I think these depths are all equal
        vec = interp1(new_z(i,j,:), squeeze(w_old(i,j,:)),zlevels,'linear','extrap');
        for k = 1:length(zlevels)
            if zlevels(k) > h(i,j)
                vec(k) = NaN;
            end
        end
        w_new(i,j,:) = vec;
    end
end 

% u_new(isnan(u_new)) = zuFill;
% v_new(isnan(v_new)) = zvFill;
% w_new(isnan(w_new)) = zwFill;
% figure
% pcolor(u_new(:,:,1)')
% 
% figure
% pcolor(u_new(:,:,4)')
% 
% figure
% pcolor(u_new(:,:,15)')
% 
% figure
% pcolor(u_new(:,:,25)')

ncwrite(nest_name, 'zu', u_new); 
ncwrite(nest_name, 'zv', v_new);
ncwrite(nest_name, 'zw', w_new);
ncwrite(nest_name, 'Depth', zlevels');
ncwriteatt(nest_name, '/', 'Edited', 'Sigma to Z');
end

