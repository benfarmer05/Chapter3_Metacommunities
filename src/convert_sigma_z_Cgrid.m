function [u_new,v_new,w_new,t_new,s_new] = convert_sigma_z(u,v,w,t,s,theta_s,theta_b,hc,h,zeta,zlevels,layers)
% %testing purposes
% u = uvel;
% v = vvel;
% w = wvel;
% t = temp;
% s = salt;
% zeta = ssh;
% %testing purposes

cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023'
numlevels = length(zlevels);

%% interpolate u,v,w,t,s at discrete depths on each of their respective grids

% % NOTE - may want to take this out actually, since it may be easier to
% % simply let set_depth think there is land where there isn't and then apply
% % the NaNs further below. interp1 can handle those NaNs ?
% %apply NaN landmask to bathymetry before passing to set_depth
% h(isnan(w(:,:,32))) = nan;

%convert convert sigma-level to z-level at every lon/lat for u,v,w,t,s
u_zlev = set_depth(2,4,theta_s,theta_b,hc,layers,3,h,zeta,0); %'3' indicates u-point grid
v_zlev = set_depth(2,4,theta_s,theta_b,hc,layers,4,h,zeta,0); %v-point grid
rho_zlev = set_depth(2,4,theta_s,theta_b,hc,layers,1,h,zeta,0); %rho-point grid

%flip the depth dimension so that it is ordered surface to bottom
u_zlev = u_zlev(:,:,layers:-1:1);
v_zlev = v_zlev(:,:,layers:-1:1);
rho_zlev = rho_zlev(:,:,layers:-1:1);

% NOTE - maybe come back to this. this will always avoid "making up" data
% past the last known point, but it also limits to what extent the output
% can mimic real values. I think letting the interpolation cook is really
% the best outcome, it always produces reasonable things from what I can
% tell. but check with the roms2cms_test script!
%set the first depth level to zero, to avoid extrapolating to the surface
u_zlev(:,:,1)=0;
v_zlev(:,:,1)=0;
rho_zlev(:,:,1)=0;

%%
% NOTE: I don't think doing interp2 really works  -
%   because it would be more about interpolating from one grid to another,
%   when in reality I don't want to change anything about the grid
%   resolution. interp1 with the current setup makes the most sense! Also,
%   Akima method is very similar to spline but I think spline is slightly
%   more elegant and I'll just stick with it. the interpolation is running
%   great - now just time to check up on the CMS and fill value issue! oh,
%   and check on the extrapolating to surface thing...update: I did, and
%   honestly "extrapolating" here makes sense

%%

% %
% figure(1);imagescn(lon_u',lat_u',u_zlev(:,:,2)');cb = colorbar; ylabel(cb,'depth (m)');cmocean deep;axis equal; %clim([-.6 .6]);
% figure(1);imagescn(lon_u',lat_u',u_zlev(:,:,10)');cb = colorbar; ylabel(cb,'depth (m)');cmocean deep;axis equal; %clim([-.6 .6]);
% figure(1);imagescn(lon_u',lat_u',u_zlev(:,:,32)');cb = colorbar; ylabel(cb,'depth (m)');cmocean deep;axis equal; %clim([-.6 .6]);
% %
% 
% %
% plot3(lon_rho, lat_rho, -h);zlim([-100 0])
% h(isnan(h))
% h(h==0)
% h(h<0)
% h(h>0)
% min(min(h))
% max(max(h))
% h(h==10)
% %
% 
% %
% figure(2);imagescn(lon_u',lat_u',u(:,:,4)');cb = colorbar; ylabel(cb,'depth (m)');cmocean deep;axis equal; %clim([-.6 .6]);
% figure(2);imagescn(lon_u',lat_u',u(:,:,10)');cb = colorbar; ylabel(cb,'depth (m)');cmocean deep;axis equal; %clim([-.6 .6]);
% figure(2);imagescn(lon_u',lat_u',u(:,:,32)');cb = colorbar; ylabel(cb,'depth (m)');cmocean deep;axis equal; %clim([-.6 .6]);
% %
% 
% 
% %
% x = 0:pi/4:2*pi; 
% v = sin(x);
% xq = 0:pi/16:2*pi;
% 
% figure
% vq1 = interp1(x,v,xq);
% plot(x,v,'o',xq,vq1,':.');
% xlim([0 2*pi]);
% title('(Default) Linear Interpolation');
% %
% 
% %
% for i = 1:200000 %200
%     dims = size(u_zlev);
%     row = randi([1 dims(1)]); %these are grid rows (oriented in latitudinal direction)
%     col = randi([1 dims(2)]); %these are grid columns (oriented in longitudinal direction)
%     x = abs(squeeze(u_zlev(row,col,:))); %the known depths
%     v = squeeze(u(row,col,:)); %the known values of interest (current speed, temp, salinity) at those depths
%     xq = zlevels; %the depths I want to interpolate to
% 
%     lonnies = lon_rho(row,col);
%     latties = lat_rho(row,col);
% 
%     if ~isnan(v(1))
%         vq1 = interp1(x,v,xq, 'linear'); %the interpolated values of interest
%         vq2 = interp1(x,v,xq, 'spline');
%         vq3 = makima(x,v,xq);
%     end
% 
%     % figure(1);plot(x,v,'o',xq,vq1,':.');xlim([0 200])
%     % figure(2);plot(x,v,'o',xq,vq2,':.');xlim([0 200])
%     % figure(3);plot(x,v,'o',xq,vq3,':.');xlim([0 200])
% 
%     % figure(1);plot3(lon_rho, lat_rho, -h); hold on
%     % plot3(lonnies,latties,-10,'rp','MarkerSize',15,'MarkerFaceColor','r')
%     % figure(2);imagescn(lon_u',lat_u',u(:,:,4)');cb = colorbar; ylabel(cb,'depth (m)');cmocean deep;axis equal; hold on %clim([-.6 .6]);
%     % plot(lonnies,latties,'rp','MarkerSize',15,'MarkerFaceColor','r')
% end

%%

%u
% j=86; i=152; %(close to land)
% j = 57 / 182 and i = 62 / 313
% %
% dims = size(u_zlev);
% i = randi([1 dims(1)]); %these are grid rows (oriented in latitudinal direction)
% j = randi([1 dims(2)]); %these are grid columns (oriented in longitudinal direction)
u_new = zeros(size(u_zlev,1), size(u_zlev,2), numlevels);
u_new(u_new == 0) = NaN;
for i = 1:size(u_new,1) %row, 189 ;; 74
    for j = 1:size(u_new,2) %column, 124 ;; 62
        vec = interp1(abs(squeeze(u_zlev(i,j,:))), squeeze(u(i,j,:)), zlevels, 'linear'); %could try akima/spline/cubic fits. but linear seems solid and can ignore NaNs on land easily
        for k = 1:numlevels
            if zlevels(k) > h(i,j)
                vec(k) = NaN;
            end
        end
        u_new(i,j,:) = vec;
    end
end

% %
% figure(1);imagescn(lon_u',lat_u',u_zlev(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean deep
% figure(1);imagescn(lon_u',lat_u',u_zlev(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean deep
% figure(1);imagescn(lon_u',lat_u',u_zlev(:,:,32)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean deep
% %
% 
% %
% figure(1);imagescn(lon_u',lat_u',u(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;clim([-.6 .6]);axis equal
% figure(1);imagescn(lon_u',lat_u',u(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;clim([-.6 .6]);axis equal
% figure(1);imagescn(lon_u',lat_u',u(:,:,32)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;clim([-.6 .6]);axis equal
% %
% 
% %
% figure(1);imagescn(lon_u',lat_u',u_new(:,:,1)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean curl;clim([-1.3 1.3]);axis equal
% figure(1);imagescn(lon_u',lat_u',u_new(:,:,20)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;clim([-.6 .6]);axis equal
% figure(1);imagescn(lon_u',lat_u',u_new(:,:,32)');cb = colorbar; ylabel(cb,'mean velocity (m/s)');cmocean speed;clim([-.6 .6]);axis equal
% %

% % %TEST
% % % plot(x,v,'o',xq,vq1,':.');
% % figure(1);plot(abs(squeeze(u_zlev(i,j,:))),squeeze(u(i,j,:)),'o',zlevels,vec,':.');
% % figure(2);plot(abs(squeeze(u_zlev(i,j,:))),squeeze(u(i,j,:)),'o',zlevels,vec2,':.');
% % title('(Default) Linear Interpolation');
% % %TEST

%v
v_new = zeros(size(v_zlev,1), size(v_zlev,2), numlevels);
v_new(v_new == 0) = NaN;
for i = 1:size(v_new,1)
    for j = 1:size(v_new,2)
        vec = interp1(abs(squeeze(v_zlev(i,j,:))), squeeze(v(i,j,:)), zlevels, 'linear');
        for k = 1:numlevels
            if zlevels(k) > h(i,j)
                vec(k) = NaN;
            end
        end
        v_new(i,j,:) = vec;
    end
end

%w,t,s [rho dimensions]
w_new = zeros(size(rho_zlev,1), size(rho_zlev,2), numlevels);
w_new(w_new == 0) = NaN;
t_new = zeros(size(rho_zlev,1), size(rho_zlev,2), numlevels);
t_new(t_new == 0) = NaN;
s_new = zeros(size(rho_zlev,1), size(rho_zlev,2), numlevels);
s_new(s_new == 0) = NaN;
for i = 1:size(w_new,1)
    for j = 1:size(w_new,2)
        vec_w = interp1(abs(squeeze(rho_zlev(i,j,:))), squeeze(w(i,j,:)), zlevels, 'linear');
        vec_t = interp1(abs(squeeze(rho_zlev(i,j,:))), squeeze(t(i,j,:)), zlevels, 'linear');
        vec_s = interp1(abs(squeeze(rho_zlev(i,j,:))), squeeze(s(i,j,:)), zlevels, 'linear');
        for k = 1:numlevels
            if zlevels(k) > h(i,j)
                vec_w(k) = NaN;
                vec_t(k) = NaN;
                vec_s(k) = NaN;
            end
        end
        w_new(i,j,:) = vec_w;
        t_new(i,j,:) = vec_t;
        s_new(i,j,:) = vec_s;
    end
end

end
