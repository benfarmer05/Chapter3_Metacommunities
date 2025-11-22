function [u_new,v_new,w_new,t_new,s_new] = convert_sigma_z(u,v,w,t,s,theta_s,theta_b,hc,h,zeta,zlevels,layers)
% %testing purposes
% u = uvel;
% v = vvel;
% w = wvel;
% t = temp;
% s = salt;
% zeta = ssh;
% %testing purposes
% 

% cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023'
numlevels = length(zlevels);

%% interpolate u,v,w,t,s at discrete depths on each of their respective grids

%convert convert sigma-level to z-level at every lon/lat for u,v,w,t,s
u_zlev = set_depth(2,4,theta_s,theta_b,hc,layers,3,h,zeta,0); %'3' indicates u-point grid
v_zlev = set_depth(2,4,theta_s,theta_b,hc,layers,4,h,zeta,0); %v-point grid
rho_zlev = set_depth(2,4,theta_s,theta_b,hc,layers,1,h,zeta,0); %rho-point grid

%flip the depth dimension so that it is ordered surface to bottom
u_zlev = u_zlev(:,:,layers:-1:1);
v_zlev = v_zlev(:,:,layers:-1:1);
rho_zlev = rho_zlev(:,:,layers:-1:1);

%set the first depth level to zero, to avoid extrapolating to the surface
u_zlev(:,:,1)=0;
v_zlev(:,:,1)=0;
rho_zlev(:,:,1)=0;

%%

% %
% figure(1);imagescn(lon_u',lat_u',u_zlev(:,:,3)');cb = colorbar; ylabel(cb,'depth (m)');cmocean deep;axis equal; %clim([-.6 .6]);
% figure(1);imagescn(lon_u',lat_u',u_zlev(:,:,10)');cb = colorbar; ylabel(cb,'depth (m)');cmocean deep;axis equal; %clim([-.6 .6]);
% figure(1);imagescn(lon_u',lat_u',u_zlev(:,:,72)');cb = colorbar; ylabel(cb,'depth (m)');cmocean deep;axis equal; %clim([-.6 .6]);
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
% figure(2);imagescn(lon_u',lat_u',u(:,:,1)');cb = colorbar; ylabel(cb,'depth (m)');cmocean deep;axis equal; clim([-.6 .6]);
% figure(2);imagescn(lon_u',lat_u',u(:,:,10)');cb = colorbar; ylabel(cb,'depth (m)');cmocean deep;axis equal; clim([-.6 .6]);
% figure(3);imagescn(lon_u',lat_u',u(:,:,72)');cb = colorbar; ylabel(cb,'depth (m)');cmocean deep;axis equal; clim([-.6 .6]);
% %
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
%
% for i = 1:1 %200
%     dims = size(u_zlev); %dims = size(u_zlev);rho_zlev
% 
%     % %option to limit random draws to shallow STT/STJ reefy area of domain
%     % row = randi([100 350]); %these are grid rows (oriented in latitudinal direction)
%     % col = randi([300 500]); %these are grid columns (oriented in longitudinal direction)
% 
%     %option to draw over entire domain
%     row = randi([1 dims(1)]); %these are grid rows (oriented in latitudinal direction)
%     col = randi([1 dims(2)]); %these are grid columns (oriented in longitudinal direction)
% 
%     x = abs(squeeze(u_zlev(row,col,:))); %the known depths. x = abs(squeeze(u_zlev(row,col,:)));rho_zlev
%     v = squeeze(u(row,col,:)); %the known values of interest (current speed, temp, salinity) at those depths. v = squeeze(t(row,col,:));
%     xq = single(zlevels); %the depths I want to interpolate to
% 
%     lonnies = single(lon_rho(row,col));
%     latties = single(lat_rho(row,col));
% 
%     if ~isnan(v(1))
%         vq1 = interp1(x,v,xq, 'linear'); %the interpolated values of interest
%         vq2 = interp1(x,v,xq, 'spline');
%         vq3 = makima(x,v,xq);
%     end
% 
%     %is this silly when we could simply cut the interpolaion at the max of
%     %the bathymetry at that lat/lon? I don't think so, because the
%     %interpolation works best with layers of identical length.
%     for k = 1:length(zlevels)
%         if zlevels(k) > h(row,col)
%             vq1(k) = NaN;
%             vq2(k) = NaN;
%             vq3(k) = NaN;
%         end
%     end  
% 
%     figure(1);plot(x,v,'o',xq,vq1,':.');xlabel('Depth (m)');ylabel('Temperature (C)');title('Linear');%xlim([0 20])
%     % figure(2);plot(x,v,'o',xq,vq2,':.');xlabel('Depth (m)');ylabel('Temperature (C)');xlim([0 250])
%     % figure(2);plot(x,v,'o',xq,vq3,':.');xlabel('Depth (m)');ylabel('Temperature (C)');title('Akima')%;xlim([0 250])
% 
%     figure(3);imagescn(lon_rho',lat_rho',t(:,:,1)');cb = colorbar; ylabel(cb,'Temperature (C)');cmocean thermal;axis equal; hold on %clim([-.6 .6]);
%     plot(lonnies,latties,'rp','MarkerSize',15,'MarkerFaceColor','r')
%     % figure(1);plot3(lon_rho, lat_rho, -h); hold on
%     % plot3(lonnies,latties,-10,'rp','MarkerSize',15,'MarkerFaceColor','r')
% end

%%
%u
% j=86; i=152; %(close to land)
% j = 57 / 182 and i = 62 / 313
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
%         vec2 = interp1(abs(squeeze(u_zlev(i,j,:))), squeeze(u(i,j,:)),zlevels,'linear','extrap');
%         for k = 1:length(zlevels)
%             if zlevels(k) > h(i,j) 
%                 vec2(k) = NaN;
%             end
%         end
        u_new(i,j,:) = vec;
    end
end

% %TEST
% % plot(x,v,'o',xq,vq1,':.');
% figure(1);plot(abs(squeeze(u_zlev(i,j,:))),squeeze(u(i,j,:)),'o',zlevels,vec,':.');
% figure(2);plot(abs(squeeze(u_zlev(i,j,:))),squeeze(u(i,j,:)),'o',zlevels,vec2,':.');
% title('(Default) Linear Interpolation');
% %TEST

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

%% regrid u,v to the rho dimensions (w/t/s variables are already set up this way). so, converting from Arakawa C-grid to Arakawa A-grid
% this is bicubic interpolation by default

u_Agrid = zeros(size(rho_zlev,1), size(rho_zlev,2), numlevels);
v_Agrid = zeros(size(rho_zlev,1), size(rho_zlev,2), numlevels);

for i = 1:numlevels
    u_Agrid(:,:,i) = nanimresize(u_new(:,:,i), [size(rho_zlev,1), size(rho_zlev,2)]); %interpolate u_new to make it 536x716 (rho grid)
    v_Agrid(:,:,i) = nanimresize(v_new(:,:,i), [size(rho_zlev,1), size(rho_zlev,2)]);
end

u_new = u_Agrid;
v_new = v_Agrid;

end
