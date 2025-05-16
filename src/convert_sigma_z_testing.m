function [u_new,v_new,w_new,t_new,s_new] = convert_sigma_z(u,v,w,t,s,h,udepths,vdepths,rhodepths,zlevels,scriptDrive) %layers, theta_s, theta_b, hc, zeta, warner, cmsname, fullname
% %testing purposes
% u = uvel;
% v = vvel;
% w = wvel;
% t = temp;
% s = salt;
% udepths = u_zlev;
% vdepths = v_zlev;
% rhodepths = rho_zlev;
% % zeta = ssh;
% % warner = warning_log;
% % cmsname = cms_fname;
% % fullname = fullFileName;
% %testing purposes

cd(scriptDrive)
numlevels = length(zlevels);

%% interpolate u,v,w,t,s at discrete depths on each of their respective grids

% %
% figure(1);imagescn(lon_u',lat_u',udepths(:,:,3)');cb = colorbar; ylabel(cb,'depth (m)');cmocean deep;axis equal; %clim([-.6 .6]);
% figure(1);imagescn(lon_u',lat_u',udepths(:,:,10)');cb = colorbar; ylabel(cb,'depth (m)');cmocean deep;axis equal; %clim([-.6 .6]);
% figure(1);imagescn(lon_u',lat_u',udepths(:,:,72)');cb = colorbar; ylabel(cb,'depth (m)');cmocean deep;axis equal; %clim([-.6 .6]);
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
%     dims = size(udepths); %dims = size(udepths);rhodepths
% 
%     % %option to limit random draws to shallow STT/STJ reefy area of domain
%     % row = randi([100 350]); %these are grid rows (oriented in latitudinal direction)
%     % col = randi([300 500]); %these are grid columns (oriented in longitudinal direction)
% 
%     %option to draw over entire domain
%     row = randi([1 dims(1)]); %these are grid rows (oriented in latitudinal direction)
%     col = randi([1 dims(2)]); %these are grid columns (oriented in longitudinal direction)
% 
%     x = abs(squeeze(udepths(row,col,:))); %the known depths. x = abs(squeeze(udepths(row,col,:)));rhodepths
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

% %testing
% % functionalI = 229; functionalJ = 429; %shallow
% functionalI = 152; functionalJ = 86; %deep
% 
% functionalX = abs(squeeze(udepths(functionalI,functionalJ,:)));
% functionalY = squeeze(u(functionalI,functionalJ,:));
% functionalTest = interp1(functionalX, functionalY, zlevels, 'linear');
% 
% functionalX_nozero = abs(squeeze(udepths_nozero(functionalI,functionalJ,:)));
% 
% zlevelsdiff = zlevels;
% zlevelsdiff(1) = functionalX_nozero(1);
% functionalTest_zlev = interp1(functionalX, functionalY, zlevelsdiff, 'linear');
% functionalTest_nozero_zlev = interp1(functionalX_nozero, functionalY, zlevelsdiff, 'linear');
% functionalTest_nozero = interp1(functionalX_nozero, functionalY, zlevels, 'linear');
% %testing

u_new = zeros(size(udepths,1), size(udepths,2), numlevels);
u_new(u_new == 0) = NaN;
for i = 1:size(u_new,1) %row, 189 ;; 74
    for j = 1:size(u_new,2) %column, 124 ;; 62

        xinterp = abs(squeeze(udepths(i,j,:)));

        % if ~issorted(xinterp) || any(xinterp(2:end) == 0) %check whether the data is not sorted properly from surface to bottom, or if there are errant zeros (potentially due to how USCROMS handles waves and how 'set_depth' above interacts with that)
        %     % xinterp = smooth(xinterp); %minor change to smooth out the values to more closely reflect reality. greatest effect on indices 1-4
        % 
        %     xsmooth = xinterp(1:6);
        %     xsmooth = smooth(xsmooth);
        %     % xinterp2 = xinterp;
        %     xinterp = [xsmooth; xinterp(7:end)];
        %     % xinterp2 = [xsmooth; xinterp2(7:end)];
        % end
        % 
        % %check for any identical values and apply a minuscule shift to the
        % % value located at the greater index. this avoids issues with
        % % interpolation
        % duplicates = diff(xinterp) == 0; %find where consecutive elements are identical
        % xinterp(find(duplicates) + 1) = xinterp(find(duplicates) + 1) + eps(xinterp(find(duplicates) + 1)); %apply shift
        % % duplicates = diff(xinterp2) == 0;
        % % xinterp2(find(duplicates) + 1) = xinterp2(find(duplicates) + 1) + eps(xinterp2(find(duplicates) + 1));

        vec = interp1(xinterp, squeeze(u(i,j,:)), zlevels, 'linear');

        for k = 1:numlevels
            if zlevels(k) > h(i,j)
                vec(k) = NaN;
            end
        end
        u_new(i,j,:) = vec;
    end
end

% %TEST
% % plot(x,v,'o',xq,vq1,':.');
% figure(1);plot(abs(squeeze(udepths(i,j,:))),squeeze(u(i,j,:)),'o',zlevels,vec,':.');
% figure(2);plot(abs(squeeze(udepths(i,j,:))),squeeze(u(i,j,:)),'o',zlevels,vec2,':.');
% title('(Default) Linear Interpolation');
% %TEST

% %testing
% functionalI = 319; functionalJ = 414;
% nonfunctionalI = 320; nonfunctionalJ = 414; %344, 248
% 
% functionalX = abs(squeeze(vdepths(functionalI,functionalJ,:)));
% functionalY = squeeze(v(functionalI,functionalJ,:));
% 
% nonfunctionalX = abs(squeeze(vdepths(nonfunctionalI,nonfunctionalJ,:)));
% nonfunctionalY = squeeze(v(nonfunctionalI,nonfunctionalJ,:));
% 
% functionalTest = interp1(functionalX, functionalY, zlevels, 'linear'); functionalTest
% 
% % if ~issorted(nonfunctionalX)
% %     xsmooth = nonfunctionalX(1:6);
% %     xsmooth = smooth(xsmooth);
% %     test = [xsmooth; nonfunctionalX(7:end)];
% % end
% 
% nonfunctionalTest = interp1(nonfunctionalX, nonfunctionalY, zlevels, 'linear'); nonfunctionalTest
% %testing

%v
v_new = zeros(size(vdepths,1), size(vdepths,2), numlevels);
v_new(v_new == 0) = NaN;
for i = 1:size(v_new,1)
    for j = 1:size(v_new,2)

        xinterp = abs(squeeze(vdepths(i,j,:)));

        % if ~issorted(xinterp) || any(xinterp(2:end) == 0) %check whether the data is not sorted properly from surface to bottom, or if there are errant zeros (potentially due to how USCROMS handles waves and how 'set_depth' above interacts with that)
        %     % xinterp = smooth(xinterp); %minor change to smooth out the values to more closely reflect reality. greatest effect on indices 1-4
        % 
        %     xsmooth = xinterp(1:6);
        %     xsmooth = smooth(xsmooth);
        %     % xinterp2 = xinterp;
        %     xinterp = [xsmooth; xinterp(7:end)];
        %     % xinterp2 = [xsmooth; xinterp2(7:end)];
        % end
        % 
        % %check for any identical values and apply a minuscule shift to the
        % % value located at the greater index. this avoids issues with
        % % interpolation
        % duplicates = diff(xinterp) == 0; %find where consecutive elements are identical
        % xinterp(find(duplicates) + 1) = xinterp(find(duplicates) + 1) + eps(xinterp(find(duplicates) + 1)); %apply shift
        % % duplicates = diff(xinterp2) == 0;
        % % xinterp2(find(duplicates) + 1) = xinterp2(find(duplicates) + 1) + eps(xinterp2(find(duplicates) + 1));

        vec = interp1(xinterp, squeeze(v(i,j,:)), zlevels, 'linear');

        for k = 1:numlevels
            if zlevels(k) > h(i,j)
                vec(k) = NaN;
            end
        end
        v_new(i,j,:) = vec;
    end
end

%w,t,s [rho dimensions]
w_new = zeros(size(rhodepths,1), size(rhodepths,2), numlevels);
w_new(w_new == 0) = NaN;
t_new = zeros(size(rhodepths,1), size(rhodepths,2), numlevels);
t_new(t_new == 0) = NaN;
s_new = zeros(size(rhodepths,1), size(rhodepths,2), numlevels);
s_new(s_new == 0) = NaN;

% %testing
% functionalI = 229; functionalJ = 429;
% nonfunctionalI = 230; nonfunctionalJ = 429; %344, 248
% 
% functionalX = abs(squeeze(rhodepths(functionalI,functionalJ,:)));
% functionalY = squeeze(t(functionalI,functionalJ,:));
% 
% nonfunctionalX = abs(squeeze(rhodepths(nonfunctionalI,nonfunctionalJ,:)));
% nonfunctionalY = squeeze(t(nonfunctionalI,nonfunctionalJ,:));
% 
% functionalTest = interp1(functionalX, functionalY, zlevels, 'linear'); functionalTest
% 
% nonfunctionalTest = interp1(nonfunctionalX, nonfunctionalY, zlevels, 'linear'); nonfunctionalTest
% 
% 
% 
% functionalI_zeta = 229; functionalJ_zeta = 429;
% nonfunctionalI_zeta = 230; nonfunctionalJ_zeta = 429; %344, 248
% 
% functionalX_zeta = abs(squeeze(rhodepths_zeta(functionalI_zeta,functionalJ_zeta,:)));
% functionalY_zeta = squeeze(t(functionalI_zeta,functionalJ_zeta,:));
% 
% nonfunctionalX_zeta = abs(squeeze(rhodepths_zeta(nonfunctionalI_zeta,nonfunctionalJ_zeta,:)));
% nonfunctionalY_zeta = squeeze(t(nonfunctionalI_zeta,nonfunctionalJ_zeta,:));
% 
% functionalTest_zeta = interp1(functionalX_zeta, functionalY_zeta, zlevels, 'linear'); functionalTest_zeta
% 
% nonfunctionalTest_zeta = interp1(nonfunctionalX_zeta, nonfunctionalY_zeta, zlevels, 'linear'); nonfunctionalTest_zeta
% 
% %testing

% STOPPING POINT - bunch of fucked stuff here

for i = 1:size(w_new,1)
    for j = 1:size(w_new,2)

        xinterp = abs(squeeze(rhodepths(i,j,:)));

        % % % version where the 2nd/4th index bug is caught first
        % % %
        % % %catch small minor 'set_depth' bug where the 3rd index in the depth
        % % %values for a given lat/lon in the u, v, or rho-grid(s) is
        % % %occasionally '0', which otherwise throws off the interpolation.
        % % %also, cases where the 3rd index is strangely a shallower depth
        % % %than the 2nd index, or identical to the 2nd
        % % if xinterp(2) == xinterp(4) || xinterp(2) == 0 %another bug where there can be identical values in the 2nd and 4th indices
        % %     xinterp = smooth(xinterp); %minor change to smooth out the values to more closely reflect reality. greatest effect on indices 1-4
        % % elseif xinterp(3) == 0 || xinterp(3) < xinterp(2) || xinterp(3) == xinterp(2)
        % %         xinterp(3) = (xinterp(2) + xinterp(4)) / 2; %simple average of the depth values preceding and following the incorrect '0'
        % % end
        % 
        % if ~issorted(xinterp) %check whether the data is not sorted properly from surface to bottom (potentially due to how USCROMS handles waves and how 'set_depth' above interacts with that)
        %     % xinterp = smooth(xinterp); %minor change to smooth out the values to more closely reflect reality. greatest effect on indices 1-4
        % 
        %     xinterp2 = xinterp;
        %     xinterp2 = sort(xinterp2);
        % 
        %     xsmooth = xinterp2(1:6);
        %     xsmooth = smooth(xsmooth);
        %     % xinterp2 = xinterp;
        %     xinterp2 = [xsmooth; xinterp2(7:end)];
        %     % xinterp2 = [xsmooth; xinterp2(7:end)];
        % end
        % 
        % %check for any identical values and apply a minuscule shift to the
        % % value located at the greater index. this avoids issues with
        % % interpolation
        % duplicates = diff(xinterp) == 0; %find where consecutive elements are identical
        % xinterp(find(duplicates) + 1) = xinterp(find(duplicates) + 1) + eps(xinterp(find(duplicates) + 1)); %apply shift
        % % duplicates = diff(xinterp2) == 0;
        % % xinterp2(find(duplicates) + 1) = xinterp2(find(duplicates) + 1) + eps(xinterp2(find(duplicates) + 1));
        % 
        % % %something I was testing to properly sort 'xinterp' by depth if it
        % % % isn't already after the checks above. testing showed that this
        % % % made no difference in the final product, at least in shallow
        % % % zones, and only very small differences in the deep zones (where
        % % % there aren't really these kinds of errors anyways generally)
        % % if ~issorted(xinterp2)
        % %     xinterp2 = sort(xinterp2);
        % %     % xinterp2 = xinterp2(randperm(length(xinterp2)));
        % % end
        % 
        % % xinterp2 = smooth(xinterp2);
        % % xinterp2(1:5) = xinterp2(randperm(5));
        
        vec_w = interp1(xinterp, squeeze(w(i,j,:)), zlevels, 'linear');
        vec_t = interp1(xinterp, squeeze(t(i,j,:)), zlevels, 'linear');
        vec_s = interp1(xinterp, squeeze(s(i,j,:)), zlevels, 'linear');

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


% % TEST for cutting off edges
% u_Agrid = zeros(size(udepths,1), size(vdepths,2), numlevels);
% v_Agrid = zeros(size(udepths,1), size(vdepths,2), numlevels);
% rho_Agrid = zeros(size(udepths,1), size(vdepths,2), numlevels);
% 
% u_test = u_new;
% u_test(:, 716, :) = [];

% % TEST for interp2
% % note - to avoid any kind of extrapolation outside the boundaries of each
% % grid, I think it makes sense to take one of the smaller grids (let's say
% % u-grid), and chop off one slice from its dimension with an extra slice
% % (in this case, going from 716 to 715). now, there is
% % a template grid which is smaller than both the rho-grid and v-grid. from
% % here, we can interpolate from rho- and v- into the newly chopped u-grid.
% u_test = u_new;
% u_test(:, 716, :) = []; %take away first row of second dimension from u-grid
% 
% [X,Y] = meshgrid(-3:3);
% V = peaks(X,Y);
% figure;surf(X,Y,V);title('Original Sampling');
% [Xq,Yq] = meshgrid(-3:0.25:3);
% Vq = interp2(X,Y,V,Xq,Yq);
% figure;surf(Xq,Yq,Vq);title('Linear Interpolation Using Finer Grid');

% scatteredinterpolant test
u_Agrid = zeros(size(rhodepths,1), size(rhodepths,2), numlevels);
v_Agrid = zeros(size(rhodepths,1), size(rhodepths,2), numlevels);
for k = 1:size(u_new, 3) % Loop over depth levels
    % Extract data for current depth level
    u_slice = u_new(:,:,k);
    F = scatteredInterpolant(lon_u(:), lat_u(:), u_slice(:), 'natural', 'boundary');
    % Evaluate interpolant at rho-grid points
    u_Agrid(:,:,k) = F(lon_rho, lat_rho);
end

for k = 1:size(v_new, 3) % Loop over depth levels
    % Extract data for current depth level
    v_slice = v_new(:,:,k);
    F = scatteredInterpolant(lon_v(:), lat_v(:), v_slice(:), 'natural', 'boundary');
    % Evaluate interpolant at rho-grid points
    v_Agrid(:,:,k) = F(lon_rho, lat_rho);
end

% figure(1);imagescn(lon_rho',lat_rho',u_Agrid(:,:,k)');axis equal;xlim([-64.95, -64.5]);ylim([17.6, 17.9]);title('griddata method')
% figure(2);imagescn(lon_rho',lat_rho',v_Agrid(:,:,k)');axis equal;xlim([-64.95, -64.5]);ylim([17.6, 17.9]);title('griddata method')

% % nanimresize test
% u_Agrid = zeros(size(rhodepths,1), size(rhodepths,2), numlevels);
% v_Agrid = zeros(size(rhodepths,1), size(rhodepths,2), numlevels);
% 
% for i = 1:numlevels
%     u_Agrid(:,:,i) = nanimresize(u_new(:,:,i), [size(rhodepths,1), size(rhodepths,2)]); %interpolate u_new (535x716) to rho-grid (536x716)
%     v_Agrid(:,:,i) = nanimresize(v_new(:,:,i), [size(rhodepths,1), size(rhodepths,2)]); %interpolate v_new (536x715) to rho-grid (536x716)
% end
% 
% figure(3);imagescn(lon_rho',lat_rho',u_Agrid(:,:,1)');axis equal;xlim([-64.95, -64.5]);ylim([17.6, 17.9]);title('nanimresize method')
% figure(4);imagescn(lon_rho',lat_rho',v_Agrid(:,:,1)');axis equal;xlim([-64.95, -64.5]);ylim([17.6, 17.9]);title('nanimresize method')

% 
% % uv2rho test
% u_Agrid2 = uv2rho(u_new, 'u'); %interpolate u_new (535x716) to rho-grid (536x716)
% v_Agrid2 = uv2rho(v_new, 'v'); %interpolate v_new (536x715) to rho-grid (536x716)
% 
% figure(3);imagescn(lon_rho',lat_rho',u_Agrid2(:,:,1)');axis equal;xlim([-64.95, -64.5]);ylim([17.6, 17.9]);title('uv2rho method')

u_new = u_Agrid;
v_new = v_Agrid;

end


% other options:
%griddata
%scatteredinterpolant
%interp2
%uv2rho