
clear;clc

myDir = pwd;
files = dir(fullfile(myDir,'traj*'));
outputrearrange = struct();

%% rearrange to structure. can export as .mat
%       alternative option here would be saving a structure for every ASCII
%       output - or in our case, one structure/ASCII per processor (so up
%       to probably 640 or so processors depending on size of the HPC
%       simulation).
% Note - this gets HUGE too, and may need to be broken up by variable or
% structured in some different kind of way.
tic
% %testing
% i = 1; l = 1; p = 1;
for i = 1:length(files)
    
    baseFileName = files(i).name;
    traj = readmatrix(baseFileName);
    % traj = 'traj_file_15'; % Can loop this
    % traj = readmatrix(traj);
    
    locs = unique(traj(:,1)); % Would help to have the vector of poly IDs/release lines from the release file here
    % outputrearrange = struct();

    for l = 1:length(locs)
    
        nextindex = length(outputrearrange)+1;

        subtraj = traj(traj(:,1) == locs(l),:); %subset traj file by location

        % NOTE - very important to subtract one here. GIS treats the first
        % unique ID as ZERO, not ONE. The CMS treats the first unique ID as
        % ONE, not ZERO. So, when you're going from CMS back to GIS, need
        % to translate back to the correct ID #.
        outputrearrange(nextindex).loc = locs(l)-1;
        partcls = unique(subtraj(:,2));
        times = unique(subtraj(:,3));
        nanmatrix = nan(length(times),length(partcls));
        % newloc = nanmatrix;
        % newpart = nanmatrix;
        newtime = nanmatrix;
        newlon = nanmatrix;
        newlat = nanmatrix;
        newdep = nanmatrix;
        newdist = nanmatrix;
        newcode = nanmatrix;
        newdate = nanmatrix;
    
        for p = 1:length(partcls)
            subp = subtraj(subtraj(:,2) == partcls(p),:); %subset again by particle
            [~, I] = sort(subp(:,3));
            subp = subp(I,:);
            % newloc(1:size(subp,1),p) = subp(:,1);
            % newpart(1:size(subp,1),p) = subp(:,2);
            newtime(1:size(subp,1),p) = subp(:,3);
            newlon(1:size(subp,1),p) = subp(:,4)-360;
            newlat(1:size(subp,1),p) = subp(:,5);
            newdep(1:size(subp,1),p) = subp(:,6);
            newdist(1:size(subp,1),p) = subp(:,7);
            newcode(1:size(subp,1),p) = subp(:,8);
            newdate(1:size(subp,1),p) = subp(:,9);
        end
    

        % Note - I think this needs to get edited so there is an indexing
        % of the current structure row, which is then referenced to add all
        % the material

        % if l==1
        %     nextindex = 1;
        % 
        % else
        %     nextindex = length(outputrearrange)+1;
        % end

        % nextindex = length(outputrearrange)+1;
        
        % outputrearrange(nextindex).loc = newloc;
        % outputrearrange(nextindex).part = newpart;
        outputrearrange(nextindex).time = newtime;
        outputrearrange(nextindex).lon = newlon;
        outputrearrange(nextindex).lat = newlat;
        outputrearrange(nextindex).dep = newdep;
        outputrearrange(nextindex).dist = newdist;
        outputrearrange(nextindex).code = newcode;
        outputrearrange(nextindex).date = newdate;
    end
    % save(strcat(trajlist(k).folder,'/',traj,'.mat'),'outputrearrange');
end
toc

outputrearrange(:,1) = [];
whos outputrearrange

%% PLOT
%bathymetry
ROMS_examplenest = 'croco_his.01227.nc';
ncdisp(ROMS_examplenest)
ROMS_grid = 'vigrid.nc';
lons = ncread(ROMS_grid, 'lon_rho');
lats = ncread(ROMS_grid, 'lat_rho');
mask = ncread(ROMS_grid, 'mask_rho');
depths = ncread(ROMS_grid, 'h');
depthcut1 = depths;
depthcut1(depthcut1==0) = NaN;
depthcut1(mask == 0) = NaN;
cb = colorbar();hold on
ylabel(cb,'Depth','FontSize',13,'Rotation',270)
imagescn(lons', lats', -depthcut1');set(gca, 'ydir', 'normal');clim([-250 -1]);axis equal; hold on
title('5-day simulation of SCTLD dispersal','FontSize',13)

% pcolor(lons, lats, -depthcut1);set(gca, 'ydir', 'normal');colorbar;clim([-250 -1]);axis equal; hold on

% DAN SUGGESTION FOR PLOTTING THAT DOESN'T MAKE IT LOOK LIKE THE LANDMASK
% IS ALL WRONG:
%
% [15:46] Daniel Holstein

% I think I figured this out once... I think you want to use boundary()
% 
% [15:47] Daniel Holstein

% with a shrink factor of 1

%
% DAN SUGGESTION FOR PLOTTING THAT DOESN'T MAKE IT LOOK LIKE THE LANDMASK


%release points
% [relpnts] = readmatrix('points.xlsx');
[relpnts] = readmatrix('points_650_none-on-land.csv');
plot(relpnts(:,11),relpnts(:,12), '*r', 'MarkerSize', 1)

% gif('myfile.gif')
% gif('myfile.gif','DelayTime',1/24000,'resolution',400)
gif('myfile.gif','DelayTime',1/24000)


% %trajectories - by lat/lon/particle
% for i = 1:length(outputrearrange) %1:100:length(outputrearrange) 
%     for j = 1:size(outputrearrange(i).lon, 2)
%         plot(outputrearrange(i).lon(:,j),outputrearrange(i).lat(:,j));hold on
%     end
%     % gif
% end

%rearrange again for better GIF plotting - quick way to get this done
gif_outputrearrange = outputrearrange;
for i = 2:length(gif_outputrearrange)/2
    gif_outputrearrange(1).lon = horzcat(gif_outputrearrange(1).lon, gif_outputrearrange(i).lon);
    gif_outputrearrange(1).lat = horzcat(gif_outputrearrange(1).lat, gif_outputrearrange(i).lat);
end

for i = (length(gif_outputrearrange)/2+2):length(gif_outputrearrange)
    gif_outputrearrange(2).lon = horzcat(gif_outputrearrange(2).lon, gif_outputrearrange(i).lon);
    gif_outputrearrange(2).lat = horzcat(gif_outputrearrange(2).lat, gif_outputrearrange(i).lat);
end
gif_outputrearrange = gif_outputrearrange(:,1:2);

%%

% %trajectories - timeseries
for i = 1:length(gif_outputrearrange)
    for j = 1:10:size(gif_outputrearrange(i).lon, 2) %goes up to 1000 (#/particles)
        % for k = 1:size(gif_outputrearrange(i).lon, 1) %goes up to 106 (#/time steps)
            plot(gif_outputrearrange(i).lon(:,j),gif_outputrearrange(i).lat(:,j));hold on
        % end
        gif
    end
end

% 
% for i = 1:length(gif_outputrearrange)
%     % for j = 1:size(gif_outputrearrange(i).lon, 2) %goes up to 1000 (#/particles)
%         for k = 1:size(gif_outputrearrange(i).lon, 1) %goes up to 106 (#/time steps)
%             plot(gif_outputrearrange(i).lon(k,:),gif_outputrearrange(i).lat(k,:));hold on
%             gif
%         end
%     % end
% end

% for i = 1:length(gif_outputrearrange)
%     for j = 1:size(gif_outputrearrange(i).lon, 1) %goes up to 106 (#/time steps)
%         for k = 1:size(gif_outputrearrange(i).lon, 2) %goes up to 1000 (#/particles)
%             plot(gif_outputrearrange(i).lon(j,k),gif_outputrearrange(i).lat(j,k));hold on
%         end
%         gif
%     end
% end