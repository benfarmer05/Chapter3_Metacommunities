

% % THIS SEEMS RIDICULOUS
% traj = zeros(1,9);
% myDir = pwd;
% files = dir(fullfile(myDir,'traj*')); %create an index for every netcdf file in directory
% for i = 1:length(files) %loop through each file
%     baseFileName = files(i).name;
%     currtraj = readmatrix(baseFileName);
%     traj = vertcat(traj,currtraj);
% end
% traj(1,:) = [];

%ASCII output in CMS format:
% Location  Particle  Time    Longitude  Latitude  Depth  Distance  Exit_code  Release_date
% %Julian date conversion:
% datetime(508440,'convertfrom','juliandate')

clc;clear
cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed/traj_file_troubleshoot/

%can't figure out how it makes sense to do this...these output files get
%HUGE very quickly!
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
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed/Testing'
ROMS_examplenest = 'croco_his.01227.nc';
ncdisp(ROMS_examplenest)
ROMS_grid = 'vigrid.nc';
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023'
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
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/Release';     
% [relpnts] = readmatrix('points.xlsx');
[relpnts] = readmatrix('points_650_none-on-land.csv');
plot(relpnts(:,11),relpnts(:,12), '*r', 'MarkerSize', 1)

gif('myfile.gif')
% gif('myfile.gif','DelayTime',1/2400,'resolution',400)

% %trajectories - by lat/lon/particle
% for i = 1:length(outputrearrange) %1:100:length(outputrearrange) 
%     for j = 1:size(outputrearrange(i).lon, 2)
%         plot(outputrearrange(i).lon(:,j),outputrearrange(i).lat(:,j));hold on
%     end
%     % gif
% end

%rearrange again for better GIF plotting - extremely terrible way to do
%this!
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
    for j = 1:size(gif_outputrearrange(i).lon, 2) %goes up to 1000 (#/particles)
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





%% identify bad points

badptlist = zeros(1,1);
% i = 1; j = 1; %testing
for i = 1:length(outputrearrange)
    currID = outputrearrange(i).loc;

    %this is simply throwing away an entire lat/lon if ANY of its particles
    %end up exiting the model domain (-1) or started on land.
    for j = 1:size(outputrearrange(i).lon, 2)
        currcode = outputrearrange(i).code(:,j);
     
        if ismember(-1, currcode) || ismember(-2, currcode) 
            nextindex = length(badptlist)+1;
            badptlist(nextindex) = currID;
            break
        end
    end
end
badptlist(1) = [];
badptlist = badptlist';

% preparing shapefiles for bad release points to export for GIS!
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/Release';
[relpnts] = readmatrix('points_650_none-on-land.csv');

% any release point for which a particle leaves 
badboys = zeros(length(badptlist),3);
badboys(:,1) = badptlist(:,1);
badrelpts = relpnts(ismember(relpnts(:,9),badptlist),:);
badrelpts(:,[1 2 3 4 5 6]) = [];

% write to a csv, for table joining in GIS to ID bad points and
% rescue/remove
writematrix(badrelpts, 'badrelpts.csv', 'delimiter', '\t');
% writematrix(badboys2, 'badpoints2_FINAL.csv', 'delimiter', '\t');



%

% Define the output shapefile name (without extension)
shapefileName = 'badrelpts';

% Create a geospatial structure array
S = struct('Geometry', 'Point', 'X', num2cell(badrelpts(:, 4)), 'Y', num2cell(badrelpts(:, 5)), 'ID', num2cell(badrelpts(:, 3)));

% Write the shapefile
shapewrite(S, shapefileName);

% Inform the user that the shapefile has been created
fprintf('Shapefile "%s.shp" has been created.\n', shapefileName);
% Export the same data to a CSV file
csvFileName = 'badrelpts.csv';
csvHeader = {'unique_ID', 'Longitude', 'Latitude'};
dataWithHeaders = [csvHeader; num2cell(badrelpts)];

% Write to CSV
cell2csv(csvFileName, dataWithHeaders, ',');

%

% plotting the baddies here for reference! these can be huge operations if you have a full release, careful!
%bathymetry
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed/Testing'
ROMS_examplenest = 'croco_his.01227.nc';
ncdisp(ROMS_examplenest)
ROMS_grid = 'vigrid.nc';
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023'
lons = ncread(ROMS_grid, 'lon_rho');
lats = ncread(ROMS_grid, 'lat_rho');
mask = ncread(ROMS_grid, 'mask_rho');
depths = ncread(ROMS_grid, 'h');
depthcut1 = depths;
depthcut1(depthcut1==0) = NaN;
depthcut1(mask == 0) = NaN;
% imagescn(lons', lats', -depthcut1');set(gca, 'ydir', 'normal');colorbar;clim([-250 -1]);axis equal; hold on
pcolor(lons, lats, -depthcut1);set(gca, 'ydir', 'normal');colorbar;clim([-250 -1]);axis equal; hold on

%release points
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/Release';     
% [relpnts] = readmatrix('points.xlsx');
[relpnts] = readmatrix('points_650_none-on-land.csv');
plot(relpnts(:,10),relpnts(:,11), '*k', 'MarkerSize', 1)

plot(badrelpts(:,4), badrelpts(:,5),'ob'); %show blue for -1 exit code
% plot(badboys2(:,2)-360, badboys2(:,3),'or'); %show red for -2 exit code

%% CMS make_mtx.m

%% parameters

% con_file = CMS connectivity output
con_file=load('con_file_1');
% n = number of polygons
n=200;


%% compute a square connectivity matrix(ixj)

settle =zeros(n,n);

% Build the sparse matrix
for i = 1:n
    
    %recruitment polygon in column #2
    pol =con_file(con_file(:,2)==i,:);
    if isempty(pol)
        a = zeros(n,1);
    end
    
    for j=1:n
        %source polygon in column #1
        source=pol(pol(:,1)==j,1);
        a (j,:) =numel(source);
    end
    
    settle(:,i)=a;
end

%% plot the matrix

% number of particles transiting form node i to node j
clf;
[xi,yi]=meshgrid(1:1:n, 1:1:n);
pcolor(xi,yi,settle);hold on
axis square
ylabel('Source Node');
xlabel('Receiving Node');
colorbar

%% CMS draw_traj_ascii.m

%% parameters

clear;clc
cd /Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed/traj_file_troubleshoot
examplenest='nest_1_20190221000000.nc';
ncdisp(examplenest)
traj_filename='traj_file_15';

%% read in the nest data
 
lon = ncread(examplenest, 'Longitude');
lat = ncread(examplenest, 'Latitude');
uvel = ncread(examplenest, 'zu');
vvel = ncread(examplenest, 'zv');

%% draw the land

mask = squeeze(uvel(:,:,1));
mask(mask>100) = 0;
mask(mask<100 & mask~=0) = 1;

%draw the land and water
figure;contourf(lon,lat,mask',[0.5 0.5],'linestyle','none');axis equal
% figure();imagescn(lon', lat', mask');set(gca, 'ydir', 'normal');hold on;colorbar;axis equal

%color of the land in rgb color model
colormap([0.5 0.6 0.7])

%% read in the trajectory data

%open trajectory file
data=load(traj_filename);
%sort the particles
[tmp,I] = sort(data(:,2));
data = data(I,:);
[tmp,I] = sort(data(:,1));
data = data(I,:);
%get data from particle file
id=data(:,1);
time=data(:,3);
lon=data(:,4);
lat=data(:,5);
depth=data(:,6);
%get the start and end of each trajectory
starts = find(time==0);
ends=[starts(2:end)-1; length(lon)];
%number of trajectories
num_traj=length(starts);
%different color for each trajectory
colors=jet(num_traj);

%% plot all the trajectories

hold on;

for i=1:num_traj

   nums=starts(i):ends(i);
 
   %plot(lon(nums), lat(nums),'color',colors(i,:));
   plot(lon(nums), lat(nums),'color','red');axis equal
 
end

hold off;

%print title
title(['Plotted ',num2str(num_traj),' trajectories']);

%% plot particle speed

num_colors = 100;

%set land velocities to NaN
uvel(uvel>100)= NaN;
vvel(vvel>100)= NaN;

%calculate speed
speed = sqrt(uvel .* uvel + vvel .* vvel);
%calculate min and max value of speed
minval = min(min(speed(:,:,1,1)));
maxval = max(max(speed(:,:,1,1)));
%draw speed
contourlevels = linspace(minval, maxval,num_colors);
colormap(jet(num_colors));
contourf(lon,lat,speed(:,:,1,1)',contourlevels,'linestyle','none');
colorbar;
axis equal;

% figure();imagescn(lon', lat',speed(:,:,1,1)');cb = colorbar; ylabel(cb,'velocity (m/s)');cmocean speed;axis equal

