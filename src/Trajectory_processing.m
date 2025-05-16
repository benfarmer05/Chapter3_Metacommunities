%% bring in trajectories
cd '/Users/benja/OneDrive - Louisiana State University/Seascape Lab/Temporary_FTP/...'

trajecty = 'traj_file_007.nc';
ncdisp(trajecty)
coder = ncread(trajecty, 'exitcode');

% plot one trajectory
lons = ncread(trajecty,'lon');
lats = ncread(trajecty,'lat');
plot(lons-360,lats, '-'); hold on
% code = ncread(trajecty,'exitcode');
% disp(code) %take a look at the exit codes. 0 is good; everything else is ass my dude

%% identifying bad points
baddies1 = [];
baddies2 = [];
Lons = [];
Lats = [];
Depths = [];

% plot all the trajectories! also identify bad exit codes in relation to
% their release line (i.e. release point ID from GIS)
myDir = uigetdir; %gets directory
myFiles = dir(fullfile(myDir, '*.nc'));
hold on
for k = 1:length(myFiles)
  baseFileName = myFiles(k).name;
  fullFileName = fullfile(myDir, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  code = ncread(fullFileName,'exitcode');
  
  %identify bad exit codes and depths
  location = ncread(fullFileName, 'location');
  depth = ncread(fullFileName, 'depth');
  bad1 = code==-1; %particle left the model area
  bad2 = code==-2; %particle was too close to land, didn't move at all prob
  bad1locs = location(bad1);
  bad2locs = location(bad2);
  baddies1 = [baddies1; bad1locs];
  baddies2 = [baddies2; bad2locs]; %kind of a slow way to do this, but works in this case! basically adding ON TOP (or below, in terms of additional rows) of itself each loop
  disp(code)
  disp(baseFileName)
  
  %plot trajectories
  lon = ncread(fullFileName, 'lon');
  lat = ncread(fullFileName, 'lat');
  Lons = [Lons lon]; %concatenate the lons, lats, and depths
  Lats = [Lats lat];
  Depths = [Depths depth];
  plot(lon-360, lat, 'r')
end
%hold off

% preparing shapefiles for bad release points to export for GIS!
Releasy = load('MAR_test_baddies2_ReleaseFile.txt');
Releasy = load('L.analis_rel_2014.txt');

% any release point for which a particle leaves 
badboys1 = zeros(length(baddies1),3);
badboys1(:,1) = baddies1(:,1);
badboys1(:,2) = Releasy(baddies1,2); %exitcode -1 long
badboys1(:,3) = Releasy(baddies1,3); %exitcode -1 lat

badboys2 = zeros(length(baddies2),3);
badboys2(:,1) = baddies2(:,1);
badboys2(:,2) = Releasy(baddies2,2); %exitcode -1 long
badboys2(:,3) = Releasy(baddies2,3); %exitcode -1 lat

% write to a csv, for table joining in GIS to ID bad points and
% rescue/remove
cd '/Users/benja/MATLAB-Drive/Rare/Post_processing'
writematrix(badboys1, 'badpoints1_FINAL.csv', 'delimiter', '\t');
writematrix(badboys2, 'badpoints2_FINAL.csv', 'delimiter', '\t');

% plotting the baddies here for reference! these can be huge operations if you have a full release, careful!
plot(Lons-360,Lats,'Color',[.5 .5 .5 .3]);
hold on
axis equal
plot(Releasy(:,2)-360,Releasy(:,3),'*k');
plot(badboys1(:,2)-360, badboys1(:,3),'ob'); %show blue for -1 exit code
plot(badboys2(:,2)-360, badboys2(:,3),'or'); %show red for -2 exit code

plot(-Depths(:,baddies1)) %could try plot3 to look at this in 3d!!

%% more things to do...
% NEXT
% try and figure out how to export these "bad" points as a shapefile
% most bad points are coastal. eyeball and shift points by hand
% prioritize the polygons which have a good amount of water, and can easily
% take up a new point!
% basically get rid of blue points outside of domain; figure out red ones
% which are salvageable
% if I make any movements of red points, then re-run this stuff in MATLAB
% and see how things look

% %% put over hydro
% cd '/Users/benja/MATLAB-Drive/USVI_hydrodynamics';
% ncdisp('nest_1_20070101000000.nc') % Display the attributes of the nest file
% nester = 'nest_1_20100101000000.nc';
% 
% % Load in data from hydrodynamics nest file
% lons = ncread(nester, 'Longitude'); % Extract the longitudes to a variable, lon
% lats = ncread(nester, 'Latitude'); % Extract the lats
% depth = ncread(nester, 'h'); % Extract the bathymetry
% zu = ncread(nester, 'zu'); % Extract the northward velocity
% zu(zu==0) = NaN; % Make land super obvious! there can't be any velocity on land; this is a nifty workaround
% depth(depth==0) = NaN;
% 
% % Look at the where the land is
% imagesc(lons, lats, zu(:,:,1)'); %be careful 
% set(gca, 'ydir', 'normal') % Should be able to tell here how low res this is (prob very)
% colorbar
% hold on
% 
% imagesc(lons, lats, depth');
% set(gca, 'ydir', 'normal');
% hold on
% colorbar
% caxis([1 100])

% %% plotting over landmask
% cd '/Users/benja/OneDrive - Louisiana State University/Seascape Lab/GIS/Dan_USVI_maps';
% 
% LRL = shaperead('Low_Res_Landmask.shp');
% infoLRL = shapeinfo('Low_Res_Landmask.shp');
% projLRL = infoLRL.CoordinateReferenceSystem;
% 
% for i = 1:length(LRL)
%     [polylat,polylon] = projinv(projLRL,[LRL(i).X], [LRL(i).Y]);
%     LRL(i).lat = polylat;
%     LRL(i).lon = polylon;
%     plot(LRL(i).lon, LRL(i).lat,'Linewidth',3); hold on
% end
% 
% 
% HRL = shaperead('High_Res_Landmask.shp');
% infoHRL = shapeinfo('High_Res_Landmask.shp');
% 
% for i = 1:length(LRL)
%     [polylat,polylon] = projinv(projLRL,[LRL(i).X], [LRL(i).Y]);
%     LRL(i).lat = polylat;
%     LRL(i).lon = polylon;
%     plot(LRL(i).lon, LRL(i).lat,'Linewidth',2); hold on
% end
% 

% %% check the #/particles per processor and sundry
% cd '/Users/benja/MATLAB-Drive/Rare/Release_file';
% 
% reader = load('OCHY_5_2015rel.txt');
% 
% nodenum = 5; %arbitrary - change this to meet computing requirements
% procnum = nodenum*20; %HPC automatically assigns 20 processors to each node
% length(reader) %shows how many release lines we have
% particletot = sum(reader(:,5)); %how many particles we have
% particleproc = particletot/procnum;
% format long g %stop exponential notation
% disp(particleproc) %shows how many particles we are giving each processor
% 
% lobsterparticleproc = round(particleproc);
% lobstersimtime = 16934400; %196 days, in seconds
% lobsteroutput = 43200; %12 hours; frequency of output. internal time step is 2 hours (7200 s)
% lobstersteps = lobstersimtime/lobsteroutput;
% disp(lobstersteps) %number of step instances throughout pelagic larval duration (PLD)
% disp(lobsterparticleproc)
% disp(lobstersteps*lobsterparticleproc)
% 
% parrotparticleproc = round(particleproc);
% parrotsimtime = 4924800; %57 days
% parrotoutput = 21600; %6 hours. internal time step is 2 hours
% parrotsteps = parrotsimtime/parrotoutput;
% disp(parrotsteps)
% disp(parrotparticleproc)
% disp(parrotparticleproc*parrotsteps)
% 
% muttonparticleproc = round(particleproc);
% muttonsimtime = 3456000; %40 days
% muttonoutput = 21600; %6 hours. internal 2 h
% muttonsteps = muttonsimtime/muttonoutput;
% disp(muttonsteps)
% disp(muttonparticleproc)
% disp(muttonsteps*muttonparticleproc)
% 
% yellowparticleproc = round(particleproc);
% yellowsimtime = 4060800; %47 days
% yellowoutput = 21600; %6 hours. internal 2 h
% yellowsteps = yellowsimtime/yellowoutput;
% disp(yellowsteps)
% disp(yellowparticleproc)
% disp(yellowsteps*yellowparticleproc)
