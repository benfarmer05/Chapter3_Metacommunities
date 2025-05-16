%% Ben Farmer output 03/15/2021
cd '/Users/benja/OneDrive - Louisiana State University/Seascape Lab/Temporary_FTP/Traj_USVI_2k'
trajecty = 'traj_file_001.nc';
ncdisp(trajecty)

% output frequency is 3600 seconds, or 1 hour. internal time step is half
% that

myDir = pwd; %get directory. can do pwd to get current dir, or uigetdir to choose
myFiles = dir(fullfile(myDir,'traj*'));

% Get unique release dates. could make code below more efficient by
% pre-allocating...would require knowing how many release dates you have
% (e.g. once a month in a year, so 12).
allRelDats = [];
for i = 1:length(myFiles)
    baseFileName = myFiles(i).name;
    fullFileName = fullfile(myDir, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    reldat = ncread(fullFileName,'releasedate');
    reldats = unique(reldat);
    allRelDats = [allRelDats;
                  reldats];
end
uniqueRelDats = unique(allRelDats);

% plot randomtraj. can be intensive!!
for i = 1:length(myFiles)
    relFile = myFiles(i).name;
    lon = ncread(relFile,'lon');
    lat = ncread(relFile,'lat');
    plotrand = randsample(size(lon,2),1000);
    plot(lon(:,plotrand),lat(:,plotrand));hold on;
end
axis equal

%% To do
% 1. Have to decide on a decay curve
% 2. Have to run inpolygons - FOR EACH UNIQUE RELEASE DAY
% 3. Apply decay curve to inpolygons results for EACH PARTICLE's probability
% of reaching other polygons
% 4. Scale so that each column sums to 1
% 5. Accumulate all probabilities from i to j in each unique release day
% (into square conn matrix)

%% 1. Create virulence curve
cn = 0; % Competency time in seconds
mn = 10*24*60*60; % Max PLD time in seconds

spind = [cn; mn];

tstep = 1800 % This is the timestep FROM THE CMS! - needs to be updated
sect = spind(1):tstep:spind(2); 
xt = 1:length(sect);
yt = [];
    for i = 1:length(xt);
        k = log(.01)/((spind(2)-spind(1))/tstep);
        yt(i) = 1*(1+k).^xt(i); 
    end
xt = 1:(length(xt) + spind(1)/tstep); % add time to competency
yt = [repelem(0,(spind(1)/tstep)) yt];
h = plot(xt,yt); hold on; 

xticks([0:5*24*60*60/tstep:50*24*60*60/tstep])
xticklabels([0:5:50])
xlabel('Days of dispersal')
print('VirulenceCurve','-dpdf')

%% 2. Run inpolygons
% This script will cycle through  traj files and begin to analyze them
% for connectivity. 

% Things you need:
% Habitat polygons (reefpolys), with 1st column long, and 2nd column lat
% Release locations/release files
% Virulency decay curve
% inpolygons.m

load ReleaseFile.txt
uniquerels = unique(ReleaseFile(:,1));
relpnts = ReleaseFile(1:length(uniquerels),1:3);

reefpolys = shaperead('reefpolys.shp');
info = shapeinfo('reefpolys.shp');
proj = info.CoordinateReferenceSystem;
for i = 1:length(reefpolys)
    [polylat,polylon] = projinv(proj,[reefpolys(i).X], [reefpolys(i).Y]);
    reefpolys(i).lat = polylat;
    reefpolys(i).lon = polylon;
end

Xs=[];
Ys=[];
IDs = [];
for i = 1:length(reefpolys)
    X = [reefpolys(i).lon];
    Y = [reefpolys(i).lat];
    ID = repmat(i,1,length(X));
    
    Xs = [Xs X];
    Ys = [Ys Y];
    IDs = [IDs ID];
end

reefpolymat = [Xs;
               Ys;
               IDs];
    

% % Prepare and preallocate matrices 
filenum = length(myFiles);
In = cell(1,filenum);
Ind_d = cell(1,filenum);
Locs = cell(1,filenum);
Rels = cell(1,filenum);
Lons = cell(1,filenum);
Lats = cell(1,filenum); 

for k = 1%filenum
    relFile = myFiles(k).name;
    lon = ncread(relFile,'lon');
    lat = ncread(relFile,'lat');
    reldat = ncread(relFile,'releasedate');
    loc = ncread(relFile,'location'); %LOCATION IS THE RELEASE LINE
    
    Lons{k} = lon;
    Lats{k} = lat;
    Locs{k} = loc;
    Rels{k} = reldat;
    
    maxl = size(lon,2);
    brks = 1:1000:maxl; % Only do 1000 at a time
    brks(end + 1) = maxl;
    for m = 1:length(brks)-1
        [poly_in, poly_ind] = inpolygons(lon(:,brks(m):brks(m+1)-1),lat(:,brks(m):brks(m+1)-1),reefpolymat(1,:)+360,reefpolymat(2,:)); % NOTE! +360 is used here because CMS outputs longitudes in 0-360 - CHANGE IF NECESSARY!
        
        In{k}(:,brks(m):brks(m+1)-1) = poly_in;
        Ind(:,brks(m):brks(m+1)-1) = poly_ind;
        Ind_d{k} = cell2mat(Ind);
    end
    clear poly_in poly_ind lon lat loc Ind
end
save WIP -v7.3   

%% Next step
% 1.	For each particle during its competency period, the timesteps during 
% which it is inside a polygon, and that polygon?s ID, are determined. 
% Thus, the most information-dense matrix (Ind) is sized t x n, n = the # 
% of particles released, and t = the number of timesteps in the competency
% period. (above)

% 2.	Next, each column of this matrix is multiplied by a vector, y(t), 
% which contains the competency/mortality curve. Thus, the maximum number 
% in this matrix (yt_mat) is 1. This matrix is still t x n.

% Note that for computational reasons, these matrices are often broken into 
% smaller chunks and stored in cells.

yt_mat = cell(1,filenum);
for k = 1 %:filenum
    yt_mat{k} = zeros(size(In{k}));
    for i=1:size(In{k},2)
        yt_mat{k}(:,i) = In{k}(:,i) .* yt'; % For each trajectory, apply the competency/mortality curve.
    end
end

save yt_mat yt_mat -v7.3
% clear In

% 3.	Next, we search this large matrix to find connections i,j for each 
% particle?s trajectory. This matrix is j x n, where j is the number of 
% polygons (or possible connections). We accumulate the probabilities from 
% yt_mat that fit the condition of connecting i and j for each particle. 
% Finally, the columns of this large matrix (ConnLocMat) are normalized to 
% sum to 1.
% a.	To accumulate probabilities, you must ask the opposite question ? 
% what is the probability of a connection NOT occurring?
% 
% 1 ? prod(1 ? pct)
% 
% Where pct is the probability of a single connection i to j along a single 
% particle?s trajectory. There may be multiple instances of an i,j 
% connection during a single trajectory (t).

% Allocate
ConnLocMat = cell(1,filenum);
for k = 1:filenum
    ConnLocMat{k} = zeros(length(relpnts),length(yt_mat{k}));
end

% Fill. May take hours-days to complete this loop.
for k = 1 %:cellnum
    tic
    fprintf('%d ', k);
    for i = 1:size(yt_mat{k},2)
        if rem(i,10000) ==0 % Just a printed output to see that status of the loop
            fprintf('%d ',i);
            fprintf('\n');
        end
        if sum(yt_mat{k}(:,i)) == 0 % If there are no connections in the trajectory, skip it.
            continue
        else
            for j = 1:length(relpnts)
                ConnLocMat{k}(j,i) = 1-prod(1-[yt_mat{k}(Ind_d{k}(:,i) == j,i)]); % Accumulate the probabilities for the particle's unique connections
            end
        end
        ConnLocMat{k}(:,i) = ConnLocMat{k}(:,i)/sum(ConnLocMat{k}(:,i)); % Normalize to 1. "Each particle has a probability to end up somewhere"
    end
    toc
end

tic
save ConnLocMat ConnLocMat -v7.3
toc

% 4.	We summarize the large matrix ConnLocMat into an i x j matrix, 
% ConnAccum (or ConnAccumTot). Referencing the release file, we again 
% accumulate connectivity probabilities from each trajectory going i to j. 
% Once complete, the matrix is again normalized along rows (along source 
% nodes) to sum to 1. ?The probability of arrival from point i".
% a.	To accumulate probabilities, you must ask the opposite question ? 
% what is the probability of a connection NOT occurring?

% 1 ? prod(1 ? pcm)

% Where pcm is a single particle?s (migrant?s) accumulated connection 
% probability i to j, and there may be multiple particles (m) making the same connection.

% Allocate
ConnMats = cell(1,filenum);
for k = 1:8
    ConnMats{k} = zeros(5177,5177);
end

% Fill. May take hours to days.
for k = 1:cellnum
    tic
    
    fprintf('%d ', k); 
    for i = 1:5177
        Loclin = polyrep14(Locs{k}) == i; % NOT SURE THIS IS RIGHT!!!!!!! SHould it be: find(polyrep14(Locs{k}) == i);??
        if rem(i,100) == 0 % Just an output to keep track of the loop
            fprintf('%d ',i);
            fprintf('\n');
        end
        for j = 1:5177
            ConnMats{k}(i,j) = 1-prod(1-[ConnLocMat{k}(j,Loclin)]); % Accumulate all of the i,j connections across particle trajectories
        end
    end
    toc
end


%% 19 April: Envelopes

% Start disease at 1725 Flat Cay (remember this is 1 removed from actual
% attribute table ID because BVI point was removed)

% since we're starting at 1725, this is the first release (this is just
% January). Gets more complicated if you're looking at other months, since
% at that point the lines start iterating and February might be like
% 3000-something (well...to be precise it would be 1725*2. Can easily code
% this)

% NOTE -  all this was changed a few
% times by dan so need to go back and adjust accordingly
for k = 1:length(Locs)
    if ismember(startloc, cell2mat(Locs(k)))
        trajfile = k;
        partcls = find(cell2mat(Locs(k)) == startloc);
        Lons_start = Lons(k);
    else
        continue
    end
end
