function [bigstruct, connection_struct, comp_count_matrix, scaled_count_matrix,...
    norm_comp_matrix, norm_scaled_matrix] = process_year_damselfish(year)
base_dir = '/Users/hannah/Library/CloudStorage/Box-Box/CMS/Fish_Model';
traj_dir = fullfile(base_dir,'2025_runs/', num2str(year));
out_dir = fullfile(base_dir, '2025_runs/', num2str(year));
shapefile_dir = fullfile(base_dir,'FGB_Shapefiles/');
scale_dir = fullfile(base_dir, 'Damsel/Release Scaling Damselfish/');

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

bigstruct = fullfile(traj_dir, sprintf('traj_all%d.mat', year));
load(bigstruct);

% Load trajectories
% trajlist = dir(fullfile(traj_dir,'traj*'));
% bigstruct_0 = struct();
% for i = 1:length(trajlist)
%     filename = fullfile(traj_dir, sprintf('traj_file_%03d.nc', i));
%     if isfile(filename)
%         % read from temporary file
%         bigstruct_0(i).time = ncread(filename, 'time');
%         bigstruct_0(i).location = ncread(filename,'location');
%         bigstruct_0(i).lon = ncread(filename, 'lon');
%         bigstruct_0(i).lat = ncread(filename, 'lat');
%         bigstruct_0(i).depth = ncread(filename, 'depth');
%         bigstruct_0(i).distance = ncread(filename, 'distance');
%         bigstruct_0(i).exitcode = ncread(filename, 'exitcode');
%         bigstruct_0(i).releasedate = ncread(filename, 'releasedate');
% 
%         delete(filename);
%     end
% end
% 
% bigstruct = struct();
% 
% % Flatten to structure
% row_no = 0;
% for i = 1:length(bigstruct_0)
%     locs_in_file = unique(bigstruct_0(i).location);
%     for j = 1:length(locs_in_file)
%         row_no = row_no + 1;
%         r_index = bigstruct_0(i).location == locs_in_file(j);
%         bigstruct(row_no).time = bigstruct_0(i).time;
%         bigstruct(row_no).location = bigstruct_0(i).location(r_index);
%         bigstruct(row_no).lon = bigstruct_0(i).lon(:,r_index);
%         bigstruct(row_no).lat = bigstruct_0(i).lat(:,r_index);
%         bigstruct(row_no).depth = bigstruct_0(i).depth(:,r_index);
%         bigstruct(row_no).distance = bigstruct_0(i).distance(:,r_index);
%         bigstruct(row_no).exitcode = bigstruct_0(i).exitcode(r_index);
%         bigstruct(row_no).releasedate = bigstruct_0(i).releasedate(r_index);
%     end
% end
% 
% clear bigstruct_0;
% 
% % Add reef information
% reefs = readtable(fullfile(shapefile_dir, 'FGB.csv'));
% reef_names = {'Stetson','WFGB','Horseshoe','EFGB','Bright','Geyer','McGrail','MacNeil','Sonnier','Rankin',...
%         '28 Fathom','Elvers','Bouma','Rezak','Sidner','Parker','Alderice'};
% vectors = cell(1,length(reef_names));
% 
% for k = 1:length(reef_names)
%     reef_rows = strcmp(reefs.Name, reef_names{k});
%     vectors{k} = reefs(reef_rows,:);
% end
% 
% tol = 1.55e-5;
% 
% for i = 1:length(bigstruct)
%     bigstruct(i).release_reef = '';
% 
%     % select the first element in the lat/lon cells of bigstruct (release
%     % location)
%     release_lat = bigstruct(i).lat(1);
%     release_lon = bigstruct(i).lon(1);
% 
%     for j = 1:length(vectors)
%         reef_lat = table2array(vectors{j}(:,4));
%         reef_lon = table2array(vectors{j}(:,5));
% 
%         lat_match = abs(release_lat - reef_lat) < tol;
%         lon_match = abs(release_lon - reef_lon) < tol;
% 
%         both_match = lat_match & lon_match;
% 
%         if any(both_match)
%             bigstruct(i).release_reef = reef_names{j};
%             break; % exit when match is found 
%         end
%     end
% end
% 
% Competency - this doesn't change with year
mpld = 40*24*60*60;
comp = 5*24*60*60;
spind = [comp,mpld]; 
t = 10800; % 3 hr timestep
sect = spind(1):t:spind(2);
xt = 1:length(sect);
k = log(.05) / ((spind(2)-spind(1))/t);
yt = 1 * (1+k).^xt;
xt = 1:(length(xt)+spind(1)/t);
yt = [repelem(0,(spind(1)/t)) yt];
comp_curve = yt(:);

% Build connection_struct with in_polygons
fgb = shaperead(fullfile(shapefile_dir, 'FGB_Habitat_FINAL.shp'));
fgb(17) = []; % drop Bryant
for i = 1:length(fgb)
    fgb(i).X = fgb(i).X + 360;
end

% Want to order them from west to east (left to right)
meanX = zeros(length(fgb),1);
for i = 1:length(fgb)
    meanX(i) = mean(fgb(i).X(~isnan(fgb(i).X)));
end

[~, sortIdx] = sort(meanX);
fgb = fgb(sortIdx);

% Build index 
fgbpolyvector = [];
for i = 1:length(fgb)
    Xs = fgb(i).X';
    Ys = fgb(i).Y';
    Is = repmat(i, length(Xs),1);
    Is(end) = NaN;
    fgbpolyvector = [fgbpolyvector; [Xs Ys Is]];
end

connection_struct = struct('release',{},'logic_in',{},'location_in',{});
for i = 1:length(bigstruct)
    [poly_in, poly_ind] = inpolygons(bigstruct(i).lon,bigstruct(i).lat,fgbpolyvector(:,1),fgbpolyvector(:,2));
    connection_struct(i).release = i;
    connection_struct(i).logic_in = poly_in;
    connection_struct(i).location_in = poly_ind;
    connection_struct(i).competency = comp_curve;
    connection_struct(i).comp_probability = comp_curve .* poly_in;
    connection_struct(i).Repro_Scale = NaN; % preallocating 
end

% Competency Probability
% for i = 1:length(connection_struct)
%     curr_logic = connection_struct(i).logic_in;
%     curr_comp = connection_struct(i).competency;
%     connection_struct(i).comp_probability = curr_comp .* curr_logic;
% end

% Scaling
scale_file = fullfile(scale_dir, sprintf('larva_%d_full.mat', year));

if isfile(scale_file)
    s = load(scale_file);

    vars = fieldnames(s);
    data = [];
    for v = 1:numel(vars)
        candidate = s.(vars{v});
        if (istable(candidate) || ismatrix(candidate)) && size(candidate,1) > 0
            data = candidate;
            break;
        end
    end
    if isempty(data)
        warning('No usable variable found in %s', scale_file);
    end
else
    warning('No scaling file for year %d: %s', year, scale_file);
    data = [];
end

% convert regular to julian dates
if ~isempty(data)
    % need numeric or array data
    if istable(data)
        data = table2array(data);
    end

    nrows = size(data,1);
    julian = zeros(nrows,1);

    for j = 1:nrows
        doy = data(j,1);
        date_vec = datetime(year, 1, 1) + days(doy - 1);
        julian(j) = juliandate(date_vec);
    end
    data(:,end+1) = round(julian);
end

% Assign reproduction scale thru matching release julian dates
if ~isempty(data)
    for i = 1:length(bigstruct)
        reldate = bigstruct(i).releasedate(1);
        release_julian = round(reldate);

        % find match in 'data'
        idx = find(data(:,end) == release_julian, 1);
        if ~isempty(idx)
            % column 4 has the repro scale
            connection_struct(i).Repro_Scale = data(idx, 4);
            %bigstruct(i).Repro_Scale = data(idx,4);
        else
            connection_struct(i).Repro_Scale = NaN;
            %bigstruct(i).Repro_Scale = NaN;
        end
    end
else
    % no data so field stays NaN
end

%save(fullfile(out_dir, sprintf('traj_all%d.mat',year)), 'bigstruct','-v7.3');

% Connection matrices
n_release = numel(connection_struct);
n_reef = 17;

comp_count_matrix = zeros(n_reef,n_reef);

for r = 1:n_release
    locations = connection_struct(r).location_in; % cell array ._.
    comp_logic = connection_struct(r).comp_probability;

    source = locations(1,:);

    for p = 1:size(locations,2)
        src = source{p};
        if src == 0, continue; end
        
        j_targets = locations(2:end,p);
        mask = comp_logic(2:end,p) > 0;
        j_targets = [j_targets{mask}];

        if ~isempty(j_targets)
            j_targets = j_targets([true; diff(j_targets(:)) ~= 0]); % ignores consecutive connections (just taking the first one)
        end
        
        for j = 1:numel(j_targets)
            tgt = j_targets(j);
            comp_count_matrix(src,tgt) = comp_count_matrix(src,tgt) + 1;
        end
    end
end

row_sum = sum(comp_count_matrix,2);
norm_comp_matrix = comp_count_matrix ./ row_sum;
norm_comp_matrix(isnan(norm_comp_matrix)) = 0;

save(fullfile(out_dir, sprintf('comp_count_matrix_%d.mat',year)), 'comp_count_matrix','-v7.3');
save(fullfile(out_dir, sprintf('norm_comp_matrix_%d.mat',year)), 'norm_comp_matrix','-v7.3');

% Now add repro scaling into the mix
for r = 1:length(connection_struct)
    curr_release = connection_struct(r);
    curr_prob = curr_release.comp_probability;
    curr_repro = curr_release.Repro_Scale;

    scaled_prob = curr_prob .* curr_repro;

    connection_struct(r).scaled_prob = scaled_prob; 
end

% Square matrix
scaled_count_matrix = zeros(n_reef, n_reef);

for r = 1:n_release
    locations = connection_struct(r).location_in; % cell array so use {}
    scaled_logic = connection_struct(r).scaled_prob;

    source = locations(1,:);

    % find where comp > 0 for each particle
    for p = 1:size(locations,2)
        src = source{p};
        if src == 0, continue; end

        j_targets = locations(2:end,p);
        mask = scaled_logic(2:end,p) > 0;
        j_targets = [j_targets{mask}];

        % % this collapses duplicate observations and only counts consecutive
        % % connections as 1 total connection (when in a row of timesteps)
        if ~isempty(j_targets)
            j_targets = j_targets([true; diff(j_targets(:)) ~= 0]);
        end

        for j = 1:numel(j_targets)
            tgt = j_targets(j);
            scaled_count_matrix(src,tgt) = scaled_count_matrix(src,tgt) + 1;
        end
    end
end

% Normalize
row_sum = sum(scaled_count_matrix, 2);
norm_scaled_matrix = scaled_count_matrix ./ row_sum;
norm_scaled_matrix(isnan(norm_scaled_matrix)) = 0;

save(fullfile(out_dir, sprintf('scaled_count_matrix_%d.mat',year)), 'scaled_count_matrix','-v7.3');
save(fullfile(out_dir, sprintf('norm_scaled_matrix_%d.mat',year)), 'norm_scaled_matrix','-v7.3');

end


