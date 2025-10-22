%%% OPTIMIZED PASSING-THROUGH CONNECTIVITY WITH PREPROCESSING %%%




% NOTE - this could still have an issue with indexing, but I think I fixed
% that part....need to make sure 'time' is properly included as well
%
% NOTE ALSO - check big commented out version below all this, which may be
% working (not sure??) and also save other stuff that's needed like depth
% and exit code!





clear; clc;

%% ========== CONFIGURATION ==========
DO_PREPROCESSING = true;  % Set to false to load existing preprocessed data
YEAR = 2019;
QUARTER = 1;  % 1, 2, 3, or 4

fprintf('=== REEF CONNECTIVITY ANALYSIS ===\n');
fprintf('Year: %d, Quarter: %d\n', YEAR, QUARTER);
fprintf('Preprocessing mode: %s\n', string(DO_PREPROCESSING));

%% Setup paths
projectPath = matlab.project.rootProject().RootFolder;
dataPath = fullfile(projectPath, 'data');
outputPath = fullfile('D:\Dissertation\CMS_traj\output');

% Input/output directory based on quarter
quarter_name = sprintf('Q%d_%d', QUARTER, YEAR);
tempPath = fullfile('D:\Dissertation\CMS_traj', quarter_name);

% Output directory for preprocessed data
preprocessPath = fullfile(outputPath, 'CMS_traj', quarter_name);
if ~exist(preprocessPath, 'dir')
    mkdir(preprocessPath);
end

%% Load reef centroids and build polygons
centroids = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
unique_IDs = centroids(:,1);
Xs = [centroids(:,8) centroids(:,10) centroids(:,12) centroids(:,14) centroids(:,8)];
Ys = [centroids(:,9) centroids(:,11) centroids(:,13) centroids(:,15) centroids(:,9)];
n_locations = size(centroids,1);

fprintf('Loaded %d reef polygons.\n', n_locations);

%% Load release file
releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019_Q1.txt'));
release_reef_IDs = releasefile(:,1);

%% ========== PREPROCESSING OR LOADING ==========

if DO_PREPROCESSING
    %% PREPROCESSING MODE: Generate event files for all unique events
    fprintf('\n=== PREPROCESSING MODE ===\n');
    fprintf('This will process all files and create separate data files for each release event.\n');
    
    %% Step 1: Identify all unique release events (or load if exists)
    unique_events_file = fullfile(preprocessPath, 'unique_events.mat');
    
    if exist(unique_events_file, 'file')
        fprintf('\nLoading existing unique release events from: %s\n', unique_events_file);
        load(unique_events_file, 'file_release_dates', 'unique_release_dates', 'trajlist');
        n_events = length(unique_release_dates);
        fprintf('Found %d unique release events (loaded from file).\n', n_events);
    else
        fprintf('\nStep 1: Scanning trajectory files for release events...\n');
        trajlist = dir(fullfile(tempPath, 'traj*.nc'));
        fprintf('Found %d trajectory files.\n', numel(trajlist));
        
        file_release_dates = zeros(numel(trajlist), 1);
        
        fprintf('Reading release dates from files...\n');
        tic;
        for i = 1:numel(trajlist)
            if mod(i, 50) == 0
                elapsed = toc;
                rate = i / elapsed;
                remaining = (numel(trajlist) - i) / rate;
                fprintf('  Progress: %d/%d files (%.1f%%) - %.1f files/sec - ETA: %.1f sec\n', ...
                    i, numel(trajlist), 100*i/numel(trajlist), rate, remaining);
            end
            
            filename = fullfile(tempPath, trajlist(i).name);
            ncid = -1;
            try
                ncid = netcdf.open(filename, 'NC_NOWRITE');
                varid = netcdf.inqVarID(ncid, 'releasedate');
                releasedate = netcdf.getVar(ncid, varid, 0, 1, 'double');
                file_release_dates(i) = releasedate;
                netcdf.close(ncid);
                ncid = -1;
            catch ME
                warning('Could not read %s: %s', trajlist(i).name, ME.message);
                file_release_dates(i) = NaN;
                if ncid ~= -1
                    try netcdf.close(ncid); catch; end
                end
            end
        end
        elapsed_total = toc;
        fprintf('Completed reading %d files in %.2f seconds (%.1f files/sec)\n', ...
            numel(trajlist), elapsed_total, numel(trajlist)/elapsed_total);
        
        file_release_dates = file_release_dates(~isnan(file_release_dates));
        unique_release_dates = unique(file_release_dates);
        n_events = length(unique_release_dates);
        
        fprintf('\nFound %d unique release events.\n', n_events);
        
        % Save unique events info
        fprintf('Saving unique release events to: %s\n', unique_events_file);
        save(unique_events_file, 'file_release_dates', 'unique_release_dates', 'trajlist', '-v7.3');
    end
    
    %% Step 2: Process each unique event and save to separate files
    fprintf('\nStep 2: Processing each release event and saving data...\n');
    
    for event_idx = 1:n_events
        target_date = unique_release_dates(event_idx);
        calendar_date = datetime(target_date, 'ConvertFrom', 'juliandate');
        
        % Create descriptive filename based on date
        date_str = datestr(calendar_date, 'yyyy_mmmdd_HHMMSS');
        event_filename = fullfile(preprocessPath, sprintf('event_%s.mat', date_str));
        
        % Skip if this event file already exists
        if exist(event_filename, 'file')
            fprintf('\n--- Event %d/%d already exists, skipping ---\n', event_idx, n_events);
            fprintf('  File: %s\n', event_filename);
            continue;
        end
        
        fprintf('\n--- Processing Event %d/%d ---\n', event_idx, n_events);
        fprintf('  Julian date: %.2f\n', target_date);
        fprintf('  Calendar date: %s\n', datestr(calendar_date, 'dd-mmm-yyyy HH:MM:SS'));
        
        % Identify files containing this event
        files_with_event = abs(file_release_dates - target_date) < 0.01;
        relevant_files = find(files_with_event);
        n_relevant_files = numel(relevant_files);
        
        fprintf('  Found in %d files\n', n_relevant_files);
        
        % Collect particles from relevant files
        lon_cells = cell(n_relevant_files, 1);
        lat_cells = cell(n_relevant_files, 1);
        src_cells = cell(n_relevant_files, 1);
        
        tic;
        for idx = 1:n_relevant_files
            i = relevant_files(idx);
            filename = fullfile(tempPath, trajlist(i).name);
            
            try
                lon = ncread(filename, 'lon');
                lat = ncread(filename, 'lat');
                % location = ncread(filename, 'location') + 1;
                location = ncread(filename, 'location');
                releasedate = ncread(filename, 'releasedate');
                
                event_mask = abs(releasedate - target_date) < 0.01;
                if ~any(event_mask)
                    continue;
                end
                
                lon(lon > 180) = lon(lon > 180) - 360;
                
                % Get particles for this event (KEEP AS 2D!)
                lon_event = lon(:, event_mask);      % [timesteps × particles]
                lat_event = lat(:, event_mask);      % [timesteps × particles]
                loc_event = location(event_mask);    % [1 × particles]
                
                % Map to source indices
                reef_IDs = release_reef_IDs(loc_event);
                [~, src_idx] = ismember(reef_IDs, unique_IDs);
                
                % Store 2D matrices (NOT FLATTENED!)
                lon_cells{idx} = lon_event;
                lat_cells{idx} = lat_event;
                src_cells{idx} = src_idx';  % Row vector of source indices
                
            catch ME
                warning('Error processing %s: %s', trajlist(i).name, ME.message);
            end
        end
        
        % Concatenate horizontally across particles (not vertically!)
        all_lon = horzcat(lon_cells{:});           % [timesteps × total_particles]
        all_lat = horzcat(lat_cells{:});           % [timesteps × total_particles]
        all_source_idx = horzcat(src_cells{:});    % [1 × total_particles]
        
        elapsed = toc;
        n_timesteps = size(all_lon, 1);
        n_particles = size(all_lon, 2);
        fprintf('  Collected %d particles × %d timesteps = %d positions in %.1f seconds\n', ...
            n_particles, n_timesteps, numel(all_lon), elapsed);
        
        % Organize into bigstruct-style array by source location
        fprintf('  Organizing into struct array by source location...\n');
        unique_sources = unique(all_source_idx);
        n_sources = length(unique_sources);
        
        event_data = struct('source_reef_id', {}, 'source_reef_idx', {}, ...
                           'lon', {}, 'lat', {}, 'n_particles', {}, 'n_timesteps', {});
        
        for s = 1:n_sources
            src_id = unique_sources(s);
            
            % Get particles from this source
            particle_mask = all_source_idx == src_id;
            
            event_data(s).source_reef_idx = src_id;
            event_data(s).source_reef_id = unique_IDs(src_id);
            event_data(s).lon = all_lon(:, particle_mask);      % [timesteps × particles]
            event_data(s).lat = all_lat(:, particle_mask);      % [timesteps × particles]
            event_data(s).n_particles = sum(particle_mask);
            event_data(s).n_timesteps = n_timesteps;
            event_data(s).releasedate = target_date;
            event_data(s).calendar_date = calendar_date;
        end
        
        fprintf('  Organized into %d source locations (bigstruct-style)\n', n_sources);
        
        % Save metadata separately
        event_metadata = struct();
        event_metadata.event_idx = event_idx;
        event_metadata.target_date = target_date;
        event_metadata.calendar_date = calendar_date;
        event_metadata.n_sources = n_sources;
        event_metadata.n_particles_total = n_particles;
        event_metadata.n_timesteps = n_timesteps;
        
        save(event_filename, 'event_data', 'event_metadata', '-v7.3');
        fprintf('  Saved to: %s\n', event_filename);
    end
    
    % Save/update metadata file
    metadata = struct();
    metadata.year = YEAR;
    metadata.quarter = QUARTER;
    metadata.n_events = n_events;
    metadata.unique_release_dates = unique_release_dates;
    metadata.n_locations = n_locations;
    metadata.preprocessing_date = datetime('now');
    
    % Add list of event files for reference
    event_files = dir(fullfile(preprocessPath, 'event_*.mat'));
    metadata.event_files = {event_files.name}';
    
    metadata_file = fullfile(preprocessPath, 'metadata.mat');
    save(metadata_file, 'metadata', '-v7.3');
    
    fprintf('\n=== PREPROCESSING COMPLETE ===\n');
    fprintf('Processed %d events and saved to: %s\n', n_events, preprocessPath);
    fprintf('Event files:\n');
    for i = 1:length(metadata.event_files)
        fprintf('  %s\n', metadata.event_files{i});
    end
    fprintf('\nSet DO_PREPROCESSING = false to run analysis on these events.\n');
    
    return;  % Exit after preprocessing
    
else
    %% ANALYSIS MODE: Load preprocessed data
    fprintf('\n=== ANALYSIS MODE ===\n');
    
    % Load metadata
    metadata_file = fullfile(preprocessPath, 'metadata.mat');
    if ~exist(metadata_file, 'file')
        error('Metadata file not found. Run with DO_PREPROCESSING = true first.');
    end
    
    load(metadata_file, 'metadata');
    n_events = metadata.n_events;
    
    fprintf('Found %d preprocessed events in: %s\n', n_events, preprocessPath);
    
    % Find all event files
    event_files = dir(fullfile(preprocessPath, 'event_*.mat'));
    if isempty(event_files)
        error('No event files found. Run with DO_PREPROCESSING = true first.');
    end
    
    % Pick one event to analyze
    event_idx = randi(length(event_files));
    event_filename = fullfile(preprocessPath, event_files(event_idx).name);
    
    fprintf('\nLoading Event %d/%d: %s\n', event_idx, length(event_files), event_files(event_idx).name);
    load(event_filename, 'event_data', 'event_metadata');
    
    fprintf('  Julian date: %.2f\n', event_metadata.target_date);
    fprintf('  Calendar date: %s\n', datestr(event_metadata.calendar_date, 'dd-mmm-yyyy HH:MM:SS'));
    fprintf('  Source locations: %d\n', event_metadata.n_sources);
    fprintf('  Total particles: %d × %d timesteps\n', ...
        event_metadata.n_particles_total, event_metadata.n_timesteps);
    
    target_date = event_metadata.target_date;
    calendar_date = event_metadata.calendar_date;
    
    % Flatten bigstruct-style array for polygon checking
    fprintf('\nFlattening bigstruct for connectivity analysis...\n');
    fprintf('  event_data is a struct array with %d elements (one per source reef)\n', length(event_data));
    
    % Pre-allocate for speed
    total_positions = event_metadata.n_particles_total * event_metadata.n_timesteps;
    all_lon = zeros(total_positions, 1);
    all_lat = zeros(total_positions, 1);
    all_source_idx = zeros(total_positions, 1);
    
    pos_idx = 1;
    for s = 1:length(event_data)
        % Get data from this source
        lon_src = event_data(s).lon(:);  % Flatten [timesteps × particles]
        lat_src = event_data(s).lat(:);
        n_pos = length(lon_src);
        
        % Store in flat arrays
        all_lon(pos_idx:pos_idx+n_pos-1) = lon_src;
        all_lat(pos_idx:pos_idx+n_pos-1) = lat_src;
        all_source_idx(pos_idx:pos_idx+n_pos-1) = event_data(s).source_reef_idx;
        
        pos_idx = pos_idx + n_pos;
    end
    
    fprintf('  Flattened %d source reefs into %d total positions\n', ...
        length(event_data), total_positions);
end

%% ========== ANALYSIS CONTINUES FROM HERE ==========
fprintf('\n=== STARTING CONNECTIVITY ANALYSIS ===\n');

%% Remove NaN positions
valid = ~isnan(all_lon) & ~isnan(all_lat) & all_source_idx > 0;
lon_all = all_lon(valid);
lat_all = all_lat(valid);
source_all = all_source_idx(valid);

fprintf('Valid positions after filtering: %d\n', numel(lon_all));

%% OPTIMIZATION 1: Pre-compute bounding boxes for quick rejection
fprintf('\nPre-computing polygon bounding boxes...\n');
bbox_xmin = min(Xs(:,1:4), [], 2);
bbox_xmax = max(Xs(:,1:4), [], 2);
bbox_ymin = min(Ys(:,1:4), [], 2);
bbox_ymax = max(Ys(:,1:4), [], 2);

%% OPTIMIZATION 2: Data already prepared
fprintf('Data already flattened during collection.\n');
fprintf('Total valid positions to check: %d\n', numel(lon_all));

%% Check for parallel pool
pool = gcp('nocreate');
if isempty(pool)
    fprintf('Starting parallel pool...\n');
    pool = parpool('local');
else
    fprintf('Using existing parallel pool with %d workers.\n', pool.NumWorkers);
end

%% Check if INPOLY2 is available (with helper functions)
use_inpoly = false;
if exist('inpoly2', 'file') == 2
    try
        test_result = inpoly2([0 0], [0 0; 1 0; 1 1; 0 1]);
        use_inpoly = true;
        fprintf('INPOLY2 detected and working - will use fast algorithm.\n');
    catch
        fprintf('INPOLY2 found but missing helper functions (inpoly2_mat.m).\n');
        fprintf('Download complete package: https://github.com/dengwirda/inpoly\n');
        fprintf('Falling back to standard inpolygon.\n');
    end
else
    fprintf('INPOLY2 not found - using standard inpolygon.\n');
    fprintf('For faster performance, download: https://github.com/dengwirda/inpoly\n');
end

%% OPTIMIZATION 3: Parallel processing with bounding box pre-filtering
fprintf('\nProcessing with optimized parallel method...\n');
tic;

% Process in chunks to manage memory
chunk_size = 50000;
n_chunks = ceil(numel(lon_all) / chunk_size);

% Pre-allocate cell array for parallel results
chunk_connections = cell(n_chunks, 1);

% Copy variables for parfor (ensures proper broadcasting)
Xs_par = Xs;
Ys_par = Ys;
bbox_xmin_par = bbox_xmin;
bbox_xmax_par = bbox_xmax;
bbox_ymin_par = bbox_ymin;
bbox_ymax_par = bbox_ymax;
n_locations_par = n_locations;
use_inpoly_par = use_inpoly;

fprintf('Starting parallel processing of %d chunks...\n', n_chunks);

% PARALLEL LOOP over chunks
parfor chunk = 1:n_chunks
    % Get chunk indices
    idx_start = (chunk-1)*chunk_size + 1;
    idx_end = min(chunk*chunk_size, numel(lon_all));
    
    % Slice data for this chunk
    lon_chunk = lon_all(idx_start:idx_end);
    lat_chunk = lat_all(idx_start:idx_end);
    src_chunk = source_all(idx_start:idx_end);
    
    n_pts = numel(lon_chunk);
    dest_chunk = zeros(n_pts, 1);
    
    % Loop over polygons with bounding box pre-filtering
    for j = 1:n_locations_par
        % FAST: Bounding box check (vectorized)
        in_bbox = lon_chunk >= bbox_xmin_par(j) & lon_chunk <= bbox_xmax_par(j) & ...
                  lat_chunk >= bbox_ymin_par(j) & lat_chunk <= bbox_ymax_par(j);
        
        if ~any(in_bbox)
            continue;  % Skip this polygon entirely
        end
        
        % Only test points that passed bounding box
        candidates = find(in_bbox);
        
        if use_inpoly_par
            % INPOLY2: Fast algorithm with O((N+M)*log(N)) complexity
            in_poly = inpoly2([lon_chunk(candidates), lat_chunk(candidates)], ...
                             [Xs_par(j,1:4)', Ys_par(j,1:4)']);
        else
            % Standard MATLAB inpolygon
            in_poly = inpolygon(lon_chunk(candidates), lat_chunk(candidates), ...
                               Xs_par(j,:), Ys_par(j,:));
        end
        
        % Mark destinations
        dest_chunk(candidates(in_poly)) = j;
    end
    
    % Store connections as [source, destination] pairs
    valid_conn = dest_chunk > 0;
    chunk_connections{chunk} = [src_chunk(valid_conn), dest_chunk(valid_conn)];
end

fprintf('Parallel processing complete. Aggregating results...\n');
ConnMatrix = zeros(n_locations, n_locations);

for chunk = 1:n_chunks
    connections = chunk_connections{chunk};
    for i = 1:size(connections, 1)
        ConnMatrix(connections(i,1), connections(i,2)) = ...
            ConnMatrix(connections(i,1), connections(i,2)) + 1;
    end
end

elapsed = toc;
fprintf('Processing complete in %.1f seconds.\n', elapsed);

%% Normalize connectivity matrix (row normalization)
fprintf('\nNormalizing connectivity matrix...\n');

% Raw counts matrix
ConnMatrix_raw = ConnMatrix;

% Row normalization: each row sums to 1
% Interpretation: fraction of particle-time from source i spent in destination j
row_sums = sum(ConnMatrix, 2);
row_sums(row_sums == 0) = 1;  % Avoid division by zero for sources with no connections

ConnMatrix_normalized = ConnMatrix ./ row_sums;

fprintf('Normalization complete.\n');
fprintf('  Each row now sums to 1 (or 0 if no connections from that source).\n');

%% Summary
fprintf('\n=== RESULTS ===\n');
fprintf('RAW COUNTS:\n');
fprintf('  Total connections recorded: %d\n', nnz(ConnMatrix_raw));

[src_indices, dst_indices] = find(ConnMatrix_raw);
active_sources = unique(src_indices);
active_dests = unique(dst_indices);

fprintf('  Active sources: %d\n', numel(active_sources));
fprintf('  Active destinations: %d\n', numel(active_dests));
fprintf('  Self-recruitment counts: %d\n', sum(diag(ConnMatrix_raw)));

fprintf('\nNORMALIZED MATRIX:\n');
fprintf('  Non-zero entries: %d\n', nnz(ConnMatrix_normalized));
fprintf('  Self-recruitment probabilities range: [%.4f, %.4f]\n', ...
    min(diag(ConnMatrix_normalized(active_sources, active_sources))), ...
    max(diag(ConnMatrix_normalized(active_sources, active_sources))));

%% Network metrics (using raw counts)
out_degree = sum(ConnMatrix_raw, 2);
in_degree = sum(ConnMatrix_raw, 1)';
total_connectivity = out_degree + in_degree;

fprintf('\nTop 5 most connected locations:\n');
[~, sorted_idx] = sort(total_connectivity, 'descend');
for i = 1:min(5, nnz(total_connectivity))
    idx = sorted_idx(i);
    fprintf('  Location %d (ID=%d): Out=%d, In=%d, Total=%d\n', ...
        idx, unique_IDs(idx), out_degree(idx), in_degree(idx), ...
        total_connectivity(idx));
end

%% Visualizations
% 1. Connectivity matrix heatmap (NORMALIZED)
figure('Position', [100 100 800 600]);
active_locs = unique([active_sources; active_dests]);
sub_conn = ConnMatrix_normalized(active_locs, active_locs);
imagesc(sub_conn);
colorbar;
caxis([0 1]);  % Probabilities from 0 to 1
xlabel('Destination Polygon');
ylabel('Source Polygon');
title(sprintf('Normalized Connectivity - %s', datestr(calendar_date, 'dd-mmm-yyyy')));
axis square;

% 2. Raw counts heatmap for comparison
figure('Position', [920 100 800 600]);
sub_conn_raw = ConnMatrix_raw(active_locs, active_locs);
imagesc(log10(sub_conn_raw + 1));
colorbar;
xlabel('Destination Polygon');
ylabel('Source Polygon');
title(sprintf('Raw Connectivity (log10) - %s', datestr(calendar_date, 'dd-mmm-yyyy')));
axis square;

% 3. Spatial plot (using raw counts for node size)
figure('Position', [100 750 1000 800]);
hold on;

for i = 1:numel(active_locs)
    idx = active_locs(i);
    plot(Xs(idx,:), Ys(idx,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
end

centroids_x = mean(Xs(:,1:4), 2);
centroids_y = mean(Ys(:,1:4), 2);

conn_normalized = total_connectivity(active_locs);
conn_normalized = conn_normalized / max(conn_normalized) * 200 + 20;

scatter(centroids_x(active_locs), centroids_y(active_locs), ...
    conn_normalized, total_connectivity(active_locs), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

colorbar;
xlabel('Longitude');
ylabel('Latitude');
title(sprintf('Reef Connectivity - Raw Counts - %s', datestr(calendar_date, 'dd-mmm-yyyy')));
axis equal tight;
grid on;

% 4. Spatial plot (using normalized probabilities)
figure('Position', [920 750 1000 800]);
hold on;

for i = 1:numel(active_locs)
    idx = active_locs(i);
    plot(Xs(idx,:), Ys(idx,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
end

% Calculate normalized connectivity metrics
out_degree_norm = sum(ConnMatrix_normalized, 2);  % Sum of outgoing probabilities
in_degree_norm = sum(ConnMatrix_normalized, 1)';  % Sum of incoming probabilities
total_connectivity_norm = out_degree_norm + in_degree_norm;

conn_size_norm = total_connectivity_norm(active_locs);
conn_size_norm = conn_size_norm / max(conn_size_norm) * 200 + 20;

scatter(centroids_x(active_locs), centroids_y(active_locs), ...
    conn_size_norm, total_connectivity_norm(active_locs), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

colorbar;
xlabel('Longitude');
ylabel('Latitude');
title(sprintf('Reef Connectivity - Normalized - %s', datestr(calendar_date, 'dd-mmm-yyyy')));
axis equal tight;
grid on;

fprintf('\n=== COMPLETE ===\n');
fprintf('Optimizations applied:\n');
fprintf('  1. Bounding box pre-filtering\n');
fprintf('  2. Vectorized trajectory data\n');
fprintf('  3. Parallel processing (%d workers)\n', pool.NumWorkers);
if use_inpoly
    fprintf('  4. INPOLY2 fast algorithm (O((N+M)*log(N)))\n');
else
    fprintf('  4. Standard inpolygon (consider installing INPOLY for more speed)\n');
end
















% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%% OPTIMIZED PASSING-THROUGH CONNECTIVITY WITH PREPROCESSING %%%
% % Updated to save ALL available variables including time, depth, distance, exitcode
% 
% 
% 
% 
% % NOTE - this could still have an issue with indexing, but I think I fixed
% % that part....need to make sure 'time' is properly included as well
% 
% 
% 
% 
% 
% clear; clc;
% 
% %% ========== CONFIGURATION ==========
% DO_PREPROCESSING = true;  % Set to false to load existing preprocessed data
% YEAR = 2019;
% QUARTER = 1;  % 1, 2, 3, or 4
% 
% fprintf('=== REEF CONNECTIVITY ANALYSIS ===\n');
% fprintf('Year: %d, Quarter: %d\n', YEAR, QUARTER);
% fprintf('Preprocessing mode: %s\n', string(DO_PREPROCESSING));
% 
% %% Setup paths
% projectPath = matlab.project.rootProject().RootFolder;
% dataPath = fullfile(projectPath, 'data');
% outputPath = fullfile('D:\Dissertation\CMS_traj\output');
% 
% % Input/output directory based on quarter
% quarter_name = sprintf('Q%d_%d', QUARTER, YEAR);
% tempPath = fullfile('D:\Dissertation\CMS_traj', quarter_name);
% 
% % Output directory for preprocessed data
% preprocessPath = fullfile(outputPath, 'CMS_traj', quarter_name);
% if ~exist(preprocessPath, 'dir')
%     mkdir(preprocessPath);
% end
% 
% %% Load reef centroids and build polygons
% centroids = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
% unique_IDs = centroids(:,1);
% Xs = [centroids(:,8) centroids(:,10) centroids(:,12) centroids(:,14) centroids(:,8)];
% Ys = [centroids(:,9) centroids(:,11) centroids(:,13) centroids(:,15) centroids(:,9)];
% n_locations = size(centroids,1);
% 
% fprintf('Loaded %d reef polygons.\n', n_locations);
% 
% %% Load release file
% releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019_Q1.txt'));
% release_reef_IDs = releasefile(:,1);
% 
% %% ========== PREPROCESSING OR LOADING ==========
% 
% if DO_PREPROCESSING
%     %% PREPROCESSING MODE: Generate event files for all unique events
%     fprintf('\n=== PREPROCESSING MODE ===\n');
%     fprintf('This will process all files and create separate data files for each release event.\n');
% 
%     %% Step 1: Identify all unique release events (or load if exists)
%     unique_events_file = fullfile(preprocessPath, 'unique_events.mat');
% 
%     if exist(unique_events_file, 'file')
%         fprintf('\nLoading existing unique release events from: %s\n', unique_events_file);
%         load(unique_events_file, 'file_release_dates', 'unique_release_dates', 'trajlist');
%         n_events = length(unique_release_dates);
%         fprintf('Found %d unique release events (loaded from file).\n', n_events);
%     else
%         fprintf('\nStep 1: Scanning trajectory files for release events...\n');
%         trajlist = dir(fullfile(tempPath, 'traj*.nc'));
%         fprintf('Found %d trajectory files.\n', numel(trajlist));
% 
%         file_release_dates = zeros(numel(trajlist), 1);
% 
%         fprintf('Reading release dates from files...\n');
%         tic;
%         for i = 1:numel(trajlist)
%             if mod(i, 50) == 0
%                 elapsed = toc;
%                 rate = i / elapsed;
%                 remaining = (numel(trajlist) - i) / rate;
%                 fprintf('  Progress: %d/%d files (%.1f%%) - %.1f files/sec - ETA: %.1f sec\n', ...
%                     i, numel(trajlist), 100*i/numel(trajlist), rate, remaining);
%             end
% 
%             filename = fullfile(tempPath, trajlist(i).name);
%             ncid = -1;
%             try
%                 ncid = netcdf.open(filename, 'NC_NOWRITE');
%                 varid = netcdf.inqVarID(ncid, 'releasedate');
%                 releasedate = netcdf.getVar(ncid, varid, 0, 1, 'double');
%                 file_release_dates(i) = releasedate;
%                 netcdf.close(ncid);
%                 ncid = -1;
%             catch ME
%                 warning('Could not read %s: %s', trajlist(i).name, ME.message);
%                 file_release_dates(i) = NaN;
%                 if ncid ~= -1
%                     try netcdf.close(ncid); catch; end
%                 end
%             end
%         end
%         elapsed_total = toc;
%         fprintf('Completed reading %d files in %.2f seconds (%.1f files/sec)\n', ...
%             numel(trajlist), elapsed_total, numel(trajlist)/elapsed_total);
% 
%         file_release_dates = file_release_dates(~isnan(file_release_dates));
%         unique_release_dates = unique(file_release_dates);
%         n_events = length(unique_release_dates);
% 
%         fprintf('\nFound %d unique release events.\n', n_events);
% 
%         % Save unique events info
%         fprintf('Saving unique release events to: %s\n', unique_events_file);
%         save(unique_events_file, 'file_release_dates', 'unique_release_dates', 'trajlist', '-v7.3');
%     end
% 
%     %% Step 2: Process each unique event and save to separate files
%     fprintf('\nStep 2: Processing each release event and saving data...\n');
% 
%     for event_idx = 1:n_events
%         target_date = unique_release_dates(event_idx);
%         calendar_date = datetime(target_date, 'ConvertFrom', 'juliandate');
% 
%         % Create descriptive filename based on date
%         date_str = datestr(calendar_date, 'yyyy_mmmdd_HHMMSS');
%         event_filename = fullfile(preprocessPath, sprintf('event_%s.mat', date_str));
% 
%         % Skip if this event file already exists
%         if exist(event_filename, 'file')
%             fprintf('\n--- Event %d/%d already exists, skipping ---\n', event_idx, n_events);
%             fprintf('  File: %s\n', event_filename);
%             continue;
%         end
% 
%         fprintf('\n--- Processing Event %d/%d ---\n', event_idx, n_events);
%         fprintf('  Julian date: %.2f\n', target_date);
%         fprintf('  Calendar date: %s\n', datestr(calendar_date, 'dd-mmm-yyyy HH:MM:SS'));
% 
%         % Identify files containing this event
%         files_with_event = abs(file_release_dates - target_date) < 0.01;
%         relevant_files = find(files_with_event);
%         n_relevant_files = numel(relevant_files);
% 
%         fprintf('  Found in %d files\n', n_relevant_files);
% 
%         % Collect particles from relevant files - NOW INCLUDING ALL VARIABLES
%         lon_cells = cell(n_relevant_files, 1);
%         lat_cells = cell(n_relevant_files, 1);
%         depth_cells = cell(n_relevant_files, 1);
%         distance_cells = cell(n_relevant_files, 1);
%         exitcode_cells = cell(n_relevant_files, 1);
%         src_cells = cell(n_relevant_files, 1);
% 
%         % Time is shared across all particles, so we'll store it once
%         time_array = [];
% 
%         tic;
%         for idx = 1:n_relevant_files
%             i = relevant_files(idx);
%             filename = fullfile(tempPath, trajlist(i).name);
% 
%             try
%                 % Read all variables
%                 lon = ncread(filename, 'lon');
%                 lat = ncread(filename, 'lat');
%                 location = ncread(filename, 'location');
%                 releasedate = ncread(filename, 'releasedate');
% 
%                 % Read time array (only once since it's the same for all particles)
%                 if isempty(time_array)
%                     time_array = ncread(filename, 'time');
%                 end
% 
%                 event_mask = abs(releasedate - target_date) < 0.01;
%                 if ~any(event_mask)
%                     continue;
%                 end
% 
%                 % Safety check: ensure event_mask doesn't exceed array bounds
%                 n_particles = size(lon, 2);
%                 if length(event_mask) > n_particles
%                     warning('event_mask length (%d) exceeds lon particles (%d) in %s - truncating', ...
%                         length(event_mask), n_particles, trajlist(i).name);
%                     event_mask = event_mask(1:n_particles);
%                 end
% 
%                 lon(lon > 180) = lon(lon > 180) - 360;
% 
%                 % Read additional variables - with dimension checking
%                 try
%                     depth = ncread(filename, 'depth');
%                     if size(depth, 2) == n_particles
%                         depth_event = depth(:, event_mask);
%                     else
%                         depth_event = nan(size(lon, 1), sum(event_mask));
%                     end
%                 catch
%                     depth_event = nan(size(lon, 1), sum(event_mask));
%                 end
% 
%                 try
%                     distance = ncread(filename, 'distance');
%                     if size(distance, 2) == n_particles
%                         distance_event = distance(:, event_mask);
%                     else
%                         distance_event = nan(size(lon, 1), sum(event_mask));
%                     end
%                 catch
%                     distance_event = nan(size(lon, 1), sum(event_mask));
%                 end
% 
%                 try
%                     exitcode = ncread(filename, 'exitcode');
%                     if size(exitcode, 2) == n_particles
%                         exitcode_event = exitcode(:, event_mask);
%                     else
%                         exitcode_event = nan(size(lon, 1), sum(event_mask));
%                     end
%                 catch
%                     exitcode_event = nan(size(lon, 1), sum(event_mask));
%                 end
% 
%                 % Get particles for this event (KEEP AS 2D!)
%                 lon_event = lon(:, event_mask);         % [timesteps × particles]
%                 lat_event = lat(:, event_mask);         % [timesteps × particles]
%                 loc_event = location(event_mask);       % [1 × particles]
% 
%                 % Map to source indices
%                 reef_IDs = release_reef_IDs(loc_event);
%                 [~, src_idx] = ismember(reef_IDs, unique_IDs);
% 
%                 % Store 2D matrices (NOT FLATTENED!)
%                 lon_cells{idx} = lon_event;
%                 lat_cells{idx} = lat_event;
%                 depth_cells{idx} = depth_event;
%                 distance_cells{idx} = distance_event;
%                 exitcode_cells{idx} = exitcode_event;
%                 src_cells{idx} = src_idx';  % Row vector of source indices
% 
%             catch ME
%                 warning('Error processing %s: %s', trajlist(i).name, ME.message);
%             end
%         end
% 
%         % Concatenate horizontally across particles (not vertically!)
%         all_lon = horzcat(lon_cells{:});           % [timesteps × total_particles]
%         all_lat = horzcat(lat_cells{:});           % [timesteps × total_particles]
%         all_depth = horzcat(depth_cells{:});       % [timesteps × total_particles]
%         all_distance = horzcat(distance_cells{:}); % [timesteps × total_particles]
%         all_exitcode = horzcat(exitcode_cells{:}); % [timesteps × total_particles]
%         all_source_idx = horzcat(src_cells{:});    % [1 × total_particles]
% 
%         elapsed = toc;
%         n_timesteps = size(all_lon, 1);
%         n_particles = size(all_lon, 2);
%         fprintf('  Collected %d particles × %d timesteps = %d positions in %.1f seconds\n', ...
%             n_particles, n_timesteps, numel(all_lon), elapsed);
% 
%         % Organize into bigstruct-style array by source location
%         fprintf('  Organizing into struct array by source location...\n');
%         unique_sources = unique(all_source_idx);
%         n_sources = length(unique_sources);
% 
%         event_data = struct('source_reef_id', {}, 'source_reef_idx', {}, ...
%                            'lon', {}, 'lat', {}, 'depth', {}, 'distance', {}, ...
%                            'exitcode', {}, 'time', {}, ...
%                            'n_particles', {}, 'n_timesteps', {});
% 
%         for s = 1:n_sources
%             src_id = unique_sources(s);
% 
%             % Get particles from this source
%             particle_mask = all_source_idx == src_id;
% 
%             event_data(s).source_reef_idx = src_id;
%             event_data(s).source_reef_id = unique_IDs(src_id);
%             event_data(s).lon = all_lon(:, particle_mask);          % [timesteps × particles]
%             event_data(s).lat = all_lat(:, particle_mask);          % [timesteps × particles]
%             event_data(s).depth = all_depth(:, particle_mask);      % [timesteps × particles]
%             event_data(s).distance = all_distance(:, particle_mask); % [timesteps × particles]
%             event_data(s).exitcode = all_exitcode(:, particle_mask); % [timesteps × particles]
%             event_data(s).time = time_array;                        % [timesteps × 1] - same for all
%             event_data(s).n_particles = sum(particle_mask);
%             event_data(s).n_timesteps = n_timesteps;
%             event_data(s).releasedate = target_date;
%             event_data(s).calendar_date = calendar_date;
%         end
% 
%         fprintf('  Organized into %d source locations (bigstruct-style)\n', n_sources);
%         fprintf('  Variables saved: lon, lat, depth, distance, exitcode, time\n');
% 
%         % Save metadata separately
%         event_metadata = struct();
%         event_metadata.event_idx = event_idx;
%         event_metadata.target_date = target_date;
%         event_metadata.calendar_date = calendar_date;
%         event_metadata.n_sources = n_sources;
%         event_metadata.n_particles_total = n_particles;
%         event_metadata.n_timesteps = n_timesteps;
%         event_metadata.time_array = time_array;  % Store time in metadata too
%         event_metadata.time_units = 'seconds since release';
%         event_metadata.variables_saved = {'lon', 'lat', 'depth', 'distance', 'exitcode', 'time'};
% 
%         save(event_filename, 'event_data', 'event_metadata', '-v7.3');
%         fprintf('  Saved to: %s\n', event_filename);
%     end
% 
%     % Save/update metadata file
%     metadata = struct();
%     metadata.year = YEAR;
%     metadata.quarter = QUARTER;
%     metadata.n_events = n_events;
%     metadata.unique_release_dates = unique_release_dates;
%     metadata.n_locations = n_locations;
%     metadata.preprocessing_date = datetime('now');
%     metadata.variables_saved = {'lon', 'lat', 'depth', 'distance', 'exitcode', 'time'};
% 
%     % Add list of event files for reference
%     event_files = dir(fullfile(preprocessPath, 'event_*.mat'));
%     metadata.event_files = {event_files.name}';
% 
%     metadata_file = fullfile(preprocessPath, 'metadata.mat');
%     save(metadata_file, 'metadata', '-v7.3');
% 
%     fprintf('\n=== PREPROCESSING COMPLETE ===\n');
%     fprintf('Processed %d events and saved to: %s\n', n_events, preprocessPath);
%     fprintf('Variables saved: lon, lat, depth, distance, exitcode, time\n');
%     fprintf('Event files:\n');
%     for i = 1:length(metadata.event_files)
%         fprintf('  %s\n', metadata.event_files{i});
%     end
%     fprintf('\nSet DO_PREPROCESSING = false to run analysis on these events.\n');
% 
%     return;  % Exit after preprocessing
% 
% else
%     %% ANALYSIS MODE: Load preprocessed data
%     fprintf('\n=== ANALYSIS MODE ===\n');
% 
%     % Load metadata
%     metadata_file = fullfile(preprocessPath, 'metadata.mat');
%     if ~exist(metadata_file, 'file')
%         error('Metadata file not found. Run with DO_PREPROCESSING = true first.');
%     end
% 
%     load(metadata_file, 'metadata');
%     n_events = metadata.n_events;
% 
%     fprintf('Found %d preprocessed events in: %s\n', n_events, preprocessPath);
%     fprintf('Variables available: %s\n', strjoin(metadata.variables_saved, ', '));
% 
%     % Find all event files
%     event_files = dir(fullfile(preprocessPath, 'event_*.mat'));
%     if isempty(event_files)
%         error('No event files found. Run with DO_PREPROCESSING = true first.');
%     end
% 
%     % Pick one event to analyze
%     event_idx = randi(length(event_files));
%     event_filename = fullfile(preprocessPath, event_files(event_idx).name);
% 
%     fprintf('\nLoading Event %d/%d: %s\n', event_idx, length(event_files), event_files(event_idx).name);
%     load(event_filename, 'event_data', 'event_metadata');
% 
%     fprintf('  Julian date: %.2f\n', event_metadata.target_date);
%     fprintf('  Calendar date: %s\n', datestr(event_metadata.calendar_date, 'dd-mmm-yyyy HH:MM:SS'));
%     fprintf('  Source locations: %d\n', event_metadata.n_sources);
%     fprintf('  Total particles: %d × %d timesteps\n', ...
%         event_metadata.n_particles_total, event_metadata.n_timesteps);
%     fprintf('  Time range: %d to %d seconds\n', ...
%         min(event_metadata.time_array), max(event_metadata.time_array));
% 
%     target_date = event_metadata.target_date;
%     calendar_date = event_metadata.calendar_date;
% 
%     % Flatten bigstruct-style array for polygon checking
%     fprintf('\nFlattening bigstruct for connectivity analysis...\n');
%     fprintf('  event_data is a struct array with %d elements (one per source reef)\n', length(event_data));
% 
%     % Pre-allocate for speed
%     total_positions = event_metadata.n_particles_total * event_metadata.n_timesteps;
%     all_lon = zeros(total_positions, 1);
%     all_lat = zeros(total_positions, 1);
%     all_source_idx = zeros(total_positions, 1);
% 
%     pos_idx = 1;
%     for s = 1:length(event_data)
%         % Get data from this source
%         lon_src = event_data(s).lon(:);  % Flatten [timesteps × particles]
%         lat_src = event_data(s).lat(:);
%         n_pos = length(lon_src);
% 
%         % Store in flat arrays
%         all_lon(pos_idx:pos_idx+n_pos-1) = lon_src;
%         all_lat(pos_idx:pos_idx+n_pos-1) = lat_src;
%         all_source_idx(pos_idx:pos_idx+n_pos-1) = event_data(s).source_reef_idx;
% 
%         pos_idx = pos_idx + n_pos;
%     end
% 
%     fprintf('  Flattened %d source reefs into %d total positions\n', ...
%         length(event_data), total_positions);
% end
% 
% %% ========== ANALYSIS CONTINUES FROM HERE ==========
% fprintf('\n=== STARTING CONNECTIVITY ANALYSIS ===\n');
% 
% %% Remove NaN positions
% valid = ~isnan(all_lon) & ~isnan(all_lat) & all_source_idx > 0;
% lon_all = all_lon(valid);
% lat_all = all_lat(valid);
% source_all = all_source_idx(valid);
% 
% fprintf('Valid positions after filtering: %d\n', numel(lon_all));
% 
% %% OPTIMIZATION 1: Pre-compute bounding boxes for quick rejection
% fprintf('\nPre-computing polygon bounding boxes...\n');
% bbox_xmin = min(Xs(:,1:4), [], 2);
% bbox_xmax = max(Xs(:,1:4), [], 2);
% bbox_ymin = min(Ys(:,1:4), [], 2);
% bbox_ymax = max(Ys(:,1:4), [], 2);
% 
% %% OPTIMIZATION 2: Data already prepared
% fprintf('Data already flattened during collection.\n');
% fprintf('Total valid positions to check: %d\n', numel(lon_all));
% 
% %% Check for parallel pool
% pool = gcp('nocreate');
% if isempty(pool)
%     fprintf('Starting parallel pool...\n');
%     pool = parpool('local');
% else
%     fprintf('Using existing parallel pool with %d workers.\n', pool.NumWorkers);
% end
% 
% %% Check if INPOLY2 is available (with helper functions)
% use_inpoly = false;
% if exist('inpoly2', 'file') == 2
%     try
%         test_result = inpoly2([0 0], [0 0; 1 0; 1 1; 0 1]);
%         use_inpoly = true;
%         fprintf('INPOLY2 detected and working - will use fast algorithm.\n');
%     catch
%         fprintf('INPOLY2 found but missing helper functions (inpoly2_mat.m).\n');
%         fprintf('Download complete package: https://github.com/dengwirda/inpoly\n');
%         fprintf('Falling back to standard inpolygon.\n');
%     end
% else
%     fprintf('INPOLY2 not found - using standard inpolygon.\n');
%     fprintf('For faster performance, download: https://github.com/dengwirda/inpoly\n');
% end
% 
% %% OPTIMIZATION 3: Parallel processing with bounding box pre-filtering
% fprintf('\nProcessing with optimized parallel method...\n');
% tic;
% 
% % Process in chunks to manage memory
% chunk_size = 50000;
% n_chunks = ceil(numel(lon_all) / chunk_size);
% 
% % Pre-allocate cell array for parallel results
% chunk_connections = cell(n_chunks, 1);
% 
% % Copy variables for parfor (ensures proper broadcasting)
% Xs_par = Xs;
% Ys_par = Ys;
% bbox_xmin_par = bbox_xmin;
% bbox_xmax_par = bbox_xmax;
% bbox_ymin_par = bbox_ymin;
% bbox_ymax_par = bbox_ymax;
% n_locations_par = n_locations;
% use_inpoly_par = use_inpoly;
% 
% fprintf('Starting parallel processing of %d chunks...\n', n_chunks);
% 
% % PARALLEL LOOP over chunks
% parfor chunk = 1:n_chunks
%     % Get chunk indices
%     idx_start = (chunk-1)*chunk_size + 1;
%     idx_end = min(chunk*chunk_size, numel(lon_all));
% 
%     % Slice data for this chunk
%     lon_chunk = lon_all(idx_start:idx_end);
%     lat_chunk = lat_all(idx_start:idx_end);
%     src_chunk = source_all(idx_start:idx_end);
% 
%     n_pts = numel(lon_chunk);
%     dest_chunk = zeros(n_pts, 1);
% 
%     % Loop over polygons with bounding box pre-filtering
%     for j = 1:n_locations_par
%         % FAST: Bounding box check (vectorized)
%         in_bbox = lon_chunk >= bbox_xmin_par(j) & lon_chunk <= bbox_xmax_par(j) & ...
%                   lat_chunk >= bbox_ymin_par(j) & lat_chunk <= bbox_ymax_par(j);
% 
%         if ~any(in_bbox)
%             continue;  % Skip this polygon entirely
%         end
% 
%         % Only test points that passed bounding box
%         candidates = find(in_bbox);
% 
%         if use_inpoly_par
%             % INPOLY2: Fast algorithm with O((N+M)*log(N)) complexity
%             in_poly = inpoly2([lon_chunk(candidates), lat_chunk(candidates)], ...
%                              [Xs_par(j,1:4)', Ys_par(j,1:4)']);
%         else
%             % Standard MATLAB inpolygon
%             in_poly = inpolygon(lon_chunk(candidates), lat_chunk(candidates), ...
%                                Xs_par(j,:), Ys_par(j,:));
%         end
% 
%         % Mark destinations
%         dest_chunk(candidates(in_poly)) = j;
%     end
% 
%     % Store connections as [source, destination] pairs
%     valid_conn = dest_chunk > 0;
%     chunk_connections{chunk} = [src_chunk(valid_conn), dest_chunk(valid_conn)];
% end
% 
% fprintf('Parallel processing complete. Aggregating results...\n');
% ConnMatrix = zeros(n_locations, n_locations);
% 
% for chunk = 1:n_chunks
%     connections = chunk_connections{chunk};
%     for i = 1:size(connections, 1)
%         ConnMatrix(connections(i,1), connections(i,2)) = ...
%             ConnMatrix(connections(i,1), connections(i,2)) + 1;
%     end
% end
% 
% elapsed = toc;
% fprintf('Processing complete in %.1f seconds.\n', elapsed);
% 
% %% Normalize connectivity matrix (row normalization)
% fprintf('\nNormalizing connectivity matrix...\n');
% 
% % Raw counts matrix
% ConnMatrix_raw = ConnMatrix;
% 
% % Row normalization: each row sums to 1
% % Interpretation: fraction of particle-time from source i spent in destination j
% row_sums = sum(ConnMatrix, 2);
% row_sums(row_sums == 0) = 1;  % Avoid division by zero for sources with no connections
% 
% ConnMatrix_normalized = ConnMatrix ./ row_sums;
% 
% fprintf('Normalization complete.\n');
% fprintf('  Each row now sums to 1 (or 0 if no connections from that source).\n');
% 
% %% Summary
% fprintf('\n=== RESULTS ===\n');
% fprintf('RAW COUNTS:\n');
% fprintf('  Total connections recorded: %d\n', nnz(ConnMatrix_raw));
% 
% [src_indices, dst_indices] = find(ConnMatrix_raw);
% active_sources = unique(src_indices);
% active_dests = unique(dst_indices);
% 
% fprintf('  Active sources: %d\n', numel(active_sources));
% fprintf('  Active destinations: %d\n', numel(active_dests));
% fprintf('  Self-recruitment counts: %d\n', sum(diag(ConnMatrix_raw)));
% 
% fprintf('\nNORMALIZED MATRIX:\n');
% fprintf('  Non-zero entries: %d\n', nnz(ConnMatrix_normalized));
% fprintf('  Self-recruitment probabilities range: [%.4f, %.4f]\n', ...
%     min(diag(ConnMatrix_normalized(active_sources, active_sources))), ...
%     max(diag(ConnMatrix_normalized(active_sources, active_sources))));
% 
% %% Network metrics (using raw counts)
% out_degree = sum(ConnMatrix_raw, 2);
% in_degree = sum(ConnMatrix_raw, 1)';
% total_connectivity = out_degree + in_degree;
% 
% fprintf('\nTop 5 most connected locations:\n');
% [~, sorted_idx] = sort(total_connectivity, 'descend');
% for i = 1:min(5, nnz(total_connectivity))
%     idx = sorted_idx(i);
%     fprintf('  Location %d (ID=%d): Out=%d, In=%d, Total=%d\n', ...
%         idx, unique_IDs(idx), out_degree(idx), in_degree(idx), ...
%         total_connectivity(idx));
% end
% 
% %% Visualizations
% % 1. Connectivity matrix heatmap (NORMALIZED)
% figure('Position', [100 100 800 600]);
% active_locs = unique([active_sources; active_dests]);
% sub_conn = ConnMatrix_normalized(active_locs, active_locs);
% imagesc(sub_conn);
% colorbar;
% caxis([0 1]);  % Probabilities from 0 to 1
% xlabel('Destination Polygon');
% ylabel('Source Polygon');
% title(sprintf('Normalized Connectivity - %s', datestr(calendar_date, 'dd-mmm-yyyy')));
% axis square;
% 
% % 2. Raw counts heatmap for comparison
% figure('Position', [920 100 800 600]);
% sub_conn_raw = ConnMatrix_raw(active_locs, active_locs);
% imagesc(log10(sub_conn_raw + 1));
% colorbar;
% xlabel('Destination Polygon');
% ylabel('Source Polygon');
% title(sprintf('Raw Connectivity (log10) - %s', datestr(calendar_date, 'dd-mmm-yyyy')));
% axis square;
% 
% % 3. Spatial plot (using raw counts for node size)
% figure('Position', [100 750 1000 800]);
% hold on;
% 
% for i = 1:numel(active_locs)
%     idx = active_locs(i);
%     plot(Xs(idx,:), Ys(idx,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
% end
% 
% centroids_x = mean(Xs(:,1:4), 2);
% centroids_y = mean(Ys(:,1:4), 2);
% 
% conn_normalized = total_connectivity(active_locs);
% conn_normalized = conn_normalized / max(conn_normalized) * 200 + 20;
% 
% scatter(centroids_x(active_locs), centroids_y(active_locs), ...
%     conn_normalized, total_connectivity(active_locs), 'filled', ...
%     'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% 
% colorbar;
% xlabel('Longitude');
% ylabel('Latitude');
% title(sprintf('Reef Connectivity - Raw Counts - %s', datestr(calendar_date, 'dd-mmm-yyyy')));
% axis equal tight;
% grid on;
% 
% % 4. Spatial plot (using normalized probabilities)
% figure('Position', [920 750 1000 800]);
% hold on;
% 
% for i = 1:numel(active_locs)
%     idx = active_locs(i);
%     plot(Xs(idx,:), Ys(idx,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
% end
% 
% % Calculate normalized connectivity metrics
% out_degree_norm = sum(ConnMatrix_normalized, 2);  % Sum of outgoing probabilities
% in_degree_norm = sum(ConnMatrix_normalized, 1)';  % Sum of incoming probabilities
% total_connectivity_norm = out_degree_norm + in_degree_norm;
% 
% conn_size_norm = total_connectivity_norm(active_locs);
% conn_size_norm = conn_size_norm / max(conn_size_norm) * 200 + 20;
% 
% scatter(centroids_x(active_locs), centroids_y(active_locs), ...
%     conn_size_norm, total_connectivity_norm(active_locs), 'filled', ...
%     'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% 
% colorbar;
% xlabel('Longitude');
% ylabel('Latitude');
% title(sprintf('Reef Connectivity - Normalized - %s', datestr(calendar_date, 'dd-mmm-yyyy')));
% axis equal tight;
% grid on;
% 
% fprintf('\n=== COMPLETE ===\n');
% fprintf('Optimizations applied:\n');
% fprintf('  1. Bounding box pre-filtering\n');
% fprintf('  2. Vectorized trajectory data\n');
% fprintf('  3. Parallel processing (%d workers)\n', pool.NumWorkers);
% if use_inpoly
%     fprintf('  4. INPOLY2 fast algorithm (O((N+M)*log(N)))\n');
% else
%     fprintf('  4. Standard inpolygon (consider installing INPOLY for more speed)\n');
% end
























%% maybe even better ?? but slow as shit. this is one that has produced 5 events so far





%%% MEMORY-OPTIMIZED REEF CONNECTIVITY PREPROCESSING %%%
clear; clc;

%% ========== CONFIGURATION ==========
DO_PREPROCESSING = true;
YEAR = 2019;
QUARTER = 1;

fprintf('=== REEF CONNECTIVITY ANALYSIS ===\n');
fprintf('Year: %d, Quarter: %d | Preprocessing: %s\n', YEAR, QUARTER, string(DO_PREPROCESSING));

%% Setup paths
projectPath = matlab.project.rootProject().RootFolder;
dataPath = fullfile(projectPath, 'data');
outputPath = fullfile('D:\Dissertation\CMS_traj\output');
quarter_name = sprintf('Q%d_%d', QUARTER, YEAR);
tempPath = fullfile('D:\Dissertation\CMS_traj', quarter_name);
preprocessPath = fullfile(outputPath, 'CMS_traj', quarter_name);
if ~exist(preprocessPath, 'dir'), mkdir(preprocessPath); end

%% Load reef data
centroids = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
unique_IDs = centroids(:,1);
Xs = [centroids(:,8) centroids(:,10) centroids(:,12) centroids(:,14) centroids(:,8)];
Ys = [centroids(:,9) centroids(:,11) centroids(:,13) centroids(:,15) centroids(:,9)];
n_locations = size(centroids,1);
fprintf('Loaded %d reef polygons.\n', n_locations);

releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019_Q1.txt'));
release_reef_IDs = releasefile(:,1);

%% ========== PREPROCESSING OR ANALYSIS ==========
if DO_PREPROCESSING
    fprintf('\n=== PREPROCESSING MODE ===\n');
    
    %% Step 1: Get unique release events
    unique_events_file = fullfile(preprocessPath, 'unique_events.mat');
    
    if exist(unique_events_file, 'file')
        fprintf('Loading existing unique events...\n');
        load(unique_events_file, 'file_release_dates', 'unique_release_dates', 'trajlist');
    else
        fprintf('Scanning trajectory files for release events...\n');
        trajlist = dir(fullfile(tempPath, 'traj*.nc'));
        fprintf('Found %d trajectory files.\n', numel(trajlist));
        
        file_release_dates = zeros(numel(trajlist), 1);
        tic;
        for i = 1:numel(trajlist)
            if mod(i, 50) == 0
                fprintf('  Progress: %d/%d (%.1f%%) - %.1f files/sec\n', ...
                    i, numel(trajlist), 100*i/numel(trajlist), i/toc);
            end
            
            ncid = -1;
            try
                ncid = netcdf.open(fullfile(tempPath, trajlist(i).name), 'NC_NOWRITE');
                varid = netcdf.inqVarID(ncid, 'releasedate');
                file_release_dates(i) = netcdf.getVar(ncid, varid, 0, 1, 'double');
                netcdf.close(ncid); ncid = -1;
            catch ME
                warning('Could not read %s: %s', trajlist(i).name, ME.message);
                file_release_dates(i) = NaN;
                if ncid ~= -1, try netcdf.close(ncid); catch; end; end
            end
        end
        fprintf('Completed in %.2f seconds\n', toc);
        
        file_release_dates = file_release_dates(~isnan(file_release_dates));
        unique_release_dates = unique(file_release_dates);
        save(unique_events_file, 'file_release_dates', 'unique_release_dates', 'trajlist', '-v7.3');
    end
    
    n_events = length(unique_release_dates);
    fprintf('Found %d unique release events.\n', n_events);
    
    %% Step 2: Process each event (MEMORY OPTIMIZED)
    fprintf('\n=== Processing events with direct-to-struct approach ===\n');
    
    for event_idx = 1:n_events
        target_date = unique_release_dates(event_idx);
        calendar_date = datetime(target_date, 'ConvertFrom', 'juliandate');
        date_str = datestr(calendar_date, 'yyyy_mmmdd_HHMMSS');
        event_filename = fullfile(preprocessPath, sprintf('event_%s.mat', date_str));
        
        if exist(event_filename, 'file')
            fprintf('\n--- Event %d/%d exists, skipping ---\n', event_idx, n_events);
            continue;
        end
        
        fprintf('\n--- Event %d/%d: %s ---\n', event_idx, n_events, datestr(calendar_date, 'dd-mmm-yyyy HH:MM:SS'));
        
        % Find relevant files
        relevant_files = find(abs(file_release_dates - target_date) < 0.01);
        fprintf('  Found in %d files\n', numel(relevant_files));
        
        % Initialize Map for direct accumulation (MEMORY OPTIMIZATION)
        time_array = [];
        source_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
        
        tic;
        for idx = 1:numel(relevant_files)
            file_idx = relevant_files(idx);
            filename = fullfile(tempPath, trajlist(file_idx).name);
            
            try
                % Read essential variables
                lon = ncread(filename, 'lon');
                lat = ncread(filename, 'lat');
                location = ncread(filename, 'location');
                releasedate = ncread(filename, 'releasedate');
                
                if isempty(time_array)
                    time_array = ncread(filename, 'time');
                end
                
                event_mask = abs(releasedate - target_date) < 0.01;
                if ~any(event_mask), continue; end
                
                n_particles = size(lon, 2);
                if length(event_mask) > n_particles
                    event_mask = event_mask(1:n_particles);
                end
                
                lon(lon > 180) = lon(lon > 180) - 360;
                lon_event = lon(:, event_mask);
                lat_event = lat(:, event_mask);
                loc_event = location(event_mask);
                
                % Read optional variables with dimension checking
                try
                    depth = ncread(filename, 'depth');
                    if size(depth, 2) == n_particles
                        depth_event = depth(:, event_mask);
                    else
                        depth_event = nan(size(lon_event));
                    end
                catch
                    depth_event = nan(size(lon_event));
                end
                
                try
                    distance = ncread(filename, 'distance');
                    if size(distance, 2) == n_particles
                        distance_event = distance(:, event_mask);
                    else
                        distance_event = nan(size(lon_event));
                    end
                catch
                    distance_event = nan(size(lon_event));
                end
                
                try
                    exitcode = ncread(filename, 'exitcode');
                    if size(exitcode, 2) == n_particles
                        exitcode_event = exitcode(:, event_mask);
                    else
                        exitcode_event = nan(size(lon_event));
                    end
                catch
                    exitcode_event = nan(size(lon_event));
                end
                
                % Map particles to source reef indices
                reef_IDs = release_reef_IDs(loc_event);
                [~, src_indices] = ismember(reef_IDs, unique_IDs);
                
                % Add each particle directly to its source bin
                for p = 1:length(src_indices)
                    source_id = int32(src_indices(p));
                    
                    if isKey(source_map, source_id)
                        bin_data = source_map(source_id);
                        bin_data.lon = [bin_data.lon, lon_event(:, p)];
                        bin_data.lat = [bin_data.lat, lat_event(:, p)];
                        bin_data.depth = [bin_data.depth, depth_event(:, p)];
                        bin_data.distance = [bin_data.distance, distance_event(:, p)];
                        bin_data.exitcode = [bin_data.exitcode, exitcode_event(:, p)];
                        bin_data.n_particles = bin_data.n_particles + 1;
                        source_map(source_id) = bin_data;
                    else
                        bin_data = struct();
                        bin_data.source_reef_idx = double(source_id);
                        bin_data.source_reef_id = unique_IDs(double(source_id));
                        bin_data.lon = lon_event(:, p);
                        bin_data.lat = lat_event(:, p);
                        bin_data.depth = depth_event(:, p);
                        bin_data.distance = distance_event(:, p);
                        bin_data.exitcode = exitcode_event(:, p);
                        bin_data.n_particles = 1;
                        source_map(source_id) = bin_data;
                    end
                end
                
            catch ME
                warning('Error processing %s: %s', trajlist(file_idx).name, ME.message);
            end
        end
        
        % Convert Map to struct array
        source_keys = cell2mat(source_map.keys());
        n_sources = length(source_keys);
        event_data = repmat(struct('source_reef_id', [], 'source_reef_idx', [], ...
                                   'lon', [], 'lat', [], 'depth', [], 'distance', [], ...
                                   'exitcode', [], 'time', [], 'n_particles', [], 'n_timesteps', []), ...
                           n_sources, 1);
        
        for s = 1:n_sources
            bin_data = source_map(source_keys(s));
            event_data(s).source_reef_idx = bin_data.source_reef_idx;
            event_data(s).source_reef_id = bin_data.source_reef_id;
            event_data(s).lon = bin_data.lon;
            event_data(s).lat = bin_data.lat;
            event_data(s).depth = bin_data.depth;
            event_data(s).distance = bin_data.distance;
            event_data(s).exitcode = bin_data.exitcode;
            event_data(s).time = time_array;
            event_data(s).n_particles = bin_data.n_particles;
            event_data(s).n_timesteps = length(time_array);
            event_data(s).releasedate = target_date;
            event_data(s).calendar_date = calendar_date;
        end
        
        fprintf('  Organized %d source locations in %.1f sec\n', n_sources, toc);
        
        % Save with metadata
        event_metadata = struct('event_idx', event_idx, 'target_date', target_date, ...
                               'calendar_date', calendar_date, 'n_sources', n_sources, ...
                               'n_particles_total', sum([event_data.n_particles]), ...
                               'n_timesteps', length(time_array), 'time_array', time_array, ...
                               'time_units', 'seconds since release', ...
                               'variables_saved', {{'lon', 'lat', 'depth', 'distance', 'exitcode', 'time'}});
        
        save(event_filename, 'event_data', 'event_metadata', '-v7.3');
        fprintf('  Saved: %s\n', event_filename);
    end
    
    % Save overall metadata
    metadata = struct('year', YEAR, 'quarter', QUARTER, 'n_events', n_events, ...
                     'unique_release_dates', unique_release_dates, 'n_locations', n_locations, ...
                     'preprocessing_date', datetime('now'), ...
                     'variables_saved', {{'lon', 'lat', 'depth', 'distance', 'exitcode', 'time'}});
    
    event_files = dir(fullfile(preprocessPath, 'event_*.mat'));
    metadata.event_files = {event_files.name}';
    save(fullfile(preprocessPath, 'metadata.mat'), 'metadata', '-v7.3');
    
    fprintf('\n=== PREPROCESSING COMPLETE ===\n');
    fprintf('Processed %d events | Saved to: %s\n', n_events, preprocessPath);
    fprintf('Set DO_PREPROCESSING = false to run analysis.\n');
    return;
    
else
    %% ANALYSIS MODE
    fprintf('\n=== ANALYSIS MODE ===\n');
    
    metadata_file = fullfile(preprocessPath, 'metadata.mat');
    if ~exist(metadata_file, 'file')
        error('Metadata file not found. Run with DO_PREPROCESSING = true first.');
    end
    
    load(metadata_file, 'metadata');
    event_files = dir(fullfile(preprocessPath, 'event_*.mat'));
    
    fprintf('Found %d preprocessed events\n', length(event_files));
    fprintf('Variables: %s\n', strjoin(metadata.variables_saved, ', '));
    
    % Load random event
    event_idx = randi(length(event_files));
    event_filename = fullfile(preprocessPath, event_files(event_idx).name);
    
    fprintf('\nLoading Event %d/%d: %s\n', event_idx, length(event_files), event_files(event_idx).name);
    load(event_filename, 'event_data', 'event_metadata');
    
    fprintf('  Date: %s\n', datestr(event_metadata.calendar_date, 'dd-mmm-yyyy HH:MM:SS'));
    fprintf('  Sources: %d | Particles: %d | Timesteps: %d\n', ...
        event_metadata.n_sources, event_metadata.n_particles_total, event_metadata.n_timesteps);
    
    % Flatten for connectivity analysis
    fprintf('\nFlattening for connectivity analysis...\n');
    total_positions = event_metadata.n_particles_total * event_metadata.n_timesteps;
    all_lon = zeros(total_positions, 1);
    all_lat = zeros(total_positions, 1);
    all_source_idx = zeros(total_positions, 1);
    
    pos_idx = 1;
    for s = 1:length(event_data)
        n_pos = numel(event_data(s).lon);
        all_lon(pos_idx:pos_idx+n_pos-1) = event_data(s).lon(:);
        all_lat(pos_idx:pos_idx+n_pos-1) = event_data(s).lat(:);
        all_source_idx(pos_idx:pos_idx+n_pos-1) = event_data(s).source_reef_idx;
        pos_idx = pos_idx + n_pos;
    end
    
    fprintf('  Flattened %d positions from %d sources\n', total_positions, length(event_data));
end