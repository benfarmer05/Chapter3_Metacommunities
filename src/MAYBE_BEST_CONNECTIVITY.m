%% ========== STREAMING CONNECTIVITY ANALYSIS (FULLY INDEPENDENT) ==========
% Processes trajectory NetCDF files directly without any preprocessing
% Memory-efficient: never loads full events into memory
clear; clc;

%% Configuration
YEAR = 2019;
QUARTER = 1;
EVENT_TO_ANALYZE = 'all';  % Which event: number (1, 2, 3...) or 'all'

fprintf('=== STREAMING CONNECTIVITY ANALYSIS ===\n');

%% Setup paths
projectPath = matlab.project.rootProject().RootFolder;
dataPath = fullfile(projectPath, 'data');
quarter_name = sprintf('Q%d_%d', QUARTER, YEAR);
tempPath = fullfile('D:\Dissertation\CMS_traj', quarter_name);
outputPath = fullfile('D:\Dissertation\CMS_traj\output', 'CMS_traj', quarter_name);

if ~exist(outputPath, 'dir'), mkdir(outputPath); end

%% Load reef geometry data
fprintf('Loading reef polygons...\n');
centroids = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
unique_IDs = centroids(:,1);
Xs = [centroids(:,8) centroids(:,10) centroids(:,12) centroids(:,14) centroids(:,8)];
Ys = [centroids(:,9) centroids(:,11) centroids(:,13) centroids(:,15) centroids(:,9)];
n_locations = size(centroids,1);
fprintf('Loaded %d reef polygons.\n', n_locations);

% Pre-compute bounding boxes for optimization
bbox_xmin = min(Xs(:,1:4), [], 2);
bbox_xmax = max(Xs(:,1:4), [], 2);
bbox_ymin = min(Ys(:,1:4), [], 2);
bbox_ymax = max(Ys(:,1:4), [], 2);

%% Load release file mapping
releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019_Q1.txt'));
release_reef_IDs = releasefile(:,1);

%% STEP 1: Scan trajectory files (with caching)
fprintf('\n=== IDENTIFYING RELEASE EVENTS ===\n');

% Check for cached scan results
scan_cache_file = fullfile(outputPath, 'event_scan_cache.mat');

if exist(scan_cache_file, 'file')
    fprintf('Loading cached event scan results...\n');
    load(scan_cache_file, 'file_release_dates', 'unique_release_dates', 'trajlist', 'scan_date');
    n_files = length(trajlist);
    n_events = length(unique_release_dates);
    fprintf('Loaded from cache (scanned on %s)\n', datestr(scan_date));
    fprintf('  %d trajectory files\n', n_files);
    fprintf('  %d unique events\n', n_events);
else
    fprintf('No cache found. Scanning trajectory files...\n');
    
    trajlist = dir(fullfile(tempPath, 'traj*.nc'));
    n_files = length(trajlist);
    fprintf('Found %d trajectory files.\n', n_files);
    
    if n_files == 0
        error('No trajectory files found in %s', tempPath);
    end
    
    fprintf('Scanning for unique release events...\n');
    file_release_dates = zeros(n_files, 1);
    tic;
    
    for i = 1:n_files
        if mod(i, 50) == 0
            fprintf('  Scanned %d/%d files (%.1f%%) - %.1f files/sec\n', ...
                i, n_files, 100*i/n_files, i/toc);
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
    
    fprintf('Scan complete in %.1f seconds\n', toc);
    
    % Remove failed reads and find unique events
    file_release_dates = file_release_dates(~isnan(file_release_dates));
    unique_release_dates = unique(file_release_dates);
    n_events = length(unique_release_dates);
    
    fprintf('Found %d unique release events.\n', n_events);
    
    % Save cache for future runs
    scan_date = datetime('now');
    save(scan_cache_file, 'file_release_dates', 'unique_release_dates', 'trajlist', 'scan_date', '-v7.3');
    fprintf('Cached results saved to: %s\n', scan_cache_file);
end

drawnow;

%% STEP 2: Determine which events to process
if ischar(EVENT_TO_ANALYZE) && strcmpi(EVENT_TO_ANALYZE, 'all')
    events_to_process = 1:n_events;
    fprintf('Will process all %d events\n', n_events);
else
    events_to_process = EVENT_TO_ANALYZE;
    if max(events_to_process) > n_events
        error('Requested event %d but only %d events available', max(events_to_process), n_events);
    end
    fprintf('Will process %d event(s)\n', length(events_to_process));
end

%% Check for INPOLY2
use_inpoly = false;
if exist('inpoly2', 'file') == 2
    try
        test_result = inpoly2([0 0], [0 0; 1 0; 1 1; 0 1]);
        use_inpoly = true;
        fprintf('Using INPOLY2 fast algorithm.\n');
    catch
        fprintf('Using standard inpolygon.\n');
    end
else
    fprintf('Using standard inpolygon.\n');
end

%% Check parallel pool
pool = gcp('nocreate');
if isempty(pool)
    fprintf('Starting parallel pool...\n');
    pool = parpool('local');
else
    fprintf('Using existing parallel pool with %d workers.\n', pool.NumWorkers);
end

%% STEP 3: Process each event
for event_idx = events_to_process
    fprintf('\n========================================\n');
    fprintf('PROCESSING EVENT %d/%d\n', event_idx, n_events);
    fprintf('========================================\n');
    fprintf('Press Ctrl+C to cancel\n');
    drawnow;
    
    target_date = unique_release_dates(event_idx);
    calendar_date = datetime(target_date, 'ConvertFrom', 'juliandate');
    fprintf('Date: %s\n', datestr(calendar_date, 'dd-mmm-yyyy HH:MM:SS'));
    
    % Find trajectory files for THIS event
    relevant_file_indices = find(abs(file_release_dates - target_date) < 0.01);
    n_relevant_files = length(relevant_file_indices);
    fprintf('This event spans %d trajectory files\n', n_relevant_files);
    
    if n_relevant_files == 0
        warning('No files found for this event. Skipping.');
        continue;
    end
    
    % Initialize connectivity matrix
    ConnMatrix = zeros(n_locations, n_locations);
    
    %% Process each relevant file (STREAMING)
    fprintf('\nProcessing trajectory files...\n');
    fprintf('Note: Once parallel processing starts, cannot be interrupted.\n');
    drawnow;
    
    overall_timer = tic;
    
    for file_counter = 1:n_relevant_files
        if mod(file_counter, 5) == 0 || file_counter == n_relevant_files
            fprintf('  File %d/%d (%.1f%%) - %.1f sec elapsed\n', ...
                file_counter, n_relevant_files, 100*file_counter/n_relevant_files, toc(overall_timer));
            drawnow;
        end
        
        file_idx = relevant_file_indices(file_counter);
        filename = fullfile(tempPath, trajlist(file_idx).name);
        
        try
            % Read only necessary data from NetCDF
            lon = ncread(filename, 'lon');
            lat = ncread(filename, 'lat');
            location = ncread(filename, 'location');
            releasedate = ncread(filename, 'releasedate');
            
            % Filter to this specific event (double-check)
            event_mask = abs(releasedate - target_date) < 0.01;
            if ~any(event_mask), continue; end
            
            n_particles = size(lon, 2);
            if length(event_mask) > n_particles
                event_mask = event_mask(1:n_particles);
            end
            
            % Extract event particles
            lon(lon > 180) = lon(lon > 180) - 360;
            lon_event = lon(:, event_mask);
            lat_event = lat(:, event_mask);
            loc_event = location(event_mask);
            
            % Map particles to source reef indices
            reef_IDs = release_reef_IDs(loc_event);
            [~, src_indices] = ismember(reef_IDs, unique_IDs);
            
            % Flatten trajectories for this file
            [n_timesteps, n_part] = size(lon_event);
            lon_flat = lon_event(:);
            lat_flat = lat_event(:);
            src_flat = repelem(src_indices, n_timesteps);
            
            % Remove invalid positions
            valid = ~isnan(lon_flat) & ~isnan(lat_flat) & src_flat > 0;
            lon_flat = lon_flat(valid);
            lat_flat = lat_flat(valid);
            src_flat = src_flat(valid);
            
            if isempty(lon_flat)
                clear lon lat location releasedate lon_event lat_event loc_event
                clear lon_flat lat_flat src_flat
                continue;
            end
            
            % Process in chunks with parallel computing
            chunk_size = 50000;
            n_chunks = ceil(length(lon_flat) / chunk_size);
            
            % Copy variables for parfor
            Xs_par = Xs;
            Ys_par = Ys;
            bbox_xmin_par = bbox_xmin;
            bbox_xmax_par = bbox_xmax;
            bbox_ymin_par = bbox_ymin;
            bbox_ymax_par = bbox_ymax;
            n_locations_par = n_locations;
            use_inpoly_par = use_inpoly;
            
            chunk_results = cell(n_chunks, 1);
            
            parfor chunk = 1:n_chunks
                idx_start = (chunk-1)*chunk_size + 1;
                idx_end = min(chunk*chunk_size, length(lon_flat));
                
                lon_chunk = lon_flat(idx_start:idx_end);
                lat_chunk = lat_flat(idx_start:idx_end);
                src_chunk = src_flat(idx_start:idx_end);
                
                dest_chunk = zeros(length(lon_chunk), 1);
                
                % Test each polygon with bounding box pre-filtering
                for j = 1:n_locations_par
                    in_bbox = lon_chunk >= bbox_xmin_par(j) & lon_chunk <= bbox_xmax_par(j) & ...
                              lat_chunk >= bbox_ymin_par(j) & lat_chunk <= bbox_ymax_par(j);
                    
                    if ~any(in_bbox), continue; end
                    
                    candidates = find(in_bbox);
                    
                    if use_inpoly_par
                        in_poly = inpoly2([lon_chunk(candidates), lat_chunk(candidates)], ...
                                         [Xs_par(j,1:4)', Ys_par(j,1:4)']);
                    else
                        in_poly = inpolygon(lon_chunk(candidates), lat_chunk(candidates), ...
                                           Xs_par(j,:), Ys_par(j,:));
                    end
                    
                    dest_chunk(candidates(in_poly)) = j;
                end
                
                % Store connections
                valid_conn = dest_chunk > 0;
                chunk_results{chunk} = [src_chunk(valid_conn), dest_chunk(valid_conn)];
            end
            
            % Accumulate connections into connectivity matrix
            for chunk = 1:n_chunks
                connections = chunk_results{chunk};
                for i = 1:size(connections, 1)
                    ConnMatrix(connections(i,1), connections(i,2)) = ...
                        ConnMatrix(connections(i,1), connections(i,2)) + 1;
                end
            end
            
            % Clear file data before next iteration
            clear lon lat location releasedate lon_event lat_event loc_event
            clear lon_flat lat_flat src_flat chunk_results connections
            
        catch ME
            warning('Error processing %s: %s', trajlist(file_idx).name, ME.message);
        end
    end
    
    elapsed = toc(overall_timer);
    fprintf('Connectivity computed in %.1f seconds (%.1f minutes)\n', elapsed, elapsed/60);
    drawnow;
    
    %% Normalize connectivity matrix
    ConnMatrix_raw = ConnMatrix;
    row_sums = sum(ConnMatrix, 2);
    row_sums(row_sums == 0) = 1;
    ConnMatrix_normalized = ConnMatrix ./ row_sums;
    
    %% Results and statistics
    fprintf('\n=== CONNECTIVITY RESULTS ===\n');
    fprintf('Total connections: %d\n', nnz(ConnMatrix_raw));
    
    [src_indices, dst_indices] = find(ConnMatrix_raw);
    active_sources = unique(src_indices);
    active_dests = unique(dst_indices);
    
    fprintf('Active sources: %d\n', length(active_sources));
    fprintf('Active destinations: %d\n', length(active_dests));
    fprintf('Self-recruitment events: %d\n', sum(diag(ConnMatrix_raw)));
    
    % Network metrics
    out_degree = sum(ConnMatrix_raw, 2);
    in_degree = sum(ConnMatrix_raw, 1)';
    total_connectivity = out_degree + in_degree;
    
    fprintf('\nTop 5 most connected reefs:\n');
    [~, sorted_idx] = sort(total_connectivity, 'descend');
    for i = 1:min(5, nnz(total_connectivity))
        idx = sorted_idx(i);
        fprintf('  Reef %d (ID=%d): Out=%d, In=%d, Total=%d\n', ...
            idx, unique_IDs(idx), out_degree(idx), in_degree(idx), total_connectivity(idx));
    end
    
    %% Save results
    results_file = fullfile(outputPath, sprintf('connectivity_%s.mat', ...
                            datestr(calendar_date, 'yyyy_mmmdd_HHMMSS')));
    
    connectivity_results = struct();
    connectivity_results.ConnMatrix_raw = ConnMatrix_raw;
    connectivity_results.ConnMatrix_normalized = ConnMatrix_normalized;
    connectivity_results.calendar_date = calendar_date;
    connectivity_results.target_date = target_date;
    connectivity_results.processing_time_sec = elapsed;
    connectivity_results.n_files_processed = n_relevant_files;
    connectivity_results.out_degree = out_degree;
    connectivity_results.in_degree = in_degree;
    connectivity_results.total_connectivity = total_connectivity;
    connectivity_results.method = 'streaming_independent';
    connectivity_results.year = YEAR;
    connectivity_results.quarter = QUARTER;
    
    save(results_file, 'connectivity_results', '-v7.3');
    fprintf('\nResults saved: %s\n', results_file);
    drawnow;
    
    %% Generate visualizations
    fprintf('Generating visualizations...\n');
    drawnow;
    
    active_locs = unique([active_sources; active_dests]);
    centroids_x = mean(Xs(:,1:4), 2);
    centroids_y = mean(Ys(:,1:4), 2);
    
    % Normalized connectivity heatmap
    fig1 = figure('Position', [100 100 800 600], 'Visible', 'off');
    sub_conn = ConnMatrix_normalized(active_locs, active_locs);
    imagesc(sub_conn);
    colorbar;
    caxis([0 1]);
    xlabel('Destination Reef');
    ylabel('Source Reef');
    title(sprintf('Normalized Connectivity - %s', datestr(calendar_date, 'dd-mmm-yyyy')));
    axis square;
    saveas(fig1, fullfile(outputPath, sprintf('conn_norm_%s.png', ...
                          datestr(calendar_date, 'yyyy_mmmdd_HHMMSS'))));
    close(fig1);
    
    % Raw counts heatmap
    fig2 = figure('Position', [100 100 800 600], 'Visible', 'off');
    sub_conn_raw = ConnMatrix_raw(active_locs, active_locs);
    imagesc(log10(sub_conn_raw + 1));
    colorbar;
    xlabel('Destination Reef');
    ylabel('Source Reef');
    title(sprintf('Raw Connectivity (log10) - %s', datestr(calendar_date, 'dd-mmm-yyyy')));
    axis square;
    saveas(fig2, fullfile(outputPath, sprintf('conn_raw_%s.png', ...
                          datestr(calendar_date, 'yyyy_mmmdd_HHMMSS'))));
    close(fig2);
    
    % Spatial connectivity map
    fig3 = figure('Position', [100 100 1000 800], 'Visible', 'off');
    hold on;
    
    for i = 1:length(active_locs)
        idx = active_locs(i);
        plot(Xs(idx,:), Ys(idx,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    end
    
    conn_size = total_connectivity(active_locs);
    conn_size = conn_size / max(conn_size) * 200 + 20;
    
    scatter(centroids_x(active_locs), centroids_y(active_locs), ...
        conn_size, total_connectivity(active_locs), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    
    colorbar;
    xlabel('Longitude');
    ylabel('Latitude');
    title(sprintf('Reef Connectivity - %s', datestr(calendar_date, 'dd-mmm-yyyy')));
    axis equal tight;
    grid on;
    saveas(fig3, fullfile(outputPath, sprintf('conn_spatial_%s.png', ...
                          datestr(calendar_date, 'yyyy_mmmdd_HHMMSS'))));
    close(fig3);
    
    fprintf('Figures saved.\n');
    drawnow;
    
    %% Clear memory before next event
    fprintf('Clearing memory...\n');
    clear ConnMatrix ConnMatrix_raw ConnMatrix_normalized connectivity_results
    clear sub_conn sub_conn_raw conn_size active_locs
    java.lang.System.gc();
    fprintf('Memory cleared.\n');
end

fprintf('\n=== ANALYSIS COMPLETE ===\n');
fprintf('Processed %d event(s) using streaming approach\n', length(events_to_process));
fprintf('Total trajectory files in quarter: %d\n', n_files);
fprintf('Unique events found: %d\n', n_events);




















% %% investigate
% 
% 
% 
% %% VISUAL VERIFICATION: Sample Trajectories Over Reef Grid
% % Simple trajectory plotting to verify connectivity makes sense
% clear; clc;
% 
% %% Configuration
% YEAR = 2019;
% QUARTER = 1;
% EVENT_TO_CHECK = 1;
% N_TRAJECTORIES = 100;  % How many trajectories to plot total
% 
% %% Setup paths
% projectPath = matlab.project.rootProject().RootFolder;
% dataPath = fullfile(projectPath, 'data');
% quarter_name = sprintf('Q%d_%d', QUARTER, YEAR);
% tempPath = fullfile('D:\Dissertation\CMS_traj', quarter_name);
% outputPath = fullfile('D:\Dissertation\CMS_traj\output', 'CMS_traj', quarter_name);
% 
% fprintf('=== TRAJECTORY VISUALIZATION ===\n');
% 
% %% Load reef polygons
% fprintf('Loading reef polygons...\n');
% centroids = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
% Xs = [centroids(:,8) centroids(:,10) centroids(:,12) centroids(:,14) centroids(:,8)];
% Ys = [centroids(:,9) centroids(:,11) centroids(:,13) centroids(:,15) centroids(:,9)];
% 
% %% Load event scan cache
% scan_cache = fullfile(outputPath, 'event_scan_cache.mat');
% if ~exist(scan_cache, 'file')
%     error('Event scan cache not found. Run streaming connectivity script first.');
% end
% 
% load(scan_cache, 'file_release_dates', 'unique_release_dates', 'trajlist');
% 
% % Get target date for this event
% target_date = unique_release_dates(EVENT_TO_CHECK);
% calendar_date = datetime(target_date, 'ConvertFrom', 'juliandate');
% fprintf('Event: %s\n', datestr(calendar_date, 'dd-mmm-yyyy HH:MM:SS'));
% 
% % Find relevant files
% relevant_files = find(abs(file_release_dates - target_date) < 0.01);
% fprintf('Event spans %d trajectory files\n', length(relevant_files));
% 
% %% Sample trajectories
% fprintf('Sampling %d trajectories...\n', N_TRAJECTORIES);
% 
% % Pre-allocate
% selected_data = struct('lon', [], 'lat', []);
% trajectories_collected = 0;
% 
% % Random seed for reproducibility
% rng('shuffle');
% seed = randi(10000);
% rng(seed);
% fprintf('Random seed: %d\n', seed);
% 
% % Sample from relevant files
% n_files_to_sample = min(5, length(relevant_files));
% sampled_file_indices = randsample(length(relevant_files), n_files_to_sample, false);
% 
% for i = sampled_file_indices'
%     if trajectories_collected >= N_TRAJECTORIES, break; end
% 
%     file_idx = relevant_files(i);
%     filename = fullfile(tempPath, trajlist(file_idx).name);
% 
%     try
%         location = ncread(filename, 'location');
%         lon = ncread(filename, 'lon');
%         lat = ncread(filename, 'lat');
%         releasedate = ncread(filename, 'releasedate');
% 
%         % Filter to this event
%         event_mask = abs(releasedate - target_date) < 0.01;
%         if ~any(event_mask), continue; end
% 
%         % Adjust mask size if needed
%         n_particles = size(lon, 2);
%         if length(event_mask) > n_particles
%             event_mask = event_mask(1:n_particles);
%         end
% 
%         % Convert longitude if needed
%         lon(lon > 180) = lon(lon > 180) - 360;
% 
%         % Get particles from this event
%         lon_event = lon(:, event_mask);
%         lat_event = lat(:, event_mask);
%         location_event = location(event_mask);
% 
%         % Sample unique locations
%         unique_locs = unique(location_event);
%         n_take = min(length(unique_locs), N_TRAJECTORIES - trajectories_collected);
%         sampled_locs = randsample(unique_locs, n_take, false);
% 
%         for loc = sampled_locs'
%             if trajectories_collected >= N_TRAJECTORIES, break; end
% 
%             loc_idx = find(location_event == loc, 1);
%             trajectories_collected = trajectories_collected + 1;
%             selected_data(trajectories_collected).lon = lon_event(:, loc_idx);
%             selected_data(trajectories_collected).lat = lat_event(:, loc_idx);
%         end
%     catch ME
%         warning('Error reading file: %s', ME.message);
%         continue;
%     end
% end
% 
% selected_data = selected_data(1:trajectories_collected);
% fprintf('Collected %d trajectories\n', trajectories_collected);
% 
% %% Create figure
% figure('Position', [100 100 1400 1000]);
% hold on;
% 
% % Plot reef polygons
% fprintf('Drawing reef polygons...\n');
% for i = 1:size(Xs, 1)
%     plot(Xs(i,:), Ys(i,:), 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
% end
% 
% % Plot trajectories
% fprintf('Plotting trajectories...\n');
% colors = lines(10);
% start_lons = zeros(trajectories_collected, 1);
% start_lats = zeros(trajectories_collected, 1);
% 
% for i = 1:trajectories_collected
%     % Plot trajectory
%     plot(selected_data(i).lon, selected_data(i).lat, ...
%         'Color', colors(mod(i-1,10)+1,:), 'LineWidth', 0.5);
% 
%     % Store starting position
%     valid_idx = find(~isnan(selected_data(i).lon) & ~isnan(selected_data(i).lat), 1);
%     if ~isempty(valid_idx)
%         start_lons(i) = selected_data(i).lon(valid_idx);
%         start_lats(i) = selected_data(i).lat(valid_idx);
%     end
% end
% 
% % Plot starting locations as red stars
% valid_starts = ~isnan(start_lons);
% plot(start_lons(valid_starts), start_lats(valid_starts), 'r*', ...
%     'MarkerSize', 8, 'LineWidth', 1.5);
% 
% % Finalize plot
% xlabel('Longitude', 'FontSize', 14, 'FontWeight', 'bold');
% ylabel('Latitude', 'FontSize', 14, 'FontWeight', 'bold');
% title(sprintf('Particle Trajectories - %s\n%d trajectories (seed: %d)', ...
%     datestr(calendar_date, 'dd-mmm-yyyy'), trajectories_collected, seed), ...
%     'FontSize', 16, 'FontWeight', 'bold');
% axis equal tight;
% grid on;
% box on;
% 
% fprintf('\n=== VISUALIZATION COMPLETE ===\n');
% fprintf('Red stars = starting locations\n');
% fprintf('Colored lines = particle trajectories\n');
% 
% %% Load and visualize connectivity matrices
% fprintf('\n=== LOADING CONNECTIVITY RESULTS ===\n');
% 
% % Find connectivity file for this event
% conn_files = dir(fullfile(outputPath, 'connectivity_*.mat'));
% if isempty(conn_files)
%     warning('No connectivity files found. Skipping matrix visualization.');
%     return;
% end
% 
% % Sort by date and get matching event
% conn_dates = zeros(length(conn_files), 1);
% for i = 1:length(conn_files)
%     temp = load(fullfile(outputPath, conn_files(i).name), 'connectivity_results');
%     conn_dates(i) = temp.connectivity_results.target_date;
% end
% [~, sort_idx] = sort(conn_dates);
% conn_files = conn_files(sort_idx);
% 
% % Load the connectivity for this event
% conn_file = fullfile(outputPath, conn_files(EVENT_TO_CHECK).name);
% fprintf('Loading: %s\n', conn_files(EVENT_TO_CHECK).name);
% load(conn_file, 'connectivity_results');
% 
% ConnMatrix_raw = connectivity_results.ConnMatrix_raw;
% ConnMatrix_normalized = connectivity_results.ConnMatrix_normalized;
% 
% fprintf('Matrix size: %d x %d\n', size(ConnMatrix_raw, 1), size(ConnMatrix_raw, 2));
% fprintf('Non-zero connections: %d\n', nnz(ConnMatrix_raw));
% 
% %% Calculate alternative normalization: by total particle-timesteps
% % This shows: probability that any single particle-timestep from source i lands in dest j
% % Denominator = n_particles * n_timesteps for each source
% 
% % We need to know how many particle-timesteps each source had
% % This requires going back to the event data
% fprintf('Calculating particle-timestep normalization...\n');
% 
% % Check if we have n_files_processed in results
% if isfield(connectivity_results, 'n_files_processed')
%     fprintf('Note: Using row-sum approximation for particle-timesteps\n');
%     % Approximate: assume each source had roughly the same number of particle-timesteps
%     % Total particle-timesteps per source ≈ total connections / avg connectivity
%     total_positions_per_source = sum(ConnMatrix_raw, 2);  % Total hits per source
% 
%     % For sources with no hits, use a small value to avoid division by zero
%     total_positions_per_source(total_positions_per_source == 0) = 1;
% 
%     % Create normalization by dividing by total positions
%     % But this is actually the same as row normalization...
%     % We need the ACTUAL total particle-timesteps, not just hits
% 
%     % Better approach: assume uniform particle release
%     % Total particle-timesteps = connections / connection_rate
%     % For now, use a fixed estimate based on typical simulation
%     assumed_timesteps = 1345;  % typical from your data
%     assumed_particles_per_source = 100;  % typical
%     total_particle_timesteps = assumed_timesteps * assumed_particles_per_source;
% 
%     ConnMatrix_by_particle_time = ConnMatrix_raw / total_particle_timesteps;
% 
%     fprintf('  Using assumed values: %d timesteps × %d particles = %d particle-timesteps per source\n', ...
%         assumed_timesteps, assumed_particles_per_source, total_particle_timesteps);
% else
%     % If we don't have the info, use row normalization as placeholder
%     ConnMatrix_by_particle_time = ConnMatrix_normalized;
% end
% 
% %% Plot connectivity matrices
% figure('Position', [100 100 1600 700]);
% 
% % Raw connectivity (log scale)
% subplot(1, 2, 1);
% imagesc(log10(ConnMatrix_raw + 1));
% colorbar;
% xlabel('Destination Reef', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Source Reef', 'FontSize', 12, 'FontWeight', 'bold');
% title(sprintf('Raw Connectivity (log10)\n%s', ...
%     datestr(calendar_date, 'dd-mmm-yyyy')), ...
%     'FontSize', 14, 'FontWeight', 'bold');
% axis square;
% 
% % Normalized connectivity
% subplot(1, 2, 2);
% imagesc(ConnMatrix_normalized);
% colorbar;
% caxis([0 1]);
% xlabel('Destination Reef', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Source Reef', 'FontSize', 12, 'FontWeight', 'bold');
% title(sprintf('Normalized Connectivity\n%s', ...
%     datestr(calendar_date, 'dd-mmm-yyyy')), ...
%     'FontSize', 14, 'FontWeight', 'bold');
% axis square;
% 
% sgtitle('Connectivity Matrices', 'FontSize', 16, 'FontWeight', 'bold');
% 
% %% Plot histograms of random row
% % Pick a random source reef that has connections
% sources_with_connections = find(sum(ConnMatrix_raw, 2) > 0);
% if isempty(sources_with_connections)
%     warning('No sources with connections found.');
%     return;
% end
% 
% random_source = sources_with_connections(randi(length(sources_with_connections)));
% fprintf('\nAnalyzing random source reef: %d\n', random_source);
% 
% % Get this source's connections
% raw_row = ConnMatrix_raw(random_source, :);
% norm_row = ConnMatrix_normalized(random_source, :);
% 
% % Only look at non-zero values for histograms
% raw_nonzero = raw_row(raw_row > 0);
% norm_nonzero = norm_row(norm_row > 0);
% 
% fprintf('  Raw connections: %d destinations\n', length(raw_nonzero));
% fprintf('  Raw range: [%d, %d]\n', min(raw_nonzero), max(raw_nonzero));
% fprintf('  Normalized range: [%.4f, %.4f]\n', min(norm_nonzero), max(norm_nonzero));
% 
% figure('Position', [100 100 1600 700]);
% 
% % Histogram of raw values
% subplot(1, 2, 1);
% histogram(raw_nonzero, 50, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'k');
% xlabel('Connection Count', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Frequency', 'FontSize', 12, 'FontWeight', 'bold');
% title(sprintf('Raw Connectivity Distribution\nSource Reef %d → All Destinations', random_source), ...
%     'FontSize', 14, 'FontWeight', 'bold');
% grid on;
% 
% % Histogram of normalized values
% subplot(1, 2, 2);
% histogram(norm_nonzero, 50, 'FaceColor', [0.9 0.5 0.3], 'EdgeColor', 'k');
% xlabel('Connection Probability', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Frequency', 'FontSize', 12, 'FontWeight', 'bold');
% title(sprintf('Normalized Connectivity Distribution\nSource Reef %d → All Destinations', random_source), ...
%     'FontSize', 14, 'FontWeight', 'bold');
% grid on;
% 
% sgtitle(sprintf('Connection Value Distributions (Non-zero only, n=%d destinations)', ...
%     length(raw_nonzero)), 'FontSize', 16, 'FontWeight', 'bold');
% 
% fprintf('\n=== ALL VISUALIZATIONS COMPLETE ===\n');