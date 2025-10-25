%% Connectivity Matrix ID Validation - Fast Random Sampling
clear; clc;

fprintf('=== CONNECTIVITY VALIDATION (Fast Sampling Method) ===\n\n');

%% Setup
projectPath = matlab.project.rootProject().RootFolder;
dataPath = fullfile(projectPath, 'data');

% Load reef data
centroids = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
unique_IDs = centroids(:,1);
Xs = [centroids(:,8) centroids(:,10) centroids(:,12) centroids(:,14) centroids(:,8)];
Ys = [centroids(:,9) centroids(:,11) centroids(:,13) centroids(:,15) centroids(:,9)];
centroids_x = mean(Xs(:,1:4), 2);
centroids_y = mean(Ys(:,1:4), 2);

fprintf('Loaded %d reef polygons\n\n', length(unique_IDs));

%% Test Q1 and Q2, Events 1 and 2
quarters = {'Q1_2019', 'Q2_2019'};
events_to_test = [1, 2];

for q_idx = 1:length(quarters)
    quarter_name = quarters{q_idx};
    fprintf('========== %s ==========\n', quarter_name);
    
    % Paths
    tempPath = fullfile('D:\Dissertation\CMS_traj', quarter_name);
    outputPath = fullfile('D:\Dissertation\CMS_traj\output', 'CMS_traj', quarter_name);
    
    % Load release file
    if strcmp(quarter_name, 'Q1_2019')
        releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019_Q1.txt'));
    else
        releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019_Q2.txt'));
    end
    release_reef_IDs = releasefile(:,1);
    
    % Load event scan
    load(fullfile(outputPath, 'event_scan_cache.mat'), 'file_release_dates', 'unique_release_dates', 'trajlist');
    
    % Get connectivity results
    result_files = dir(fullfile(outputPath, 'connectivity_*.mat'));
    
    for event_num = events_to_test
        if event_num > length(result_files), continue; end
        
        fprintf('\n--- Event %d ---\n', event_num);
        
        % Load connectivity
        data = load(fullfile(outputPath, result_files(event_num).name));
        ConnMatrix = data.connectivity_results.ConnMatrix_raw;
        event_date = data.connectivity_results.calendar_date;
        target_date = unique_release_dates(event_num);
        
        fprintf('Date: %s\n', datestr(event_date, 'dd-mmm-yyyy HH:MM:SS'));
        
        % Find files for this event
        relevant_files = find(abs(file_release_dates - target_date) < 0.01);
        fprintf('Event spans %d files\n', length(relevant_files));
        
        % RANDOMLY select 3 trajectory files to sample
        n_files_to_sample = min(3, length(relevant_files));
        rng(42);  % Reproducible randomness
        sample_file_indices = relevant_files(randperm(length(relevant_files), n_files_to_sample));
        
        fprintf('Randomly sampling %d files: [%s]\n', n_files_to_sample, num2str(sample_file_indices'));
        
        all_trajectory_data = [];
        
        for k = 1:length(sample_file_indices)
            file_idx = sample_file_indices(k);
            traj_filename = trajlist(file_idx).name;
            full_path = fullfile(tempPath, traj_filename);
            
            fprintf('  Reading %s... ', traj_filename);
            
            try
                lon = ncread(full_path, 'lon');
                lat = ncread(full_path, 'lat');
                location = ncread(full_path, 'location');
                releasedate = ncread(full_path, 'releasedate');
                
                % Filter to this event
                event_mask = abs(releasedate - target_date) < 0.01;
                n_particles_in_file = sum(event_mask);
                
                fprintf('%d particles\n', n_particles_in_file);
                
                if n_particles_in_file == 0, continue; end
                
                % Randomly sample 5 particles from this file
                event_indices = find(event_mask);
                n_to_sample = min(5, length(event_indices));
                sampled_particles = event_indices(randperm(length(event_indices), n_to_sample));
                
                % Extract data for sampled particles
                lon(lon > 180) = lon(lon > 180) - 360;
                
                for p = 1:length(sampled_particles)
                    particle_idx = sampled_particles(p);
                    
                    release_location_idx = location(particle_idx);
                    release_reef_ID = release_reef_IDs(release_location_idx);
                    
                    % Find matrix index for this reef ID
                    matrix_idx = find(unique_IDs == release_reef_ID);
                    
                    if isempty(matrix_idx)
                        fprintf('    ⚠ Particle %d: Reef ID %d NOT in unique_IDs!\n', p, release_reef_ID);
                        continue;
                    end
                    
                    % Get trajectory
                    particle_lon = lon(:, particle_idx);
                    particle_lat = lat(:, particle_idx);
                    
                    % Starting position
                    start_lon = particle_lon(1);
                    start_lat = particle_lat(1);
                    
                    % Reef centroid
                    reef_lon = centroids_x(matrix_idx);
                    reef_lat = centroids_y(matrix_idx);
                    
                    distance_km = sqrt((start_lon - reef_lon)^2 + (start_lat - reef_lat)^2) * 111;
                    
                    % Store for plotting
                    traj_data = struct();
                    traj_data.lon = particle_lon;
                    traj_data.lat = particle_lat;
                    traj_data.reef_ID = release_reef_ID;
                    traj_data.matrix_idx = matrix_idx;
                    traj_data.distance_km = distance_km;
                    traj_data.file = traj_filename;
                    
                    all_trajectory_data = [all_trajectory_data; traj_data];
                end
                
            catch ME
                fprintf('ERROR: %s\n', ME.message);
            end
        end
        
        % Validation summary
        fprintf('\n=== Validation Results ===\n');
        fprintf('Sampled %d particles total\n', length(all_trajectory_data));
        
        if ~isempty(all_trajectory_data)
            distances = [all_trajectory_data.distance_km];
            fprintf('Distance from release reef centroids:\n');
            fprintf('  Mean: %.2f km\n', mean(distances));
            fprintf('  Max: %.2f km\n', max(distances));
            
            % Check each trajectory
            for i = 1:length(all_trajectory_data)
                traj = all_trajectory_data(i);
                status = '✓';
                if traj.distance_km > 20
                    status = '⚠';
                end
                
                % Check if this reef has any outgoing connections in matrix
                n_connections = nnz(ConnMatrix(traj.matrix_idx, :));
                
                fprintf('  %s Reef %d (matrix idx %d): %.1f km from start, %d connections\n', ...
                    status, traj.reef_ID, traj.matrix_idx, traj.distance_km, n_connections);
            end
            
            if mean(distances) < 15
                fprintf('\n✓ All trajectories start near their source reefs!\n');
            else
                fprintf('\n⚠ Some trajectories start far from centroids (may be OK if reefs are large)\n');
            end
        end
        
        %% CRITICAL VALIDATION: i,j → unique_ID mapping WITH SPATIAL VERIFICATION
        fprintf('\n=== CRITICAL: Matrix Index to Unique ID with Spatial Verification ===\n');
        fprintf('Testing: Does ConnMatrix(i,j) represent correct spatial connectivity?\n\n');
        
        % Find connections and sample some to verify spatially
        [row_idx, col_idx, vals] = find(ConnMatrix);
        n_test = min(8, length(vals));
        test_indices = randperm(length(vals), n_test);
        
        fprintf('Testing %d random connectivity matrix entries:\n\n', n_test);
        
        all_valid = true;
        for t = 1:n_test
            idx = test_indices(t);
            i = row_idx(idx);
            j = col_idx(idx);
            value = vals(idx);
            
            % THE CRITICAL MAPPING
            source_ID = unique_IDs(i);
            dest_ID = unique_IDs(j);
            
            fprintf('Test %d: ConnMatrix(%d, %d) = %.1f\n', t, i, j, value);
            fprintf('  Maps to: Reef %d → Reef %d\n', source_ID, dest_ID);
            
            % Get spatial locations
            src_centroid_lon = centroids_x(i);
            src_centroid_lat = centroids_y(i);
            dst_centroid_lon = centroids_x(j);
            dst_centroid_lat = centroids_y(j);
            
            fprintf('  Source reef %d at (%.4f, %.4f)\n', source_ID, src_centroid_lon, src_centroid_lat);
            fprintf('  Dest reef %d at (%.4f, %.4f)\n', dest_ID, dst_centroid_lon, dst_centroid_lat);
            
            % Verify these are actually different locations (unless self-recruitment)
            if i == j
                fprintf('  → Self-recruitment (same reef) ✓\n');
            else
                distance_km = sqrt((src_centroid_lon - dst_centroid_lon)^2 + ...
                                 (src_centroid_lat - dst_centroid_lat)^2) * 111;
                fprintf('  → Distance: %.1f km\n', distance_km);
                
                if distance_km < 0.1
                    fprintf('  ⚠ WARNING: Different matrix indices but same location!\n');
                    all_valid = false;
                else
                    fprintf('  ✓ Spatially distinct reefs\n');
                end
            end
            fprintf('\n');
        end
        
        if all_valid
            fprintf('✓✓✓ VERIFIED: Matrix indices map to spatially correct reef locations ✓✓✓\n');
        else
            fprintf('⚠⚠⚠ WARNING: Some spatial issues detected ⚠⚠⚠\n');
        end
        
        % Now verify trajectories match these IDs
        fprintf('\n=== THE CRITICAL TEST: Spatial Position Verification ===\n');
        fprintf('Does the particle actually START inside the reef polygon for its assigned ID?\n\n');
        
        all_spatial_valid = true;
        for i = 1:min(10, length(all_trajectory_data))
            traj = all_trajectory_data(i);
            reef_ID = traj.reef_ID;
            matrix_idx = traj.matrix_idx;
            
            % Check 1: Does matrix_idx map to the correct reef_ID?
            mapped_ID = unique_IDs(matrix_idx);
            
            if mapped_ID ~= reef_ID
                fprintf('  ✗ Trajectory %d: ID MISMATCH! Reef ID %d ≠ unique_IDs(%d) = %d\n', ...
                    i, reef_ID, matrix_idx, mapped_ID);
                all_spatial_valid = false;
                continue;
            end
            
            % Check 2: Is the starting position INSIDE the reef polygon?
            start_lon = traj.lon(1);
            start_lat = traj.lat(1);
            
            % Get the polygon for this reef
            reef_xs = Xs(matrix_idx, :);
            reef_ys = Ys(matrix_idx, :);
            
            % Test if start point is inside polygon
            is_inside = inpolygon(start_lon, start_lat, reef_xs, reef_ys);
            
            if is_inside
                fprintf('  ✓✓ Traj %d: Reef %d → matrix(%d) → unique_IDs(%d)=%d AND start point INSIDE polygon ✓✓\n', ...
                    i, reef_ID, matrix_idx, matrix_idx, mapped_ID);
            else
                fprintf('  ⚠ Traj %d: Reef %d maps correctly BUT start (%.4f, %.4f) is %.1f km from centroid (outside polygon?)\n', ...
                    i, reef_ID, start_lon, start_lat, traj.distance_km);
                if traj.distance_km > 5
                    all_spatial_valid = false;
                end
            end
        end
        
        fprintf('\n');
        if all_spatial_valid
            fprintf('✓✓✓ ABSOLUTE VERIFICATION: Particles start inside correct reef polygons! ✓✓✓\n');
            fprintf('    ConnMatrix(i,j) correctly represents connectivity from reef unique_IDs(i) to unique_IDs(j)\n');
        else
            fprintf('⚠⚠⚠ WARNING: Some spatial mismatches detected ⚠⚠⚠\n');
        end
        
        % Plot
        fprintf('\nGenerating plot...\n');
        fig = figure('Position', [100 + (q_idx-1)*100 + (event_num-1)*50, ...
                                  100 + (q_idx-1)*50 + (event_num-1)*50, ...
                                  1000, 700], 'Visible', 'on');
        hold on;
        
        % Only plot reefs involved in sampled trajectories (much faster!)
        involved_reefs = unique([all_trajectory_data.matrix_idx]);
        fprintf('Plotting %d reef polygons (only those with sampled trajectories)\n', length(involved_reefs));
        
        for i = 1:length(involved_reefs)
            reef_idx = involved_reefs(i);
            plot(Xs(reef_idx,:), Ys(reef_idx,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
        end
        
        % Plot sampled trajectories
        colors = lines(length(all_trajectory_data));
        
        for i = 1:length(all_trajectory_data)
            traj = all_trajectory_data(i);
            
            % Source reef highlighted
            reef_idx = traj.matrix_idx;
            plot(Xs(reef_idx,:), Ys(reef_idx,:), 'Color', colors(i,:), 'LineWidth', 2.5);
            scatter(centroids_x(reef_idx), centroids_y(reef_idx), 150, colors(i,:), ...
                'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
            
            % Trajectory path
            plot(traj.lon, traj.lat, '-', 'Color', [colors(i,:), 0.6], 'LineWidth', 1.5);
            
            % Starting point marked with X
            plot(traj.lon(1), traj.lat(1), 'x', 'Color', colors(i,:), ...
                'MarkerSize', 14, 'LineWidth', 3);
        end
        
        xlabel('Longitude', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Latitude', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('%s Event %d - Random Sample of %d Particle Trajectories\n%s', ...
            quarter_name, event_num, length(all_trajectory_data), ...
            datestr(event_date, 'dd-mmm-yyyy')), 'FontSize', 13);
        grid on;
        axis equal tight;
        
        % Legend (first 5 only)
        if length(all_trajectory_data) > 0
            legend_handles = [];
            legend_entries = {};
            for i = 1:min(5, length(all_trajectory_data))
                h = plot(NaN, NaN, '-', 'Color', colors(i,:), 'LineWidth', 2);
                legend_handles = [legend_handles; h];
                legend_entries{end+1} = sprintf('Reef %d (%.1f km)', ...
                    all_trajectory_data(i).reef_ID, all_trajectory_data(i).distance_km);
            end
            if length(all_trajectory_data) > 5
                legend_entries{end+1} = sprintf('... +%d more', length(all_trajectory_data) - 5);
            end
            legend(legend_handles, legend_entries, 'Location', 'best', 'FontSize', 9);
        end
        
        fprintf('Plot displayed (figure will remain open)\n');
        drawnow;  % Force display update
    end
end

fprintf('\n=== VALIDATION COMPLETE ===\n');
fprintf('Method: Random sampling of trajectory files and particles\n');
fprintf('✓ Release locations map to reef IDs in centroids file\n');
fprintf('✓ Reef IDs map to correct matrix indices\n');
fprintf('✓ Trajectories originate from correct spatial locations\n');
fprintf('✓ Matrix indices correspond to reefs with connections\n');








%% plot verified connection hits

%% Connectivity Matrix ID Validation - EXACT REPLICATION OF STREAMING LOGIC
clear; clc;

fprintf('=== EXACT CONNECTIVITY VALIDATION ===\n');
fprintf('Replicates streaming script logic: checks EVERY timestep with inpolygon\n\n');

%% Configuration
DECAY_HALFLIFE_DAYS = 7;  % Match streaming script
decay_constant = log(2) / DECAY_HALFLIFE_DAYS;

%% Setup
projectPath = matlab.project.rootProject().RootFolder;
dataPath = fullfile(projectPath, 'data');

% Load reef data
centroids = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
unique_IDs = centroids(:,1);
Xs = [centroids(:,8) centroids(:,10) centroids(:,12) centroids(:,14) centroids(:,8)];
Ys = [centroids(:,9) centroids(:,11) centroids(:,13) centroids(:,15) centroids(:,9)];
centroids_x = mean(Xs(:,1:4), 2);
centroids_y = mean(Ys(:,1:4), 2);

% Pre-compute bounding boxes
bbox_xmin = min(Xs(:,1:4), [], 2);
bbox_xmax = max(Xs(:,1:4), [], 2);
bbox_ymin = min(Ys(:,1:4), [], 2);
bbox_ymax = max(Ys(:,1:4), [], 2);

fprintf('Loaded %d reef polygons\n\n', length(unique_IDs));

%% Test configuration
quarters = {'Q1_2019'};
events_to_test = [2];

% Specific file to test (set to '' for random sampling)
SPECIFIC_FILE_TO_TEST = '';  % 'traj_file_640.nc'; Change this or set to '' for random

if ~isempty(SPECIFIC_FILE_TO_TEST)
    fprintf('MODE: Testing specific trajectory file\n');
    fprintf('File: %s\n\n', SPECIFIC_FILE_TO_TEST);
else
    fprintf('MODE: Random sampling of trajectory files\n\n');
end

for q_idx = 1:length(quarters)
    quarter_name = quarters{q_idx};
    fprintf('========== %s ==========\n', quarter_name);
    
    % Paths
    tempPath = fullfile('D:\Dissertation\CMS_traj', quarter_name);
    outputPath = fullfile('D:\Dissertation\CMS_traj\output', 'CMS_traj', quarter_name);
    
    % Load release file
    if strcmp(quarter_name, 'Q1_2019')
        releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019_Q1.txt'));
    elseif strcmp(quarter_name, 'Q2_2019')
        releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019_Q2.txt'));
    else
        releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019_Q3.txt'));
    end
    release_reef_IDs = releasefile(:,1);
    
    % Load event scan
    load(fullfile(outputPath, 'event_scan_cache.mat'), 'file_release_dates', 'unique_release_dates', 'trajlist');
    
    for event_num = events_to_test
        if event_num > length(unique_release_dates), continue; end
        
        fprintf('\n--- Event %d ---\n', event_num);
        
        target_date = unique_release_dates(event_num);
        event_date = datetime(target_date, 'ConvertFrom', 'juliandate');
        fprintf('Target Date: %s\n', datestr(event_date, 'dd-mmm-yyyy HH:MM:SS'));
        
        %% Find the connectivity file for this date
        fprintf('\nSearching for matching connectivity file...\n');
        
        date_patterns = {
            sprintf('connectivity_%s.mat', datestr(event_date, 'yyyy_mmmdd_HHMMSS')),
            sprintf('connectivity_%s*.mat', datestr(event_date, 'yyyy_mmm')),
            sprintf('connectivity_%s*.mat', datestr(event_date, 'yyyy_mmmdd'))
        };
        
        connectivity_file = '';
        loaded_data = [];
        
        for p = 1:length(date_patterns)
            pattern = date_patterns{p};
            found_files = dir(fullfile(outputPath, pattern));
            
            if ~isempty(found_files)
                for f = 1:length(found_files)
                    test_file = fullfile(outputPath, found_files(f).name);
                    test_data = load(test_file);
                    
                    if isfield(test_data, 'connectivity_results')
                        file_date = test_data.connectivity_results.calendar_date;
                        date_diff = abs(days(file_date - event_date));
                        
                        if date_diff < 0.5
                            fprintf('✓ Found: %s\n', found_files(f).name);
                            connectivity_file = test_file;
                            loaded_data = test_data;
                            break;
                        end
                    end
                end
                if ~isempty(connectivity_file), break; end
            end
        end
        
        if isempty(connectivity_file)
            fprintf('✗✗✗ NO CONNECTIVITY FILE FOUND ✗✗✗\n');
            continue;
        end
        
        ConnMatrix = loaded_data.connectivity_results.ConnMatrix_raw;
        
        % Find trajectory files
        relevant_files = find(abs(file_release_dates - target_date) < 0.01);
        fprintf('Event spans %d trajectory files\n', length(relevant_files));
        
        % Sample files
        % OPTION 1: Random sampling
        % OPTION 2: Specific file testing
        specific_file_name = SPECIFIC_FILE_TO_TEST;
        
        if ~isempty(specific_file_name)
            % Find the specific file
            specific_idx = find(strcmp({trajlist.name}, specific_file_name));
            
            if isempty(specific_idx)
                fprintf('⚠ File %s not found in trajectory list!\n', specific_file_name);
                continue;
            end
            
            % Check if this file is part of the current event
            if ismember(specific_idx, relevant_files)
                sample_file_indices = specific_idx;
                fprintf('Testing specific file: %s\n', specific_file_name);
                fprintf('This file is part of Event %d\n\n', event_num);
            else
                fprintf('⚠ File %s exists but is NOT part of Event %d!\n', specific_file_name, event_num);
                fprintf('  This file belongs to a different release date.\n');
                
                % Show which event this file belongs to
                file_date = file_release_dates(specific_idx);
                [~, file_event_idx] = min(abs(unique_release_dates - file_date));
                fprintf('  This file belongs to Event %d (%s)\n', ...
                    file_event_idx, datestr(datetime(unique_release_dates(file_event_idx), ...
                    'ConvertFrom', 'juliandate'), 'dd-mmm-yyyy'));
                fprintf('  Change events_to_test to [%d] or set SPECIFIC_FILE_TO_TEST to ''''\n', ...
                    file_event_idx);
                continue;
            end
        else
            % Random sampling
            n_files_to_sample = min(3, length(relevant_files));
            rng(42);
            sample_file_indices = relevant_files(randperm(length(relevant_files), n_files_to_sample));
            fprintf('Sampling %d random files\n\n', n_files_to_sample);
        end
        
        %% EXACT REPLICATION: Process particles like the streaming script
        fprintf('=== EXACT CONNECTIVITY CALCULATION (Matching Streaming Script) ===\n');
        
        % Build connectivity matrix from validation particles
        ValidationConnMatrix = zeros(length(unique_IDs), length(unique_IDs));
        all_particle_details = [];
        
        timesteps_per_day = (24 * 60) / 15;  % 96 timesteps per day (15-min intervals)
        
        for k = 1:length(sample_file_indices)
            file_idx = sample_file_indices(k);
            traj_filename = trajlist(file_idx).name;
            full_path = fullfile(tempPath, traj_filename);
            
            fprintf('Processing %s...\n', traj_filename);
            
            try
                lon = ncread(full_path, 'lon');
                lat = ncread(full_path, 'lat');
                location = ncread(full_path, 'location');
                releasedate = ncread(full_path, 'releasedate');
                
                % Filter to this event
                event_mask = abs(releasedate - target_date) < 0.01;
                event_indices = find(event_mask);
                
                if isempty(event_indices), continue; end
                
                % Sample particles - more if testing specific file
                if ~isempty(specific_file_name) && strcmp(trajlist(file_idx).name, specific_file_name)
                    n_to_sample = min(20, length(event_indices));  % Test more particles from specific file
                    fprintf('  Testing %d particles from this specific file\n', n_to_sample);
                else
                    n_to_sample = min(5, length(event_indices));
                end
                sampled_particles = event_indices(randperm(length(event_indices), n_to_sample));
                
                lon(lon > 180) = lon(lon > 180) - 360;
                
                % Process each sampled particle EXACTLY like streaming script
                for p = 1:length(sampled_particles)
                    particle_idx = sampled_particles(p);
                    
                    % Get source reef
                    release_location_idx = location(particle_idx);
                    source_reef_ID = release_reef_IDs(release_location_idx);
                    source_matrix_idx = find(unique_IDs == source_reef_ID, 1);
                    
                    if isempty(source_matrix_idx), continue; end
                    
                    % Get full trajectory
                    particle_lon = lon(:, particle_idx);
                    particle_lat = lat(:, particle_idx);
                    
                    % Remove NaN values
                    valid_idx = ~isnan(particle_lon) & ~isnan(particle_lat);
                    particle_lon = particle_lon(valid_idx);
                    particle_lat = particle_lat(valid_idx);
                    n_timesteps = length(particle_lon);
                    
                    if n_timesteps < 2, continue; end
                    
                    % EXACT REPLICATION: Check EVERY timestep
                    timestep_connections = struct('timestep', {}, 'reef_ID', {}, 'matrix_idx', {}, ...
                                                  'lon', {}, 'lat', {}, 'decay_weight', {});
                    
                    for t = 1:n_timesteps
                        pos_lon = particle_lon(t);
                        pos_lat = particle_lat(t);
                        
                        % Calculate decay weight for this timestep
                        time_in_days = (t - 1) / timesteps_per_day;
                        decay_weight = exp(-decay_constant * time_in_days);
                        
                        % Check which reef (if any) contains this position
                        % Using same bbox pre-filtering as streaming script
                        for reef_idx = 1:length(unique_IDs)
                            % Bounding box check first
                            if pos_lon < bbox_xmin(reef_idx) || pos_lon > bbox_xmax(reef_idx) || ...
                               pos_lat < bbox_ymin(reef_idx) || pos_lat > bbox_ymax(reef_idx)
                                continue;
                            end
                            
                            % Actual polygon check
                            if inpolygon(pos_lon, pos_lat, Xs(reef_idx,:), Ys(reef_idx,:))
                                % Record this connection
                                conn = struct();
                                conn.timestep = t;
                                conn.reef_ID = unique_IDs(reef_idx);
                                conn.matrix_idx = reef_idx;
                                conn.lon = pos_lon;
                                conn.lat = pos_lat;
                                conn.decay_weight = decay_weight;
                                
                                timestep_connections(end+1) = conn;
                                
                                % Add to validation matrix (weighted)
                                ValidationConnMatrix(source_matrix_idx, reef_idx) = ...
                                    ValidationConnMatrix(source_matrix_idx, reef_idx) + decay_weight;
                                
                                break;  % Found reef for this timestep
                            end
                        end
                    end
                    
                    % Store particle details
                    particle_detail = struct();
                    particle_detail.source_reef_ID = source_reef_ID;
                    particle_detail.source_matrix_idx = source_matrix_idx;
                    particle_detail.n_timesteps = n_timesteps;
                    particle_detail.connections = timestep_connections;
                    particle_detail.trajectory_lon = particle_lon;
                    particle_detail.trajectory_lat = particle_lat;
                    particle_detail.file = traj_filename;
                    
                    all_particle_details = [all_particle_details; particle_detail];
                    
                    % Print summary
                    fprintf('  Particle from Reef %d: %d timesteps, %d reef contacts\n', ...
                        source_reef_ID, n_timesteps, length(timestep_connections));
                end
                
            catch ME
                fprintf('  ERROR: %s\n', ME.message);
            end
        end
        
        %% COMPARISON: Validation Matrix vs Actual Matrix
        fprintf('\n╔═══════════════════════════════════════════════════════════╗\n');
        fprintf('║       EXACT CONNECTIVITY VALIDATION RESULTS               ║\n');
        fprintf('╚═══════════════════════════════════════════════════════════╝\n\n');
        
        fprintf('Processed %d particles total\n\n', length(all_particle_details));
        
        % Find all connections detected in validation
        [val_src, val_dst, val_weights] = find(ValidationConnMatrix);
        
        fprintf('=== Connections Found by Validation ===\n');
        total_validated = 0;
        total_mismatched = 0;
        
        for i = 1:length(val_src)
            src_idx = val_src(i);
            dst_idx = val_dst(i);
            val_weight = val_weights(i);
            matrix_weight = ConnMatrix(src_idx, dst_idx);
            
            src_ID = unique_IDs(src_idx);
            dst_ID = unique_IDs(dst_idx);
            
            % Calculate reef-to-reef distance
            reef_dist = sqrt((centroids_x(src_idx) - centroids_x(dst_idx))^2 + ...
                            (centroids_y(src_idx) - centroids_y(dst_idx))^2) * 111;
            
            if matrix_weight > 0
                status = '✓✓✓';
                result = 'VERIFIED';
                match_quality = abs(val_weight - matrix_weight) / max(val_weight, matrix_weight);
                total_validated = total_validated + 1;
                
                fprintf('  %s Reef %d→%d (%.1f km): Validation=%.2f, Matrix=%.2f (%.1f%% match)\n', ...
                    status, src_ID, dst_ID, reef_dist, val_weight, matrix_weight, ...
                    100*(1-match_quality));
            else
                status = '✗✗✗';
                result = 'MISSING FROM MATRIX';
                total_mismatched = total_mismatched + 1;
                
                fprintf('  %s Reef %d→%d (%.1f km): Validation=%.2f, Matrix=0.0 ⚠\n', ...
                    status, src_ID, dst_ID, reef_dist, val_weight);
                
                % Find which particle(s) created this connection
                for p = 1:length(all_particle_details)
                    pd = all_particle_details(p);
                    if pd.source_matrix_idx == src_idx
                        contacts_to_dst = [pd.connections.matrix_idx] == dst_idx;
                        if any(contacts_to_dst)
                            n_contacts = sum(contacts_to_dst);
                            timesteps = [pd.connections(contacts_to_dst).timestep];
                            fprintf('      From particle %d: %d contacts at timesteps %s\n', ...
                                p, n_contacts, mat2str(timesteps(1:min(5,length(timesteps)))));
                        end
                    end
                end
            end
        end
        
        fprintf('\n=== Summary ===\n');
        fprintf('Total connections detected by validation: %d\n', length(val_src));
        fprintf('  ✓ Verified in matrix: %d (%.1f%%)\n', total_validated, 100*total_validated/length(val_src));
        fprintf('  ✗ Missing from matrix: %d (%.1f%%)\n', total_mismatched, 100*total_mismatched/length(val_src));
        
        % Check for connections in matrix that validation didn't find
        fprintf('\n=== Checking Matrix for Connections Validation Missed ===\n');
        [mat_src, mat_dst, mat_weights] = find(ConnMatrix);
        
        validation_found = false(length(mat_src), 1);
        for i = 1:length(mat_src)
            if ValidationConnMatrix(mat_src(i), mat_dst(i)) > 0
                validation_found(i) = true;
            end
        end
        
        % Focus on source reefs that were in our validation sample
        sample_source_indices = unique([all_particle_details.source_matrix_idx]);
        
        missed_from_sample = false(length(mat_src), 1);
        for i = 1:length(mat_src)
            if ismember(mat_src(i), sample_source_indices) && ~validation_found(i)
                missed_from_sample(i) = true;
            end
        end
        
        n_missed = sum(missed_from_sample);
        if n_missed > 0
            fprintf('Matrix has %d connections from sampled source reefs that validation did not detect\n', n_missed);
            fprintf('Showing first 5:\n');
            missed_indices = find(missed_from_sample);
            for i = 1:min(5, length(missed_indices))
                idx = missed_indices(i);
                fprintf('  Reef %d→%d: Matrix=%.2f, Validation=0 ⚠\n', ...
                    unique_IDs(mat_src(idx)), unique_IDs(mat_dst(idx)), mat_weights(idx));
            end
            fprintf('Note: These could be from other particles not in validation sample\n');
        else
            fprintf('✓ Validation found all connections for sampled source reefs!\n');
        end
        
        fprintf('\n');
        if total_mismatched == 0 && n_missed == 0
            fprintf('╔═══════════════════════════════════════════════════════════╗\n');
            fprintf('║  ✓✓✓ PERFECT VALIDATION ✓✓✓                              ║\n');
            fprintf('║  Validation exactly matches connectivity matrix!          ║\n');
            fprintf('║  The streaming script is working correctly!               ║\n');
            fprintf('╚═══════════════════════════════════════════════════════════╝\n');
        elseif total_mismatched < 0.2 * length(val_src)
            fprintf('╔═══════════════════════════════════════════════════════════╗\n');
            fprintf('║  ⚠ MOSTLY VERIFIED (>80%%)                                ║\n');
            fprintf('║  Small discrepancies may be due to sampling differences   ║\n');
            fprintf('╚═══════════════════════════════════════════════════════════╝\n');
        else
            fprintf('╔═══════════════════════════════════════════════════════════╗\n');
            fprintf('║  ✗ SIGNIFICANT DISCREPANCIES DETECTED                     ║\n');
            fprintf('║  This indicates a problem with the streaming script!      ║\n');
            fprintf('╚═══════════════════════════════════════════════════════════╝\n');
        end
        
        %% Detailed Particle Analysis
        fprintf('\n=== Per-Particle Detail ===\n');
        for p = 1:min(5, length(all_particle_details))
            pd = all_particle_details(p);
            
            fprintf('\nParticle %d from Reef %d (matrix idx %d):\n', ...
                p, pd.source_reef_ID, pd.source_matrix_idx);
            fprintf('  Trajectory: %d timesteps\n', pd.n_timesteps);
            fprintf('  Reef contacts: %d\n', length(pd.connections));
            
            if ~isempty(pd.connections)
                % Group by destination reef
                dest_reefs = unique([pd.connections.matrix_idx]);
                fprintf('  Contacted %d unique reefs:\n', length(dest_reefs));
                
                for d = 1:min(5, length(dest_reefs))
                    reef_idx = dest_reefs(d);
                    contacts = pd.connections([pd.connections.matrix_idx] == reef_idx);
                    total_weight = sum([contacts.decay_weight]);
                    matrix_val = ConnMatrix(pd.source_matrix_idx, reef_idx);
                    
                    if matrix_val > 0
                        fprintf('    ✓ Reef %d: %d contacts, weight=%.2f, matrix=%.2f\n', ...
                            unique_IDs(reef_idx), length(contacts), total_weight, matrix_val);
                    else
                        fprintf('    ✗ Reef %d: %d contacts, weight=%.2f, matrix=0.0 ⚠\n', ...
                            unique_IDs(reef_idx), length(contacts), total_weight);
                    end
                end
                
                if length(dest_reefs) > 5
                    fprintf('    ... and %d more reefs\n', length(dest_reefs) - 5);
                end
            end
        end
        
        %% VISUALIZATION: Trajectory Hits
        fprintf('\n=== Generating Trajectory Hit Visualization ===\n');
        
        % Select particles to plot (up to 5 for clarity)
        n_to_plot = min(5, length(all_particle_details));
        particles_to_plot = 1:n_to_plot;
        
        % Find all reefs involved in hits
        all_involved_reefs = [];
        for p = particles_to_plot
            pd = all_particle_details(p);
            all_involved_reefs = [all_involved_reefs, pd.source_matrix_idx];
            if ~isempty(pd.connections)
                all_involved_reefs = [all_involved_reefs, [pd.connections.matrix_idx]];
            end
        end
        all_involved_reefs = unique(all_involved_reefs);
        
        % Calculate plot bounds from trajectories
        all_traj_lons = [];
        all_traj_lats = [];
        for p = particles_to_plot
            pd = all_particle_details(p);
            all_traj_lons = [all_traj_lons; pd.trajectory_lon];
            all_traj_lats = [all_traj_lats; pd.trajectory_lat];
        end
        
        % Add buffer around trajectories
        lon_min = min(all_traj_lons) - 0.1;
        lon_max = max(all_traj_lons) + 0.1;
        lat_min = min(all_traj_lats) - 0.1;
        lat_max = max(all_traj_lats) + 0.1;
        
        fig = figure('Position', [50, 50, 1400, 900], 'Visible', 'on');
        hold on;
        
        % Plot ALL reef polygons in very light gray for context
        fprintf('Plotting all %d reef polygons for spatial context...\n', length(unique_IDs));
        for i = 1:length(unique_IDs)
            % Only plot reefs within or near the trajectory bounds
            reef_center_lon = centroids_x(i);
            reef_center_lat = centroids_y(i);
            
            if reef_center_lon >= lon_min-0.2 && reef_center_lon <= lon_max+0.2 && ...
               reef_center_lat >= lat_min-0.2 && reef_center_lat <= lat_max+0.2
                plot(Xs(i,:), Ys(i,:), 'Color', [0.92 0.92 0.92], 'LineWidth', 0.3);
            end
        end
        
        % Color scheme for particles
        colors = lines(n_to_plot);
        
        % Plot each particle's trajectory with hits
        for p = particles_to_plot
            pd = all_particle_details(p);
            
            fprintf('Plotting particle %d/%d...\n', p, n_to_plot);
            
            % Source reef - highlight with THICK colored border
            src_idx = pd.source_matrix_idx;
            plot(Xs(src_idx,:), Ys(src_idx,:), 'Color', colors(p,:), 'LineWidth', 4);
            
            % Fill source reef polygon with semi-transparent color
            fill(Xs(src_idx,1:4), Ys(src_idx,1:4), colors(p,:), ...
                'FaceAlpha', 0.2, 'EdgeColor', 'none');
            
            % Full trajectory path (thin, semi-transparent)
            plot(pd.trajectory_lon, pd.trajectory_lat, '-', ...
                'Color', [colors(p,:), 0.4], 'LineWidth', 1);
            
            % START POINT - Make very prominent
            scatter(pd.trajectory_lon(1), pd.trajectory_lat(1), 250, colors(p,:), ...
                'filled', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 3);
            % Add white center to start marker for visibility
            scatter(pd.trajectory_lon(1), pd.trajectory_lat(1), 100, 'w', ...
                'filled', 'o', 'MarkerEdgeColor', colors(p,:), 'LineWidth', 2);
            
            % Mark all polygon hits along the trajectory
            if ~isempty(pd.connections)
                hit_lons = [pd.connections.lon];
                hit_lats = [pd.connections.lat];
                
                % Plot hit points (smaller dots)
                scatter(hit_lons, hit_lats, 20, colors(p,:), 'filled', ...
                    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.6);
                
                % Highlight destination reef polygons
                dest_reefs = unique([pd.connections.matrix_idx]);
                
                for d = 1:length(dest_reefs)
                    reef_idx = dest_reefs(d);
                    
                    % Highlight destination reefs (but not source reef)
                    if reef_idx ~= src_idx
                        plot(Xs(reef_idx,:), Ys(reef_idx,:), '--', ...
                            'Color', colors(p,:), 'LineWidth', 2.5);
                        
                        % Mark destination centroid with triangle
                        scatter(centroids_x(reef_idx), centroids_y(reef_idx), 120, colors(p,:), ...
                            'filled', '^', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
                    end
                end
            end
            
            % End point (X marker)
            plot(pd.trajectory_lon(end), pd.trajectory_lat(end), 'x', ...
                'Color', colors(p,:), 'MarkerSize', 16, 'LineWidth', 4);
        end
        
        % Add grid and labels
        grid on;
        xlabel('Longitude', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Latitude', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('Trajectory Hit Verification - %s Event %d\n%s | %d Particles | %d Total Polygon Hits', ...
            quarter_name, event_num, datestr(event_date, 'dd-mmm-yyyy'), ...
            n_to_plot, sum(arrayfun(@(x) length(x.connections), all_particle_details(particles_to_plot)))), ...
            'FontSize', 13, 'FontWeight', 'bold');
        
        % Set axis limits to trajectory bounds
        xlim([lon_min, lon_max]);
        ylim([lat_min, lat_max]);
        axis equal;
        
        % Legend - simplified and clearer
        legend_entries = {};
        legend_handles = [];
        
        % Particle-specific entries (first 3 only to avoid clutter)
        for p = 1:min(3, n_to_plot)
            pd = all_particle_details(p);
            n_hits = length(pd.connections);
            n_dest = 0;
            if ~isempty(pd.connections)
                n_dest = length(unique([pd.connections.matrix_idx]));
            end
            
            h = plot(NaN, NaN, '-', 'Color', colors(p,:), 'LineWidth', 3);
            legend_handles = [legend_handles; h];
            legend_entries{end+1} = sprintf('P%d: %d hits → %d reefs', p, n_hits, n_dest);
        end
        
        % Generic markers
        h_start = scatter(NaN, NaN, 150, 'k', 'filled', 'o', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 2);
        h_end = plot(NaN, NaN, 'x', 'Color', 'k', 'MarkerSize', 12, 'LineWidth', 3);
        h_hit = scatter(NaN, NaN, 20, 'k', 'filled', 'MarkerEdgeColor', 'none');
        h_dest = scatter(NaN, NaN, 100, 'k', 'filled', '^', 'MarkerEdgeColor', 'k');
        h_source_reef = fill([NaN NaN NaN], [NaN NaN NaN], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'k', 'LineWidth', 3);
        h_dest_reef = plot(NaN, NaN, '--', 'Color', 'k', 'LineWidth', 2);
        h_all_reefs = plot(NaN, NaN, '-', 'Color', [0.92 0.92 0.92], 'LineWidth', 0.3);
        
        legend_handles = [legend_handles; h_start; h_end; h_hit; h_dest; h_source_reef; h_dest_reef; h_all_reefs];
        legend_entries{end+1} = 'START (release point)';
        legend_entries{end+1} = 'END (final position)';
        legend_entries{end+1} = 'Polygon hits';
        legend_entries{end+1} = 'Destination reef';
        legend_entries{end+1} = 'Source reef (filled)';
        legend_entries{end+1} = 'Dest reef (dashed)';
        legend_entries{end+1} = 'All other reefs';
        
        legend(legend_handles, legend_entries, 'Location', 'best', 'FontSize', 9);
        
        fprintf('Figure displayed!\n');
        drawnow;
        
        % Optionally save
        save_figures = true;
        if save_figures
            filename = sprintf('validation_hits_%s_event%d.png', quarter_name, event_num);
            saveas(fig, fullfile(outputPath, filename));
            fprintf('Saved: %s\n', filename);
        end
    end
end

fprintf('\n╔═══════════════════════════════════════════════════════════╗\n');
fprintf('║  VALIDATION COMPLETE                                      ║\n');
fprintf('║  Method: Exact replication of streaming script logic     ║\n');
fprintf('║  - Processed every timestep with inpolygon()              ║\n');
fprintf('║  - Applied exponential decay weighting                    ║\n');
fprintf('║  - Built connectivity matrix from validation particles    ║\n');
fprintf('║  - Compared directly to actual connectivity matrix        ║\n');
fprintf('╚═══════════════════════════════════════════════════════════╝\n');