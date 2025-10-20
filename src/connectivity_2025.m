clear; clc

%% Initialize paths
projectPath = matlab.project.rootProject().RootFolder;
tempPath = fullfile(projectPath, 'temp');
dataPath = fullfile(projectPath, 'data');
outputPath = fullfile(projectPath, 'output');

fprintf('=== CONNECTIVITY ANALYSIS (Dobbelaere 2020 approach) ===\n');

%% =================================================================
%% 1. LOAD RELEASE POINTS (vertices/centroids)
%% =================================================================
fprintf('\n1. Loading release points...\n');
relpoints = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));

% Columns: unique_ID, centroid_lon, centroid_lat, coral metrics, v1-v4 vertices
unique_IDs = relpoints(:, 1);
centroid_lons = relpoints(:, 2);
centroid_lats = relpoints(:, 3);

n_locations = length(unique_IDs);
fprintf('Loaded %d release locations from CSV\n', n_locations);

% Load release file to map CMS line numbers to unique IDs
fprintf('Loading release file...\n');
releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019.txt'));
release_unique_IDs = releasefile(:, 1);  % First column is unique_ID
n_release_lines = length(release_unique_IDs);
fprintf('Release file has %d lines\n', n_release_lines);

% Create mapping: CMS location (line number) -> index in relpoints
% CMS location indexing: check if it's 0-based or 1-based
cms_to_relpoints_idx = zeros(n_release_lines, 1);
for i = 1:n_release_lines
    rel_ID = release_unique_IDs(i);
    idx = find(unique_IDs == rel_ID, 1);
    if ~isempty(idx)
        cms_to_relpoints_idx(i) = idx;
    end
end
fprintf('Mapped %d release file entries to centroid data\n', sum(cms_to_relpoints_idx > 0));

%% =================================================================
%% 2. CREATE DISEASE VIABILITY/DECAY CURVE
%% =================================================================
fprintf('\n2. Creating disease viability curve...\n');

% Get timestep from first trajectory file
trajlist = dir(fullfile(tempPath, 'traj*.nc'));
if isempty(trajlist)
    error('No trajectory files found in %s', tempPath);
end

sample_file = fullfile(tempPath, trajlist(1).name);
time_data = ncread(sample_file, 'time');
timestep_seconds = mean(diff(time_data));

% Disease parameters
competency_start_days = 0;      % Days until disease becomes viable
max_dur_days = 14;              % Maximum particle survival duration (days)

% Convert to timesteps
comp_start_steps = round(competency_start_days * 24 * 3600 / timestep_seconds);
max_dur_steps = round(max_dur_days * 24 * 3600 / timestep_seconds);

% Determine actual trajectory length
n_timesteps = length(time_data);
fprintf('Trajectory has %d timesteps (%.1f days)\n', n_timesteps, n_timesteps * timestep_seconds / 86400);

% Create exponential decay curve for disease viability
viability_window = min(comp_start_steps + max_dur_steps, n_timesteps);
decay_curve = zeros(n_timesteps, 1);

if viability_window > comp_start_steps
    k = log(0.01) / (viability_window - comp_start_steps);  % Decay to 1% at max duration
    active_steps = (comp_start_steps + 1):viability_window;
    decay_curve(active_steps) = exp(k * (0:(length(active_steps)-1)));
end

fprintf('Disease viability: steps %d to %d (%.1f to %.1f days)\n', ...
    comp_start_steps, viability_window, ...
    comp_start_steps * timestep_seconds / 86400, ...
    viability_window * timestep_seconds / 86400);

%% =================================================================
%% 3. PROCESS TRAJECTORY FILES TO BUILD CONNECTIVITY MATRIX
%% =================================================================
fprintf('\n3. Building connectivity matrix...\n');

% Check CMS location indexing (0-based or 1-based)
fprintf('Checking CMS location indexing...\n');
sample_location = ncread(sample_file, 'location');
min_loc = min(sample_location);
max_loc = max(sample_location);
fprintf('  CMS location range: %d to %d\n', min_loc, max_loc);

if min_loc == 0
    fprintf('  CMS uses 0-based indexing (location 0 = line 1 in release file)\n');
    cms_offset = 1;  % Add 1 to convert CMS location to MATLAB line number
else
    fprintf('  CMS uses 1-based indexing (location 1 = line 1 in release file)\n');
    cms_offset = 0;  % No offset needed
end

% Initialize connectivity matrix
ConnMatrix = zeros(n_locations, n_locations);

% Diagnostics
total_particles_processed = 0;
particles_per_source = zeros(n_locations, 1);
hits_per_dest = zeros(n_locations, 1);

fprintf('Processing %d trajectory files...\n', length(trajlist));

tic;
for file_idx = 1:length(trajlist)
    if mod(file_idx, 5) == 0
        fprintf('  File %d/%d (%.1f%%)...\n', file_idx, length(trajlist), 100*file_idx/length(trajlist));
    end
    
    filename = fullfile(tempPath, trajlist(file_idx).name);
    
    try
        % Read trajectory data
        lon = ncread(filename, 'lon');
        lat = ncread(filename, 'lat');
        location = ncread(filename, 'location');
        
        % Adjust longitude if needed (CMS outputs 0-360)
        if max(lon(:)) > 180
            lon = mod(lon + 180, 360) - 180;
        end
        
        n_particles = size(lon, 2);
        
        % Process particles in chunks
        chunk_size = 500;
        for chunk_start = 1:chunk_size:n_particles
            chunk_end = min(chunk_start + chunk_size - 1, n_particles);
            
            chunk_lon = lon(:, chunk_start:chunk_end);
            chunk_lat = lat(:, chunk_start:chunk_end);
            chunk_loc = location(chunk_start:chunk_end);
            
            % Process each particle
            for p = 1:(chunk_end - chunk_start + 1)
                cms_location = chunk_loc(p);
                
                % Convert CMS location to release file line number
                release_line = cms_location + cms_offset;
                
                % Check bounds
                if release_line < 1 || release_line > n_release_lines
                    continue;
                end
                
                % Get index in relpoints array
                source_idx = cms_to_relpoints_idx(release_line);
                if source_idx == 0
                    continue;  % No mapping found
                end
                
                total_particles_processed = total_particles_processed + 1;
                particles_per_source(source_idx) = particles_per_source(source_idx) + 1;
                
                particle_lon = chunk_lon(:, p);
                particle_lat = chunk_lat(:, p);
                
                % Remove NaN positions
                valid_idx = ~isnan(particle_lon) & ~isnan(particle_lat);
                if ~any(valid_idx)
                    continue;
                end
                
                valid_lon = particle_lon(valid_idx);
                valid_lat = particle_lat(valid_idx);
                valid_times = find(valid_idx);
                
                % Check each destination polygon
                for dest_idx = 1:n_locations
                    % Get vertices from CSV columns 8-15
                    poly_x = [relpoints(dest_idx, 8), relpoints(dest_idx, 10), ...
                              relpoints(dest_idx, 12), relpoints(dest_idx, 14)];
                    poly_y = [relpoints(dest_idx, 9), relpoints(dest_idx, 11), ...
                              relpoints(dest_idx, 13), relpoints(dest_idx, 15)];
                    
                    % Vectorized point-in-polygon check
                    in_poly = inpolygon(valid_lon, valid_lat, poly_x, poly_y);
                    
                    if any(in_poly)
                        hits_per_dest(dest_idx) = hits_per_dest(dest_idx) + 1;
                        
                        % Get timesteps where particle was in this polygon
                        hit_times = valid_times(in_poly);
                        hit_probs = decay_curve(hit_times);
                        
                        % Accumulate probability: 1 - prod(1 - p_i)
                        connection_prob = 1 - prod(1 - hit_probs);
                        
                        if connection_prob > 0
                            ConnMatrix(source_idx, dest_idx) = ...
                                1 - (1 - ConnMatrix(source_idx, dest_idx)) * (1 - connection_prob);
                        end
                    end
                end
            end
        end
        
    catch ME
        warning('Error processing file %s: %s', trajlist(file_idx).name, ME.message);
        continue;
    end
end

processing_time = toc;
fprintf('Connectivity matrix built in %.1f seconds\n', processing_time);

% Print diagnostics
fprintf('\n=== DIAGNOSTICS ===\n');
fprintf('Total particles processed: %d\n', total_particles_processed);
fprintf('Sources with particles: %d / %d\n', sum(particles_per_source > 0), n_locations);
fprintf('Destinations with hits: %d / %d\n', sum(hits_per_dest > 0), n_locations);
fprintf('Non-zero connections: %d / %d (%.2f%%)\n', ...
    nnz(ConnMatrix), numel(ConnMatrix), 100*nnz(ConnMatrix)/numel(ConnMatrix));
fprintf('Source IDs present in trajectories: [%d ... %d]\n', ...
    min(unique_IDs(particles_per_source > 0)), max(unique_IDs(particles_per_source > 0)));

%% =================================================================
%% 4. NORMALIZE AND ANALYZE CONNECTIVITY MATRIX
%% =================================================================
fprintf('\n4. Analyzing connectivity...\n');

% Normalize rows to sum to 1
row_sums = sum(ConnMatrix, 2);
ConnMatrix_normalized = ConnMatrix ./ row_sums;
ConnMatrix_normalized(isnan(ConnMatrix_normalized)) = 0;

% Calculate Weighted Connectivity Length (WCL)
WCL = zeros(n_locations, 1);
for i = 1:n_locations
    if row_sums(i) > 0
        distances = zeros(n_locations, 1);
        for j = 1:n_locations
            % Haversine distance formula - INLINE, NO FUNCTION CALL
            lat1_rad = centroid_lats(i) * pi / 180;
            lat2_rad = centroid_lats(j) * pi / 180;
            lon1_rad = centroid_lons(i) * pi / 180;
            lon2_rad = centroid_lons(j) * pi / 180;
            
            dlat = lat2_rad - lat1_rad;
            dlon = lon2_rad - lon1_rad;
            
            a = sin(dlat/2)^2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon/2)^2;
            c = 2 * atan2(sqrt(a), sqrt(1-a));
            distances(j) = 6371 * c;  % Earth radius = 6371 km
        end
        WCL(i) = sum(ConnMatrix_normalized(i, :) .* distances') / sum(ConnMatrix_normalized(i, :));
    end
end

% Self-recruitment
self_recruitment = diag(ConnMatrix_normalized);

% Fraction exchanged
fraction_exchanged = 1 - self_recruitment;

% Summary statistics
fprintf('Mean WCL: %.2f km (SD: %.2f km)\n', mean(WCL(WCL>0)), std(WCL(WCL>0)));
fprintf('Mean self-recruitment: %.3f\n', mean(self_recruitment));
fprintf('Mean fraction exchanged: %.3f\n', mean(fraction_exchanged));

%% =================================================================
%% 5. VISUALIZE CONNECTIVITY MATRIX
%% =================================================================
fprintf('\n5. Creating visualizations...\n');

% Plot 1: Connectivity matrix heatmap
figure('Position', [100 100 1000 800]);

subplot(2,2,1);
imagesc(log10(ConnMatrix_normalized + 1e-10));
colorbar;
title('Connectivity Matrix (log10 scale)');
xlabel('Destination Location');
ylabel('Source Location');
axis square;

% Plot 2: Connectivity metrics
subplot(2,2,2);
histogram(WCL(WCL > 0), 30);
xlabel('Weighted Connectivity Length (km)');
ylabel('Frequency');
title(sprintf('Mean WCL: %.1f km', mean(WCL(WCL>0))));

subplot(2,2,3);
histogram(self_recruitment, 30);
xlabel('Self-recruitment');
ylabel('Frequency');
title(sprintf('Mean: %.3f', mean(self_recruitment)));

subplot(2,2,4);
histogram(fraction_exchanged, 30);
xlabel('Fraction Exchanged');
ylabel('Frequency');
title(sprintf('Mean: %.3f', mean(fraction_exchanged)));

sgtitle(sprintf('Connectivity Analysis - %d locations, %.1f day trajectories', ...
    n_locations, n_timesteps * timestep_seconds / 86400));

% Plot 3: Spatial connectivity map
figure('Position', [100 100 1200 600]);

subplot(1,2,1);
scatter(centroid_lons, centroid_lats, 50, WCL, 'filled');
colorbar;
title('Weighted Connectivity Length');
xlabel('Longitude');
ylabel('Latitude');
axis equal tight;

subplot(1,2,2);
scatter(centroid_lons, centroid_lats, 50, self_recruitment, 'filled');
colorbar;
title('Self-recruitment');
xlabel('Longitude');
ylabel('Latitude');
axis equal tight;

%% =================================================================
%% 6. SAVE RESULTS
%% =================================================================
fprintf('\n6. Saving results...\n');

% Save connectivity matrix and metrics
results = struct();
results.ConnMatrix = ConnMatrix;
results.ConnMatrix_normalized = ConnMatrix_normalized;
results.unique_IDs = unique_IDs;
results.centroid_lons = centroid_lons;
results.centroid_lats = centroid_lats;
results.WCL = WCL;
results.self_recruitment = self_recruitment;
results.fraction_exchanged = fraction_exchanged;
results.disease_viability_curve = decay_curve;
results.timestep_seconds = timestep_seconds;

save(fullfile(outputPath, 'connectivity_results.mat'), 'results');

% Export to CSV
output_table = table(unique_IDs, centroid_lons, centroid_lats, WCL, ...
                     self_recruitment, fraction_exchanged, ...
                     'VariableNames', {'unique_ID', 'lon', 'lat', 'WCL_km', ...
                     'self_recruitment', 'fraction_exchanged'});
writetable(output_table, fullfile(outputPath, 'connectivity_metrics.csv'));

fprintf('\nAnalysis complete!\n');
fprintf('Results saved to: %s\n', outputPath);