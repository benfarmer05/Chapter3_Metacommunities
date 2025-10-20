clear; clc

%% Initialize paths
projectPath = matlab.project.rootProject().RootFolder;
tempPath = fullfile(projectPath, 'temp');
dataPath = fullfile(projectPath, 'data');
outputPath = fullfile(projectPath, 'output');

fprintf('=== CONNECTIVITY ANALYSIS (Optimized) ===\n');

%% =================================================================
%% 1. LOAD RELEASE POINTS AND MAPPING
%% =================================================================
fprintf('\n1. Loading release points...\n');
relpoints = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));

unique_IDs = relpoints(:, 1);
centroid_lons = relpoints(:, 2);
centroid_lats = relpoints(:, 3);
n_locations = length(unique_IDs);
fprintf('Loaded %d locations from CSV\n', n_locations);

% Load release file
fprintf('Loading release file...\n');
releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019.txt'));
release_unique_IDs = releasefile(:, 1);
n_release_lines = length(release_unique_IDs);

% Create mapping: release file line → relpoints index
cms_to_relpoints_idx = zeros(n_release_lines, 1);
for i = 1:n_release_lines
    idx = find(unique_IDs == release_unique_IDs(i), 1);
    if ~isempty(idx)
        cms_to_relpoints_idx(i) = idx;
    end
end
fprintf('Mapped %d release entries\n', sum(cms_to_relpoints_idx > 0));

%% =================================================================
%% 2. PREPARE POLYGON DATA (vectorized format)
%% =================================================================
fprintf('\n2. Preparing polygon data...\n');

% Create polygon arrays in format needed by inpolygons
polyXs = [];
polyYs = [];
polyIDs = [];

for i = 1:n_locations
    % Get 4 vertices and close polygon
    vx = [relpoints(i, 8), relpoints(i, 10), relpoints(i, 12), relpoints(i, 14), relpoints(i, 8)];
    vy = [relpoints(i, 9), relpoints(i, 11), relpoints(i, 13), relpoints(i, 15), relpoints(i, 9)];
    
    polyXs = [polyXs vx];
    polyYs = [polyYs vy];
    polyIDs = [polyIDs repmat(i, 1, length(vx))];
end
fprintf('Created polygon arrays for %d locations\n', n_locations);

%% =================================================================
%% 3. CREATE DISEASE VIABILITY CURVE
%% =================================================================
fprintf('\n3. Creating disease viability curve...\n');

trajlist = dir(fullfile(tempPath, 'traj*.nc'));
if isempty(trajlist)
    error('No trajectory files found');
end

sample_file = fullfile(tempPath, trajlist(1).name);
time_data = ncread(sample_file, 'time');
timestep_seconds = mean(diff(time_data));
n_timesteps = length(time_data);

% Disease parameters
competency_start_days = 0;
max_dur_days = 14;

comp_start_steps = round(competency_start_days * 24 * 3600 / timestep_seconds);
max_dur_steps = round(max_dur_days * 24 * 3600 / timestep_seconds);

viability_window = min(comp_start_steps + max_dur_steps, n_timesteps);
decay_curve = zeros(n_timesteps, 1);

if viability_window > comp_start_steps
    k = log(0.01) / (viability_window - comp_start_steps);
    active_steps = (comp_start_steps + 1):viability_window;
    decay_curve(active_steps) = exp(k * (0:(length(active_steps)-1)));
end

fprintf('Viability: %.1f to %.1f days\n', ...
    comp_start_steps * timestep_seconds / 86400, ...
    viability_window * timestep_seconds / 86400);

% Detect CMS indexing
sample_location = ncread(sample_file, 'location');
cms_offset = (min(sample_location) == 0) * 1;
if cms_offset == 0
    fprintf('CMS indexing: 1-based (offset=0)\n');
else
    fprintf('CMS indexing: 0-based (offset=1)\n');
end

%% =================================================================
%% 4. BUILD CONNECTIVITY MATRIX (VECTORIZED - SERIAL ONLY)
%% =================================================================
fprintf('\n4. Building connectivity matrix...\n');
fprintf('Processing %d files...\n', length(trajlist));

ConnMatrix = zeros(n_locations, n_locations);

tic;
for file_idx = 1:length(trajlist)
    if mod(file_idx, 10) == 0
        fprintf('  File %d/%d (%.1f%%)...\n', file_idx, length(trajlist), 100*file_idx/length(trajlist));
    end
    
    filename = fullfile(tempPath, trajlist(file_idx).name);
    
    try
        lon = ncread(filename, 'lon');
        lat = ncread(filename, 'lat');
        location = ncread(filename, 'location');
        
        if max(lon(:)) > 180
            lon = mod(lon + 180, 360) - 180;
        end
        
        n_particles = size(lon, 2);
        
        % Process each particle
        for p = 1:n_particles
            release_line = location(p) + cms_offset;
            if release_line < 1 || release_line > length(cms_to_relpoints_idx)
                continue;
            end
            
            source_idx = cms_to_relpoints_idx(release_line);
            if source_idx == 0
                continue;
            end
            
            particle_lon = lon(:, p);
            particle_lat = lat(:, p);
            
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
                % Get vertices for this polygon
                poly_x = [relpoints(dest_idx, 8), relpoints(dest_idx, 10), ...
                          relpoints(dest_idx, 12), relpoints(dest_idx, 14)];
                poly_y = [relpoints(dest_idx, 9), relpoints(dest_idx, 11), ...
                          relpoints(dest_idx, 13), relpoints(dest_idx, 15)];
                
                % Check all timesteps at once using built-in inpolygon
                in_poly = inpolygon(valid_lon, valid_lat, poly_x, poly_y);
                
                if any(in_poly)
                    % Get timesteps when particle was in this polygon
                    hit_times = valid_times(in_poly);
                    hit_probs = decay_curve(hit_times);
                    
                    % Accumulate probability
                    connection_prob = 1 - prod(1 - hit_probs);
                    
                    if connection_prob > 0
                        ConnMatrix(source_idx, dest_idx) = ...
                            1 - (1 - ConnMatrix(source_idx, dest_idx)) * (1 - connection_prob);
                    end
                end
            end
        end
    catch ME
        warning('Error processing file %s: %s', trajlist(file_idx).name, ME.message);
    end
end

processing_time = toc;
fprintf('Processing complete in %.1f seconds (%.1f minutes)\n', processing_time, processing_time/60);

fprintf('Non-zero connections: %d / %d (%.3f%%)\n', ...
    nnz(ConnMatrix), numel(ConnMatrix), 100*nnz(ConnMatrix)/numel(ConnMatrix));

%% =================================================================
%% 5. ANALYZE CONNECTIVITY MATRIX
%% =================================================================
fprintf('\n5. Analyzing connectivity...\n');

% Normalize
row_sums = sum(ConnMatrix, 2);
ConnMatrix_normalized = ConnMatrix ./ row_sums;
ConnMatrix_normalized(isnan(ConnMatrix_normalized)) = 0;

% Calculate WCL
WCL = zeros(n_locations, 1);
for i = 1:n_locations
    if row_sums(i) > 0
        distances = zeros(n_locations, 1);
        for j = 1:n_locations
            lat1 = centroid_lats(i) * pi / 180;
            lat2 = centroid_lats(j) * pi / 180;
            lon1 = centroid_lons(i) * pi / 180;
            lon2 = centroid_lons(j) * pi / 180;
            
            dlat = lat2 - lat1;
            dlon = lon2 - lon1;
            
            a = sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2;
            c = 2 * atan2(sqrt(a), sqrt(1-a));
            distances(j) = 6371 * c;
        end
        WCL(i) = sum(ConnMatrix_normalized(i, :) .* distances') / sum(ConnMatrix_normalized(i, :));
    end
end

self_recruitment = diag(ConnMatrix_normalized);
fraction_exchanged = 1 - self_recruitment;

fprintf('Mean WCL: %.2f km (SD: %.2f km)\n', mean(WCL(WCL>0)), std(WCL(WCL>0)));
fprintf('Mean self-recruitment: %.3f\n', mean(self_recruitment));
fprintf('Mean fraction exchanged: %.3f\n', mean(fraction_exchanged));

%% =================================================================
%% 6. VISUALIZE
%% =================================================================
fprintf('\n6. Creating visualizations...\n');

figure('Position', [100 100 1000 800]);

subplot(2,2,1);
imagesc(log10(ConnMatrix_normalized + 1e-10));
colorbar;
title('Connectivity Matrix (log10)');
xlabel('Destination');
ylabel('Source');
axis square;

subplot(2,2,2);
histogram(WCL(WCL > 0), 30);
xlabel('WCL (km)');
ylabel('Frequency');
title(sprintf('Mean: %.1f km', mean(WCL(WCL>0))));

subplot(2,2,3);
histogram(self_recruitment, 30);
xlabel('Self-recruitment');
title(sprintf('Mean: %.3f', mean(self_recruitment)));

subplot(2,2,4);
histogram(fraction_exchanged, 30);
xlabel('Fraction Exchanged');
title(sprintf('Mean: %.3f', mean(fraction_exchanged)));

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
%% 7. SAVE RESULTS
%% =================================================================
fprintf('\n7. Saving results...\n');

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

output_table = table(unique_IDs, centroid_lons, centroid_lats, WCL, ...
                     self_recruitment, fraction_exchanged, ...
                     'VariableNames', {'unique_ID', 'lon', 'lat', 'WCL_km', ...
                     'self_recruitment', 'fraction_exchanged'});
writetable(output_table, fullfile(outputPath, 'connectivity_metrics.csv'));

fprintf('\nAnalysis complete!\n');
fprintf('Results saved to: %s\n', outputPath);