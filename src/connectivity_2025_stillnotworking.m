clear; clc

%% Initialize paths
projectPath = matlab.project.rootProject().RootFolder;
tempPath = fullfile(projectPath, 'temp');
dataPath = fullfile(projectPath, 'data');
outputPath = fullfile(projectPath, 'output');

fprintf('=== CONNECTIVITY ANALYSIS BY RELEASE EVENT ===\n');

%% =================================================================
%% 1. LOAD DATA
%% =================================================================
fprintf('\n1. Loading data...\n');

% Load centroids
relpoints = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
unique_IDs = relpoints(:, 1);
centroid_lons = relpoints(:, 2);
centroid_lats = relpoints(:, 3);
n_locations = length(unique_IDs);

% Load release file
releasefile = readmatrix(fullfile(tempPath, 'ReleaseFile_USVI_2019.txt'));
release_unique_IDs = releasefile(:, 1);

% Create mapping
cms_to_relpoints_idx = zeros(length(release_unique_IDs), 1);
for i = 1:length(release_unique_IDs)
    idx = find(unique_IDs == release_unique_IDs(i), 1);
    if ~isempty(idx)
        cms_to_relpoints_idx(i) = idx;
    end
end

fprintf('  %d locations\n', n_locations);

%% =================================================================
%% 2. PREPARE POLYGON DATA (like original script)
%% =================================================================
fprintf('\n2. Preparing polygons for inpolygons...\n');

% Create polygon matrix like original: [Xs; Ys; IDs]
Xs = [];
Ys = [];
IDs = [];

for i = 1:n_locations
    vx = [relpoints(i,8), relpoints(i,10), relpoints(i,12), relpoints(i,14), relpoints(i,8)];
    vy = [relpoints(i,9), relpoints(i,11), relpoints(i,13), relpoints(i,15), relpoints(i,9)];
    
    Xs = [Xs vx];
    Ys = [Ys vy];
    IDs = [IDs repmat(i, 1, length(vx))];
end

reefpolymat = [Xs; Ys; IDs];
fprintf('  Polygon matrix created\n');

%% =================================================================
%% 3. CREATE DECAY CURVE (following Dobbelaere 2020)
%% =================================================================
fprintf('\n3. Creating decay curve...\n');

trajlist = dir(fullfile(tempPath, 'traj*.nc'));
sample_file = fullfile(tempPath, trajlist(1).name);
time_data = ncread(sample_file, 'time');

% Calculate timestep ignoring NaNs
time_diffs = diff(time_data);
timestep_sec = mean(time_diffs(~isnan(time_diffs)));

% Disease parameters (following Dobbelaere 2020)
% Particle half-life approach
halflife_days = 7;  % Half-life of disease agents in days
halflife_sec = halflife_days * 24 * 3600;

% Calculate max trajectory duration we might see (14 days)
max_trajectory_days = 14;
max_steps = ceil(max_trajectory_days * 24 * 3600 / timestep_sec);

% Create decay curve based on exponential decay with half-life
% N(t) = N0 * exp(-lambda * t), where lambda = ln(2) / halflife
lambda = log(2) / (halflife_sec / timestep_sec);  % Decay constant in steps
yt = exp(-lambda * (0:(max_steps-1)))';  % Exponential decay

fprintf('  Timestep: %.1f sec (%.2f hours)\n', timestep_sec, timestep_sec/3600);
fprintf('  Half-life: %.1f days (%.1f steps)\n', halflife_days, halflife_sec/timestep_sec);
fprintf('  Max trajectory: %d steps (%.1f days)\n', max_steps, max_steps*timestep_sec/86400);
fprintf('  Decay curve: start=%.4f, after halflife=%.4f, end=%.4f\n', ...
    yt(1), yt(round(halflife_sec/timestep_sec)), yt(end));

n_timesteps = max_steps;  % Use this as our reference length

% Plot decay curve over real time
figure('Position', [100, 100, 800, 400]);
time_days = (0:(max_steps-1)) * timestep_sec / 86400;  % Convert steps to days
plot(time_days, yt, 'LineWidth', 2);
hold on;
plot(halflife_days, 0.5, 'r*', 'MarkerSize', 12, 'LineWidth', 2);  % Mark half-life
xlabel('Time (days)');
ylabel('Disease Agent Viability');
title(sprintf('Decay Curve (Half-life = %.1f days)', halflife_days));
grid on;
legend('Viability', 'Half-life point', 'Location', 'northeast');
ylim([0, 1.1]);
xlim([0, max_trajectory_days]);

% Detect CMS indexing
sample_location = ncread(sample_file, 'location');
cms_offset = (min(sample_location) == 0) * 1;

%% =================================================================
%% 4. IDENTIFY RELEASE EVENTS
%% =================================================================
fprintf('\n4. Identifying release events...\n');

all_release_dates = [];
for i = 1:length(trajlist)
    filename = fullfile(tempPath, trajlist(i).name);
    releasedate = ncread(filename, 'releasedate');
    all_release_dates = [all_release_dates; releasedate(:)];
end

unique_release_dates = unique(all_release_dates);
n_events = length(unique_release_dates);
release_datetimes = datetime(unique_release_dates, 'ConvertFrom', 'juliandate');

fprintf('  Found %d events\n', n_events);
for i = 1:n_events
    fprintf('    Event %d: %s\n', i, datestr(release_datetimes(i)));
end

%% =================================================================
%% 5. BUILD CONNECTIVITY (LIKE ORIGINAL SCRIPT)
%% =================================================================
fprintf('\n5. Building connectivity matrices...\n');

ConnMatrices = zeros(n_locations, n_locations, n_events);

total_start = tic;

for event_idx = 1:n_events
    fprintf('\n--- EVENT %d/%d: %s ---\n', event_idx, n_events, datestr(release_datetimes(event_idx)));
    event_release_date = unique_release_dates(event_idx);
    event_start = tic;
    
    % Process all files for this event
    for file_idx = 1:length(trajlist)
        if mod(file_idx, 50) == 0
            fprintf('  File %d/%d...\n', file_idx, length(trajlist));
        end
        
        filename = fullfile(tempPath, trajlist(file_idx).name);
        
        try
            % Read data
            lon = ncread(filename, 'lon');
            lat = ncread(filename, 'lat');
            location = ncread(filename, 'location');
            releasedate = ncread(filename, 'releasedate');
            
            % Filter to this event's particles
            event_mask = abs(releasedate - event_release_date) < 0.01;
            if ~any(event_mask)
                continue;
            end
            
            % Adjust longitude
            if max(lon(:)) > 180
                lon_adj = lon;
            else
                lon_adj = lon + 360;  % Match original: +360 for CMS 0-360 output
            end
            
            % Get event particles only
            lon_event = lon_adj(:, event_mask);
            lat_event = lat(:, event_mask);
            loc_event = location(event_mask);
            
            % RUN INPOLYGONS ONCE (like original)
            [poly_in, poly_ind] = inpolygons(lon_event, lat_event, reefpolymat(1,:), reefpolymat(2,:));
            
            % Apply decay curve to each particle (like original yt_mat)
            % Handle particles that may have different lengths due to early exit
            yt_mat = zeros(size(poly_in));
            for p = 1:size(poly_in, 2)
                particle_length = size(poly_in, 1);
                % Use decay curve up to particle length (or pad if particle is longer)
                if particle_length <= length(yt)
                    yt_mat(:, p) = poly_in(:, p) .* yt(1:particle_length);
                else
                    % Particle lasted longer than our decay curve - extend with zeros
                    yt_extended = [yt; zeros(particle_length - length(yt), 1)];
                    yt_mat(:, p) = poly_in(:, p) .* yt_extended;
                end
            end
            
            % ============================================================
            % BUILD CONNECTIVITY FROM THIS FILE'S PARTICLES
            % ============================================================
            % Loop through each particle that was released in this event
            for p = 1:size(yt_mat, 2)
                
                % --- Step 1: Get source location ---
                % Convert CMS location (line number) to release file line
                release_line = loc_event(p) + cms_offset;
                if release_line < 1 || release_line > length(cms_to_relpoints_idx)
                    continue;  % Invalid line number
                end
                
                % Map release file line to our relpoints index (source)
                source_idx = cms_to_relpoints_idx(release_line);
                if source_idx == 0
                    continue;  % No mapping found
                end
                
                % Skip if particle never entered any polygon
                if sum(yt_mat(:,p)) == 0
                    continue;
                end
                
                % --- Step 2: Loop through each timestep for this particle ---
                % poly_ind is a CELL ARRAY: poly_ind{timestep, particle}
                % Each cell contains the polygon ID(s) the particle was in
                particle_length = size(poly_in, 1);  % Actual length of this particle's trajectory
                
                for t = 1:particle_length
                    
                    % Check if particle was in any polygon at this timestep
                    if ~poly_in(t, p)
                        continue;  % Not in any polygon at time t
                    end
                    
                    % --- Step 3: Get destination polygon ID(s) ---
                    % poly_ind{t,p} returns the polygon ID(s) at timestep t
                    % (can be empty, single ID, or multiple IDs if overlapping)
                    dest_ids = poly_ind{t, p};
                    
                    if isempty(dest_ids)
                        continue;  % No polygon ID found (shouldn't happen if poly_in=true)
                    end
                    
                    % --- Step 4: Apply decay curve weight ---
                    % Get disease viability at this timestep
                    weight = yt(t);
                    if weight == 0
                        continue;  % Outside viability window
                    end
                    
                    % --- Step 5: Accumulate connectivity ---
                    % For each destination polygon this particle hit at time t
                    for dest_idx = dest_ids(:)'
                        % Accumulate probability: 1 - (1-old) * (1-new)
                        % This is the "probability of NOT missing both connections"
                        ConnMatrices(source_idx, dest_idx, event_idx) = ...
                            1 - (1 - ConnMatrices(source_idx, dest_idx, event_idx)) * (1 - weight);
                    end
                end
            end
            
            % Print diagnostic for first file with particles
            if file_idx == find(cellfun(@(x) any(abs(x - event_release_date) < 0.01), ...
                    arrayfun(@(f) ncread(fullfile(tempPath, f.name), 'releasedate'), trajlist, 'UniformOutput', false)), 1)
                fprintf('    First file diagnostic:\n');
                fprintf('      Particles processed: %d\n', size(yt_mat, 2));
                fprintf('      Total polygon hits: %d\n', nnz(poly_in));
                fprintf('      Current connections: %d\n', nnz(ConnMatrices(:,:,event_idx)));
            end
            
        catch ME
            warning('Error in file %s: %s', trajlist(file_idx).name, ME.message);
        end
    end
    
    event_time = toc(event_start);
    fprintf('  Event complete in %.1f min\n', event_time/60);
    fprintf('  Non-zero: %d (%.3f%%)\n', nnz(ConnMatrices(:,:,event_idx)), ...
        100*nnz(ConnMatrices(:,:,event_idx))/numel(ConnMatrices(:,:,event_idx)));
    
    elapsed = toc(total_start);
    est_remaining = (elapsed / event_idx) * (n_events - event_idx);
    fprintf('  Estimated remaining: %.1f min\n', est_remaining/60);
end

fprintf('\n=== COMPLETE in %.1f minutes ===\n', toc(total_start)/60);

%% =================================================================
%% 6. SAVE
%% =================================================================
fprintf('\n6. Saving...\n');

results = struct();
results.ConnMatrices = ConnMatrices;
results.release_dates = unique_release_dates;
results.release_datetimes = release_datetimes;
results.unique_IDs = unique_IDs;
results.centroid_lons = centroid_lons;
results.centroid_lats = centroid_lats;

save(fullfile(outputPath, 'connectivity_by_event.mat'), 'results', '-v7.3');

for i = 1:n_events
    fname = sprintf('conn_event_%02d_%s.csv', i, datestr(release_datetimes(i), 'yyyy-mm-dd'));
    writematrix(ConnMatrices(:,:,i), fullfile(outputPath, fname));
end

fprintf('Done!\n');