clear; clc

%% setup
projectPath = matlab.project.rootProject().RootFolder;
tempPath = fullfile(projectPath, 'temp');
dataPath = fullfile(projectPath, 'data');
outputPath = fullfile(projectPath, 'output');
seascapePath = fullfile(outputPath, 'seascape_SIR');

RECALCULATE_CONNECTIVITY = false;  % Set to false to load from temp/P.mat
RUN_PARAMETER_SWEEP = true;        % Set to false to load from temp/ParameterSweep_Results.mat
USE_PARALLEL = true;               % Set to true to use parallel computing (requires Parallel Computing Toolbox)

% ===== NEW FILTERING TOGGLES =====
FILTER_LOW_COVER = false;           % Remove sites with <1% total coral cover
FILTER_MISSING_GROUPS = false;       % Remove sites with 0% in any susceptibility group
FILTER_NO_TRAJECTORIES = true;      % Remove sites that never produced trajectories
% ==================================

% ===== 2018 CONNECTIVITY BACKFILL =====
USE_2018_BACKFILL = true;           % If true, use Jan 2019 connectivity for any dates before 2019
% =======================================

% DATE_RANGE = [];  % Empty = create/load all connectivity matrices
% DATE_RANGE = [datetime(2019,1,1), datetime(2019,3,31)];  % just Q1
% DATE_RANGE = [datetime(2019,1,1), datetime(2019,6,30)];  % Q1-Q2
% DATE_RANGE = [datetime(2019,1,1), datetime(2019,12,31)];  % Q1-Q4
DATE_RANGE = [datetime(2018,12,1), datetime(2019,12,31)];  % 2018-2019 example

%% Load reef data from CSV

reefDataFile = fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv');
reefData_original = readtable(reefDataFile);

fprintf('========================================\n');
fprintf('INITIAL REEF DATA\n');
fprintf('Total sites loaded: %d\n', height(reefData_original));
fprintf('========================================\n\n');

%% STAGE 1: Apply site filtering (Filters A & B - CSV-based only)

% Store original indices before filtering
reefData = reefData_original;
reefData.original_index = (1:height(reefData))';
sites_to_keep = true(height(reefData), 1);

% Filter a) Sites with <1% coral cover
if FILTER_LOW_COVER
    low_cover_mask = reefData.mean_coral_cover < 0.01;
    sites_removed_low = sum(low_cover_mask);
    sites_to_keep = sites_to_keep & ~low_cover_mask;
    
    fprintf('========================================\n');
    fprintf('FILTER A: Low Cover Sites (<1%%)\n');
    fprintf('Sites removed: %d\n', sites_removed_low);
    fprintf('Sites remaining: %d\n', sum(sites_to_keep));
    fprintf('========================================\n\n');
end

% Filter b) Sites with 0% in any susceptibility group
if FILTER_MISSING_GROUPS
    missing_groups_mask = (reefData.low_coral_cover == 0) | ...
                          (reefData.moderate_coral_cover == 0) | ...
                          (reefData.high_coral_cover == 0);
    sites_removed_missing = sum(missing_groups_mask & sites_to_keep);
    sites_to_keep = sites_to_keep & ~missing_groups_mask;
    
    fprintf('========================================\n');
    fprintf('FILTER B: Missing Susceptibility Groups\n');
    fprintf('Sites removed: %d\n', sites_removed_missing);
    fprintf('Sites remaining: %d\n', sum(sites_to_keep));
    fprintf('========================================\n\n');
end

% Apply filter to reef data
if FILTER_LOW_COVER || FILTER_MISSING_GROUPS
    stage1_indices = reefData.original_index(sites_to_keep);
    reefData = reefData(sites_to_keep, :);
    
    fprintf('========================================\n');
    fprintf('STAGE 1 FILTERING SUMMARY (A & B)\n');
    fprintf('Original sites: %d\n', length(sites_to_keep));
    fprintf('Filtered sites: %d\n', height(reefData));
    fprintf('Reduction: %.1f%%\n', 100 * (1 - height(reefData)/length(sites_to_keep)));
    fprintf('========================================\n\n');
else
    stage1_indices = (1:height(reefData))';
end

%% STAGE 1 CACHE: Create or load existing connectivity matrices (filters A & B only)

% Generate connectivity (conn_structs) cache filename based on date range
if ~isempty(DATE_RANGE)
    stage1_cacheFilename = sprintf('conn_structs_STAGE1_%s_to_%s', ...
                           string(DATE_RANGE(1), 'yyyyMMdd'), ...
                           string(DATE_RANGE(2), 'yyyyMMdd'));
else
    stage1_cacheFilename = 'conn_structs_STAGE1';  % Full dataset
end

% Add filtering suffix to cache filename
if FILTER_LOW_COVER || FILTER_MISSING_GROUPS
    filter_suffix = '';
    if FILTER_LOW_COVER; filter_suffix = [filter_suffix '_lowcov']; end
    if FILTER_MISSING_GROUPS; filter_suffix = [filter_suffix '_missing']; end
    stage1_cacheFilename = [stage1_cacheFilename filter_suffix];
end

% Add 2018 backfill suffix to cache filename
if USE_2018_BACKFILL && ~isempty(DATE_RANGE) && year(DATE_RANGE(1)) < 2019
    stage1_cacheFilename = [stage1_cacheFilename '_2018bf'];
end

stage1_cacheFilename = [stage1_cacheFilename '.mat'];
stage1_cacheFilePath = fullfile(tempPath, stage1_cacheFilename);

if RECALCULATE_CONNECTIVITY
    fprintf('========================================\n');
    fprintf('STAGE 1: Loading connectivity matrices from disk\n');
    if ~isempty(DATE_RANGE)
        fprintf('Date range filter: %s to %s\n', ...
                string(DATE_RANGE(1), 'dd-MMM-yyyy'), ...
                string(DATE_RANGE(2), 'dd-MMM-yyyy'));
    end
    if USE_2018_BACKFILL && ~isempty(DATE_RANGE) && year(DATE_RANGE(1)) < 2019
        fprintf('2018 backfill: ENABLED (using Jan 2019 connectivity)\n');
    end
    fprintf('========================================\n\n');
    
    quarters = {'Q1_2019', 'Q2_2019', 'Q3_2019', 'Q4_2019'};
    conn_structs_stage1 = struct();
    connCount = 0;

    % NOTE: decay_weights normalization - currently using fixed value
    % TODO: Properly define decay_weights vector in future revision
    norm_val = 65 * 727.7434; % Placeholder: 65*sum(decay_weights)

    % Store first 2019 matrix for 2018 backfill
    first_2019_matrix = [];
    first_2019_date = datetime(2019,1,1);

    for q = 1:length(quarters)
        fprintf('Processing quarter: %s\n', quarters{q});
        quarterPath = fullfile(outputPath, 'CMS_traj', quarters{q});
        
        if exist(quarterPath, 'dir') == 0
            warning('  --> Quarter folder not found: %s\n', quarterPath);
            continue;
        end
        
        myFiles = dir(fullfile(quarterPath, 'connectivity_*.mat'));
        fprintf('  Found %d connectivity files\n', length(myFiles));
        filesLoaded = 0;
        
        for i = 1:length(myFiles)
            filename = myFiles(i).name;
            filenam = fullfile(quarterPath, filename);
            
            tokens = regexp(filename, 'connectivity_(\d{4})_(\w{3})(\d{2})_\d+\.mat', 'tokens');
            if ~isempty(tokens)
                year_str = tokens{1}{1};
                month_str = tokens{1}{2};
                day_str = tokens{1}{3};
                Dat = datetime(str2double(year_str), month(datetime(['1-' month_str '-2000'])), str2double(day_str), 12, 0, 0);
            else
                warning('Could not parse date from filename: %s', filename);
                continue;
            end
            
            if ~isempty(DATE_RANGE)
                if Dat < DATE_RANGE(1) || Dat > DATE_RANGE(2)
                    continue;
                end
            end
            
            connCount = connCount + 1;
            filesLoaded = filesLoaded + 1;
            fprintf('  [%d] Loading: %s (Date: %s)\n', connCount, filename, string(Dat, 'dd-MMM-yyyy'));
            
            Con = load(filenam);
            Dat_actual = Con.connectivity_results.calendar_date;
            if abs(days(Dat - Dat_actual)) > 0.1
                warning('Filename date mismatch: %s vs %s', string(Dat, 'dd-MMM-yyyy'), string(Dat_actual, 'dd-MMM-yyyy'));
            end
            
            conn_structs_stage1(connCount).DY = day(Dat);
            conn_structs_stage1(connCount).MO = month(Dat);
            conn_structs_stage1(connCount).YR = year(Dat);
            conn_structs_stage1(connCount).Date = Dat;
            
            conmat = sparse(Con.connectivity_results.ConnMatrix_raw);
            
            % ===== APPLY SITE FILTERING TO CONNECTIVITY (STAGE 1: A & B ONLY) =====
            if FILTER_LOW_COVER || FILTER_MISSING_GROUPS
                conmat = conmat(stage1_indices, stage1_indices);
            end
            % =======================================================================
            
            % Remove self retention (diagonal)
            conmat = spdiags(zeros(size(conmat,1),1), 0, conmat);
            
            % Normalize
            conmat = conmat ./ norm_val;
            
            conn_structs_stage1(connCount).full = conmat;
            
            % Store first 2019 matrix for backfill
            if isempty(first_2019_matrix) && year(Dat) == 2019
                first_2019_matrix = conmat;
                first_2019_date = Dat;
            end
        end
        
        if filesLoaded > 0
            fprintf('  Quarter %s complete: %d matrices loaded\n\n', quarters{q}, filesLoaded);
        else
            fprintf('  Quarter %s: No files in date range\n\n', quarters{q});
        end
    end

    % Add 2018 backfill matrices if needed
    if USE_2018_BACKFILL && ~isempty(DATE_RANGE) && year(DATE_RANGE(1)) < 2019 && ~isempty(first_2019_matrix)
        fprintf('========================================\n');
        fprintf('Adding 2018 backfill connectivity matrices\n');
        fprintf('Using Jan 1 2019 connectivity for all 2018 dates\n');
        
        % Generate dates from DATE_RANGE(1) to Dec 31 2018, every 14 days
        backfill_start = DATE_RANGE(1);
        backfill_end = datetime(2018,12,31);
        backfill_dates = backfill_start:days(14):backfill_end;
        
        for bf_idx = 1:length(backfill_dates)
            connCount = connCount + 1;
            bf_date = backfill_dates(bf_idx);
            
            conn_structs_stage1(connCount).DY = day(bf_date);
            conn_structs_stage1(connCount).MO = month(bf_date);
            conn_structs_stage1(connCount).YR = year(bf_date);
            conn_structs_stage1(connCount).Date = bf_date;
            conn_structs_stage1(connCount).full = first_2019_matrix;
            
            fprintf('  [%d] Backfilled: %s (using Jan 1 2019 connectivity)\n', ...
                    connCount, string(bf_date, 'dd-MMM-yyyy'));
        end
        
        fprintf('Backfill complete: Added %d matrices for 2018\n', length(backfill_dates));
        fprintf('========================================\n\n');
    end

    % Sort by date
    fprintf('Sorting %d matrices by date...\n', connCount);
    [~, sortIdx] = sort([conn_structs_stage1.Date]);
    conn_structs_stage1 = conn_structs_stage1(sortIdx);

    fprintf('Saving STAGE 1 connectivity data to %s...\n', stage1_cacheFilename);
    save(stage1_cacheFilePath, 'conn_structs_stage1', 'stage1_indices', '-v7.3')
    
    fprintf('\n========================================\n');
    fprintf('STAGE 1 COMPLETE: Loaded %d connectivity matrices\n', connCount);
    fprintf('Matrix dimensions: %d x %d\n', height(reefData), height(reefData));
    if ~isempty(DATE_RANGE)
        fprintf('Date range: %s to %s\n', string(DATE_RANGE(1), 'dd-MMM-yyyy'), ...
            string(DATE_RANGE(2), 'dd-MMM-yyyy'));
    end
    fprintf('Cached to: %s\n', stage1_cacheFilename);
    fprintf('========================================\n\n');
    
else
    % Load pre-calculated connectivity matrices (Stage 1)
    fprintf('========================================\n');
    fprintf('STAGE 1: Loading cached connectivity matrices from %s...\n', stage1_cacheFilename);
    
    if exist(stage1_cacheFilePath, 'file') == 0
        error(['Stage 1 cache file not found: %s\n' ...
               'Set RECALCULATE_CONNECTIVITY = true to generate it.'], stage1_cacheFilePath);
    end
    
    load(stage1_cacheFilePath, 'conn_structs_stage1', 'stage1_indices');
    
    % Reload original unfiltered data to apply stage1_indices correctly
    reefData = reefData_original(stage1_indices, :);
    
    fprintf('STAGE 1 LOADED: %d connectivity matrices from cache\n', length(conn_structs_stage1));
    fprintf('Matrix dimensions: %d x %d\n', height(reefData), height(reefData));
    if ~isempty(DATE_RANGE)
        fprintf('Date range: %s to %s\n', ...
                string(min([conn_structs_stage1.Date]), 'dd-MMM-yyyy'), string(max([conn_structs_stage1.Date]), 'dd-MMM-yyyy'));
    end
    fprintf('========================================\n\n');
end

%% STAGE 2: Apply Filter C (trajectory-based filter) and create final connectivity

% Generate final cache filename
if ~isempty(DATE_RANGE)
    final_cacheFilename = sprintf('conn_structs_FINAL_%s_to_%s', ...
                           string(DATE_RANGE(1), 'yyyyMMdd'), ...
                           string(DATE_RANGE(2), 'yyyyMMdd'));
else
    final_cacheFilename = 'conn_structs_FINAL';
end

% Add ALL filter suffixes to final cache filename
if FILTER_LOW_COVER || FILTER_MISSING_GROUPS || FILTER_NO_TRAJECTORIES
    filter_suffix = '';
    if FILTER_LOW_COVER; filter_suffix = [filter_suffix '_lowcov']; end
    if FILTER_MISSING_GROUPS; filter_suffix = [filter_suffix '_missing']; end
    if FILTER_NO_TRAJECTORIES; filter_suffix = [filter_suffix '_notraj']; end
    final_cacheFilename = [final_cacheFilename filter_suffix];
end

if USE_2018_BACKFILL && ~isempty(DATE_RANGE) && year(DATE_RANGE(1)) < 2019
    final_cacheFilename = [final_cacheFilename '_2018bf'];
end

final_cacheFilename = [final_cacheFilename '.mat'];
final_cacheFilePath = fullfile(tempPath, final_cacheFilename);

% Determine if we need to apply trajectory filter
if RECALCULATE_CONNECTIVITY && FILTER_NO_TRAJECTORIES
    % Apply trajectory filter to stage 1 data
    fprintf('\n========================================\n');
    fprintf('STAGE 2: Applying trajectory filter (Filter C)\n');
    
    num_sites_stage1 = height(reefData);
    
    % Find reefs that ever had outgoing connections
    all_sources_ever = false(num_sites_stage1, 1);
    for i = 1:length(conn_structs_stage1)
        [src, ~] = find(conn_structs_stage1(i).full);
        all_sources_ever(unique(src)) = true;
    end
    
    inactive_sources = find(~all_sources_ever);
    fprintf('Sites removed: %d (%.1f%%)\n', length(inactive_sources), 100*length(inactive_sources)/num_sites_stage1);
    fprintf('Sites remaining: %d\n', sum(all_sources_ever));
    
    if ~isempty(inactive_sources)
        fprintf('Sample inactive reef IDs: ');
        fprintf('%d ', reefData.unique_ID(inactive_sources(1:min(10, length(inactive_sources)))));
        fprintf('...\n');
    end
    
    % Create stage 2 filter mask (relative to stage 1)
    stage2_filter_mask = all_sources_ever;
    
    % Update reef data
    reefData = reefData(stage2_filter_mask, :);
    
    % Create final indices (relative to original data)
    final_indices = stage1_indices(stage2_filter_mask);
    
    % Apply trajectory filter to all connectivity matrices
    conn_structs = struct();
    for i = 1:length(conn_structs_stage1)
        conn_structs(i).DY = conn_structs_stage1(i).DY;
        conn_structs(i).MO = conn_structs_stage1(i).MO;
        conn_structs(i).YR = conn_structs_stage1(i).YR;
        conn_structs(i).Date = conn_structs_stage1(i).Date;
        conn_structs(i).full = conn_structs_stage1(i).full(stage2_filter_mask, stage2_filter_mask);
    end
    
    fprintf('========================================\n\n');
    
    % Save final cache
    fprintf('Saving FINAL connectivity data to %s...\n', final_cacheFilename);
    save(final_cacheFilePath, 'conn_structs', 'final_indices', '-v7.3')
    
    fprintf('\n========================================\n');
    fprintf('STAGE 2 COMPLETE\n');
    fprintf('Final matrix dimensions: %d x %d\n', height(reefData), height(reefData));
    fprintf('Cached to: %s\n', final_cacheFilename);
    fprintf('========================================\n\n');
    
elseif ~RECALCULATE_CONNECTIVITY
    % Try to load final cache
    fprintf('========================================\n');
    fprintf('Attempting to load FINAL cache: %s\n', final_cacheFilename);
    
    if exist(final_cacheFilePath, 'file')
        % Final cache exists - load it
        fprintf('Final cache found! Loading...\n');
        load(final_cacheFilePath, 'conn_structs', 'final_indices');
        
        % Reload and reapply ALL filters to reef data using final_indices
        reefData = reefData_original(final_indices, :);
        
        fprintf('FINAL CACHE LOADED: %d connectivity matrices\n', length(conn_structs));
        fprintf('Matrix dimensions: %d x %d\n', height(reefData), height(reefData));
        fprintf('========================================\n\n');
    else
        % Final cache doesn't exist
        fprintf('Final cache not found.\n');
        
        if FILTER_NO_TRAJECTORIES
            % Need to create final cache from stage 1
            fprintf('Trajectory filter is enabled - creating final cache from stage 1 data...\n');
            
            num_sites_stage1 = size(conn_structs_stage1(1).full, 1);
            all_sources_ever = false(num_sites_stage1, 1);
            for i = 1:length(conn_structs_stage1)
                [src, ~] = find(conn_structs_stage1(i).full);
                all_sources_ever(unique(src)) = true;
            end
            
            stage2_filter_mask = all_sources_ever;
            
            % Update reef data (currently has stage1 filtering applied)
            reefData = reefData(stage2_filter_mask, :);
            final_indices = stage1_indices(stage2_filter_mask);
            
            % Create filtered connectivity
            conn_structs = struct();
            for i = 1:length(conn_structs_stage1)
                conn_structs(i).DY = conn_structs_stage1(i).DY;
                conn_structs(i).MO = conn_structs_stage1(i).MO;
                conn_structs(i).YR = conn_structs_stage1(i).YR;
                conn_structs(i).Date = conn_structs_stage1(i).Date;
                conn_structs(i).full = conn_structs_stage1(i).full(stage2_filter_mask, stage2_filter_mask);
            end
            
            % Save final cache
            save(final_cacheFilePath, 'conn_structs', 'final_indices', '-v7.3')
            fprintf('Final cache created and saved.\n');
        else
            % Trajectory filter disabled - use stage 1 as final
            fprintf('Trajectory filter is disabled - using stage 1 as final.\n');
            conn_structs = conn_structs_stage1;
            final_indices = stage1_indices;
        end
        fprintf('========================================\n\n');
    end
else
    % RECALCULATE=true but FILTER_NO_TRAJECTORIES=false
    % Just use stage 1 as final
    conn_structs = conn_structs_stage1;
    final_indices = stage1_indices;
    fprintf('Trajectory filter disabled - using Stage 1 as final connectivity.\n\n');
end

%% Extract final reef information (after all filtering is complete)

unique_IDs = reefData.unique_ID;
locations = [reefData.centroid_lon, reefData.centroid_lat];
num_sites = height(reefData);

fprintf('\n=== FINAL FILTERING SUMMARY ===\n');
fprintf('Original sites: %d\n', height(reefData_original));
fprintf('Final sites: %d\n', num_sites);
fprintf('Total reduction: %.1f%%\n', 100 * (1 - num_sites/height(reefData_original)));
fprintf('================================\n\n');

fprintf('\n=== CONNECTIVITY DIAGNOSTICS ===\n');
first_matrix = conn_structs(1).full;
max_val = full(max(first_matrix(:)));
mean_val = full(mean(nonzeros(first_matrix)));
num_connections = nnz(first_matrix);

fprintf('Max connectivity value: %f\n', max_val);
fprintf('Mean nonzero connectivity: %f\n', mean_val);
fprintf('Number of nonzero connections: %d\n', num_connections);
fprintf('================================\n\n');


% Define day vector for SIR function (days since Jan 1, or other defined
% reference start date)
dates = [conn_structs.Date];

% Set reference date to the earliest date in the connectivity data
% This ensures conn_days starts at 0 (or close to it) rather than having negative values
ref = min(dates);  
conn_days = days(dates - ref);

%% set state variables & parameters

% 'N' is coral cover pre-SCTLD
N_LS = reefData.low_coral_cover;
N_MS = reefData.moderate_coral_cover;
N_HS = reefData.high_coral_cover;

N_site = reefData.mean_coral_cover; %same as N_LS + N_MS + N_HS

% NOTE - 9 Nov 2025
%    The below combo I think is a good starting point for beginning to
%     "slow down" the outbreak! This is using rate-based option (see ODE
%     function), and Dobbelaere-based transmission smooth
%   - include all sites; seed at 5 Flat sites; [old removal method; no 'negative' in opts]
%   - seed_frac = 0.00001;
%   - export_thresh = 0.00001;
%   - I0 = 0.3 * 0.00001; tau = I0 / 10; % ACTUALLY THIS HAD I0 OFF!
%   - flux_scale = 1;
%   - flux_shape = -3;
%
% this below one was similar maybe?
%   - include all sites; seed at 5 Flat sites; [old removal method; no 'negative' in opts]
%   - seed_frac = 0.0001;
%   - export_thresh = 0.00002;
%   - I0 = 0.3 * 0.00001; tau = I0 / 10; % ACTUALLY THIS HAD I0 OFF!
%   - flux_scale = 1;
%   - flux_shape = -3;
%
% keep an eye on this one
%   - include only sites with non-zero cover in each susc. group; seed at 5 Flat sites; [old removal method; no 'negative' in opts]
%   - seed_frac = 0.00001;
%   - export_thresh = 0.000025;
%   - I0 = 0.3 * 0.00001; tau = I0 / 10; % ACTUALLY THIS HAD I0 OFF!
%   - flux_scale = 1;
%   - flux_shape = -3;
%
% starting to get there!! WAIT hold on it seems like I0 is doing literally
% nothing
%   - include only sites with non-zero cover in each susc. group; seed at 5 Flat sites; [old removal method; no 'negative' in opts]
%   - seed_frac = 0.00001;
%   - export_thresh = 0.000025;
%   - I0 = 0.000008; tau = I0 / 10; % ACTUALLY THIS HAD I0 OFF!
%   - flux_scale = 1;
%   - flux_shape = -3;
%
% really tamped down infection here. seems to suggest negative flux_shape
% may not be a good thing
%   - include only sites with non-zero cover in each susc. group; seed at 5 Flat sites; [old removal method; no 'negative' in opts]
%   - seed_frac = 0.00001;
%   - export_thresh = 0.000025;
%   - I0 = 0.0000003; tau = I0 / 2;
%   - flux_scale = 1;
%   - flux_shape = -3;
%
% ran this for full year and still got runaway-disease. but seems promising
%   - include only sites with non-zero cover in each susc. group; seed at 5 Flat sites; [old removal method; no 'negative' in opts]
%   - seed_frac = 0.00001;
%   - export_thresh = 0.000025;
%   - I0 = 0.0000003; tau = I0 / 2;
%   - flux_scale = 1;
%   - flux_shape = 0.001;
%
% WOO!! this one seems like a winner! [but note takes WHILE to "take off"
%   - include only sites with non-zero cover in each susc. group; seed at 5 Flat sites; [old removal method; no 'negative' in opts]
%   - seed_frac = 0.00001;
%   - export_thresh = 0.000025;
%   - I0 = 0.0000008; tau = I0 / 10;
%   - flux_scale = 1;
%   - flux_shape = 0.001;
%
% interesting test here! still takes a long time to take off but maybe more
% realistic
%   - include only sites with non-zero cover in each susc. group; seed at 5 Flat sites; [old removal method; no 'negative' in opts]
%   - seed_frac = 0.00001;
%   - export_thresh = 0.000025;
%   - I0 = 0.0000007; tau = I0 / 10;
%   - flux_scale = 1;
%   - flux_shape = 1;
%
% running the above "winner" in a more general scenario is interesting.
% appears to result in roughly the correct extension to STJ by the end of
% 2019, but vastly underestimates SCTLD losses in many of the places in
% between (including STT northeast end). could try tweaking flux_shape to
% something positive? it's possible that weak connections are getting
% discarded too easily. and maybe also adjusting I0/tau
%   - include ALL sites; seed at 5 Flat sites; [old removal method; no 'negative' in opts]
%   - seed_frac = 0.00001;
%   - export_thresh = 0.000025;
%   - I0 = 0.0000008; tau = I0 / 10;
%   - flux_scale = 1;
%   - flux_shape = 0.001;
%
% % WORKING SCENARIO 2
% % okay now getting very interesting. if you seed at more sites, doesn't
% % necessarily help a ton (though I've left it here). also, seeding in
% % December doesn't seem to make a huge dent but left here. what is
% % interesting, is leaving all sites in, and I also went in and tried
% % something to make sure removal was working correctly. and made sure
% % values can't go negative in 'opts'
% %   - include ALL sites; seed at 15+ sites around flat; [newEST removal method; use 'negative' in opts]
% seed_frac = 0.00001;
% export_thresh = 0; %null condition
% I0 = 0.00000079; tau = I0 / 10; % 0.00000079 is weirdly good ????? 822 and 83 too restrictive; 82 too allowant; 8201 and 82011 even further ?????? 80011 is okay...
% flux_scale = 1; % limits max, ranges 0:1. can be used for null condition, but is the default for shaping flux too. best to leave unchanged for now
% flux_shape = 0.001; %null condition
% %flam1
%
% this is essentially just above but a stress test to show explosion to
% more final state
%   - include ALL sites; seed at 15+ sites around flat; [new removal method; use 'negative' in opts]
%   - seed_frac = 0.00001;
%   - export_thresh = 0;
%   - I0 = 0.0000008; tau = I0 / 10; % NOTE: 0.00000101; tau = I0 / 10; is
%   less explosive in the end but still too much...and starts too slow
%   - flux_scale = 1;
%   - flux_shape = 1.5;
%
% also interesting. random jump to STX, so flux prob too high, but spread
% is alright
%   - include ALL sites; seed at 15+ sites around flat; [new removal method; use 'negative' in opts]
%   - seed_frac = 0.00001;
%   - export_thresh = 0;
%   - I0 = 0.0000017; tau = I0 / 10;
%   - flux_scale = 1;
%   - flux_shape = 3;
%
% % WORKING SCENARIO 1
% % pretty damn good
% %   - include ALL sites; seed at 15+ sites around flat; [newEST removal method; use 'negative' in opts]
% seed_frac = 0.00001;
% export_thresh = 0;
% I0 = 0.0000010455; tau = I0 / 10;
% flux_scale = 1;
% flux_shape = 1.5;
%
% % WORKING SCENARIO 3 - interesting????
% %   - include ALL sites; seed at 15+ sites around flat; [newEST removal method; use 'negative' in opts]
% seed_frac = 0.00001;
% export_thresh = 0.00002;
% I0 = 0.0000006; tau = I0 / 10;
% flux_scale = 1;
% flux_shape = -1;
% %
% WORKING SCENARIO 3
%   - include ALL sites; seed at 15+ sites around flat; [newEST removal method; use 'negative' in opts]
seed_frac = 0.000008;
export_thresh = 0.00002012;
I0 = 0.000000652; tau = I0 / 10;
flux_scale = 1;
flux_shape = -1;


% define how to seed disease at Flat Cay (or wherever chosen starting
%   location is placed). this fraction is multiplied by raw coral cover of
%   each susceptibility group
% seed_frac = 0.00001; % been using this for a lot of runs. 0.0001 also seemed to work okay
% seed_frac = 0.01;
% seed_frac = 0.0001;
% seed_frac = 0.001; % I think this is too much to start with. Flat Cay blows up unrealistically
% seed_frac = 0.00001;

% threshold of (ABSOLUTE) diseased coral cover at which a site can begin
% to export disease to other sites
%   NOTE - should also consider thresholding based off of the relative
%   cover of a site / its susceptibility groups. or outbreak
%   stage/intensity
% export_thresh = 0.0003;
% export_thresh = 0; %null condition
% export_thresh = 999;
% export_thresh = 0.00001;
% export_thresh = 0.00002;
% export_thresh = 0.00005;
% export_thresh = 0.000025;
% export_thresh = 0.00002;


% Dobbelaere-style internal transmission threshold parameters
% I0 = 0.0001;   % Infection threshold: fraction of site that must be infected 
%                % before exponential within-site spread activates
%                % (0.0001 = 0.01% of colonies; close to seed_frac)
% I0 = 0.00001;
% I0 = 0.00005;
% I0 = 0.00003; %this plus tau = I0 / 10 was as if there were no suppression
% I0 = 0; %null condition
%
% tau = 0.00002; % Transition steepness: controls sharpness of threshold
%                % (smaller = more switch-like; roughly I0/5 for smooth transition)
% tau = I0 / 5; % Conservative (smooth transition). e.g., 0.00002 if I0 = 0.0001. 
% tau = I0 / 10; % Moderate (focused transition). e.g., 0.00001 if I0 = 0.0001
% tau = I0 / 50; % Sharp (near step-function). e.g., 0.000002 if I0 = 0.0001
% tau = 1e-10; %null condition

% sensitivity tests
%
% I0 = 0; tau = 1e-10; % null condition
% I0 = 0.000003; tau = I0 / 10; % too low suppression - disease plays out solely based on connectivity
% I0 = 0.000006; tau = I0 / 5; % too low suppresion - maybe slightly suppressed compared to null though ?
% I0 = 0.001; tau = I0 / 10; %0.001 go back to if need
% I0 = 0.0000008; tau = I0 / 10;
% I0 = 0.0000001; tau = I0 / 10;
% I0 = 0.0000007; tau = I0 / 10;
% I0 = 0.000000775; tau = I0 / 10;
% I0 = 0.0000009; tau = I0 / 10;
% I0 = 0.00000101; tau = I0 / 10;
% I0 = 0.0000010455; tau = I0 / 10;
% I0 = 0.0000017; tau = I0 / 10;

% reshape parameters for controlling the contribution of upstream disease
% mass to local disease pool in each patch (site)
%   flux_scale.*(1-exp(-flux_shape.*T(:)))/(1-exp(-flux_shape));
% flux_scale = 1; % limits max, ranges 0:1. can be used for null condition, but is the default for shaping flux too. best to leave unchanged for now
% flux_scale = 0; % limits max, ranges 0:1
% flux_shape = -4; % this, with everything else (I0, tau, export_thresh) null, started producing something interesting. slightly slower outbreak
% flux_shape = 0.001; %null condition
% flux_shape = 10;
% flux_shape = -3;
% flux_shape = 1.5; %null condition
% flux_shape = -1; %null condition

% pre-define vectors for initial infected and recovered (dead) coral cover
I_LS_init = zeros(num_sites,1);
I_MS_init = zeros(num_sites,1);
I_HS_init = zeros(num_sites,1);
R_LS_init = zeros(num_sites,1);
R_MS_init = zeros(num_sites,1);
R_HS_init = zeros(num_sites,1);

% seed site IDs (surrounding Flat Cay):
% target_IDs = [29088, 29338, 29089, 29339, 29087];
target_IDs = [28838, 28839, 29089, 29339, 29338, 29337, 29087, 28587, 28840, 29090, 29091, 29341, 28836, 28856, 28835, 29085, 29086, 29591, 29341, 29091, 29340, 29090];
target_reef_IDs = find(ismember(reefData.unique_ID, target_IDs));

if isempty(target_reef_IDs)
    warning('No Flat Cay sites found in filtered data! Disease will not be seeded.');
else
    fprintf('========================================\n');
    fprintf('DISEASE SEEDING\n');
    fprintf('Flat Cay sites found: %d\n', length(target_reef_IDs));
    fprintf('Site IDs: %s\n', num2str(reefData.unique_ID(target_reef_IDs)'));
    fprintf('========================================\n\n');
    
    I_HS_init(target_reef_IDs) = seed_frac * N_HS(target_reef_IDs);
    I_MS_init(target_reef_IDs) = seed_frac * N_MS(target_reef_IDs);
    I_LS_init(target_reef_IDs) = seed_frac * N_LS(target_reef_IDs);
end

% Pre-define vectors for initial susceptible coral cover
S_LS_init = N_LS - I_LS_init;
S_MS_init = N_MS - I_MS_init;
S_HS_init = N_HS - I_HS_init;

% Pack state variables for SIR function
Y0 = [S_LS_init; S_MS_init; S_HS_init; I_LS_init; I_MS_init; I_HS_init; R_LS_init; R_MS_init; R_HS_init];

% from Nearshore
% beta (transmission rate)
b_LS = 0.03;
b_MS = 0.14;
b_HS = 2.08;

% gamma (mortality rate)
g_LS = .05;
g_MS = .55;
g_HS = 3.33;

% % from Offshore
% % beta
% b_LS = 0.22;
% b_MS = 0.28;
% b_HS = 0.48;
% 
% % gamma
% g_LS = 2.73;
% g_MS = 0.33;
% g_HS = 2.62;

%% serial test model run

% opts = odeset('OutputFcn', @odeWaitbar);
% opts = odeset('OutputFcn', @odeWaitbar, ...
%               'RelTol', 1e-8, ...     % Much tighter (default 1e-3)
%               'AbsTol', 1e-10, ...    % Much tighter (default 1e-6)
%               'NonNegative', 1:length(Y0));  % CRITICAL: Force non-negative
opts = odeset('OutputFcn', @odeWaitbar, ...
              'NonNegative', 1:length(Y0));

% Calculate tspan based on DATE_RANGE (if specified) or default to 365 days
if ~isempty(DATE_RANGE)
    simulation_days = days(DATE_RANGE(2) - DATE_RANGE(1)) + 1;
    tspan = [1 simulation_days];
else
    tspan = [1 365];  % Default to 1 year if no DATE_RANGE specified
end 

% Informational: last connectivity file will be reused through end of simulation
max_conn_day = round(days(max([conn_structs.Date]) - ref));

if tspan(end) > max_conn_day
    fprintf('\n*** CONNECTIVITY INFO ***\n');
    fprintf('Simulation runs through day %d (%s)\n', ...
            tspan(end), string(ref + days(tspan(end)), 'dd-MMM-yyyy'));
    fprintf('Last connectivity file starts day %d (%s)\n', ...
            max_conn_day, string(max([conn_structs.Date]), 'dd-MMM-yyyy'));
    fprintf('This file will be reused for days %d-%d\n', max_conn_day, tspan(end));
    fprintf('************************\n\n');
end

% Store final tspan
tspan_final = tspan;

% The ODE solver. Note that ODEfun_SCTLD_seascape is a separate .m file 
fprintf('\n========================================\n');
fprintf('STARTING ODE SOLVER\n');
fprintf('Solver: ode45\n');
fprintf('Time span: [%d %d] days\n', tspan(1), tspan(end));
fprintf('Number of sites: %d\n', num_sites);
fprintf('State vector size: %d\n', length(Y0));
fprintf('========================================\n\n');

clear Y t
tic
[t,Y] = ode45(@(t,Y) ODEfun_SCTLD_seascape(t,Y,N_LS,N_MS,N_HS,b_LS,b_MS,b_HS,g_LS,g_MS,g_HS,export_thresh,flux_scale,flux_shape,num_sites,conn_structs,conn_days,I0,tau), tspan, Y0, opts);
elapsed = toc;

fprintf('\n========================================\n');
fprintf('SOLVER COMPLETE\n');
fprintf('Elapsed time: %.1f seconds\n', elapsed);
fprintf('Time points: %d\n', length(t));
fprintf('Final simulation day: %.1f\n', t(end));
fprintf('========================================\n\n');

%% Analyze SIR output

% Extract compartments from Y and interpolate to discrete days in simulation
S_LS_output = Y(:,1:num_sites);
S_MS_output = Y(:,num_sites+1:num_sites*2);
S_HS_output = Y(:,2*num_sites+1:num_sites*3);
S_LS_output_days = interp1(t,S_LS_output,tspan(1):tspan(2));
S_MS_output_days = interp1(t,S_MS_output,tspan(1):tspan(2));
S_HS_output_days = interp1(t,S_HS_output,tspan(1):tspan(2));

I_LS_output = Y(:,3*num_sites+1:num_sites*4);
I_MS_output = Y(:,4*num_sites+1:num_sites*5);
I_HS_output = Y(:,5*num_sites+1:num_sites*6);
I_LS_output_days = interp1(t,I_LS_output,tspan(1):tspan(2));
I_MS_output_days = interp1(t,I_MS_output,tspan(1):tspan(2));
I_HS_output_days = interp1(t,I_HS_output,tspan(1):tspan(2));

R_LS_output = Y(:,6*num_sites+1:num_sites*7);
R_MS_output = Y(:,7*num_sites+1:num_sites*8);
R_HS_output = Y(:,8*num_sites+1:num_sites*9);
R_LS_output_days = interp1(t,R_LS_output,tspan(1):tspan(2));
R_MS_output_days = interp1(t,R_MS_output,tspan(1):tspan(2));
R_HS_output_days = interp1(t,R_HS_output,tspan(1):tspan(2));

% Pick seed site for analysis
% example_site = flat_cay_site_IDs(1);
% example_site = target_reef_IDs(min(2, length(target_reef_IDs)));
example_site = target_reef_IDs(randi(length(target_reef_IDs)));

% Calculate totals for mass balance
total_S_temp = S_LS_output_days(:,example_site) + S_MS_output_days(:,example_site) + S_HS_output_days(:,example_site);
total_I_temp = I_LS_output_days(:,example_site) + I_MS_output_days(:,example_site) + I_HS_output_days(:,example_site);
total_R_temp = R_LS_output_days(:,example_site) + R_MS_output_days(:,example_site) + R_HS_output_days(:,example_site);
total_mass_temp = total_S_temp + total_I_temp + total_R_temp;

fprintf('\n=== SERIAL RUN DIAGNOSTICS ===\n');
fprintf('Mass balance check at site %d:\n', unique_IDs(example_site));
fprintf('  Initial total: %.8f\n', total_mass_temp(1));
fprintf('  Final total: %.8f\n', total_mass_temp(end));
fprintf('  Mass variation: %.2e (should be ~0)\n', max(total_mass_temp) - min(total_mass_temp));
fprintf('===============================\n\n');

% Diagnostic plots
figure('Name', 'Serial Run - Post-Run Diagnostics');

subplot(2,2,1);
I_total_at_seeds_temp = I_LS_output_days(:, target_reef_IDs) + I_MS_output_days(:, target_reef_IDs) + I_HS_output_days(:, target_reef_IDs);
plot(I_total_at_seeds_temp, 'LineWidth', 1.5);
xlabel('Day'); ylabel('Infected coral cover');
title('Infection trajectory at seed sites');
legend(arrayfun(@(x) sprintf('Site %d', unique_IDs(x)), target_reef_IDs, 'UniformOutput', false), 'Location', 'best');
grid on;

subplot(2,2,2);
P_at_seeds_temp = I_LS_output_days(:,target_reef_IDs) + I_MS_output_days(:,target_reef_IDs) + I_HS_output_days(:,target_reef_IDs);
plot(P_at_seeds_temp, 'LineWidth', 1.5);
xlabel('Day'); ylabel('Disease pool (P)');
title('Disease pool at seed sites');
yline(export_thresh, 'r--', 'LineWidth', 2);
text(tspan(2)*0.7, export_thresh*1.1, 'Export threshold', 'Color', 'r');
grid on;

subplot(2,2,3);
plot(total_S_temp, 'b', 'LineWidth', 1.5); hold on;
plot(total_I_temp, 'r', 'LineWidth', 1.5);
plot(total_R_temp, 'k', 'LineWidth', 1.5);
xlabel('Day'); ylabel('Cover');
title(sprintf('Total SIR - Site %d', unique_IDs(example_site)));
legend('S','I','R');
grid on;

subplot(2,2,4);
plot(total_mass_temp, 'LineWidth', 1.5);
xlabel('Day'); ylabel('S+I+R');
title('Mass balance (should be constant)');
grid on;

% Detailed SIR by group
figure('Name', 'Serial Run - Detailed SIR by Group', 'Position', [150 150 1400 800]);

subplot(2,2,1);
plot(S_LS_output_days(:,example_site), 'b', 'LineWidth', 1.5); hold on;
plot(I_LS_output_days(:,example_site), 'r', 'LineWidth', 1.5);
plot(R_LS_output_days(:,example_site), 'k', 'LineWidth', 1.5);
xlabel('Day'); ylabel('Cover');
title(sprintf('LS Group - Site %d', unique_IDs(example_site)));
legend('S','I','R');
grid on;

subplot(2,2,2);
plot(S_MS_output_days(:,example_site), 'b', 'LineWidth', 1.5); hold on;
plot(I_MS_output_days(:,example_site), 'r', 'LineWidth', 1.5);
plot(R_MS_output_days(:,example_site), 'k', 'LineWidth', 1.5);
xlabel('Day'); ylabel('Cover');
title(sprintf('MS Group - Site %d', unique_IDs(example_site)));
legend('S','I','R');
grid on;

subplot(2,2,3);
plot(S_HS_output_days(:,example_site), 'b', 'LineWidth', 1.5); hold on;
plot(I_HS_output_days(:,example_site), 'r', 'LineWidth', 1.5);
plot(R_HS_output_days(:,example_site), 'k', 'LineWidth', 1.5);
xlabel('Day'); ylabel('Cover');
title(sprintf('HS Group - Site %d', unique_IDs(example_site)));
legend('S','I','R');
grid on;

subplot(2,2,4);
plot(total_S_temp, 'b', 'LineWidth', 1.5); hold on;
plot(total_I_temp, 'r', 'LineWidth', 1.5);
plot(total_R_temp, 'k', 'LineWidth', 1.5);
xlabel('Day'); ylabel('Cover');
title(sprintf('Total - Site %d', unique_IDs(example_site)));
legend('S','I','R');
grid on;

%% Animation of outbreak

fprintf('Creating diagnostic movie...\n');

% Extract infected, susceptible, and removed populations over time
S_total_output_days = S_LS_output_days + S_MS_output_days + S_HS_output_days;  % Total susceptible  
I_total_output_days = I_LS_output_days + I_MS_output_days + I_HS_output_days;  % Total infected
R_total_output_days = R_LS_output_days + R_MS_output_days + R_HS_output_days;  % Total removed

% Create colormaps
breakpoints = [0 .000000001 .001 .1 1]; %very very low coral cover is going to show up visually this way
Cstart = [
    0.25 .5 0.25;      % light green
    1 1 .2;            % green
    1.00 0.60 0.00;    % yellow
    0.80 0.00 0.80;    % orange
];
Cend = [
    0.00 0.60 0.00;    % green
    1.00 0.60 0.00;    % yellow
    1.00 0.00 0.00;    % orange
    0 0 0;             % red
];
cmap_I = stackedColormap(breakpoints, Cstart, Cend, 256, [1 1 1 1]);
cmap_R = stackedColormap(breakpoints, Cstart, Cend, 256, [1 1 1 1]);

% Determine how many days actually simulated
num_days = size(I_total_output_days, 1);
fprintf('Creating animation for %d days of simulation\n', num_days);

% % Create figure (no renderer specification - MATLAB handles automatically)
% f = figure('Position', [10 10 1400 1000]);
% set(f, 'nextplot', 'replacechildren'); 
% 
% % Create 2x2 layout
% T = tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');
% 
% % Initialize plots
% %
% % absolute infected coral cover
% axv1 = nexttile(T,1);
%     h1 = scatter(axv1, locations(:,1), locations(:,2), 7, I_total_output_days(1,:)', 'filled');
%     colormap(axv1, cmap_I);
%     clim([0 .01]);
%     c1 = colorbar(axv1);
%     axis equal
% 
% % relative infected coral cover
% axv2 = nexttile(T,2);
%     h2 = scatter(axv2, locations(:,1), locations(:,2), 7, I_total_output_days(1,:)'./N_site(:), 'filled');
%     colormap(axv2, cmap_I);
%     clim([0 .01]);
%     c2 = colorbar(axv2);
%     axis equal 
% 
% % absolute removed coral cover
% axv3 = nexttile(T,3);
%     h3 = scatter(axv3, locations(:,1), locations(:,2), 7, N_site(:)-S_total_output_days(1,:)', 'filled');
%     c3 = colorbar(axv3);
%     clim([0 1]);
%     colormap(axv3, cmap_R)
%     axis equal
% 
% % relative removed coral cover
% axv4 = nexttile(T,4);
%     h4 = scatter(axv4, locations(:,1), locations(:,2), 7, (N_site(:)-S_total_output_days(1,:)')./N_site(:), 'filled');
%     c4 = colorbar(axv4);
%     clim([0 1]);
%     colormap(axv4, cmap_R)
%     axis equal
% 
% % Create video writer with MPEG-4 format
% vidname = sprintf('DisVid_Diagnostic_ET%s_I0_%s_FS%.1f_FP%.2f', ...
%                   strrep(num2str(export_thresh, '%.10f'), '.', 'p'), ...
%                   strrep(num2str(I0, '%.10f'), '.', 'p'), ...
%                   flux_scale, ...
%                   flux_shape);
% v = VideoWriter(fullfile(seascapePath, vidname), 'MPEG-4');
% v.Quality = 95;  % 0-100, higher = better quality but larger file
% v.FrameRate = 5;
% open(v);
% 
% % Threshold for disease front boundary
% %   - this is saying that only once 0.1% cover within a site (patch) has
% %       been removed, is that site considered "observable" to, for example,
% %       the naked eye of strike team divers in the water within that ~0.4
% %       m2 grid square
% removed_cover_thresh = 0.001; % in proportion 0 to 1 (not 0-100%)
% 
% % Generate animation frames
% for k = 1:1:num_days
%     % Calculate date based on reference date and actual simulation start
%     Dt = ref + days(tspan(1) + k - 2);  % Uses ref from conn_days calculation
% 
%     % Find sites with dead coral above threshold (disease front)
%     observable_sick_sites = R_total_output_days(k,:) >= removed_cover_thresh;
%     if sum(observable_sick_sites) >= 3  % Need at least 3 points for boundary
%         locs_sick_sites = locations(observable_sick_sites,:);
%         try
%             bounds_sick_sites = boundary(locs_sick_sites(:,1), locs_sick_sites(:,2), 0.7);
%             draw_boundary = true;
%         catch
%             draw_boundary = false;  % Not enough points or other issue
%         end
%     else
%         draw_boundary = false;
%     end
% 
%     % Update plot 1: Disease prevalence (absolute % cover of site which is infected)
%     set(h1, 'CData', I_total_output_days(k,:)');
%     if draw_boundary
%         p1 = patch(axv1, locs_sick_sites(bounds_sick_sites,1), locs_sick_sites(bounds_sick_sites,2), [.8 .9 1], ...
%                    'EdgeColor', 'r', 'FaceAlpha', 0.2);
%     end
%     %
%     title(axv1, 'Infected cover of seafloor (%)', ...
%           sprintf('t = %s', string(Dt, 'dd-MMM-yyyy')));
%     colormap(axv1, cmap_I);
%     clim([0 .01]);
%     axis equal
% 
%     % Update plot 2: Disease prevalence (relative % cover of site which is infected)
%     set(h2, 'CData', I_total_output_days(k,:)'./N_site(:));
%     if draw_boundary
%         p2 = patch(axv2, locs_sick_sites(bounds_sick_sites,1), locs_sick_sites(bounds_sick_sites,2), [.8 .9 1], ...
%                    'EdgeColor', 'r', 'FaceAlpha', 0.2);
%     end
%     %
%     title(axv2, 'Relative cover of infected coral (%)', ...
%           sprintf('t = %s', string(Dt, 'dd-MMM-yyyy')));
%     colormap(axv2, cmap_I);
%     clim([0 .1]);
%     axis equal
% 
%     % Update plot 3: Total coral cover lost
%     set(h3, 'CData', (N_site(:)-S_total_output_days(k,:)'));
%     if draw_boundary
%         p3 = patch(axv3, locs_sick_sites(bounds_sick_sites,1), locs_sick_sites(bounds_sick_sites,2), [.8 .9 1], ...
%                    'EdgeColor', 'r', 'FaceAlpha', 0.2);
%     end
%     %
%     title(axv3, 'Removed cover of seafloor (%)', ...
%           sprintf('t = %s', string(Dt, 'dd-MMM-yyyy')));
%     clim([0 1]);
%     colormap(axv3, cmap_R)
%     axis equal
% 
%     % Update plot 4: Proportion coral cover lost
%     set(h4, 'CData', (N_site(:)-S_total_output_days(k,:)')./N_site(:));
%     if draw_boundary
%         p4 = patch(axv4, locs_sick_sites(bounds_sick_sites,1), locs_sick_sites(bounds_sick_sites,2), [.8 .9 1], ...
%                    'EdgeColor', 'r', 'FaceAlpha', 0.2);
%     end
%     %
%     title(axv4, 'Relative cover of removed coral (%)', ...
%           sprintf('t = %s', string(Dt, 'dd-MMM-yyyy')));
%     clim([0 1]);
%     colormap(axv4, cmap_R)
%     axis equal
% 
%     % Capture frame and write to video
%     F = getframe(f);
%     writeVideo(v, F);
% 
%     % Clean up patches for next frame
%     if draw_boundary
%         delete(p1)
%         delete(p2)
%         delete(p3)
%         delete(p4)
%     end
% 
%     % Progress indicator every 10 days
%     if mod(k, 10) == 0
%         fprintf('  Frame %d/%d\n', k, num_days);
%     end
% end
% 
% close(v);
% fprintf('Animation saved as: %s.mp4\n', vidname);

%% Static figure of final state

fprintf('Creating static figure of final state...\n');

% Find last valid timepoint (no NaN values)
valid_times = find(all(~isnan(I_total_output_days), 2));
if isempty(valid_times)
    warning('No valid timepoints found!');
else
    last_valid = valid_times(end);
    final_day = tspan(1) + last_valid - 1;
    final_date = datetime(2019,1,1) + (final_day - 1);
    
    fprintf('Final valid timepoint: Day %d (%s)\n', final_day, string(final_date, 'dd-MMM-yyyy'));
    
    % Create static figure
    fig_final = figure('Position', [50 50 1400 1000]);
    T_final = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    title(T_final, sprintf('Final Disease State - Day %d (%s)', final_day, string(final_date, 'dd-MMM-yyyy')), ...
          'FontSize', 16, 'FontWeight', 'bold');
    
    % Threshold for disease front boundary
    removed_cover_thresh = 0.001;
    observable_sick_sites = R_total_output_days(last_valid,:) >= removed_cover_thresh;
    
    % Try to compute boundary
    if sum(observable_sick_sites) >= 3
        locs_sick_sites = locations(observable_sick_sites,:);
        try
            bounds_sick_sites = boundary(locs_sick_sites(:,1), locs_sick_sites(:,2), 0.7);
            has_boundary = true;
        catch
            has_boundary = false;
        end
    else
        has_boundary = false;
    end
    
    % Plot 1: Disease prevalence - proportion of bottom
    nexttile(T_final, 1);
    scatter(locations(:,1), locations(:,2), 7, I_total_output_days(last_valid,:)', 'filled');
    if has_boundary
        hold on;
        patch(locs_sick_sites(bounds_sick_sites,1), locs_sick_sites(bounds_sick_sites,2), [.8 .9 1], ...
              'EdgeColor', 'r', 'FaceAlpha', 0.2, 'LineWidth', 2);
        hold off;
    end
    colormap(gca, cmap_I);
    clim([0 .01]);
    colorbar;
    title('Disease prevalence - proportion of bottom');
    axis equal;
    
    % Plot 2: Disease prevalence - proportion of living coral
    nexttile(T_final, 2);
    scatter(locations(:,1), locations(:,2), 7, I_total_output_days(last_valid,:)'./N_site(:), 'filled');
    if has_boundary
        hold on;
        patch(locs_sick_sites(bounds_sick_sites,1), locs_sick_sites(bounds_sick_sites,2), [.8 .9 1], ...
              'EdgeColor', 'r', 'FaceAlpha', 0.2, 'LineWidth', 2);
        hold off;
    end
    colormap(gca, cmap_I);
    clim([0 .1]);
    colorbar;
    title('Disease prevalence - proportion of living coral');
    axis equal;
    
    % Plot 3: Total coral cover lost
    nexttile(T_final, 3);
    scatter(locations(:,1), locations(:,2), 7, N_site(:)-S_total_output_days(last_valid,:)', 'filled');
    if has_boundary
        hold on;
        patch(locs_sick_sites(bounds_sick_sites,1), locs_sick_sites(bounds_sick_sites,2), [.8 .9 1], ...
              'EdgeColor', 'r', 'FaceAlpha', 0.2, 'LineWidth', 2);
        hold off;
    end
    colormap(gca, cmap_R);
    clim([0 1]);
    colorbar;
    title('Total coral cover lost');
    axis equal;
    
    % Plot 4: Proportion coral cover lost
    nexttile(T_final, 4);
    scatter(locations(:,1), locations(:,2), 7, (N_site(:)-S_total_output_days(last_valid,:)')./N_site(:), 'filled');
    if has_boundary
        hold on;
        patch(locs_sick_sites(bounds_sick_sites,1), locs_sick_sites(bounds_sick_sites,2), [.8 .9 1], ...
              'EdgeColor', 'r', 'FaceAlpha', 0.2, 'LineWidth', 2);
        hold off;
    end
    colormap(gca, cmap_R);
    clim([0 1]);
    colorbar;
    title('Proportion coral cover lost');
    axis equal;
    
    % Save figure
    figname = sprintf('FinalState_day%d_ET%s_I0_%s_FS%.1f_FP%.2f', ...
                  final_day, ...
                  strrep(num2str(export_thresh, '%.10f'), '.', 'p'), ...
                  strrep(num2str(I0, '%.10f'), '.', 'p'), ...
                  flux_scale, ...
                  flux_shape);    saveas(fig_final, fullfile(seascapePath, [figname '.png']));
    saveas(fig_final, fullfile(seascapePath, [figname '.fig']));
    
    fprintf('Static figure saved as: %s\n', figname);
    
    % Summary statistics
    fprintf('\n=== FINAL STATE SUMMARY ===\n');
    fprintf('Sites with any infection: %d (%.1f%%)\n', ...
            nnz(I_total_output_days(last_valid,:) > 0), 100 * nnz(I_total_output_days(last_valid, :) > 0) / num_sites);
    fprintf('Sites with >0.1%% cover lost: %d (%.1f%%)\n', ...
            nnz(R_total_output_days(last_valid, :) > 0.001), 100 * nnz(R_total_output_days(last_valid, :) > 0.001) / num_sites);
    fprintf('Total coral cover lost: %.4f (%.2f%% of initial)\n', ...
            sum(R_total_output_days(last_valid, :)), 100 * sum(R_total_output_days(last_valid, :)) / sum(N_site));
    fprintf('Max site-level cover lost: %.4f\n', max(R_total_output_days(last_valid, :)));
    fprintf('Mean cover lost (diseased sites only): %.4f\n', ...
            mean(R_total_output_days(last_valid, R_total_output_days(last_valid, :) > 0)));
    fprintf('===========================\n\n');
end


%flam3





% %% Save whole workspace
% 
% fprintf('\n========================================\n');
% fprintf('SAVING WORKSPACE\n');
% fprintf('========================================\n');
% 
% workspace_filename = fullfile(seascapePath, 'seascape_SIR_workspace.mat');
% 
% fprintf('Saving complete workspace to:\n  %s\n', workspace_filename);
% fprintf('This may take a moment for large datasets...\n');
% 
% tic;
% save(workspace_filename, '-v7.3');
% elapsed_save = toc;
% 
% % Get file size for confirmation
% file_info = dir(workspace_filename);
% file_size_mb = file_info.bytes / (1024^2);
% 
% fprintf('Workspace saved successfully!\n');
% fprintf('  File size: %.1f MB\n', file_size_mb);
% fprintf('  Save time: %.1f seconds\n', elapsed_save);
% fprintf('========================================\n\n');



%% Save essential outputs for plotting

fprintf('\n========================================\n');
fprintf('SAVING ESSENTIAL OUTPUTS\n');
fprintf('========================================\n');

% Create structure with only what's needed for plotting
outputs = struct();

% Simulation metadata
outputs.metadata.num_sites = num_sites;
outputs.metadata.num_days = num_days;
outputs.metadata.tspan = tspan_final;
outputs.metadata.ref_date = ref;
outputs.metadata.simulation_dates = ref + days(tspan_final(1):tspan_final(2));

% Site information
outputs.sites.unique_IDs = unique_IDs;
outputs.sites.locations = locations;
outputs.sites.N_site = N_site;
outputs.sites.target_reef_IDs = target_reef_IDs;

% Parameters used
outputs.params.seed_frac = seed_frac;
outputs.params.export_thresh = export_thresh;
outputs.params.I0 = I0;
outputs.params.tau = tau;
outputs.params.flux_scale = flux_scale;
outputs.params.flux_shape = flux_shape;
outputs.params.b_LS = b_LS;
outputs.params.b_MS = b_MS;
outputs.params.b_HS = b_HS;
outputs.params.g_LS = g_LS;
outputs.params.g_MS = g_MS;
outputs.params.g_HS = g_HS;

% SIR outputs (interpolated to daily values)
outputs.SIR.S_LS = S_LS_output_days;
outputs.SIR.S_MS = S_MS_output_days;
outputs.SIR.S_HS = S_HS_output_days;
outputs.SIR.I_LS = I_LS_output_days;
outputs.SIR.I_MS = I_MS_output_days;
outputs.SIR.I_HS = I_HS_output_days;
outputs.SIR.R_LS = R_LS_output_days;
outputs.SIR.R_MS = R_MS_output_days;
outputs.SIR.R_HS = R_HS_output_days;

% Computed totals (for convenience)
outputs.totals.S_total = S_total_output_days;
outputs.totals.I_total = I_total_output_days;
outputs.totals.R_total = R_total_output_days;

% Filtering info (for reproducibility)
outputs.filtering.FILTER_LOW_COVER = FILTER_LOW_COVER;
outputs.filtering.FILTER_MISSING_GROUPS = FILTER_MISSING_GROUPS;
outputs.filtering.FILTER_NO_TRAJECTORIES = FILTER_NO_TRAJECTORIES;
outputs.filtering.final_indices = final_indices;

workspace_filename = fullfile(seascapePath, 'seascape_SIR_workspace.mat');

fprintf('Saving outputs to:\n  %s\n', workspace_filename);

tic;
save(workspace_filename, 'outputs', '-v7.3');
elapsed_save = toc;

% Get file size
file_info = dir(workspace_filename);
file_size_mb = file_info.bytes / (1024^2);

fprintf('Outputs saved successfully!\n');
fprintf('  File size: %.1f MB\n', file_size_mb);
fprintf('  Save time: %.1f seconds\n', elapsed_save);
fprintf('========================================\n\n');