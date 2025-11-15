clear; clc

%% setup
projectPath = matlab.project.rootProject().RootFolder;
tempPath = fullfile(projectPath, 'temp');
dataPath = fullfile(projectPath, 'data');
outputPath = fullfile(projectPath, 'output');
seascapePath = fullfile(outputPath, 'seascape_SIR');

RECALCULATE_CONNECTIVITY = true;  % Set to false to load from temp/P.mat
RUN_PARAMETER_SWEEP = true;        % Set to false to load from temp/ParameterSweep_Results.mat
USE_PARALLEL = true;               % Set to true to use parallel computing (requires Parallel Computing Toolbox)

% ===== NEW FILTERING TOGGLES =====
FILTER_LOW_COVER = false;           % Remove sites with <1% total coral cover
FILTER_MISSING_GROUPS = true;      % Remove sites with 0% in any susceptibility group
% ==================================

% ===== 2018 CONNECTIVITY BACKFILL =====
USE_2018_BACKFILL = false;           % If true, use Jan 2019 connectivity for any dates before 2019
% =======================================

% DATE_RANGE = [];  % Empty = create/load all connectivity matrices
DATE_RANGE = [datetime(2019,1,1), datetime(2019,3,31)];  % just Q1
% DATE_RANGE = [datetime(2019,1,1), datetime(2019,6,30)];  % Q1-Q2
% DATE_RANGE = [datetime(2019,1,1), datetime(2019,12,31)];  % Q1-Q4
% DATE_RANGE = [datetime(2018,12,1), datetime(2019,12,31)];  % 2018-2019 example

%% Load reef data from CSV

reefDataFile = fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv');
reefData = readtable(reefDataFile);

fprintf('========================================\n');
fprintf('INITIAL REEF DATA\n');
fprintf('Total sites loaded: %d\n', height(reefData));
fprintf('========================================\n\n');

%% Apply site filtering

% Store original indices before filtering
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
    original_indices = reefData.original_index(sites_to_keep);
    reefData = reefData(sites_to_keep, :);
    
    fprintf('========================================\n');
    fprintf('FILTERING SUMMARY\n');
    fprintf('Original sites: %d\n', length(sites_to_keep));
    fprintf('Filtered sites: %d\n', height(reefData));
    fprintf('Reduction: %.1f%%\n', 100 * (1 - height(reefData)/length(sites_to_keep)));
    fprintf('========================================\n\n');
else
    original_indices = (1:height(reefData))';
end

% Extract reef information
unique_IDs = reefData.unique_ID;
locations = [reefData.centroid_lon, reefData.centroid_lat];
num_sites = height(reefData);

%% Create or load existing connectivity matrices

% Generate connectivity (conn_structs) cache filename based on date range
if ~isempty(DATE_RANGE)
    cacheFilename = sprintf('conn_structs_%s_to_%s', ...
                           string(DATE_RANGE(1), 'yyyyMMdd'), ...
                           string(DATE_RANGE(2), 'yyyyMMdd'));
else
    cacheFilename = 'conn_structs';  % Full dataset
end

% Add filtering suffix to cache filename
if FILTER_LOW_COVER || FILTER_MISSING_GROUPS
    filter_suffix = '';
    if FILTER_LOW_COVER; filter_suffix = [filter_suffix '_lowcov']; end
    if FILTER_MISSING_GROUPS; filter_suffix = [filter_suffix '_missing']; end
    cacheFilename = [cacheFilename filter_suffix];
end

% Add 2018 backfill suffix to cache filename
if USE_2018_BACKFILL && ~isempty(DATE_RANGE) && year(DATE_RANGE(1)) < 2019
    cacheFilename = [cacheFilename '_2018bf'];
end

cacheFilename = [cacheFilename '.mat'];
cacheFilePath = fullfile(tempPath, cacheFilename);

if RECALCULATE_CONNECTIVITY
    fprintf('========================================\n');
    fprintf('Loading connectivity matrices from disk\n');
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
    conn_structs = struct();
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
            
            conn_structs(connCount).DY = day(Dat);
            conn_structs(connCount).MO = month(Dat);
            conn_structs(connCount).YR = year(Dat);
            conn_structs(connCount).Date = Dat;
            
            conmat = sparse(Con.connectivity_results.ConnMatrix_raw);
            
            % ===== APPLY SITE FILTERING TO CONNECTIVITY =====
            if FILTER_LOW_COVER || FILTER_MISSING_GROUPS
                conmat = conmat(original_indices, original_indices);
            end
            % ================================================
            
            % Remove self retention (diagonal)
            conmat = spdiags(zeros(size(conmat,1),1), 0, conmat);
            
            % Normalize
            conmat = conmat ./ norm_val;
            
            conn_structs(connCount).full = conmat;
            
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
            
            conn_structs(connCount).DY = day(bf_date);
            conn_structs(connCount).MO = month(bf_date);
            conn_structs(connCount).YR = year(bf_date);
            conn_structs(connCount).Date = bf_date;
            conn_structs(connCount).full = first_2019_matrix;
            
            fprintf('  [%d] Backfilled: %s (using Jan 1 2019 connectivity)\n', ...
                    connCount, string(bf_date, 'dd-MMM-yyyy'));
        end
        
        fprintf('Backfill complete: Added %d matrices for 2018\n', length(backfill_dates));
        fprintf('========================================\n\n');
    end

    % Sort by date
    fprintf('Sorting %d matrices by date...\n', connCount);
    [~, sortIdx] = sort([conn_structs.Date]);
    conn_structs = conn_structs(sortIdx);

    fprintf('Saving processed connectivity data to %s...\n', cacheFilename);
    save(cacheFilePath, 'conn_structs', 'original_indices', '-v7.3')
    
    fprintf('\n========================================\n');
    fprintf('COMPLETE: Loaded %d connectivity matrices\n', connCount);
    fprintf('Matrix dimensions: %d x %d\n', num_sites, num_sites);
    if ~isempty(DATE_RANGE)
        fprintf('Date range: %s to %s\n', string(DATE_RANGE(1), 'dd-MMM-yyyy'), ...
            string(DATE_RANGE(2), 'dd-MMM-yyyy'));
    end
    fprintf('Cached to: %s\n', cacheFilename);
    fprintf('========================================\n\n');
    
else
    % Load pre-calculated connectivity matrices
    fprintf('========================================\n');
    fprintf('Loading cached connectivity matrices from %s...\n', cacheFilename);
    
    if exist(cacheFilePath, 'file') == 0
        error(['Cache file not found: %s\n' ...
               'Set RECALCULATE_CONNECTIVITY = true to generate it.'], cacheFilePath);
    end
    
    load(cacheFilePath, 'conn_structs', 'original_indices');
    fprintf('COMPLETE: Loaded %d connectivity matrices from cache\n', length(conn_structs));
    fprintf('Matrix dimensions: %d x %d\n', num_sites, num_sites);
    if ~isempty(DATE_RANGE)
        fprintf('Date range: %s to %s\n', ...
                string(min([conn_structs.Date]), 'dd-MMM-yyyy'), string(max([conn_structs.Date]), 'dd-MMM-yyyy'));
    end
    fprintf('========================================\n\n');
end

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
ref = datetime(2019,1,1);  
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
%   - include all sites; seed at 5 Flat sites
%   - seed_frac = 0.00001;
%   - export_thresh = 0.00001;
%   - I0 = 0.3 * 0.00001; tau = I0 / 10; % ACTUALLY THIS HAD I0 OFF!
%   - flux_scale = 1;
%   - flux_shape = -3;
%
% this below one was similar maybe?
%   - include all sites; seed at 5 Flat sites
%   - seed_frac = 0.0001;
%   - export_thresh = 0.00002;
%   - I0 = 0.3 * 0.00001; tau = I0 / 10; % ACTUALLY THIS HAD I0 OFF!
%   - flux_scale = 1;
%   - flux_shape = -3;
%
% keep an eye on this one
%   - include only sites with non-zero cover in each susc. group; seed at 5 Flat sites
%   - seed_frac = 0.00001;
%   - export_thresh = 0.000025;
%   - I0 = 0.3 * 0.00001; tau = I0 / 10; % ACTUALLY THIS HAD I0 OFF!
%   - flux_scale = 1;
%   - flux_shape = -3;
%
% starting to get there!! WAIT hold on it seems like I0 is doing literally
% nothing
%   - include only sites with non-zero cover in each susc. group; seed at 5 Flat sites
%   - seed_frac = 0.00001;
%   - export_thresh = 0.000025;
%   - I0 = 0.000008; tau = I0 / 10; % ACTUALLY THIS HAD I0 OFF!
%   - flux_scale = 1;
%   - flux_shape = -3;
%
% really tamped down infection here. seems to suggest negative flux_shape
% may not be a good thing
%   - include only sites with non-zero cover in each susc. group; seed at 5 Flat sites
%   - seed_frac = 0.00001;
%   - export_thresh = 0.000025;
%   - I0 = 0.0000003; tau = I0 / 2;
%   - flux_scale = 1;
%   - flux_shape = -3;
%
% ran this for full year and still got runaway-disease. but seems promising
%   - include only sites with non-zero cover in each susc. group; seed at 5 Flat sites
%   - seed_frac = 0.00001;
%   - export_thresh = 0.000025;
%   - I0 = 0.0000003; tau = I0 / 2;
%   - flux_scale = 1;
%   - flux_shape = 0.001;
%
% WOO!! this one seems like a winner! [but note takes WHILE to "take off"
%   - include only sites with non-zero cover in each susc. group; seed at 5 Flat sites
%   - seed_frac = 0.00001;
%   - export_thresh = 0.000025;
%   - I0 = 0.0000008; tau = I0 / 10;
%   - flux_scale = 1;
%   - flux_shape = 0.001;
%
% interesting test here! still takes a long time to take off but maybe more
% realistic
%   - include only sites with non-zero cover in each susc. group; seed at 5 Flat sites
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
%   - include ALL sites; seed at 5 Flat sites
%   - seed_frac = 0.00001;
%   - export_thresh = 0.000025;
%   - I0 = 0.0000008; tau = I0 / 10;
%   - flux_scale = 1;
%   - flux_shape = 0.001;




% define how to seed disease at Flat Cay (or wherever chosen starting
%   location is placed). this fraction is multiplied by raw coral cover of
%   each susceptibility group
% seed_frac = 0.00001; % been using this for a lot of runs. 0.0001 also seemed to work okay
% seed_frac = 0.01;
% seed_frac = 0.0001;
% seed_frac = 0.001; % I think this is too much to start with. Flat Cay blows up unrealistically
seed_frac = 0.00001;

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
export_thresh = 0.000025;

%flam1

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
I0 = 0.000000775; tau = I0 / 10;

% reshape parameters for controlling the contribution of upstream disease
% mass to local disease pool in each patch (site)
%   flux_scale.*(1-exp(-flux_shape.*T(:)))/(1-exp(-flux_shape));
flux_scale = 1; % limits max, ranges 0:1. can be used for null condition, but is the default for shaping flux too. best to leave unchanged for now
% flux_scale = 0; % limits max, ranges 0:1
% flux_shape = -4; % this, with everything else (I0, tau, export_thresh) null, started producing something interesting. slightly slower outbreak
flux_shape = 0.001; %null condition
% flux_shape = 10;
% flux_shape = -3;
% flux_shape = 1; %null condition

% pre-define vectors for initial infected and recovered (dead) coral cover
I_LS_init = zeros(num_sites,1);
I_MS_init = zeros(num_sites,1);
I_HS_init = zeros(num_sites,1);
R_LS_init = zeros(num_sites,1);
R_MS_init = zeros(num_sites,1);
R_HS_init = zeros(num_sites,1);

% Flat Cay site IDs (29088 is the preferred/primary location):
target_IDs = [29088, 29338, 29089, 29339, 29087];
flat_cay_site_IDs = find(ismember(reefData.unique_ID, target_IDs));

if isempty(flat_cay_site_IDs)
    warning('No Flat Cay sites found in filtered data! Disease will not be seeded.');
else
    fprintf('========================================\n');
    fprintf('DISEASE SEEDING\n');
    fprintf('Flat Cay sites found: %d\n', length(flat_cay_site_IDs));
    fprintf('Site IDs: %s\n', num2str(reefData.unique_ID(flat_cay_site_IDs)'));
    fprintf('========================================\n\n');
    
    I_HS_init(flat_cay_site_IDs) = seed_frac * N_HS(flat_cay_site_IDs);
    I_MS_init(flat_cay_site_IDs) = seed_frac * N_MS(flat_cay_site_IDs);
    I_LS_init(flat_cay_site_IDs) = seed_frac * N_LS(flat_cay_site_IDs);
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

opts = odeset('OutputFcn', @odeWaitbar);
tspan = [1 365]; 

% Informational: last connectivity file will be reused through end of simulation
max_conn_day = round(days(max([conn_structs.Date]) - datetime(2019,1,1)));

if tspan(end) > max_conn_day
    fprintf('\n*** CONNECTIVITY INFO ***\n');
    fprintf('Simulation runs through day %d (%s)\n', ...
            tspan(end), string(datetime(2019,1,1) + tspan(end), 'dd-MMM-yyyy'));
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
% [t,Y] = ode45(@(t,Y) ODEfun_SCTLD_seascape(t,Y,N_LS,N_MS,N_HS,b_LS,b_MS,b_HS,g_LS,g_MS,g_HS,export_thresh,flux_scale,flux_shape,num_sites,conn_structs,conn_days), tspan, Y0, opts);
[t,Y] = ode45(@(t,Y) ODEfun_SCTLD_seascape(t,Y,N_LS,N_MS,N_HS,b_LS,b_MS,b_HS,g_LS,g_MS,g_HS,export_thresh,flux_scale,flux_shape,num_sites,conn_structs,conn_days,I0,tau), tspan, Y0, opts);elapsed = toc;

fprintf('\n========================================\n');
fprintf('SOLVER COMPLETE\n');
fprintf('Elapsed time: %.1f seconds\n', elapsed);
fprintf('Time points: %d\n', length(t));
fprintf('Final simulation day: %.1f\n', t(end));
fprintf('========================================\n\n');

%% POST-SERIAL-RUN DIAGNOSTICS (before parameter sweep)


% Add this right after your ODE run, before any plotting:

fprintf('\n=== INCOMING FLUX DIAGNOSTICS ===\n');

% Calculate initial disease pool at seed sites
P_initial = I_LS_init + I_MS_init + I_HS_init;

% Get initial incoming flux at all sites using first connectivity matrix
initial_flux = incomingRisk_sparse(conn_structs(1).full, P_initial, export_thresh, 'rowsAreSources');

% After reshaping (like in your ODE function)
initial_flux_reshaped = flux_scale .* (1 - exp(-flux_shape .* initial_flux)) / (1 - exp(-flux_shape));

% Check at seed sites
fprintf('Seed sites (IDs: %s):\n', mat2str(unique_IDs(flat_cay_site_IDs)));
fprintf('  Raw incoming flux range: [%.6f, %.6f]\n', ...
    min(initial_flux(flat_cay_site_IDs)), max(initial_flux(flat_cay_site_IDs)));
fprintf('  Reshaped incoming flux range: [%.6f, %.6f]\n', ...
    min(initial_flux_reshaped(flat_cay_site_IDs)), max(initial_flux_reshaped(flat_cay_site_IDs)));

% Check at nearby non-seed sites (within 5 km of seed sites)
seed_locs = locations(flat_cay_site_IDs, :);
distances = pdist2(locations, seed_locs);
min_dist = min(distances, [], 2);
nearby_sites = find(min_dist < 5000 & min_dist > 0);  % 5 km radius, excluding seed sites

if ~isempty(nearby_sites)
    fprintf('\nNearby sites (within 5km, n=%d):\n', length(nearby_sites));
    fprintf('  Raw incoming flux range: [%.6f, %.6f]\n', ...
        min(initial_flux(nearby_sites)), max(initial_flux(nearby_sites)));
    fprintf('  Reshaped incoming flux range: [%.6f, %.6f]\n', ...
        min(initial_flux_reshaped(nearby_sites)), max(initial_flux_reshaped(nearby_sites)));
    fprintf('  Mean reshaped flux: %.6f\n', mean(initial_flux_reshaped(nearby_sites)));
end

% Also check during simulation - extract at day 30
day30_idx = find(t >= 30, 1, 'first');
if ~isempty(day30_idx)
    P_day30 = Y(day30_idx, 3*num_sites+1:4*num_sites) + ...
              Y(day30_idx, 4*num_sites+1:5*num_sites) + ...
              Y(day30_idx, 5*num_sites+1:6*num_sites);
    
    % Find which connectivity matrix was used at day 30
    conn_idx_day30 = max(1, find(conn_days <= 30, 1, 'last'));
    flux_day30 = incomingRisk_sparse(conn_structs(conn_idx_day30).full, P_day30', export_thresh, 'rowsAreSources');
    flux_day30_reshaped = flux_scale .* (1 - exp(-flux_shape .* flux_day30)) / (1 - exp(-flux_shape));
    
    fprintf('\nDay 30 (after some spread):\n');
    fprintf('  Sites with P > 0: %d\n', sum(P_day30 > 0));
    fprintf('  Max reshaped flux anywhere: %.6f\n', max(flux_day30_reshaped));
    fprintf('  Mean reshaped flux (sites with flux>0): %.6f\n', mean(flux_day30_reshaped(flux_day30_reshaped > 0)));
end

fprintf('================================\n\n');




%% extract SIR output for analysis

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

% Pick first seed site for analysis
example_site = flat_cay_site_IDs(1);

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
I_total_at_seeds_temp = I_LS_output_days(:, flat_cay_site_IDs) + I_MS_output_days(:, flat_cay_site_IDs) + I_HS_output_days(:, flat_cay_site_IDs);
plot(I_total_at_seeds_temp, 'LineWidth', 1.5);
xlabel('Day'); ylabel('Infected coral cover');
title('Infection trajectory at seed sites');
legend(arrayfun(@(x) sprintf('Site %d', unique_IDs(x)), flat_cay_site_IDs, 'UniformOutput', false), 'Location', 'best');
grid on;

subplot(2,2,2);
P_at_seeds_temp = I_LS_output_days(:,flat_cay_site_IDs) + I_MS_output_days(:,flat_cay_site_IDs) + I_HS_output_days(:,flat_cay_site_IDs);
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


% Check if external flux is truly zero at seed sites
P_initial = I_LS_init + I_MS_init + I_HS_init;  % Calculate initial disease pool
test_flux = incomingRisk_sparse(conn_structs(1).full, P_initial, export_thresh, 'rowsAreSources');
fprintf('Max incoming flux with export_thresh=999: %e\n', max(test_flux));
fprintf('Sites with non-zero flux: %d out of %d\n', sum(test_flux > 0), length(test_flux));


site_idx = find(reefData.unique_ID == 29089);
fprintf('Site 29089: LS=%.4f, MS=%.4f, HS=%.4f\n', N_LS(site_idx), N_MS(site_idx), N_HS(site_idx));







% Find sites with 0 in any group (in CURRENT data - no filtering applied)
sites_with_zeros = find((N_LS == 0) | (N_MS == 0) | (N_HS == 0));

fprintf('Sites with 0 in any group: %d\n', length(sites_with_zeros));

% Calculate total removed coral (R) at final timepoint for ONLY these sites
total_R_final = R_LS_output_days(end,:) + R_MS_output_days(end,:) + R_HS_output_days(end,:);
total_R_at_zero_sites = total_R_final(sites_with_zeros);

% Get top 99.5% percentile threshold among ONLY the zero-group sites
removal_threshold = prctile(total_R_at_zero_sites, 99.5);

% Find which zero-group sites have high removal
high_removal_indices = find(total_R_at_zero_sites >= removal_threshold);
high_removal_sites = sites_with_zeros(high_removal_indices);

fprintf('Zero-group sites with top 0.5%% removal: %d\n', length(high_removal_sites));

if ~isempty(high_removal_sites)
    % Pick a random one
    random_site = high_removal_sites(randi(length(high_removal_sites)));
    
    fprintf('Selected site: ID=%d, Index=%d\n', unique_IDs(random_site), random_site);
    fprintf('  LS: %.4f, MS: %.4f, HS: %.4f\n', N_LS(random_site), N_MS(random_site), N_HS(random_site));
    fprintf('  Total removed: %.4f\n', total_R_final(random_site));
    
    % Plot SIR curves
    figure('Name', 'High-Removal Site with Missing Group');
    
    subplot(2,2,1);
    plot(S_LS_output_days(:,random_site), 'b', 'LineWidth', 1.5); hold on;
    plot(I_LS_output_days(:,random_site), 'r', 'LineWidth', 1.5);
    plot(R_LS_output_days(:,random_site), 'k', 'LineWidth', 1.5);
    title(sprintf('LS - Site %d', unique_IDs(random_site)));
    legend('S','I','R'); grid on;
    
    subplot(2,2,2);
    plot(S_MS_output_days(:,random_site), 'b', 'LineWidth', 1.5); hold on;
    plot(I_MS_output_days(:,random_site), 'r', 'LineWidth', 1.5);
    plot(R_MS_output_days(:,random_site), 'k', 'LineWidth', 1.5);
    title(sprintf('MS - Site %d', unique_IDs(random_site)));
    legend('S','I','R'); grid on;
    
    subplot(2,2,3);
    plot(S_HS_output_days(:,random_site), 'b', 'LineWidth', 1.5); hold on;
    plot(I_HS_output_days(:,random_site), 'r', 'LineWidth', 1.5);
    plot(R_HS_output_days(:,random_site), 'k', 'LineWidth', 1.5);
    title(sprintf('HS - Site %d', unique_IDs(random_site)));
    legend('S','I','R'); grid on;
    
    subplot(2,2,4);
    total_S = S_LS_output_days(:,random_site) + S_MS_output_days(:,random_site) + S_HS_output_days(:,random_site);
    total_I = I_LS_output_days(:,random_site) + I_MS_output_days(:,random_site) + I_HS_output_days(:,random_site);
    total_R = R_LS_output_days(:,random_site) + R_MS_output_days(:,random_site) + R_HS_output_days(:,random_site);
    plot(total_S, 'b', 'LineWidth', 1.5); hold on;
    plot(total_I, 'r', 'LineWidth', 1.5);
    plot(total_R, 'k', 'LineWidth', 1.5);
    title(sprintf('Total - Site %d', unique_IDs(random_site)));
    legend('S','I','R'); grid on;
end
%flam2

% figure('Name', 'Sites with Missing Groups');
% 
% % Identify sites with zeros
% has_zero = (N_LS == 0) | (N_MS == 0) | (N_HS == 0);
% 
% % Plot all sites in gray
% scatter(locations(~has_zero, 1), locations(~has_zero, 2), 20, [0.7 0.7 0.7], 'filled', 'MarkerEdgeColor', 'none');
% hold on;
% 
% % Highlight sites with zeros in red
% scatter(locations(has_zero, 1), locations(has_zero, 2), 40, 'r', 'filled', 'MarkerEdgeColor', 'k');
% 
% xlabel('Longitude'); ylabel('Latitude');
% title(sprintf('Sites with Missing Groups (n=%d of %d)', sum(has_zero), num_sites));
% legend('Complete sites', 'Missing group(s)', 'Location', 'best');
% axis equal tight;
% grid on;

figure('Name', 'Site Filtering Comparison', 'Position', [100 100 1400 400]);

% Define filtering criteria
has_zero = (N_LS == 0) | (N_MS == 0) | (N_HS == 0);
low_cover = N_site < 0.01;

% Subplot 1: Both filters
subplot(1,3,1);
filtered_both = low_cover | has_zero;
scatter(locations(~filtered_both, 1), locations(~filtered_both, 2), 20, [0.7 0.7 0.7], 'filled', 'MarkerEdgeColor', 'none');
hold on;
scatter(locations(filtered_both, 1), locations(filtered_both, 2), 40, 'r', 'filled', 'MarkerEdgeColor', 'k');
xlabel('Longitude'); ylabel('Latitude');
title(sprintf('Both Filters\nRemoved: %d of %d (%.1f%%)', sum(filtered_both), num_sites, 100*sum(filtered_both)/num_sites));
legend('Kept', 'Removed', 'Location', 'best');
axis equal tight; grid on;

% Subplot 2: Just <1% cover filter
subplot(1,3,2);
scatter(locations(~low_cover, 1), locations(~low_cover, 2), 20, [0.7 0.7 0.7], 'filled', 'MarkerEdgeColor', 'none');
hold on;
scatter(locations(low_cover, 1), locations(low_cover, 2), 40, 'b', 'filled', 'MarkerEdgeColor', 'k');
xlabel('Longitude'); ylabel('Latitude');
title(sprintf('Low Cover Filter (<1%%)\nRemoved: %d of %d (%.1f%%)', sum(low_cover), num_sites, 100*sum(low_cover)/num_sites));
legend('Kept', 'Removed', 'Location', 'best');
axis equal tight; grid on;

% Subplot 3: Just missing groups filter
subplot(1,3,3);
scatter(locations(~has_zero, 1), locations(~has_zero, 2), 20, [0.7 0.7 0.7], 'filled', 'MarkerEdgeColor', 'none');
hold on;
scatter(locations(has_zero, 1), locations(has_zero, 2), 40, 'm', 'filled', 'MarkerEdgeColor', 'k');
xlabel('Longitude'); ylabel('Latitude');
title(sprintf('Missing Groups Filter\nRemoved: %d of %d (%.1f%%)', sum(has_zero), num_sites, 100*sum(has_zero)/num_sites));
legend('Kept', 'Removed', 'Location', 'best');
axis equal tight; grid on;




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

% Create figure (no renderer specification - MATLAB handles automatically)
f = figure('Position', [10 10 1400 1000]);
set(f, 'nextplot', 'replacechildren'); 

% Create 2x2 layout
T = tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Initialize plots
%
% absolute infected coral cover
axv1 = nexttile(T,1);
    h1 = scatter(axv1, locations(:,1), locations(:,2), 7, I_total_output_days(1,:)', 'filled');
    colormap(axv1, cmap_I);
    clim([0 .01]);
    c1 = colorbar(axv1);
    axis equal

% relative infected coral cover
axv2 = nexttile(T,2);
    h2 = scatter(axv2, locations(:,1), locations(:,2), 7, I_total_output_days(1,:)'./N_site(:), 'filled');
    colormap(axv2, cmap_I);
    clim([0 .01]);
    c2 = colorbar(axv2);
    axis equal 

% absolute removed coral cover
axv3 = nexttile(T,3);
    h3 = scatter(axv3, locations(:,1), locations(:,2), 7, N_site(:)-S_total_output_days(1,:)', 'filled');
    c3 = colorbar(axv3);
    clim([0 1]);
    colormap(axv3, cmap_R)
    axis equal

% relative removed coral cover
axv4 = nexttile(T,4);
    h4 = scatter(axv4, locations(:,1), locations(:,2), 7, (N_site(:)-S_total_output_days(1,:)')./N_site(:), 'filled');
    c4 = colorbar(axv4);
    clim([0 1]);
    colormap(axv4, cmap_R)
    axis equal

% Create video writer with MPEG-4 format
vidname = sprintf('DisVid_Diagnostic_thresh%d_scale%.1f_shape%.3f', ...
                  export_thresh, flux_scale, flux_shape);
v = VideoWriter(fullfile(seascapePath, vidname), 'MPEG-4');
v.Quality = 95;  % 0-100, higher = better quality but larger file
v.FrameRate = 5;
open(v);

% Threshold for disease front boundary
%   - this is saying that only once 0.1% cover within a site (patch) has
%       been removed, is that site considered "observable" to, for example,
%       the naked eye of strike team divers in the water within that ~0.4
%       m2 grid square
removed_cover_thresh = 0.001; % in proportion 0 to 1 (not 0-100%)

% Determine how many days actually simulated
num_days = size(I_total_output_days, 1);
fprintf('Creating animation for %d days of simulation\n', num_days);

% Generate animation frames
for k = 1:1:num_days
    Dt = datetime(2019,1,1) + (tspan(1) + k - 2);  % Adjust for actual tspan start

    % Find sites with dead coral above threshold (disease front)
    observable_sick_sites = R_total_output_days(k,:) >= removed_cover_thresh;
    if sum(observable_sick_sites) >= 3  % Need at least 3 points for boundary
        locs_sick_sites = locations(observable_sick_sites,:);
        try
            bounds_sick_sites = boundary(locs_sick_sites(:,1), locs_sick_sites(:,2), 0.7);
            draw_boundary = true;
        catch
            draw_boundary = false;  % Not enough points or other issue
        end
    else
        draw_boundary = false;
    end

    % Update plot 1: Disease prevalence (absolute % cover of site which is infected)
    set(h1, 'CData', I_total_output_days(k,:)');
    if draw_boundary
        p1 = patch(axv1, locs_sick_sites(bounds_sick_sites,1), locs_sick_sites(bounds_sick_sites,2), [.8 .9 1], ...
                   'EdgeColor', 'r', 'FaceAlpha', 0.2);
    end
    %
    title(axv1, 'Infected cover of seafloor (%)', ...
          sprintf('t = %s', string(Dt, 'dd-MMM-yyyy')));
    colormap(axv1, cmap_I);
    clim([0 .01]);
    axis equal

    % Update plot 2: Disease prevalence (relative % cover of site which is infected)
    set(h2, 'CData', I_total_output_days(k,:)'./N_site(:));
    if draw_boundary
        p2 = patch(axv2, locs_sick_sites(bounds_sick_sites,1), locs_sick_sites(bounds_sick_sites,2), [.8 .9 1], ...
                   'EdgeColor', 'r', 'FaceAlpha', 0.2);
    end
    %
    title(axv2, 'Relative cover of infected coral (%)', ...
          sprintf('t = %s', string(Dt, 'dd-MMM-yyyy')));
    colormap(axv2, cmap_I);
    clim([0 .1]);
    axis equal

    % Update plot 3: Total coral cover lost
    set(h3, 'CData', (N_site(:)-S_total_output_days(k,:)'));
    if draw_boundary
        p3 = patch(axv3, locs_sick_sites(bounds_sick_sites,1), locs_sick_sites(bounds_sick_sites,2), [.8 .9 1], ...
                   'EdgeColor', 'r', 'FaceAlpha', 0.2);
    end
    %
    title(axv3, 'Removed cover of seafloor (%)', ...
          sprintf('t = %s', string(Dt, 'dd-MMM-yyyy')));
    clim([0 1]);
    colormap(axv3, cmap_R)
    axis equal

    % Update plot 4: Proportion coral cover lost
    set(h4, 'CData', (N_site(:)-S_total_output_days(k,:)')./N_site(:));
    if draw_boundary
        p4 = patch(axv4, locs_sick_sites(bounds_sick_sites,1), locs_sick_sites(bounds_sick_sites,2), [.8 .9 1], ...
                   'EdgeColor', 'r', 'FaceAlpha', 0.2);
    end
    %
    title(axv4, 'Relative cover of removed coral (%)', ...
          sprintf('t = %s', string(Dt, 'dd-MMM-yyyy')));
    clim([0 1]);
    colormap(axv4, cmap_R)
    axis equal

    % Capture frame and write to video
    F = getframe(f);
    writeVideo(v, F);

    % Clean up patches for next frame
    if draw_boundary
        delete(p1)
        delete(p2)
        delete(p3)
        delete(p4)
    end

    % Progress indicator every 10 days
    if mod(k, 10) == 0
        fprintf('  Frame %d/%d\n', k, num_days);
    end
end

close(v);
fprintf('Animation saved as: %s.mp4\n', vidname);



%% Create static figure of final state

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
    figname = sprintf('FinalState_Diagnostic_day%d_thresh%d_scale%.1f_shape%.3f', ...
                      final_day, export_thresh, flux_scale, flux_shape);
    saveas(fig_final, fullfile(seascapePath, [figname '.png']));
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