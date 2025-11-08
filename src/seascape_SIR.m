
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

% DATE_RANGE = [];  % Empty = create/load all connectivity matrices
% DATE_RANGE = [datetime(2019,1,1), datetime(2019,3,31)];  % just Q1
DATE_RANGE = [datetime(2019,1,1), datetime(2019,6,30)];  % Q1-Q2

%% Load reef data from CSV

reefDataFile = fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv');
reefData = readtable(reefDataFile);

% Extract reef information
unique_IDs = reefData.unique_ID;
locations = [reefData.centroid_lon, reefData.centroid_lat];
num_sites = height(reefData);

%% Create or load existing connectivity matrices

% Generate connectivity (conn_structs) cache filename based on date range
if ~isempty(DATE_RANGE)
    cacheFilename = sprintf('conn_structs_%s_to_%s.mat', ...
                           string(DATE_RANGE(1), 'yyyyMMdd'), ...
                           string(DATE_RANGE(2), 'yyyyMMdd'));
else
    cacheFilename = 'conn_structs.mat';  % Full dataset
end
cacheFilePath = fullfile(tempPath, cacheFilename);


if RECALCULATE_CONNECTIVITY
    fprintf('========================================\n');
    fprintf('Loading connectivity matrices from disk\n');
    if ~isempty(DATE_RANGE)
        fprintf('Date range filter: %s to %s\n', ...
                string(DATE_RANGE(1), 'dd-MMM-yyyy'), ...
                string(DATE_RANGE(2), 'dd-MMM-yyyy'));
    end
    fprintf('========================================\n\n');
    
    % Pattern to match: connectivity_2019_[Date]_120000.mat
    quarters = {'Q1_2019', 'Q2_2019', 'Q3_2019', 'Q4_2019'};
    conn_structs = struct();
    connCount = 0;

    % NOTE: decay_weights normalization - currently using fixed value
    % TODO: Properly define decay_weights vector in future revision
    norm_val = 65 * 727.7434; % Placeholder: 65*sum(decay_weights)

    for q = 1:length(quarters)
        fprintf('Processing quarter: %s\n', quarters{q});
        quarterPath = fullfile(outputPath, 'CMS_traj', quarters{q});
        
        % Check if quarter folder exists
        if exist(quarterPath, 'dir') == 0
            warning('  --> Quarter folder not found: %s\n', quarterPath);
            continue;
        end
        
        % Get all connectivity files in this quarter
        myFiles = dir(fullfile(quarterPath, 'connectivity_*.mat'));
        fprintf('  Found %d connectivity files\n', length(myFiles));
        
        filesLoaded = 0;  % Track how many files actually loaded from this quarter
        
        for i = 1:length(myFiles)
            filename = myFiles(i).name;
            filenam = fullfile(quarterPath, filename);
            
            % Parse date from filename (e.g., connectivity_2019_Mar26_120000.mat)
            % Extract the date portion using regexp
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
            
            % Skip if outside date range
            if ~isempty(DATE_RANGE)
                if Dat < DATE_RANGE(1) || Dat > DATE_RANGE(2)
                    continue;  % Skip this file without loading
                end
            end
            
            % File is in range - now load it
            connCount = connCount + 1;
            filesLoaded = filesLoaded + 1;
            fprintf('  [%d] Loading: %s (Date: %s)\n', connCount, filename, string(Dat, 'dd-MMM-yyyy'));
            
            Con = load(filenam);  % Load full file now
            
            % Verify date matches (optional sanity check)
            Dat_actual = Con.connectivity_results.calendar_date;
            if abs(days(Dat - Dat_actual)) > 0.1
                warning('Filename date mismatch: %s vs %s', string(Dat, 'dd-MMM-yyyy'), string(Dat_actual, 'dd-MMM-yyyy'));

            end
            
            % Extract date information
            conn_structs(connCount).DY = day(Dat);
            conn_structs(connCount).MO = month(Dat);
            conn_structs(connCount).YR = year(Dat);
            conn_structs(connCount).Date = Dat;
            
            % Process connectivity matrix
            conmat = sparse(Con.connectivity_results.ConnMatrix_raw);
            
            % Remove self retention (diagonal)
            conmat = spdiags(zeros(size(conmat,1),1), 0, conmat);
            
            % Normalize
            conmat = conmat ./ norm_val;
            
            conn_structs(connCount).full = conmat;
        end
        
        if filesLoaded > 0
            fprintf('  Quarter %s complete: %d matrices loaded\n\n', quarters{q}, filesLoaded);
        else
            fprintf('  Quarter %s: No files in date range\n\n', quarters{q});
        end
    end

    % Sort by date
    fprintf('Sorting %d matrices by date...\n', connCount);
    [~, sortIdx] = sort([conn_structs.Date]);
    conn_structs = conn_structs(sortIdx);

    fprintf('Saving processed connectivity data to %s...\n', cacheFilename);
    save(cacheFilePath, 'conn_structs', '-v7.3')
    
    fprintf('\n========================================\n');
    fprintf('COMPLETE: Loaded %d connectivity matrices\n', connCount);
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
    
    load(cacheFilePath, 'conn_structs');
    fprintf('COMPLETE: Loaded %d connectivity matrices from cache\n', length(conn_structs));
    if ~isempty(DATE_RANGE)
        fprintf('Date range: %s to %s\n', ...
                string(min([conn_structs.Date]), 'dd-MMM-yyyy'), string(max([conn_structs.Date]), 'dd-MMM-yyyy'));
    end
    fprintf('========================================\n\n');
end

% Define day vector for the ODE - days since reference date
dates = [conn_structs.Date];
ref = datetime(2019,1,1);  
conn_days = days(dates - ref);

%% set state variables & parameters

% 'N' is coral cover pre-SCTLD
N_LS = reefData.low_coral_cover;
N_MS = reefData.moderate_coral_cover;
N_HS = reefData.high_coral_cover;

N_site = reefData.mean_coral_cover; %same as N_LS + N_MS + N_HS

% define how to seed disease at Flat Cay (or wherever chosen starting
%   location is placed)
%       NOTE - carefully consider seed value chosen here, and how many sites are
%       released from
seed_frac = 0.001;

% pre-define vectors for initial infected and recovered (dead) coral cover
I_LS_init = zeros(num_sites,1);
I_MS_init = zeros(num_sites,1);
I_HS_init = zeros(num_sites,1);
R_LS_init = zeros(num_sites,1);
R_MS_init = zeros(num_sites,1);
R_HS_init = zeros(num_sites,1);

% Flat Cay site IDs (29088 is the preferred/primary location):
flat_cay_site_IDs = [find(reefData.unique_ID==29088) find(reefData.unique_ID==29338) ...
        find(reefData.unique_ID==29089) find(reefData.unique_ID==29339) ...
        find(reefData.unique_ID==29087)];

I_HS_init(flat_cay_site_IDs) = seed_frac * N_HS(flat_cay_site_IDs);
N_HS(flat_cay_site_IDs) = N_HS(flat_cay_site_IDs) - I_HS_init(flat_cay_site_IDs);

I_MS_init(flat_cay_site_IDs) = seed_frac * N_MS(flat_cay_site_IDs);
N_MS(flat_cay_site_IDs) = N_MS(flat_cay_site_IDs) - I_MS_init(flat_cay_site_IDs);

I_LS_init(flat_cay_site_IDs) = seed_frac * N_LS(flat_cay_site_IDs);
N_LS(flat_cay_site_IDs) = N_LS(flat_cay_site_IDs) - I_LS_init(flat_cay_site_IDs);

% Pre-define vectors for initial susceptible coral cover
S_LS_init = N_LS - I_LS_init;
S_MS_init = N_MS - I_MS_init;
S_HS_init = N_HS - I_HS_init;

% Pack state variables for SIR function
Y0 = [S_LS_init; S_MS_init; S_HS_init; I_LS_init; I_MS_init; I_HS_init; R_LS_init; R_MS_init; R_HS_init];

% beta (transmission rate)
b_LS = 0.03;
b_MS = 0.14;
b_HS = 2.08;

% gamma (mortality rate)
g_LS = .05;
g_MS = .55;
g_HS = 3.33;

% threshold of (ABSOLUTE) diseased coral cover at which a site can begin
% to export disease to other sites
%   NOTE - should also consider thresholding based off of the relative
%   cover of a site / its susceptibility groups. or outbreak
%   stage/intensity
export_thresh = 0.0003; % At .001, transmission from flat_cay_site_IDs is unlikely with current parameterization. At .0005 disease immediately begins to transmit but not super duper fast...

% reshape parameters for controlling the contribution of upstream disease
% mass to local disease pool in each patch (site)
%   flux_scale.*(1-exp(-flux_shape.*T(:)))/(1-exp(-flux_shape));
flux_scale = 1; % limits max, ranges 0:1
flux_shape = -4; % determines curve shape. <0 is concave, >0 is convex; set small (.001) for no (linear) reshape. 0 yields infinity

%% serial test model run

opts = odeset('OutputFcn', @odeWaitbar);
tspan = [1 365]; 

% Adjust tspan to not exceed available connectivity data
max_conn_day = round(days(max([conn_structs.Date]) - datetime(2019,1,1)));

if tspan(end) > max_conn_day
    fprintf('\n*** TSPAN ADJUSTMENT ***\n');
    fprintf('Original tspan requested: [%d %d]\n', tspan(1), tspan(end));
    fprintf('Available connectivity data: up to day %d (date: %s)\n', ...
            max_conn_day, string(max([conn_structs.Date]), 'dd-MMM-yyyy'));
    fprintf('Adjusting tspan to: [%d %d]\n', tspan(1), max_conn_day);
    fprintf('************************\n\n');
    tspan(end) = max_conn_day;
end

% Store final tspan for use in parameter sweep
tspan_final = tspan;

% The ODE solver. Note that ODEfun_SCTLD_seascape is a separate .m file 
clear Y t
tic
[t,Y] = ode45(@(t,Y) ODEfun_SCTLD_seascape(t,Y,N_LS,N_MS,N_HS,b_LS,b_MS,b_HS,g_LS,g_MS,g_HS,export_thresh,flux_scale,flux_shape,num_sites,conn_structs,conn_days), tspan, Y0, opts);
toc

%% Parameter sweep for threshold and shape parameter optimization

threshvec = .0001:.0002:.0009;
flux_shapeVec = -4:2:4;

if RUN_PARAMETER_SWEEP
    count = 0;
    Results = struct();

    fprintf('\n========================================\n');
    fprintf('PARAMETER SWEEP: Testing %d combinations\n', length(threshvec)*length(flux_shapeVec));
    fprintf('  Thresholds: %d values from %.4f to %.4f\n', length(threshvec), min(threshvec), max(threshvec));
    fprintf('  ShapeParams: %d values from %d to %d\n', length(flux_shapeVec), min(flux_shapeVec), max(flux_shapeVec));
    if USE_PARALLEL
        fprintf('  Mode: PARALLEL (using %d workers)\n', min(6, feature('numcores')));
    else
        fprintf('  Mode: SERIAL\n');
    end
    fprintf('========================================\n\n');

    if USE_PARALLEL
        % Parallel execution
        if isempty(gcp('nocreate'))
            parpool('local', min(6, feature('numcores'))); % Use up to 6 cores
        end
        
        total_combos = length(threshvec)*length(flux_shapeVec);
        Results(total_combos).thresh = [];  % Pre-allocate struct array
        
        tic
        parfor combo = 1:total_combos
            [tv, sv] = ind2sub([length(threshvec), length(flux_shapeVec)], combo);
            
            fprintf('  [%d/%d] Running: thresh=%.4f, flux_shape=%d\n', ...
                    combo, total_combos, threshvec(tv), flux_shapeVec(sv));
            
            Y0_par = [N_LS;N_MS;N_HS;I_LS_init;I_MS_init;I_HS_init;R_LS_init;R_MS_init;R_HS_init];
            tspan_par = tspan_final;  % Use adjusted tspan
            
            [t_par,Y_par] = ode45(@(t,Y) ODEfun_SCTLD_seascape(t,Y,N_LS,N_MS,N_HS,b_LS,b_MS,b_HS,g_LS,g_MS,g_HS,threshvec(tv),flux_scale,flux_shapeVec(sv),num_sites,conn_structs,conn_days), tspan_par, Y0_par);
            
            % Extract and interpolate results
            LSHP = Y_par(:,1:num_sites);
            MSHP = Y_par(:,num_sites+1:num_sites*2);
            HSHP = Y_par(:,2*num_sites+1:num_sites*3);
            
            LSIP = Y_par(:,3*num_sites+1:num_sites*4);
            MSIP = Y_par(:,4*num_sites+1:num_sites*5);
            HSIP = Y_par(:,5*num_sites+1:num_sites*6);
            
            LSRP = Y_par(:,6*num_sites+1:num_sites*7);
            MSRP = Y_par(:,7*num_sites+1:num_sites*8);
            HSRP = Y_par(:,8*num_sites+1:num_sites*9);
            
            interpLSHP = max(0, interp1(t_par,LSHP,tspan_par(1):tspan_par(2)));
            interpMSHP = max(0, interp1(t_par,MSHP,tspan_par(1):tspan_par(2)));
            interpHSHP = max(0, interp1(t_par,HSHP,tspan_par(1):tspan_par(2)));
            
            interpLSIP = max(0, interp1(t_par,LSIP,tspan_par(1):tspan_par(2)));
            interpMSIP = max(0, interp1(t_par,MSIP,tspan_par(1):tspan_par(2)));
            interpHSIP = max(0, interp1(t_par,HSIP,tspan_par(1):tspan_par(2)));
            
            interpLSRP = max(0, interp1(t_par,LSRP,tspan_par(1):tspan_par(2)));
            interpMSRP = max(0, interp1(t_par,MSRP,tspan_par(1):tspan_par(2)));
            interpHSRP = max(0, interp1(t_par,HSRP,tspan_par(1):tspan_par(2)));
            
            Results(combo).thresh = threshvec(tv);
            Results(combo).flux_shape = flux_shapeVec(sv);
            Results(combo).LSS = interpLSHP;
            Results(combo).MSS = interpMSHP;
            Results(combo).HSS = interpHSHP;
            Results(combo).LSI = interpLSIP;
            Results(combo).MSI = interpMSIP;
            Results(combo).HSI = interpHSIP;
            Results(combo).LSR = interpLSRP;
            Results(combo).MSR = interpMSRP;
            Results(combo).HSR = interpHSRP;
        end
        total_time = toc;
        count = length(Results);  % Get actual count from Results array
        
    else
        % Serial execution (original code)
        tic
        for tv = 1:length(threshvec)
            fprintf('--- Threshold %d/%d (%.4f) ---\n', tv, length(threshvec), threshvec(tv));
            
            for sv = 1:length(flux_shapeVec)
                count = count+1;
                fprintf('  [%d/%d] Running: thresh=%.4f, flux_shape=%d ... ', ...
                        count, length(threshvec)*length(flux_shapeVec), threshvec(tv), flux_shapeVec(sv));
                
                run_start = tic;

                Y0 = [N_LS;N_MS;N_HS;I_LS_init;I_MS_init;I_HS_init;R_LS_init;R_MS_init;R_HS_init]; 
                opts  = odeset('OutputFcn', @odeWaitbar);
                tspan_sweep = tspan_final;  % Use adjusted tspan
                clear Y t
                
                [t,Y] = ode45(@(t,Y) ODEfun_SCTLD_seascape(t,Y,N_LS,N_MS,N_HS,b_LS,b_MS,b_HS,g_LS,g_MS,g_HS,threshvec(tv),flux_scale,flux_shapeVec(sv),num_sites,conn_structs,conn_days), tspan_sweep, Y0, opts);
                
                fprintf('Done (%.1f sec)\n', toc(run_start));

                LSHP = Y(:,1:num_sites);
                MSHP = Y(:,num_sites+1:num_sites*2);
                HSHP = Y(:,2*num_sites+1:num_sites*3);
                
                LSIP = Y(:,3*num_sites+1:num_sites*4);
                MSIP = Y(:,4*num_sites+1:num_sites*5);
                HSIP = Y(:,5*num_sites+1:num_sites*6);
                
                LSRP = Y(:,6*num_sites+1:num_sites*7);
                MSRP = Y(:,7*num_sites+1:num_sites*8);
                HSRP = Y(:,8*num_sites+1:num_sites*9);
                
                interpLSHP = interp1(t,LSHP,tspan_sweep(1):tspan_sweep(2));
                interpMSHP = interp1(t,MSHP,tspan_sweep(1):tspan_sweep(2));
                interpHSHP = interp1(t,HSHP,tspan_sweep(1):tspan_sweep(2));
                
                interpLSIP = interp1(t,LSIP,tspan_sweep(1):tspan_sweep(2));
                interpMSIP = interp1(t,MSIP,tspan_sweep(1):tspan_sweep(2));
                interpHSIP = interp1(t,HSIP,tspan_sweep(1):tspan_sweep(2));
                
                interpLSRP = interp1(t,LSRP,tspan_sweep(1):tspan_sweep(2));
                interpMSRP = interp1(t,MSRP,tspan_sweep(1):tspan_sweep(2));
                interpHSRP = interp1(t,HSRP,tspan_sweep(1):tspan_sweep(2));
                
                % Just in case
                interpLSHP(interpLSHP<0) = 0;
                interpMSHP(interpMSHP<0) = 0;
                interpHSHP(interpHSHP<0) = 0;
                
                interpLSIP(interpLSIP<0) = 0;
                interpMSIP(interpMSIP<0) = 0;
                interpHSIP(interpHSIP<0) = 0;
                
                interpLSRP(interpLSRP<0) = 0;
                interpMSRP(interpMSRP<0) = 0;
                interpHSRP(interpHSRP<0) = 0;

                Results(count).thresh = threshvec(tv);
                Results(count).flux_shape = flux_shapeVec(sv);
                Results(count).LSS = interpLSHP;
                Results(count).MSS = interpMSHP;
                Results(count).HSS = interpHSHP;

                Results(count).LSI = interpLSIP;
                Results(count).MSI = interpMSIP;
                Results(count).HSI = interpHSIP;

                Results(count).LSR = interpLSRP;
                Results(count).MSR = interpMSRP;
                Results(count).HSR = interpHSRP;
            end
            fprintf('\n');
        end
        total_time = toc;
        count = length(Results);
    end
            
    fprintf('========================================\n');
    fprintf('PARAMETER SWEEP COMPLETE\n');
    fprintf('  Total simulations: %d\n', count);
    fprintf('  Total time: %.1f minutes (avg: %.1f sec per simulation)\n', total_time/60, total_time/count);
    fprintf('========================================\n\n');
    
    % Save results
    fprintf('Saving parameter sweep results to temp/ParameterSweep_Results.mat...\n');
    save(fullfile(tempPath, 'ParameterSweep_Results.mat'), 'Results', 'threshvec', 'flux_shapeVec', '-v7.3');
    fprintf('Saved successfully.\n\n');
    
else
    % Load pre-calculated sweep results
    fprintf('\n========================================\n');
    fprintf('Loading cached parameter sweep results from temp/ParameterSweep_Results.mat...\n');
    
    if exist(fullfile(tempPath, 'ParameterSweep_Results.mat'), 'file') == 0
        error(['Parameter sweep cache file not found.\n' ...
               'Set RUN_PARAMETER_SWEEP = true to generate it.']);
    end
    
    load(fullfile(tempPath, 'ParameterSweep_Results.mat'), 'Results', 'threshvec', 'flux_shapeVec');
    fprintf('COMPLETE: Loaded %d parameter combinations from cache\n', length(Results));
    fprintf('========================================\n\n');
end

%% Visualize parameter sweep results - 5x5 comparison grid

% T = tiledlayout(5,5,'TileSpacing','compact','Padding','compact');
% 
% for i=1:length(Results)
%     RTIP = Results(i).LSI + Results(i).MSI + Results(i).HSI;
%     nexttile(T,i)
%     plot(RTIP)
%     title(strcat('thresh = ',num2str(Results(i).thresh),', flux_shape = ',num2str(Results(i).flux_shape)))
% end

%% Interpolate results from manual run for movie generation
% The output from the solver is not in days, but is in unequal time steps.
% You need to interpolate back to days..
LSHP = Y(:,1:num_sites);
MSHP = Y(:,num_sites+1:num_sites*2);
HSHP = Y(:,2*num_sites+1:num_sites*3);

LSIP = Y(:,3*num_sites+1:num_sites*4);
MSIP = Y(:,4*num_sites+1:num_sites*5);
HSIP = Y(:,5*num_sites+1:num_sites*6);

LSRP = Y(:,6*num_sites+1:num_sites*7);
MSRP = Y(:,7*num_sites+1:num_sites*8);
HSRP = Y(:,8*num_sites+1:num_sites*9);

interpLSHP = interp1(t,LSHP,tspan(1):tspan(2));
interpMSHP = interp1(t,MSHP,tspan(1):tspan(2));
interpHSHP = interp1(t,HSHP,tspan(1):tspan(2));

interpLSIP = interp1(t,LSIP,tspan(1):tspan(2));
interpMSIP = interp1(t,MSIP,tspan(1):tspan(2));
interpHSIP = interp1(t,HSIP,tspan(1):tspan(2));

interpLSRP = interp1(t,LSRP,tspan(1):tspan(2));
interpMSRP = interp1(t,MSRP,tspan(1):tspan(2));
interpHSRP = interp1(t,HSRP,tspan(1):tspan(2));

% Just in case
interpLSHP(interpLSHP<0) = 0;
interpMSHP(interpMSHP<0) = 0;
interpHSHP(interpHSHP<0) = 0;

interpLSIP(interpLSIP<0) = 0;
interpMSIP(interpMSIP<0) = 0;
interpHSIP(interpHSIP<0) = 0;

interpLSRP(interpLSRP<0) = 0;
interpMSRP(interpMSRP<0) = 0;
interpHSRP(interpHSRP<0) = 0;

% TIP = interpLSIP + interpMSIP + interpHSIP;
% TSP = interpLSHP + interpMSHP + interpHSHP;
% TRP = interpLSRP + interpMSRP + interpHSRP;

%% Make results into a movie - spatial visualization of disease spread

% Use Results(6) as example: thresh=0.0003, flux_shape=-4
TIP = Results(6).LSI + Results(6).MSI + Results(6).HSI;
TSP = Results(6).LSS + Results(6).MSS + Results(6).HSS;
TRP = Results(6).LSR + Results(6).MSR + Results(6).HSR;

f = figure('renderer', 'zbuffer','Position', [10 10 1400 1000]);
set(f,'nextplot','replacechildren'); 

% Create colormaps for infection (I) and recovery/dead (R)
cmap3 = colormap(flipud(autumn(round(max(max(TIP))*1000000)+1)));
cmap3(1,:) = [0 .3 1];
cmap4 = colormap(cool(round(max(N_site)*1000)+1));
cmap4(1,:) = [0 .3 1];

edges = [0 .000000001 .001 .1 1];
Cstart = [
    0.25 .5 0.25;   % light green
    1 1 .2;   % green
    1.00 0.60 0.00;   % yellow
    0.80 0.00 0.80;   % orange
];
Cend = [
    0.00 0.60 0.00;   % green
    1.00 0.60 0.00;   % yellow
    1.00 0.00 0.00;   % orange
    0 0 0;   % red
];
cmap_I = stackedColormap(edges, Cstart, Cend, 256, [1 1 1 1]);
cmap_R = stackedColormap(edges, Cstart, Cend, 256, [1 1 1 1]);

% Create 2x2 layout for movie
T = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% Initialize plots
axv1 = nexttile(T,1);
    h1 = scatter(axv1, locations(:,1),locations(:,2),7,TIP(1,:)','filled');
    colormap(axv1,cmap_I);
    clim([0 .01]);
    c1 = colorbar(axv1);
    axis equal
    
axv2 = nexttile(T,2);
    h2 = scatter(axv2, locations(:,1),locations(:,2),7,TIP(1,:)'./N_site(:),'filled');
    colormap(axv2,cmap_I);
    clim([0 .01]);
    c2 = colorbar(axv2);
    axis equal 
    
axv3 = nexttile(T,3);
    h3 = scatter(axv3, locations(:,1),locations(:,2),7,N_site(:)-TSP(1,:)','filled');
    c3 = colorbar(axv3);
    clim([0 1]);
    colormap(axv3,cmap_R)
    axis equal
    
axv4 = nexttile(T,4);
    h4 = scatter(axv4, locations(:,1),locations(:,2),7,(N_site(:)-TSP(1,:)')./N_site(:),'filled');
    c4 = colorbar(axv4);
    clim([0 1]);
    colormap(axv4,cmap_R)
    axis equal

% Create video writer
v = VideoWriter(fullfile(seascapePath, 'DisVid_Results6'));
v.FrameRate = 5;
open(v);

% Threshold for which sites to include when drawing disease front boundary
bthresh = 0.001; % .1% of bottom

% Generate movie frames
for k = 1:1:size(TIP,1)
    Dt = datetime('1-Jan-2019')+k-1;
    
    % Find sites with dead coral above threshold (disease front)
    sick = TRP(k,:) >= bthresh;
    dflocations = locations(sick,:);
    df = boundary(dflocations(:,1),dflocations(:,2),0.7);
    
    % Update plot 1: Disease prevalence - proportion of bottom
    set(h1, 'CData', TIP(k,:).');
    p1 = patch(axv1,dflocations(df,1),dflocations(df,2),[.8 .9 1],'EdgeColor','r','FaceAlpha',0.2);
    title(axv1,'Disease prevalence - proportion of bottom',strcat('t ='," ", string(Dt), 'dd-MMM-yyyy'));
    colormap(axv1,cmap_I);
    clim([0 .01]);
    c1 = colorbar(axv1);
    axis equal

    % Update plot 2: Disease prevalence - proportion of living coral
    set(h2, 'CData', TIP(k,:).'./N_site(:));
    p2 = patch(axv2, dflocations(df, 1), dflocations(df, 2), [.8 .9 1], 'EdgeColor', 'r', 'FaceAlpha', 0.2);
    title(axv2,'Disease prevalence - proportion of living coral',strcat('t ='," ", string(Dt), 'dd-MMM-yyyy'));
    colormap(axv2,cmap_I);
    clim([0 .1]);
    c2 = colorbar(axv2);
    axis equal
     
    % Update plot 3: Total coral cover lost
    set(h3, 'CData',(N_site(:)-TSP(k,:).'));
    p3 = patch(axv3,dflocations(df,1),dflocations(df,2),[.8 .9 1],'EdgeColor','r','FaceAlpha',0.2);
    title(axv3,'Total coral cover lost',strcat('t ='," ", string(Dt), 'dd-MMM-yyyy'));
    c3 = colorbar(axv3);
    clim([0 1]);
    colormap(axv3,cmap_R)
    axis equal

    % Update plot 4: Proportion coral cover lost
    set(h4, 'CData',(N_site(:)-TSP(k,:).')./N_site(:));
    p4 = patch(axv4,dflocations(df,1),dflocations(df,2),[.8 .9 1],'EdgeColor','r','FaceAlpha',0.2);
    title(axv4,'Proportion coral cover lost',strcat('t ='," ", string(Dt), 'dd-MMM-yyyy'));
    c4 = colorbar(axv4);
    clim([0 1]);
    colormap(axv4,cmap_R)
    axis equal
    
    % Capture frame and write to video
    F = getframe(f);
    writeVideo(v, F);

    % Clean up patches for next frame
    delete(p1)
    delete(p2)
    delete(p3)
    delete(p4)
end
close(v);



% %% Visualize SIR dynamics at sites with highest coral removal
% % Add this section at the end of your script
% 
% % Configuration
% numSitesToShow = 5;  % Number of top sites to visualize
% resultIndex = 6;     % Which parameter combination to visualize (change this)
% 
% fprintf('\n========================================\n');
% fprintf('Visualizing SIR dynamics for Result #%d\n', resultIndex);
% fprintf('  Threshold: %.4f\n', Results(resultIndex).thresh);
% fprintf('  ShapeParam: %d\n', Results(resultIndex).flux_shape);
% fprintf('========================================\n\n');
% 
% % Extract data for selected result
% TIP_selected = Results(resultIndex).LSI + Results(resultIndex).MSI + Results(resultIndex).HSI;
% TSP_selected = Results(resultIndex).LSS + Results(resultIndex).MSS + Results(resultIndex).HSS;
% TRP_selected = Results(resultIndex).LSR + Results(resultIndex).MSR + Results(resultIndex).HSR;
% 
% % Calculate total coral removal at each site (final - initial)
% totalRemoval = N_site' - TSP_selected(end,:);
% 
% % Find sites with highest removal
% [sortedRemoval, siteIndices] = sort(totalRemoval, 'descend');
% topSites = siteIndices(1:numSitesToShow);
% 
% % Create figure with subplots for each top site
% fig = figure('Position', [100 100 1400 900]);
% T = tiledlayout(ceil(numSitesToShow/2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% title(T, sprintf('SIR Dynamics at Top %d Most Affected Sites - Result #%d (thresh=%.4f, flux_shape=%d)', ...
%     numSitesToShow, resultIndex, Results(resultIndex).thresh, Results(resultIndex).flux_shape), ...
%     'FontSize', 14, 'FontWeight', 'bold');
% 
% % Time vector (days)
% days = 1:size(TIP_selected, 1);
% dates = datetime(2019,1,1) + days - 1;
% 
% % Plot each site
% for i = 1:numSitesToShow
%     siteIdx = topSites(i);
% 
%     nexttile(T, i);
%     hold on;
% 
%     % Extract SIR data for this site
%     S = TSP_selected(:, siteIdx);
%     I = TIP_selected(:, siteIdx);
%     R = TRP_selected(:, siteIdx);
% 
%     % Plot SIR curves
%     plot(days, S, 'b-', 'LineWidth', 2, 'DisplayName', 'Susceptible (S)');
%     plot(days, I, 'r-', 'LineWidth', 2, 'DisplayName', 'Infected (I)');
%     plot(days, R, 'k-', 'LineWidth', 2, 'DisplayName', 'Removed (R)');
% 
%     % Add initial coral cover reference line
%     yline(N_site(siteIdx), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Initial Cover', 'Alpha', 0.5);
% 
%     % Formatting
%     xlabel('Days from Jan 1, 2019');
%     ylabel('Coral Cover Proportion');
%     title(sprintf('Site #%d (ID: %d) - %.1f%% Loss', ...
%         siteIdx, unique_IDs(siteIdx), 100*totalRemoval(siteIdx)/N_site(siteIdx)));
%     legend('Location', 'best', 'FontSize', 8);
%     grid on;
%     ylim([0, max(N_site(siteIdx)*1.1, 0.01)]);
% 
%     hold off;
% end
% 
% % Add summary statistics
% fprintf('Top %d sites by coral removal:\n', numSitesToShow);
% fprintf('Rank | Site# | Reef ID | Initial Cover | Final Cover | Removal | %% Loss\n');
% fprintf('-----|-------|---------|---------------|-------------|---------|--------\n');
% for i = 1:numSitesToShow
%     siteIdx = topSites(i);
%     initialCover = N_site(siteIdx);
%     finalCover = TSP_selected(end, siteIdx);
%     removal = totalRemoval(siteIdx);
%     percentLoss = 100 * removal / initialCover;
% 
%     fprintf('%4d | %5d | %7d | %13.4f | %11.4f | %7.4f | %6.1f%%\n', ...
%         i, siteIdx, unique_IDs(siteIdx), initialCover, finalCover, removal, percentLoss);
% end
% fprintf('\n');
% 
% %% Alternative view: All sites on one plot (normalized)
% fig2 = figure('Position', [150 150 1000 600]);
% hold on;
% 
% % Color scheme for sites
% colors = lines(numSitesToShow);
% 
% for i = 1:numSitesToShow
%     siteIdx = topSites(i);
% 
%     % Normalize by initial cover for comparison
%     S_norm = TSP_selected(:, siteIdx) / N_site(siteIdx);
%     I_norm = TIP_selected(:, siteIdx) / N_site(siteIdx);
%     R_norm = TRP_selected(:, siteIdx) / N_site(siteIdx);
% 
%     % Plot with different line styles
%     plot(days, S_norm, '-', 'Color', colors(i,:), 'LineWidth', 1.5, ...
%         'DisplayName', sprintf('Site %d - S', unique_IDs(siteIdx)));
%     plot(days, I_norm, '--', 'Color', colors(i,:), 'LineWidth', 1.5, ...
%         'DisplayName', sprintf('Site %d - I', unique_IDs(siteIdx)));
%     plot(days, R_norm, ':', 'Color', colors(i,:), 'LineWidth', 2, ...
%         'DisplayName', sprintf('Site %d - R', unique_IDs(siteIdx)));
% end
% 
% xlabel('Days from Jan 1, 2019');
% ylabel('Proportion of Initial Coral Cover');
% title(sprintf('Normalized SIR Dynamics - Result #%d (thresh=%.4f, flux_shape=%d)', ...
%     resultIndex, Results(resultIndex).thresh, Results(resultIndex).flux_shape));
% legend('Location', 'eastoutside', 'FontSize', 8);
% grid on;
% ylim([0, 1.1]);
% hold off;
% 
% %% Additional plots focusing on Infected dynamics
% 
% % Figure 3: Infected curves comparison - absolute values
% fig3 = figure('Position', [200 200 1200 500]);
% subplot(1,2,1);
% hold on;
% 
% for i = 1:numSitesToShow
%     siteIdx = topSites(i);
%     I = TIP_selected(:, siteIdx);
% 
%     plot(days, I, '-', 'Color', colors(i,:), 'LineWidth', 2, ...
%         'DisplayName', sprintf('Site %d (ID: %d)', siteIdx, unique_IDs(siteIdx)));
% end
% 
% xlabel('Days from Jan 1, 2019');
% ylabel('Infected Coral Cover Proportion');
% title(sprintf('Infected Coral Time Series - Absolute Values\nResult #%d (thresh=%.4f, flux_shape=%d)', ...
%     resultIndex, Results(resultIndex).thresh, Results(resultIndex).flux_shape));
% legend('Location', 'best', 'FontSize', 9);
% grid on;
% hold off;
% 
% % Subplot 2: Normalized by initial cover
% subplot(1,2,2);
% hold on;
% 
% for i = 1:numSitesToShow
%     siteIdx = topSites(i);
%     I_norm = TIP_selected(:, siteIdx) / N_site(siteIdx);
% 
%     plot(days, I_norm, '-', 'Color', colors(i,:), 'LineWidth', 2, ...
%         'DisplayName', sprintf('Site %d (ID: %d)', siteIdx, unique_IDs(siteIdx)));
% end
% 
% xlabel('Days from Jan 1, 2019');
% ylabel('Infected / Initial Coral Cover');
% title(sprintf('Infected Coral Time Series - Normalized\nResult #%d', resultIndex));
% legend('Location', 'best', 'FontSize', 9);
% grid on;
% hold off;
% 
% % Figure 4: Infection dynamics statistics
% fig4 = figure('Position', [250 250 1200 700]);
% T2 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% title(T2, sprintf('Infection Dynamics Analysis - Result #%d', resultIndex), ...
%     'FontSize', 14, 'FontWeight', 'bold');
% 
% % Panel 1: Peak infection timing
% nexttile(T2, 1);
% hold on;
% peakTimes = zeros(numSitesToShow, 1);
% peakValues = zeros(numSitesToShow, 1);
% 
% for i = 1:numSitesToShow
%     siteIdx = topSites(i);
%     I = TIP_selected(:, siteIdx);
%     [peakVal, peakIdx] = max(I);
%     peakTimes(i) = peakIdx;
%     peakValues(i) = peakVal;
% 
%     bar(i, peakIdx, 'FaceColor', colors(i,:));
% end
% 
% xticks(1:numSitesToShow);
% xticklabels(arrayfun(@(x) sprintf('Site %d', unique_IDs(topSites(x))), 1:numSitesToShow, 'UniformOutput', false));
% xtickangle(45);
% ylabel('Day of Peak Infection');
% title('Timing of Peak Infection');
% grid on;
% hold off;
% 
% % Panel 2: Peak infection magnitude
% nexttile(T2, 2);
% hold on;
% 
% for i = 1:numSitesToShow
%     bar(i, peakValues(i), 'FaceColor', colors(i,:));
% end
% 
% xticks(1:numSitesToShow);
% xticklabels(arrayfun(@(x) sprintf('Site %d', unique_IDs(topSites(x))), 1:numSitesToShow, 'UniformOutput', false));
% xtickangle(45);
% ylabel('Peak Infected Cover Proportion');
% title('Magnitude of Peak Infection');
% grid on;
% hold off;
% 
% % Panel 3: Infection rate (derivative of I)
% nexttile(T2, 3);
% hold on;
% 
% for i = 1:numSitesToShow
%     siteIdx = topSites(i);
%     I = TIP_selected(:, siteIdx);
%     dIdt = diff(I);  % Rate of change
% 
%     plot(days(2:end), dIdt, '-', 'Color', colors(i,:), 'LineWidth', 1.5, ...
%         'DisplayName', sprintf('Site %d', unique_IDs(siteIdx)));
% end
% 
% xlabel('Days from Jan 1, 2019');
% ylabel('dI/dt (Change in Infected)');
% title('Infection Rate Over Time');
% legend('Location', 'best', 'FontSize', 8);
% grid on;
% yline(0, 'k--', 'Alpha', 0.5);
% hold off;
% 
% % Panel 4: Cumulative infection burden (area under I curve)
% nexttile(T2, 4);
% hold on;
% cumBurden = zeros(numSitesToShow, 1);
% 
% for i = 1:numSitesToShow
%     siteIdx = topSites(i);
%     I = TIP_selected(:, siteIdx);
%     cumBurden(i) = trapz(days, I);  % Integrate infection over time
% 
%     bar(i, cumBurden(i), 'FaceColor', colors(i,:));
% end
% 
% xticks(1:numSitesToShow);
% xticklabels(arrayfun(@(x) sprintf('Site %d', unique_IDs(topSites(x))), 1:numSitesToShow, 'UniformOutput', false));
% xtickangle(45);
% ylabel('Cumulative Infection Burden');
% title('Total Infection Burden (Area Under I Curve)');
% grid on;
% hold off;
% 
% % Print infection statistics
% fprintf('\n--- INFECTION DYNAMICS STATISTICS ---\n');
% fprintf('Site# | Reef ID | Peak Day | Peak I | Peak I/C0 | Cum. Burden\n');
% fprintf('------|---------|----------|--------|-----------|------------\n');
% for i = 1:numSitesToShow
%     siteIdx = topSites(i);
%     fprintf('%5d | %7d | %8d | %6.4f | %9.4f | %11.4f\n', ...
%         siteIdx, unique_IDs(siteIdx), peakTimes(i), peakValues(i), ...
%         peakValues(i)/N_site(siteIdx), cumBurden(i));
% end
% fprintf('\n');
% 
% fprintf('Visualization complete!\n');
% fprintf('TIP: Change resultIndex variable (line 6) to view different parameter combinations\n\n');
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% Local SIR Outbreak Simulator - Test different community compositions
% % This simulates disease dynamics at isolated sites (no connectivity)
% % to understand how community composition affects outbreak patterns
% % 
% % ADD THIS SECTION to your existing script after the parameter sweep section
% 
% %% Define test scenarios - EXACTLY matching the R example
% % Each row: [LS%, MS%, HS%, Total Cover]
% % Covers are in % (will be divided by 100 to get proportions)
% scenarios = [
%     3.64,  15.7,  4.97, NaN;  % Nearshore - will calculate total
%     10,    5,     2,    NaN;  % Scenario 3
%     5,     2.5,   1,    NaN;  % Scenario 4
%     2.5,   1.25,  3,    NaN;  % Scenario 5
%     0.5,   5,     1,    NaN;  % Scenario 6
%     0.25,  2.5,   0.5,  NaN;  % Scenario 7
%     0.125, 1.25,  0.25, NaN;  % Scenario 8
%     0.1,   1,     0.02, NaN;  % Scenario 9
% ];
% 
% % Calculate total cover for each scenario
% scenarios(:,4) = sum(scenarios(:,1:3), 2);
% 
% % Convert from % to proportions (0-1 scale)
% scenarios = scenarios / 100;
% 
% numScenarios = size(scenarios, 1);
% 
% %% Parameters (using same values from your main script)
% % These are already defined in your environment, but listed here for clarity
% % b_LS = 0.03;
% % b_MS = 0.14;
% % b_HS = 2.08;
% % g_LS = 0.05;
% % g_MS = 0.55;
% % g_HS = 3.33;
% 
% % Time span for local simulations
% tspan_local = [0 365];  % 365 days
% tspan_vec_local = 0:1:365;
% 
% %% Run simulations for each scenario
% LocalResults = struct();
% 
% fprintf('========================================\n');
% fprintf('Running Local Outbreak Simulations\n');
% fprintf('Testing %d community composition scenarios\n', numScenarios);
% fprintf('========================================\n\n');
% 
% for s = 1:numScenarios
%     fprintf('Scenario %d: LS=%.2f%%, MS=%.2f%%, HS=%.2f%%, Total Cover=%.2f%%\n', ...
%         s, scenarios(s,1)*100, scenarios(s,2)*100, scenarios(s,3)*100, scenarios(s,4)*100);
% 
%     % Calculate initial conditions from proportions (already 0-1 scale)
%     CLS_init = scenarios(s,1);
%     CMS_init = scenarios(s,2);
%     CHS_init = scenarios(s,3);
%     totalCover = scenarios(s,4);
% 
%     % Initial susceptible
%     SLS0_local = CLS_init;
%     SMS0_local = CMS_init;
%     SHS0_local = CHS_init;
% 
%     % Seed initial infection (0.01% of cover, matching your working example)
%     % Priority: seed HS if available, else MS, else LS
%     epsilon = 1e-6;
%     initial_infected = 0.0001;  % 0.01% as proportion
% 
%     if CHS_init > epsilon
%         IHS0_local = initial_infected;
%         IMS0_local = 0;
%         ILS0_local = 0;
%     elseif CMS_init > epsilon
%         IHS0_local = 0;
%         IMS0_local = initial_infected;
%         ILS0_local = 0;
%     else
%         IHS0_local = 0;
%         IMS0_local = 0;
%         ILS0_local = initial_infected;
%     end
% 
%     % Adjust susceptible to account for initial infection
%     SLS0_local = SLS0_local - ILS0_local;
%     SMS0_local = SMS0_local - IMS0_local;
%     SHS0_local = SHS0_local - IHS0_local;
% 
%     % Initial removed (dead)
%     RLS0_local = 0;
%     RMS0_local = 0;
%     RHS0_local = 0;
% 
%     % Pack into initial condition vector
%     Y0_local = [SLS0_local; SMS0_local; SHS0_local; ILS0_local; IMS0_local; IHS0_local; RLS0_local; RMS0_local; RHS0_local];
% 
%     % Run ODE with TRUE dummy connectivity that nullifies spatial transmission
%     % Create connectivity structure for single isolated site
%     P_dummy(1).full = sparse(1, 1, 0);  % 1x1 matrix with zero connectivity
%     P_dummy(1).Date = datetime(2019,1,1);
%     conn_days_dummy = [0];  % Single time point
% 
%     % Critical parameters to nullify spatial transmission while keeping local:
%     % - Set c = 0: This zeros out ALL external transmission (T becomes 0)
%     % - High threshold ensures no transmission via connectivity anyway
%     % - With c=0, the reshape function produces T=0 regardless of flux_shape
%     thresh_dummy = 1e10;     % Impossibly high threshold
%     c_dummy = 0;              % Zero coefficient = no external transmission
%     flux_shape_dummy = 0.001; % Doesn't matter when c=0
%     num_sites_dummy = 1;           % Single site (no spatial structure)
% 
%     % Run using ODEfun_SCTLD_seascape with connectivity nullified
%     [t_local, Y_local] = ode45(@(t,Y) ODEfun_SCTLD_seascape(t, Y, ...
%                                                              CLS_init, CMS_init, CHS_init, ...
%                                                              b_LS, b_MS, b_HS, g_LS, g_MS, g_HS, ...
%                                                              thresh_dummy, c_dummy, flux_shape_dummy, ...
%                                                              num_sites_dummy, P_dummy, conn_days_dummy), ...
%                                tspan_local, Y0_local);
% 
%     % Extract results
%     SLS_res = Y_local(:,1);
%     SMS_res = Y_local(:,2);
%     SHS_res = Y_local(:,3);
%     ILS_res = Y_local(:,4);
%     IMS_res = Y_local(:,5);
%     IHS_res = Y_local(:,6);
%     RLS_res = Y_local(:,7);
%     RMS_res = Y_local(:,8);
%     RHS_res = Y_local(:,9);
% 
%     % Interpolate to daily values
%     SLS_interp = interp1(t_local, SLS_res, tspan_vec_local);
%     SMS_interp = interp1(t_local, SMS_res, tspan_vec_local);
%     SHS_interp = interp1(t_local, SHS_res, tspan_vec_local);
%     ILS_interp = interp1(t_local, ILS_res, tspan_vec_local);
%     IMS_interp = interp1(t_local, IMS_res, tspan_vec_local);
%     IHS_interp = interp1(t_local, IHS_res, tspan_vec_local);
%     RLS_interp = interp1(t_local, RLS_res, tspan_vec_local);
%     RMS_interp = interp1(t_local, RMS_res, tspan_vec_local);
%     RHS_interp = interp1(t_local, RHS_res, tspan_vec_local);
% 
%     % Store results
%     LocalResults(s).scenario = scenarios(s,:);
%     LocalResults(s).SLS = SLS_interp;
%     LocalResults(s).SMS = SMS_interp;
%     LocalResults(s).SHS = SHS_interp;
%     LocalResults(s).ILS = ILS_interp;
%     LocalResults(s).IMS = IMS_interp;
%     LocalResults(s).IHS = IHS_interp;
%     LocalResults(s).RLS = RLS_interp;
%     LocalResults(s).RMS = RMS_interp;
%     LocalResults(s).RHS = RHS_interp;
%     LocalResults(s).totalCover = totalCover;
% 
%     % Calculate totals
%     LocalResults(s).S_total = SLS_interp + SMS_interp + SHS_interp;
%     LocalResults(s).I_total = ILS_interp + IMS_interp + IHS_interp;
%     LocalResults(s).R_total = RLS_interp + RMS_interp + RHS_interp;
% 
%     % Add this right after the ode45 call in the loop:
%     if s == 1  % Only print for first scenario to verify dummy setup
%         fprintf('\n=== VERIFYING DUMMY CONDITIONS (SCENARIO 1) ===\n');
%         fprintf('Setup ensures PURE within-site transmission:\n');
%         fprintf('  - Single site (num_sites=1): no spatial neighbors\n');
%         fprintf('  - Connectivity matrix: 1x1 sparse with value 0\n');
%         fprintf('  - c_dummy=0: External transmission T forced to zero\n');
%         fprintf('  - thresh_dummy=1e10: Impossibly high for any transmission\n\n');
% 
%         fprintf('Initial conditions:\n');
%         fprintf('  Cover: LS=%.4f, MS=%.4f, HS=%.4f (Total=%.2f%%)\n', ...
%                 CLS_init, CMS_init, CHS_init, (CLS_init+CMS_init+CHS_init)*100);
%         fprintf('  Initial infected: LS=%.6f, MS=%.6f, HS=%.6f\n', ...
%                 ILS0_local, IMS0_local, IHS0_local);
%         fprintf('  (Seeded 0.01%% = %.6f in highest susceptibility group present)\n', initial_infected);
% 
%         fprintf('\nDynamics check:\n');
%         fprintf('  Peak infected day: LS=%.1f, MS=%.1f, HS=%.1f\n', ...
%                 t_local(find(Y_local(:,4)==max(Y_local(:,4)),1)), ...
%                 t_local(find(Y_local(:,5)==max(Y_local(:,5)),1)), ...
%                 t_local(find(Y_local(:,6)==max(Y_local(:,6)),1)));
%         fprintf('  Peak values: LS=%.6f, MS=%.6f, HS=%.6f\n', ...
%                 max(Y_local(:,4)), max(Y_local(:,5)), max(Y_local(:,6)));
%         fprintf('  Final removed: LS=%.6f, MS=%.6f, HS=%.6f (Total=%.6f)\n', ...
%                 Y_local(end,7), Y_local(end,8), Y_local(end,9), ...
%                 Y_local(end,7)+Y_local(end,8)+Y_local(end,9));
%         fprintf('=========================================\n\n');
%     end
% 
%     fprintf('  Complete.\n\n');
% end
% 
% %% Visualization 1: All scenarios comparison
% fig_local1 = figure('Position', [50 50 1400 900]);
% T1 = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
% title(T1, 'Local SIR Outbreak Dynamics - Different Community Compositions', ...
%     'FontSize', 14, 'FontWeight', 'bold');
% 
% for s = 1:numScenarios
%     nexttile(T1, s);
%     hold on;
% 
%     % Plot total S, I, R
%     plot(tspan_vec_local, LocalResults(s).S_total, 'b-', 'LineWidth', 2, 'DisplayName', 'Susceptible');
%     plot(tspan_vec_local, LocalResults(s).I_total, 'r-', 'LineWidth', 2, 'DisplayName', 'Infected');
%     plot(tspan_vec_local, LocalResults(s).R_total, 'k-', 'LineWidth', 2, 'DisplayName', 'Removed');
% 
%     xlabel('Days');
%     ylabel('Coral Cover Proportion');
%     title(sprintf('LS:%.2f%% MS:%.2f%% HS:%.2f%% | Cov:%.2f%%', ...
%         scenarios(s,1)*100, scenarios(s,2)*100, scenarios(s,3)*100, scenarios(s,4)*100));
%     legend('Location', 'best', 'FontSize', 7);
%     grid on;
%     ylim([0, LocalResults(s).totalCover * 1.1]);
%     hold off;
% end
% 
% %% Visualization 2: Detailed breakdown by susceptibility class
% fig_local2 = figure('Position', [100 100 1400 900]);
% T2 = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
% title(T2, 'SIR Dynamics by Susceptibility Class', 'FontSize', 14, 'FontWeight', 'bold');
% 
% for s = 1:numScenarios
%     nexttile(T2, s);
%     hold on;
% 
%     % Plot each class separately
%     % LS group
%     plot(tspan_vec_local, LocalResults(s).SLS, 'b-', 'LineWidth', 1.5, 'DisplayName', 'LS-S');
%     plot(tspan_vec_local, LocalResults(s).ILS, 'b--', 'LineWidth', 1.5, 'DisplayName', 'LS-I');
%     plot(tspan_vec_local, LocalResults(s).RLS, 'b:', 'LineWidth', 2, 'DisplayName', 'LS-R');
% 
%     % MS group
%     plot(tspan_vec_local, LocalResults(s).SMS, 'g-', 'LineWidth', 1.5, 'DisplayName', 'MS-S');
%     plot(tspan_vec_local, LocalResults(s).IMS, 'g--', 'LineWidth', 1.5, 'DisplayName', 'MS-I');
%     plot(tspan_vec_local, LocalResults(s).RMS, 'g:', 'LineWidth', 2, 'DisplayName', 'MS-R');
% 
%     % HS group
%     plot(tspan_vec_local, LocalResults(s).SHS, 'r-', 'LineWidth', 1.5, 'DisplayName', 'HS-S');
%     plot(tspan_vec_local, LocalResults(s).IHS, 'r--', 'LineWidth', 1.5, 'DisplayName', 'HS-I');
%     plot(tspan_vec_local, LocalResults(s).RHS, 'r:', 'LineWidth', 2, 'DisplayName', 'HS-R');
% 
%     xlabel('Days');
%     ylabel('Coral Cover Proportion');
%     title(sprintf('LS:%.2f%% MS:%.2f%% HS:%.2f%% | Cov:%.2f%%', ...
%         scenarios(s,1)*100, scenarios(s,2)*100, scenarios(s,3)*100, scenarios(s,4)*100));
%     legend('Location', 'eastoutside', 'FontSize', 6);
%     grid on;
%     hold off;
% end
% 
% %% Visualization 3: Infection dynamics comparison
% fig_local3 = figure('Position', [150 150 1200 800]);
% T3 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% title(T3, 'Infection Dynamics Comparison Across Scenarios', 'FontSize', 14, 'FontWeight', 'bold');
% 
% colors_local = lines(numScenarios);
% 
% % Panel 1: Total infected over time
% nexttile(T3, 1);
% hold on;
% for s = 1:numScenarios
%     plot(tspan_vec_local, LocalResults(s).I_total, '-', 'Color', colors_local(s,:), 'LineWidth', 2, ...
%         'DisplayName', sprintf('Sc%d: L%.1f M%.1f H%.1f', s, scenarios(s,1)*100, scenarios(s,2)*100, scenarios(s,3)*100));
% end
% xlabel('Days');
% ylabel('Total Infected');
% title('Total Infection Over Time');
% legend('Location', 'best', 'FontSize', 7);
% grid on;
% hold off;
% 
% % Panel 2: Peak infection comparison
% nexttile(T3, 2);
% hold on;
% peakInfections_local = zeros(numScenarios, 1);
% peakDays_local = zeros(numScenarios, 1);
% for s = 1:numScenarios
%     [peakInfections_local(s), peakIdx] = max(LocalResults(s).I_total);
%     peakDays_local(s) = tspan_vec_local(peakIdx);
%     bar(s, peakInfections_local(s), 'FaceColor', colors_local(s,:));
% end
% xticks(1:numScenarios);
% xticklabels(arrayfun(@(x) sprintf('Sc%d', x), 1:numScenarios, 'UniformOutput', false));
% ylabel('Peak Infected Cover');
% title('Peak Infection Magnitude');
% grid on;
% hold off;
% 
% % Panel 3: Final mortality comparison
% nexttile(T3, 3);
% hold on;
% finalMortality_local = zeros(numScenarios, 1);
% for s = 1:numScenarios
%     finalMortality_local(s) = LocalResults(s).R_total(end);
%     bar(s, finalMortality_local(s), 'FaceColor', colors_local(s,:));
% end
% xticks(1:numScenarios);
% xticklabels(arrayfun(@(x) sprintf('Sc%d', x), 1:numScenarios, 'UniformOutput', false));
% ylabel('Final Removed (Dead) Cover');
% title('Total Mortality by Scenario');
% grid on;
% hold off;
% 
% % Panel 4: Percent loss
% nexttile(T3, 4);
% hold on;
% pctLoss_local = zeros(numScenarios, 1);
% for s = 1:numScenarios
%     pctLoss_local(s) = 100 * finalMortality_local(s) / LocalResults(s).totalCover;
%     bar(s, pctLoss_local(s), 'FaceColor', colors_local(s,:));
% end
% xticks(1:numScenarios);
% xticklabels(arrayfun(@(x) sprintf('Sc%d', x), 1:numScenarios, 'UniformOutput', false));
% ylabel('Percent Coral Loss (%)');
% title('Percent Cover Lost');
% grid on;
% hold off;
% 
% %% Summary statistics
% fprintf('\n========================================\n');
% fprintf('LOCAL OUTBREAK SIMULATION SUMMARY\n');
% fprintf('========================================\n');
% fprintf('Sc | LS%% | MS%% | HS%% | Cover | Peak I | Peak Day | Final R | %% Loss\n');
% fprintf('---|------|------|------|-------|--------|----------|---------|--------\n');
% for s = 1:numScenarios
%     fprintf('%2d | %4.2f | %4.2f | %4.2f | %5.2f%% | %6.4f | %8d | %7.4f | %6.1f%%\n', ...
%         s, scenarios(s,1)*100, scenarios(s,2)*100, scenarios(s,3)*100, scenarios(s,4)*100, ...
%         peakInfections_local(s), peakDays_local(s), finalMortality_local(s), pctLoss_local(s));
% end
% fprintf('========================================\n\n');