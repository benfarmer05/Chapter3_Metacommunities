clear; clc

%% Initialize paths
projectPath = matlab.project.rootProject().RootFolder;
tempPath = fullfile(projectPath, 'temp');
dataPath = fullfile(projectPath, 'data');
outputPath = fullfile(projectPath, 'output');
seascapePath = fullfile(outputPath, 'seascape_SIR');

if ~exist(seascapePath, 'dir')
    mkdir(seascapePath);
end

%% Toggle: Load connectivity matrices or use cached version
RECALCULATE_CONNECTIVITY = false;
RUN_PARAMETER_SWEEP = true;
USE_PARALLEL = true;
DATE_RANGE = [datetime(2019,1,1), datetime(2019,3,31)];

%% Load reef data from CSV
reefDataFile = fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv');
reefData = readtable(reefDataFile);

unique_IDs = reefData.unique_ID;
XY = [reefData.centroid_lon, reefData.centroid_lat];
habs = height(reefData);

CLS1 = reefData.low_coral_cover;
CMS1 = reefData.moderate_coral_cover;
CHS1 = reefData.high_coral_cover;
CC = reefData.mean_coral_cover;

%% Load connectivity matrices
if ~isempty(DATE_RANGE)
    cacheFilename = sprintf('P_%s_to_%s.mat', ...
                           datestr(DATE_RANGE(1), 'yyyymmdd'), ...
                           datestr(DATE_RANGE(2), 'yyyymmdd'));
else
    cacheFilename = 'P.mat';
end
cacheFilePath = fullfile(tempPath, cacheFilename);

if RECALCULATE_CONNECTIVITY
    fprintf('========================================\n');
    fprintf('Loading connectivity matrices from disk\n');
    if ~isempty(DATE_RANGE)
        fprintf('Date range filter: %s to %s\n', ...
                datestr(DATE_RANGE(1)), datestr(DATE_RANGE(2)));
    end
    fprintf('========================================\n\n');
    
    quarters = {'Q1_2019', 'Q2_2019', 'Q3_2019', 'Q4_2019'};
    P = struct();
    connCount = 0;
    norm_val = 65 * 727.7434;

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
            fprintf('  [%d] Loading: %s (Date: %s)\n', connCount, filename, datestr(Dat));
            
            Con = load(filenam);
            
            Dat_actual = Con.connectivity_results.calendar_date;
            if abs(days(Dat - Dat_actual)) > 0.1
                warning('Filename date mismatch: %s vs %s', datestr(Dat), datestr(Dat_actual));
            end
            
            P(connCount).DY = day(Dat);
            P(connCount).MO = month(Dat);
            P(connCount).YR = year(Dat);
            P(connCount).Date = Dat;
            
            conmat = sparse(Con.connectivity_results.ConnMatrix_raw);
            conmat = spdiags(zeros(size(conmat,1),1), 0, conmat);
            conmat = conmat ./ norm_val;
            
            P(connCount).full = conmat;
        end
        
        if filesLoaded > 0
            fprintf('  Quarter %s complete: %d matrices loaded\n\n', quarters{q}, filesLoaded);
        else
            fprintf('  Quarter %s: No files in date range\n\n', quarters{q});
        end
    end

    fprintf('Sorting %d matrices by date...\n', connCount);
    [~, sortIdx] = sort([P.Date]);
    P = P(sortIdx);

    fprintf('Saving processed connectivity data to %s...\n', cacheFilename);
    save(cacheFilePath, 'P', '-v7.3')
    
    fprintf('\n========================================\n');
    fprintf('COMPLETE: Loaded %d connectivity matrices\n', connCount);
    if ~isempty(DATE_RANGE)
        fprintf('Date range: %s to %s\n', datestr(DATE_RANGE(1)), datestr(DATE_RANGE(2)));
    end
    fprintf('Cached to: %s\n', cacheFilename);
    fprintf('========================================\n\n');
    
else
    fprintf('========================================\n');
    fprintf('Loading cached connectivity matrices from %s...\n', cacheFilename);
    
    if exist(cacheFilePath, 'file') == 0
        error(['Cache file not found: %s\n' ...
               'Set RECALCULATE_CONNECTIVITY = true to generate it.'], cacheFilePath);
    end
    
    load(cacheFilePath, 'P');
    fprintf('COMPLETE: Loaded %d connectivity matrices from cache\n', length(P));
    if ~isempty(DATE_RANGE)
        fprintf('Date range: %s to %s\n', ...
                datestr(min([P.Date])), datestr(max([P.Date])));
    end
    fprintf('========================================\n\n');
end

dates = [P.Date];
ref = datetime(2019,1,1);  
Pdays = days(dates - ref);

%% Initial conditions
SLS1 = CLS1;
SMS1 = CMS1;
SHS1 = CHS1;

ILS1 = zeros(habs,1);
IMS1 = zeros(habs,1);
IHS1 = zeros(habs,1);

Flat = [find(reefData.unique_ID==29088) find(reefData.unique_ID==29338) ...
        find(reefData.unique_ID==29089) find(reefData.unique_ID==29339) ...
        find(reefData.unique_ID==29087)];

IHS1(Flat) = .01*CHS1(Flat);
SHS1(Flat) = SHS1(Flat)-IHS1(Flat);
IMS1(Flat) = .01*CMS1(Flat);
SMS1(Flat) = SMS1(Flat)-IMS1(Flat);
ILS1(Flat) = .01*CLS1(Flat);
SLS1(Flat) = SLS1(Flat)-ILS1(Flat);

RLS1 = zeros(habs,1);
RMS1 = zeros(habs,1);
RHS1 = zeros(habs,1);

%% Parameters
bls = 0.03;
bms = 0.14;
bhs = 2.08;

kls = .05;
kms = .55;
khs = 3.33;

% UPDATED: I0_frac replaces I0 - now a fraction (e.g., 0.01 = 1%)
% I0_frac = 0.01;  % 1% of class population triggers site-wide activation
% I0_frac = 0.0001;  % 1% of class population triggers site-wide activation
I0_frac = 0;  % 1% of class population triggers site-wide activation
% tau = 0.0001;    % tau is no longer I0-dependent - set small for sharp transition
tau = 1e-10;    % tau is no longer I0-dependent - set small for sharp transition

pars_simp = [bls; bms; bhs; kls; kms; khs];

%% Parameter sweep - EXPANDED to test connectivity variants
threshvec = [0, 0.0003];
shapeParamVec = [0.001, -4];
% c = 1;
c = 100;

Y0 = [SLS1;SMS1;SHS1;ILS1;IMS1;IHS1;RLS1;RMS1;RHS1]; 
save(fullfile(tempPath, 'inits_01.mat'), 'Y0')

tspan = [1 90];

%% DIAGNOSTIC RUN: Analyze T values with sample parameters
fprintf('\n========================================\n');
fprintf('RUNNING DIAGNOSTIC TO ANALYZE T VALUES\n');
fprintf('========================================\n\n');

test_thresh = threshvec(1);
test_shape = shapeParamVec(1);

fprintf('Testing with: thresh=%.4f, shapeParam=%.3f\n', test_thresh, test_shape);
fprintf('Running single simulation to capture T statistics...\n\n');

Y0_diag = [SLS1;SMS1;SHS1;ILS1;IMS1;IHS1;RLS1;RMS1;RHS1];

opts_diag = odeset('OutputFcn', @odeWaitbar);

tic
[t_diag, Y_diag] = ode45(@(t,Y) ODEfun_SCTLD_seascape(t,Y,CLS1,CMS1,CHS1,...
                                bls,bms,bhs,kls,kms,khs,...
                                test_thresh,c,test_shape,...
                                habs,P,Pdays,I0_frac,tau), ...
                         tspan, Y0_diag, opts_diag);
diag_time = toc;

fprintf('\nDiagnostic run complete (%.1f seconds).\n', diag_time);
fprintf('========================================\n\n');

LSHP_diag = Y_diag(:,1:habs);
MSHP_diag = Y_diag(:,habs+1:habs*2);
HSHP_diag = Y_diag(:,2*habs+1:habs*3);

LSIP_diag = Y_diag(:,3*habs+1:habs*4);
MSIP_diag = Y_diag(:,4*habs+1:habs*5);
HSIP_diag = Y_diag(:,5*habs+1:habs*6);

LSRP_diag = Y_diag(:,6*habs+1:habs*7);
MSRP_diag = Y_diag(:,7*habs+1:habs*8);
HSRP_diag = Y_diag(:,8*habs+1:habs*9);

TIP_diag = LSIP_diag + MSIP_diag + HSIP_diag;
TRP_diag = LSRP_diag + MSRP_diag + HSRP_diag;

fig_diag = figure('Position', [50 50 1400 600]);
T_diag = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(T_diag, sprintf('Diagnostic Run Results (Day %d) - thresh=%.4f, shape=%.3f', ...
                      tspan(end), test_thresh, test_shape), ...
      'FontSize', 14, 'FontWeight', 'bold');

nexttile(T_diag, 1);
scatter(XY(:,1), XY(:,2), 20, TIP_diag(end,:)'*100, 'filled');
c1 = colorbar;
c1.Label.String = 'Infected Cover (%)';
colormap(gca, hot);
clim([0, max(TIP_diag(end,:)*100)]);
title(sprintf('Infected Cover (Day %d)', tspan(end)));
xlabel('Longitude');
ylabel('Latitude');
axis equal tight;

nexttile(T_diag, 2);
scatter(XY(:,1), XY(:,2), 20, TRP_diag(end,:)'*100, 'filled');
c2 = colorbar;
c2.Label.String = 'Removed Cover (%)';
colormap(gca, hot);
clim([0, max(TRP_diag(end,:)*100)]);
title(sprintf('Removed Cover (Day %d)', tspan(end)));
xlabel('Longitude');
ylabel('Latitude');
axis equal tight;

saveas(fig_diag, fullfile(seascapePath, 'Diagnostic_Spatial_Maps.png'));
fprintf('Saved diagnostic spatial maps to: %s\n', fullfile(seascapePath, 'Diagnostic_Spatial_Maps.png'));

%% DIAGNOSTIC: Check threshold activation at seeded sites
fprintf('\n========================================\n');
fprintf('THRESHOLD ACTIVATION DIAGNOSTIC\n');
fprintf('========================================\n');

fprintf('\nInitial seeding at Flat Cay sites:\n');
for f = 1:length(Flat)
    site = Flat(f);
    
    % Calculate per-class thresholds
    I0_LS = I0_frac * CLS1(site);
    I0_MS = I0_frac * CMS1(site);
    I0_HS = I0_frac * CHS1(site);
    
    % Site minimum threshold
    thresholds = [I0_LS, I0_MS, I0_HS];
    thresholds = thresholds(thresholds > 0);
    site_min_threshold = min(thresholds);
    
    % Total infected
    site_total_I = ILS1(site) + IMS1(site) + IHS1(site);
    
    % Calculate activation
    tau_site = site_min_threshold / 10;
    beta_prime_test = (bls/2) * (1 + tanh((site_total_I - site_min_threshold) / tau_site));
    
    fprintf('  Site %d:\n', site);
    fprintf('    Class thresholds: LS=%.6f MS=%.6f HS=%.6f\n', I0_LS, I0_MS, I0_HS);
    fprintf('    Site min threshold: %.6f (%.4f%%)\n', site_min_threshold, site_min_threshold*100);
    fprintf('    Total infected: %.6f (%.4f%%)\n', site_total_I, site_total_I*100);
    fprintf('    Beta_prime/beta: %.4f (%.1f%% of full strength)\n', ...
            beta_prime_test/bls, 100*beta_prime_test/bls);
    fprintf('    Cover: LS=%.4f MS=%.4f HS=%.4f Total=%.4f\n', ...
            CLS1(site), CMS1(site), CHS1(site), CC(site));
end

fprintf('\n========================================\n\n');






%% DIAGNOSTIC: Check T (connectivity) values at seeded and nearby sites
fprintf('\n========================================\n');
fprintf('CONNECTIVITY (T) DIAGNOSTIC\n');
fprintf('========================================\n\n');

% Get initial disease pool at seeded sites
DP_init = ILS1 + IMS1 + IHS1;

% Calculate T for all sites using first connectivity matrix
Pint_test = P(1).full;
T_init = incomingRisk_sparse(Pint_test, DP_init, test_thresh, 'rowsAreSources');
T_init_raw = T_init;  % Before transformation
T_init = c.*(1-exp(-test_shape.*T_init))/(1-exp(-test_shape));  % After transformation

fprintf('Seeded sites (Flat Cay) - Outgoing connectivity:\n');
for f = 1:length(Flat)
    site = Flat(f);
    % How much does this site SEND to others?
    outgoing = full(sum(Pint_test(site,:)));
    fprintf('  Site %d: DP=%.6f, sends total prob=%.6f to network\n', ...
            site, DP_init(site), outgoing);
end

fprintf('\nNearby sites - Incoming connectivity (T values):\n');
% Find sites within some radius of Flat sites
flat_coords = XY(Flat,:);
for i = 1:habs
    if ismember(i, Flat), continue; end  % Skip seeded sites
    
    dist_to_flat = min(sqrt(sum((XY(i,:) - flat_coords).^2, 2)));
    
    if dist_to_flat < 0.1  % Within ~10km (adjust as needed)
        fprintf('  Site %d (dist=%.3f): T_raw=%.8f, T_transformed=%.8f, DP_local=%.6f\n', ...
                i, dist_to_flat, T_init_raw(i), T_init(i), DP_init(i));
    end
end

fprintf('\nSummary statistics:\n');
fprintf('  T_raw: max=%.8f, mean(nonzero)=%.8f\n', ...
        max(T_init_raw), mean(T_init_raw(T_init_raw>0)));
fprintf('  T_transformed: max=%.8f, mean(nonzero)=%.8f\n', ...
        max(T_init), mean(T_init(T_init>0)));
fprintf('  Sites with T>0: %d / %d\n', sum(T_init>0), habs);

fprintf('========================================\n\n');













pause(2);  % Give time to read output before sweep starts
max_conn_day = round(days(max([P.Date]) - datetime(2019,1,1)));



if tspan(end) > max_conn_day
    fprintf('\n*** TSPAN ADJUSTMENT ***\n');
    fprintf('Original tspan requested: [%d %d]\n', tspan(1), tspan(end));
    fprintf('Available connectivity data: up to day %d (date: %s)\n', ...
            max_conn_day, datestr(max([P.Date]), 'dd-mmm-yyyy'));
    fprintf('Adjusting tspan to: [%d %d]\n', tspan(1), max_conn_day);
    fprintf('************************\n\n');
    tspan(end) = max_conn_day;
end

tspan_final = tspan;

if RUN_PARAMETER_SWEEP
    count = 0;
    Results = struct();

    fprintf('\n========================================\n');
    fprintf('PARAMETER SWEEP: Testing %d combinations\n', length(threshvec)*length(shapeParamVec));
    fprintf('  Thresholds: %d values [%s]\n', length(threshvec), num2str(threshvec));
    fprintf('  ShapeParams: %d values [%s]\n', length(shapeParamVec), num2str(shapeParamVec));
    fprintf('  I0_frac (local threshold): %.4f (%.2f%%)\n', I0_frac, I0_frac*100);
    fprintf('  tau (transition): %.6f\n', tau);
    if USE_PARALLEL
        fprintf('  Mode: PARALLEL (using %d workers)\n', min(6, feature('numcores')));
    else
        fprintf('  Mode: SERIAL\n');
    end
    fprintf('========================================\n\n');

    if USE_PARALLEL
        if isempty(gcp('nocreate'))
            parpool('local', min(6, feature('numcores')));
        end
        
        total_combos = length(threshvec)*length(shapeParamVec);
        Results(total_combos).thresh = [];
        
        tic
        parfor combo = 1:total_combos
            [tv, sv] = ind2sub([length(threshvec), length(shapeParamVec)], combo);
            
            fprintf('  [%d/%d] Running: thresh=%.4f, shapeParam=%.3f\n', ...
                    combo, total_combos, threshvec(tv), shapeParamVec(sv));
            
            Y0_par = [SLS1;SMS1;SHS1;ILS1;IMS1;IHS1;RLS1;RMS1;RHS1];
            tspan_par = tspan_final;
            
            [t_par,Y_par] = ode45(@(t,Y) ODEfun_SCTLD_seascape(t,Y,CLS1,CMS1,CHS1,bls,bms,bhs,kls,kms,khs,threshvec(tv),c,shapeParamVec(sv),habs,P,Pdays,I0_frac,tau), tspan_par, Y0_par);
            
            LSHP = Y_par(:,1:habs);
            MSHP = Y_par(:,habs+1:habs*2);
            HSHP = Y_par(:,2*habs+1:habs*3);
            
            LSIP = Y_par(:,3*habs+1:habs*4);
            MSIP = Y_par(:,4*habs+1:habs*5);
            HSIP = Y_par(:,5*habs+1:habs*6);
            
            LSRP = Y_par(:,6*habs+1:habs*7);
            MSRP = Y_par(:,7*habs+1:habs*8);
            HSRP = Y_par(:,8*habs+1:habs*9);
            
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
            Results(combo).shapeParam = shapeParamVec(sv);
            Results(combo).I0_frac = I0_frac;
            Results(combo).tau = tau;
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
        count = length(Results);
        
    else
        tic
        for tv = 1:length(threshvec)
            fprintf('--- Threshold %d/%d (%.4f) ---\n', tv, length(threshvec), threshvec(tv));
            
            for sv = 1:length(shapeParamVec)
                count = count+1;
                fprintf('  [%d/%d] Running: thresh=%.4f, shapeParam=%.3f ... ', ...
                        count, length(threshvec)*length(shapeParamVec), threshvec(tv), shapeParamVec(sv));
                
                run_start = tic;

                Y0 = [SLS1;SMS1;SHS1;ILS1;IMS1;IHS1;RLS1;RMS1;RHS1]; 
                tspan_sweep = tspan_final;
                clear Y t
                
                [t,Y] = ode45(@(t,Y) ODEfun_SCTLD_seascape(t,Y,CLS1,CMS1,CHS1,bls,bms,bhs,kls,kms,khs,threshvec(tv),c,shapeParamVec(sv),habs,P,Pdays,I0_frac,tau), tspan_sweep, Y0);
                
                fprintf('Done (%.1f sec)\n', toc(run_start));

                LSHP = Y(:,1:habs);
                MSHP = Y(:,habs+1:habs*2);
                HSHP = Y(:,2*habs+1:habs*3);
                
                LSIP = Y(:,3*habs+1:habs*4);
                MSIP = Y(:,4*habs+1:habs*5);
                HSIP = Y(:,5*habs+1:habs*6);
                
                LSRP = Y(:,6*habs+1:habs*7);
                MSRP = Y(:,7*habs+1:habs*8);
                HSRP = Y(:,8*habs+1:habs*9);
                
                interpLSHP = max(0, interp1(t,LSHP,tspan_sweep(1):tspan_sweep(2)));
                interpMSHP = max(0, interp1(t,MSHP,tspan_sweep(1):tspan_sweep(2)));
                interpHSHP = max(0, interp1(t,HSHP,tspan_sweep(1):tspan_sweep(2)));
                
                interpLSIP = max(0, interp1(t,LSIP,tspan_sweep(1):tspan_sweep(2)));
                interpMSIP = max(0, interp1(t,MSIP,tspan_sweep(1):tspan_sweep(2)));
                interpHSIP = max(0, interp1(t,HSIP,tspan_sweep(1):tspan_sweep(2)));
                
                interpLSRP = max(0, interp1(t,LSRP,tspan_sweep(1):tspan_sweep(2)));
                interpMSRP = max(0, interp1(t,MSRP,tspan_sweep(1):tspan_sweep(2)));
                interpHSRP = max(0, interp1(t,HSRP,tspan_sweep(1):tspan_sweep(2)));

                Results(count).thresh = threshvec(tv);
                Results(count).shapeParam = shapeParamVec(sv);
                Results(count).I0_frac = I0_frac;
                Results(count).tau = tau;
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
    
    fprintf('Saving parameter sweep results to temp/ParameterSweep_Results_Enhanced.mat...\n');
    save(fullfile(tempPath, 'ParameterSweep_Results_Enhanced.mat'), 'Results', 'threshvec', 'shapeParamVec', 'I0_frac', 'tau', '-v7.3');
    fprintf('Saved successfully.\n\n');
    
    fprintf('Generating spatial maps for all parameter combinations...\n');
    fig_sweep_spatial = figure('Position', [50 50 1800 1000]);
    T_sweep = tiledlayout(length(threshvec), length(shapeParamVec)*2, ...
                          'TileSpacing', 'compact', 'Padding', 'compact');
    title(T_sweep, sprintf('Parameter Sweep Spatial Results (Day %d)', tspan_final(end)), ...
          'FontSize', 14, 'FontWeight', 'bold');
    
    max_I_global = 0;
    max_R_global = 0;
    for i = 1:length(Results)
        TIP = Results(i).LSI + Results(i).MSI + Results(i).HSI;
        TRP = Results(i).LSR + Results(i).MSR + Results(i).HSR;
        max_I_global = max(max_I_global, max(TIP(end,:)));
        max_R_global = max(max_R_global, max(TRP(end,:)));
    end
    
    for i = 1:length(Results)
        TIP = Results(i).LSI + Results(i).MSI + Results(i).HSI;
        TRP = Results(i).LSR + Results(i).MSR + Results(i).HSR;
        
        nexttile(T_sweep, (i-1)*2 + 1);
        scatter(XY(:,1), XY(:,2), 8, TIP(end,:)'*100, 'filled');
        colormap(gca, hot);
        clim([0, max_I_global*100]);
        c_i = colorbar;
        c_i.Label.String = 'I (%)';
        c_i.FontSize = 7;
        title(sprintf('Infected: T=%.4f S=%.2f', Results(i).thresh, Results(i).shapeParam), ...
              'FontSize', 8);
        axis equal tight off;
        
        nexttile(T_sweep, (i-1)*2 + 2);
        scatter(XY(:,1), XY(:,2), 8, TRP(end,:)'*100, 'filled');
        colormap(gca, hot);
        clim([0, max_R_global*100]);
        c_r = colorbar;
        c_r.Label.String = 'R (%)';
        c_r.FontSize = 7;
        title(sprintf('Removed: T=%.4f S=%.2f', Results(i).thresh, Results(i).shapeParam), ...
              'FontSize', 8);
        axis equal tight off;
    end
    
    saveas(fig_sweep_spatial, fullfile(seascapePath, 'ParameterSweep_Spatial_Maps.png'));
    fprintf('Saved parameter sweep spatial maps to: %s\n\n', ...
            fullfile(seascapePath, 'ParameterSweep_Spatial_Maps.png'));
    
else
    fprintf('\n========================================\n');
    fprintf('Loading cached parameter sweep results from temp/ParameterSweep_Results_Enhanced.mat...\n');
    
    if exist(fullfile(tempPath, 'ParameterSweep_Results_Enhanced.mat'), 'file') == 0
        error(['Parameter sweep cache file not found.\n' ...
               'Set RUN_PARAMETER_SWEEP = true to generate it.']);
    end
    
    load(fullfile(tempPath, 'ParameterSweep_Results_Enhanced.mat'), 'Results', 'threshvec', 'shapeParamVec', 'I0_frac', 'tau');
    fprintf('COMPLETE: Loaded %d parameter combinations from cache\n', length(Results));
    fprintf('========================================\n\n');
end

%% Analyze and visualize results
fprintf('\n========================================\n');
fprintf('ANALYZING RESULTS\n');
fprintf('========================================\n\n');

for i = 1:length(Results)
    TIP = Results(i).LSI + Results(i).MSI + Results(i).HSI;
    TRP = Results(i).LSR + Results(i).MSR + Results(i).HSR;
    
    infected_sites = find(max(TIP, [], 1) > 0);
    Results(i).n_infected = length(infected_sites);
    Results(i).infected_sites = infected_sites;
    
    Results(i).total_mortality = sum(TRP(end,:));
    Results(i).mean_mortality = mean(TRP(end,:));
    
    fprintf('Combo %d: thresh=%.4f, shape=%.3f\n', i, Results(i).thresh, Results(i).shapeParam);
    fprintf('  Sites infected: %d\n', Results(i).n_infected);
    fprintf('  Total mortality: %.4f\n', Results(i).total_mortality);
    fprintf('  Mean mortality: %.6f\n\n', Results(i).mean_mortality);
end

%% Plot multi-group SIR curves for selected sites
combo_idx = 1;
TIP = Results(combo_idx).LSI + Results(combo_idx).MSI + Results(combo_idx).HSI;
TRP = Results(combo_idx).LSR + Results(combo_idx).MSR + Results(combo_idx).HSR;

significant_sites = find(TRP(end,:) > 0.001);
[~, sort_idx] = sort(TRP(end,significant_sites), 'descend');
top_sites = significant_sites(sort_idx(1:min(9, length(significant_sites))));

if ~isempty(top_sites)
    fig1 = figure('Position', [100 100 1600 1000]);
    T1 = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(T1, sprintf('Multi-Group SIR Dynamics - Combo %d (thresh=%.4f, shape=%.3f)', ...
                      combo_idx, Results(combo_idx).thresh, Results(combo_idx).shapeParam), ...
          'FontSize', 14, 'FontWeight', 'bold');
    
    days_vec = tspan_final(1):tspan_final(2);
    
    for idx = 1:length(top_sites)
        site = top_sites(idx);
        
        nexttile(T1, idx);
        hold on;
        
        plot(days_vec, Results(combo_idx).LSI(:,site)*100, 'b-', 'LineWidth', 1.5, 'DisplayName', 'LS Infected');
        plot(days_vec, Results(combo_idx).MSI(:,site)*100, 'Color', [1 0.5 0], 'LineWidth', 1.5, 'DisplayName', 'MS Infected');
        plot(days_vec, Results(combo_idx).HSI(:,site)*100, 'r-', 'LineWidth', 1.5, 'DisplayName', 'HS Infected');
        
        yline(CC(site)*100, 'k--', 'LineWidth', 1, 'DisplayName', 'Total Cover', 'Alpha', 0.3);
        
        xlabel('Days');
        ylabel('Infected Cover (%)');
        title(sprintf('Site %d\nLS:%.2f%% MS:%.2f%% HS:%.2f%%', ...
                      site, CLS1(site)*100, CMS1(site)*100, CHS1(site)*100), 'FontSize', 9);
        legend('Location', 'best', 'FontSize', 7);
        grid on;
        hold off;
    end
    
    saveas(fig1, fullfile(seascapePath, sprintf('MultiGroup_SIR_Combo%d.png', combo_idx)));
    
    fig2 = figure('Position', [150 150 1600 1000]);
    T2 = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(T2, sprintf('Combined SIR Dynamics - Combo %d (thresh=%.4f, shape=%.3f)', ...
                      combo_idx, Results(combo_idx).thresh, Results(combo_idx).shapeParam), ...
          'FontSize', 14, 'FontWeight', 'bold');
    
    for idx = 1:length(top_sites)
        site = top_sites(idx);
        
        nexttile(T2, idx);
        hold on;
        
        S_total = (Results(combo_idx).LSS(:,site) + Results(combo_idx).MSS(:,site) + Results(combo_idx).HSS(:,site)) * 100;
        I_total = (Results(combo_idx).LSI(:,site) + Results(combo_idx).MSI(:,site) + Results(combo_idx).HSI(:,site)) * 100;
        R_total = (Results(combo_idx).LSR(:,site) + Results(combo_idx).MSR(:,site) + Results(combo_idx).HSR(:,site)) * 100;
        
        plot(days_vec, S_total, 'b-', 'LineWidth', 2, 'DisplayName', 'Susceptible');
        plot(days_vec, I_total, 'r-', 'LineWidth', 2, 'DisplayName', 'Infected');
        plot(days_vec, R_total, 'k-', 'LineWidth', 2, 'DisplayName', 'Removed');
        
        xlabel('Days');
        ylabel('Cover (%)');
        title(sprintf('Site %d - Total Cover: %.2f%%', site, CC(site)*100), 'FontSize', 9);
        legend('Location', 'best', 'FontSize', 8);
        grid on;
        hold off;
    end
    
    saveas(fig2, fullfile(seascapePath, sprintf('Combined_SIR_Combo%d.png', combo_idx)));
    
    fprintf('Generated SIR curve plots for %d sites\n', length(top_sites));
else
    fprintf('No sites with significant infection found for combo %d\n', combo_idx);
end

fprintf('\n========================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('Figures saved to: %s\n', seascapePath);
fprintf('========================================\n\n');