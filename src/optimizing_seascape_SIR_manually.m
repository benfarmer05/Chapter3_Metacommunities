%% Connectivity Parameter Exploration Test - ENHANCED
% Tests different threshold and shapeParam combinations with REAL connectivity
% to understand how these parameters affect disease spread dynamics
%
% Prerequisites: Must have already loaded P, reefData, and set initial conditions
% This script assumes the main script environment is already loaded

fprintf('\n========================================\n');
fprintf('CONNECTIVITY PARAMETER EXPLORATION TEST\n');
fprintf('========================================\n\n');

%% Test Parameters
% Focus on a narrow range around your current values to see sensitivity

% Threshold: minimum infected cover for site to transmit (absolute cover units)
% Your current: 0.0003
thresh_test = [0.0001, 0.0003, 0.0005, 0.001];

% ShapeParam: controls transformation of flux T to infection probability
% Negative = concave (diminishing returns), Positive = convex (accelerating)
% Your current: -4
shapeParam_test = [-10, -4, -1, 0.001, 1, 4];

% Keep c constant at 1 (max external infection probability)
c_test = 1;

fprintf('Testing %d parameter combinations:\n', length(thresh_test) * length(shapeParam_test));
fprintf('  Thresholds: [%s]\n', num2str(thresh_test));
fprintf('  ShapeParams: [%s]\n', num2str(shapeParam_test));
fprintf('  Time span: %d to %d days\n', tspan_final(1), tspan_final(2));
fprintf('\n');

%% Setup for test runs
% Use same initial conditions as main script
Y0_test = [SLS1; SMS1; SHS1; ILS1; IMS1; IHS1; RLS1; RMS1; RHS1];

% Parallel execution toggle
USE_PARALLEL = true;  % Set to false for serial execution

% NEW: Toggle to exclude sites with zero cover in ANY susceptibility class
EXCLUDE_INCOMPLETE_SITES = false;  % Set to true to exclude sites missing LS, MS, or HS

if EXCLUDE_INCOMPLETE_SITES
    % Identify sites that have all three susceptibility classes present
    sites_with_all_classes = (CLS1 > 0) & (CMS1 > 0) & (CHS1 > 0);
    n_complete_sites = sum(sites_with_all_classes);
    n_excluded_sites = habs - n_complete_sites;
    
    fprintf('\n*** INCOMPLETE SITE EXCLUSION ENABLED ***\n');
    fprintf('Sites with all classes (LS, MS, HS > 0): %d\n', n_complete_sites);
    fprintf('Sites excluded (missing one or more classes): %d\n', n_excluded_sites);
    fprintf('  - Missing LS: %d sites\n', sum(CLS1 == 0));
    fprintf('  - Missing MS: %d sites\n', sum(CMS1 == 0));
    fprintf('  - Missing HS: %d sites\n', sum(CHS1 == 0));
    fprintf('*****************************************\n\n');
    
    % Zero out initial conditions for incomplete sites
    % This effectively removes them from the simulation
    zero_mask = ~sites_with_all_classes;
    
    Y0_test_orig = Y0_test;  % Store original
    
    % Zero out S, I, R for all classes at incomplete sites
    Y0_test(zero_mask) = 0;                           % LS_S
    Y0_test(habs + find(zero_mask)) = 0;              % MS_S
    Y0_test(2*habs + find(zero_mask)) = 0;            % HS_S
    Y0_test(3*habs + find(zero_mask)) = 0;            % LS_I
    Y0_test(4*habs + find(zero_mask)) = 0;            % MS_I
    Y0_test(5*habs + find(zero_mask)) = 0;            % HS_I
    Y0_test(6*habs + find(zero_mask)) = 0;            % LS_R
    Y0_test(7*habs + find(zero_mask)) = 0;            % MS_R
    Y0_test(8*habs + find(zero_mask)) = 0;            % HS_R
    
else
    fprintf('\n*** Using ALL sites (including those with incomplete cover) ***\n\n');
    sites_with_all_classes = true(habs, 1);  % All sites included
end

%% Run parameter combinations
total_combos = length(thresh_test) * length(shapeParam_test);

if USE_PARALLEL
    fprintf('Running simulations in PARALLEL mode...\n');
    fprintf('  Using up to %d workers\n\n', min(6, feature('numcores')));
    
    % Start parallel pool if not already running
    if isempty(gcp('nocreate'))
        parpool('local', min(6, feature('numcores')));
    end
    
    % Pre-allocate struct array
    ConnTest(total_combos).thresh = [];
    
    tic
    parfor combo = 1:total_combos
        [ti, si] = ind2sub([length(thresh_test), length(shapeParam_test)], combo);
        
        curr_thresh = thresh_test(ti);
        curr_shape = shapeParam_test(si);
        
        fprintf('[%d/%d] Running: thresh=%.4f, shapeParam=%.1f\n', ...
                combo, total_combos, curr_thresh, curr_shape);
        
        % Run ODE with current parameters
        [t_test, Y_test] = ode45(@(t,Y) ODEfun_SCTLD_seascape(t, Y, ...
                                                               CLS1, CMS1, CHS1, ...
                                                               bls, bms, bhs, ...
                                                               kls, kms, khs, ...
                                                               curr_thresh, c_test, curr_shape, ...
                                                               habs, P, Pdays), ...
                                 tspan_final, Y0_test);
        
        % Extract results
        LS_S = Y_test(:, 1:habs);
        MS_S = Y_test(:, habs+1:habs*2);
        HS_S = Y_test(:, 2*habs+1:habs*3);
        
        LS_I = Y_test(:, 3*habs+1:habs*4);
        MS_I = Y_test(:, 4*habs+1:habs*5);
        HS_I = Y_test(:, 5*habs+1:habs*6);
        
        LS_R = Y_test(:, 6*habs+1:habs*7);
        MS_R = Y_test(:, 7*habs+1:habs*8);
        HS_R = Y_test(:, 8*habs+1:habs*9);
        
        % Interpolate to daily values
        days_vec = tspan_final(1):tspan_final(2);
        
        LS_S_interp = max(0, interp1(t_test, LS_S, days_vec));
        MS_S_interp = max(0, interp1(t_test, MS_S, days_vec));
        HS_S_interp = max(0, interp1(t_test, HS_S, days_vec));
        
        LS_I_interp = max(0, interp1(t_test, LS_I, days_vec));
        MS_I_interp = max(0, interp1(t_test, MS_I, days_vec));
        HS_I_interp = max(0, interp1(t_test, HS_I, days_vec));
        
        LS_R_interp = max(0, interp1(t_test, LS_R, days_vec));
        MS_R_interp = max(0, interp1(t_test, MS_R, days_vec));
        HS_R_interp = max(0, interp1(t_test, HS_R, days_vec));
        
        % Store results
        ConnTest(combo).thresh = curr_thresh;
        ConnTest(combo).shapeParam = curr_shape;
        ConnTest(combo).c = c_test;
        
        ConnTest(combo).LS_S = LS_S_interp;
        ConnTest(combo).MS_S = MS_S_interp;
        ConnTest(combo).HS_S = HS_S_interp;
        
        ConnTest(combo).LS_I = LS_I_interp;
        ConnTest(combo).MS_I = MS_I_interp;
        ConnTest(combo).HS_I = HS_I_interp;
        
        ConnTest(combo).LS_R = LS_R_interp;
        ConnTest(combo).MS_R = MS_R_interp;
        ConnTest(combo).HS_R = HS_R_interp;
        
        % Calculate totals
        ConnTest(combo).Total_I = LS_I_interp + MS_I_interp + HS_I_interp;
        ConnTest(combo).Total_S = LS_S_interp + MS_S_interp + HS_S_interp;
        ConnTest(combo).Total_R = LS_R_interp + MS_R_interp + HS_R_interp;
        
        % Summary statistics
        ConnTest(combo).total_infected_ever = sum(ConnTest(combo).Total_R(end,:) > 0);
        ConnTest(combo).mean_final_mortality = mean(ConnTest(combo).Total_R(end,:));
        ConnTest(combo).max_infected_ever = max(max(ConnTest(combo).Total_I));
        
        [~, peak_day_idx] = max(sum(ConnTest(combo).Total_I, 2));
        ConnTest(combo).peak_infection_day = days_vec(peak_day_idx);
        
        % Track which sites were included
        ConnTest(combo).sites_included = sites_with_all_classes;
    end
    total_time = toc;
    
    fprintf('\n========================================\n');
    fprintf('PARALLEL execution complete!\n');
    fprintf('  Total simulations: %d\n', total_combos);
    fprintf('  Total time: %.1f minutes (avg: %.1f sec per simulation)\n', ...
            total_time/60, total_time/total_combos);
    fprintf('========================================\n\n');
    
else
    % Serial execution
    fprintf('Running simulations in SERIAL mode...\n\n');
    
    ConnTest = struct();
    test_count = 0;
    
    tic
    for ti = 1:length(thresh_test)
        fprintf('--- Threshold %d/%d (%.4f) ---\n', ti, length(thresh_test), thresh_test(ti));
        
        for si = 1:length(shapeParam_test)
            test_count = test_count + 1;
            
            curr_thresh = thresh_test(ti);
            curr_shape = shapeParam_test(si);
            
            fprintf('  [%d/%d] Running: thresh=%.4f, shapeParam=%.1f ... ', ...
                    test_count, total_combos, curr_thresh, curr_shape);
            
            run_tic = tic;
            
            % Run ODE with current parameters
            [t_test, Y_test] = ode45(@(t,Y) ODEfun_SCTLD_seascape(t, Y, ...
                                                                   CLS1, CMS1, CHS1, ...
                                                                   bls, bms, bhs, ...
                                                                   kls, kms, khs, ...
                                                                   curr_thresh, c_test, curr_shape, ...
                                                                   habs, P, Pdays), ...
                                     tspan_final, Y0_test);
            
            fprintf('Done (%.1f sec)\n', toc(run_tic));
            
            % Extract results
            LS_S = Y_test(:, 1:habs);
            MS_S = Y_test(:, habs+1:habs*2);
            HS_S = Y_test(:, 2*habs+1:habs*3);
            
            LS_I = Y_test(:, 3*habs+1:habs*4);
            MS_I = Y_test(:, 4*habs+1:habs*5);
            HS_I = Y_test(:, 5*habs+1:habs*6);
            
            LS_R = Y_test(:, 6*habs+1:habs*7);
            MS_R = Y_test(:, 7*habs+1:habs*8);
            HS_R = Y_test(:, 8*habs+1:habs*9);
            
            % Interpolate to daily values
            days_vec = tspan_final(1):tspan_final(2);
            
            LS_S_interp = max(0, interp1(t_test, LS_S, days_vec));
            MS_S_interp = max(0, interp1(t_test, MS_S, days_vec));
            HS_S_interp = max(0, interp1(t_test, HS_S, days_vec));
            
            LS_I_interp = max(0, interp1(t_test, LS_I, days_vec));
            MS_I_interp = max(0, interp1(t_test, MS_I, days_vec));
            HS_I_interp = max(0, interp1(t_test, HS_I, days_vec));
            
            LS_R_interp = max(0, interp1(t_test, LS_R, days_vec));
            MS_R_interp = max(0, interp1(t_test, MS_R, days_vec));
            HS_R_interp = max(0, interp1(t_test, HS_R, days_vec));
            
            % Store results
            ConnTest(test_count).thresh = curr_thresh;
            ConnTest(test_count).shapeParam = curr_shape;
            ConnTest(test_count).c = c_test;
            
            ConnTest(test_count).LS_S = LS_S_interp;
            ConnTest(test_count).MS_S = MS_S_interp;
            ConnTest(test_count).HS_S = HS_S_interp;
            
            ConnTest(test_count).LS_I = LS_I_interp;
            ConnTest(test_count).MS_I = MS_I_interp;
            ConnTest(test_count).HS_I = HS_I_interp;
            
            ConnTest(test_count).LS_R = LS_R_interp;
            ConnTest(test_count).MS_R = MS_R_interp;
            ConnTest(test_count).HS_R = HS_R_interp;
            
            % Calculate totals
            ConnTest(test_count).Total_I = LS_I_interp + MS_I_interp + HS_I_interp;
            ConnTest(test_count).Total_S = LS_S_interp + MS_S_interp + HS_S_interp;
            ConnTest(test_count).Total_R = LS_R_interp + MS_R_interp + HS_R_interp;
            
            % Summary statistics
            ConnTest(test_count).total_infected_ever = sum(ConnTest(test_count).Total_R(end,:) > 0);
            ConnTest(test_count).mean_final_mortality = mean(ConnTest(test_count).Total_R(end,:));
            ConnTest(test_count).max_infected_ever = max(max(ConnTest(test_count).Total_I));
            
            [~, peak_day_idx] = max(sum(ConnTest(test_count).Total_I, 2));
            ConnTest(test_count).peak_infection_day = days_vec(peak_day_idx);
            
            % Track which sites were included
            ConnTest(test_count).sites_included = sites_with_all_classes;
        end
        fprintf('\n');
    end
    total_time = toc;
    
    fprintf('========================================\n');
    fprintf('SERIAL execution complete!\n');
    fprintf('  Total simulations: %d\n', test_count);
    fprintf('  Total time: %.1f minutes (avg: %.1f sec per simulation)\n', ...
            total_time/60, total_time/test_count);
    fprintf('========================================\n\n');
end

%% Analysis: Summary Statistics Table

fprintf('SUMMARY STATISTICS\n');
fprintf('==================\n');
if EXCLUDE_INCOMPLETE_SITES
    fprintf('NOTE: Statistics below exclude %d sites with incomplete cover\n', n_excluded_sites);
    fprintf('      (only sites with LS, MS, AND HS > 0 are included)\n');
end
fprintf('\n');
fprintf('Combo | Thresh  | Shape | Sites  | Peak Day | Peak I  | Final R\n');
fprintf('      |         | Param | Infect |          | (mean)  | (mean)\n');
fprintf('------|---------|-------|--------|----------|---------|--------\n');

for i = 1:length(ConnTest)
    fprintf('%5d | %.5f | %5.1f | %6d | %8d | %.5f | %.5f\n', ...
            i, ...
            ConnTest(i).thresh, ...
            ConnTest(i).shapeParam, ...
            ConnTest(i).total_infected_ever, ...
            ConnTest(i).peak_infection_day, ...
            ConnTest(i).max_infected_ever, ...
            ConnTest(i).mean_final_mortality);
end
fprintf('\n');

%% Visualization 1: Spatial spread patterns for each combination

fig1 = figure('Position', [50 50 1800 1200]);
T1 = tiledlayout(length(thresh_test), length(shapeParam_test), ...
                 'TileSpacing', 'compact', 'Padding', 'compact');
title(T1, 'Total Infection Over Time - Summed Across All Sites', ...
      'FontSize', 14, 'FontWeight', 'bold');

for i = 1:length(ConnTest)
    nexttile(T1, i);
    
    total_I_sum = sum(ConnTest(i).Total_I, 2);  % Sum across all sites
    
    plot(days_vec, total_I_sum, 'r-', 'LineWidth', 2);
    
    xlabel('Days');
    ylabel('Total Infected Cover');
    title(sprintf('T=%.4f, S=%.1f\n%d sites, peak day %d', ...
                  ConnTest(i).thresh, ConnTest(i).shapeParam, ...
                  ConnTest(i).total_infected_ever, ConnTest(i).peak_infection_day), ...
          'FontSize', 8);
    grid on;
end

%% Visualization 2: Spatial maps at peak infection day

fig2 = figure('Position', [100 100 1800 1200]);
T2 = tiledlayout(length(thresh_test), length(shapeParam_test), ...
                 'TileSpacing', 'compact', 'Padding', 'compact');
title(T2, 'Spatial Distribution at Peak Infection Day', ...
      'FontSize', 14, 'FontWeight', 'bold');

for i = 1:length(ConnTest)
    nexttile(T2, i);
    
    peak_day = ConnTest(i).peak_infection_day;
    peak_idx = find(days_vec == peak_day);
    
    infected_at_peak = ConnTest(i).Total_I(peak_idx, :);
    
    scatter(XY(:,1), XY(:,2), 20, infected_at_peak, 'filled');
    colormap(gca, hot);
    colorbar;
    clim([0, max(infected_at_peak)]);
    
    title(sprintf('T=%.4f, S=%.1f (Day %d)', ...
                  ConnTest(i).thresh, ConnTest(i).shapeParam, peak_day), ...
          'FontSize', 8);
    axis equal tight;
end

%% Visualization 3: Final mortality comparison

fig3 = figure('Position', [150 150 1400 600]);
T3 = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(T3, 'Final Disease Impact Comparison', 'FontSize', 14, 'FontWeight', 'bold');

% Panel 1: Number of sites infected
nexttile(T3, 1);
sites_infected = [ConnTest.total_infected_ever];
bar(1:length(ConnTest), sites_infected);
xlabel('Parameter Combination');
ylabel('Number of Sites Infected');
title(sprintf('Sites Experiencing Mortality (Total sites: %d)', habs));
xticks(1:length(ConnTest));
xticklabels(arrayfun(@(x) sprintf('%d', x), 1:length(ConnTest), 'UniformOutput', false));
xtickangle(45);
grid on;

% Panel 2: Mean final mortality across all sites
nexttile(T3, 2);
mean_mortality = [ConnTest.mean_final_mortality];
bar(1:length(ConnTest), mean_mortality);
xlabel('Parameter Combination');
ylabel('Mean Final Removed Cover');
title('Average Mortality Per Site');
xticks(1:length(ConnTest));
xticklabels(arrayfun(@(x) sprintf('%d', x), 1:length(ConnTest), 'UniformOutput', false));
xtickangle(45);
grid on;

%% Visualization 4: Thresh vs ShapeParam heatmaps

fig4 = figure('Position', [200 200 1400 500]);
T4 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(T4, 'Parameter Space Exploration: Thresh vs ShapeParam', ...
      'FontSize', 14, 'FontWeight', 'bold');

% Reshape results into matrices for heatmap
sites_infected_mat = reshape([ConnTest.total_infected_ever], ...
                              length(shapeParam_test), length(thresh_test))';
peak_day_mat = reshape([ConnTest.peak_infection_day], ...
                        length(shapeParam_test), length(thresh_test))';
mean_mortality_mat = reshape([ConnTest.mean_final_mortality], ...
                              length(shapeParam_test), length(thresh_test))';

% Heatmap 1: Sites infected
nexttile(T4, 1);
imagesc(shapeParam_test, thresh_test, sites_infected_mat);
colorbar;
xlabel('ShapeParam');
ylabel('Threshold');
title('Number of Sites Infected');
set(gca, 'YDir', 'normal');
colormap(gca, parula);

% Heatmap 2: Peak day
nexttile(T4, 2);
imagesc(shapeParam_test, thresh_test, peak_day_mat);
colorbar;
xlabel('ShapeParam');
ylabel('Threshold');
title('Peak Infection Day');
set(gca, 'YDir', 'normal');
colormap(gca, parula);

% Heatmap 3: Mean mortality
nexttile(T4, 3);
imagesc(shapeParam_test, thresh_test, mean_mortality_mat);
colorbar;
xlabel('ShapeParam');
ylabel('Threshold');
title('Mean Final Mortality');
set(gca, 'YDir', 'normal');
colormap(gca, hot);

%% NEW: Visualization 5: Spatial maps of percent coral cover removed

fig5 = figure('Position', [250 250 1800 1200]);
T5 = tiledlayout(length(thresh_test), length(shapeParam_test), ...
                 'TileSpacing', 'compact', 'Padding', 'compact');
if EXCLUDE_INCOMPLETE_SITES
    title(T5, 'Percent Coral Cover Removed at Each Site (Final State) - Complete Sites Only', ...
          'FontSize', 14, 'FontWeight', 'bold');
else
    title(T5, 'Percent Coral Cover Removed at Each Site (Final State)', ...
          'FontSize', 14, 'FontWeight', 'bold');
end

for i = 1:length(ConnTest)
    nexttile(T5, i);
    
    % Calculate % coral cover removed at each site
    final_removed = ConnTest(i).Total_R(end, :);
    pct_removed = 100 * final_removed ./ CC';
    
    % Handle sites with no initial cover
    pct_removed(CC == 0) = 0;
    
    % Calculate total removed cover across all sites
    if EXCLUDE_INCOMPLETE_SITES
        total_removed = sum(final_removed(sites_with_all_classes));
    else
        total_removed = sum(final_removed);
    end
    
    scatter(XY(:,1), XY(:,2), 20, pct_removed, 'filled');
    
    % If excluding incomplete sites, overlay excluded sites in gray
    if EXCLUDE_INCOMPLETE_SITES
        hold on;
        excluded_idx = ~sites_with_all_classes;
        scatter(XY(excluded_idx,1), XY(excluded_idx,2), 15, [0.7 0.7 0.7], 'x', 'LineWidth', 1);
        hold off;
    end
    
    colormap(gca, hot);
    c = colorbar;
    c.Label.String = '% Removed';
    clim([0, 100]);  % 0-100% removed
    
    title(sprintf('T=%.4f, S=%.1f\nTotal: %.3f', ...
                  ConnTest(i).thresh, ConnTest(i).shapeParam, total_removed), ...
          'FontSize', 8);
    axis equal tight;
    xlabel('Longitude');
    ylabel('Latitude');
end

%% NEW: Visualization 6: Individual site outbreak dynamics
% Select one parameter combination to examine in detail
% Use combination closest to your working parameters (thresh=0.0003, shapeParam=-4)
target_combo = find([ConnTest.thresh] == 0.0003 & [ConnTest.shapeParam] == -4);
if isempty(target_combo)
    % If exact match not found, use first combination
    target_combo = 1;
    fprintf('Warning: Target combination (thresh=0.0003, shapeParam=-4) not found. Using combo %d instead.\n', target_combo);
end

fprintf('\n========================================\n');
fprintf('Analyzing individual site dynamics for:\n');
fprintf('  Combination %d: thresh=%.4f, shapeParam=%.1f\n', ...
        target_combo, ConnTest(target_combo).thresh, ConnTest(target_combo).shapeParam);
fprintf('========================================\n\n');

% Calculate final loss for each site
final_loss = ConnTest(target_combo).Total_R(end, :)';
pct_loss_per_site = 100 * final_loss ./ CC;
pct_loss_per_site(CC == 0) = 0;  % Sites with no coral

% Sort sites by loss
[sorted_loss, loss_idx] = sort(final_loss, 'descend');

% Select sites to plot:
% - Top 5 highest loss
% - 3 sites with intermediate loss (around 33rd, 50th, 67th percentile)
% - Note: Only include sites with actual loss > 0 AND (if toggle enabled) complete cover

if EXCLUDE_INCOMPLETE_SITES
    sites_with_loss = find(final_loss > 0 & sites_with_all_classes);
else
    sites_with_loss = find(final_loss > 0);
end
n_sites_with_loss = length(sites_with_loss);

if n_sites_with_loss >= 8
    % High loss sites (top 5) - but only from eligible sites
    eligible_high_loss = loss_idx(ismember(loss_idx, sites_with_loss));
    high_loss_sites = eligible_high_loss(1:min(5, length(eligible_high_loss)));
    
    % Intermediate loss sites
    pct_33 = round(n_sites_with_loss * 0.33);
    pct_50 = round(n_sites_with_loss * 0.50);
    pct_67 = round(n_sites_with_loss * 0.67);
    
    % Find eligible sites at these percentiles
    eligible_loss_idx = loss_idx(ismember(loss_idx, sites_with_loss));
    mid_loss_sites = eligible_loss_idx([pct_33, pct_50, pct_67]);
    
    selected_sites = [high_loss_sites; mid_loss_sites];
else
    % If fewer sites, just take what we have
    eligible_loss_idx = loss_idx(ismember(loss_idx, sites_with_loss));
    selected_sites = eligible_loss_idx(1:min(8, length(eligible_loss_idx)));
end

n_selected = length(selected_sites);

fig6 = figure('Position', [300 300 1600 1000]);
T6 = tiledlayout(ceil(n_selected/3), 3, 'TileSpacing', 'compact', 'Padding', 'compact');
if EXCLUDE_INCOMPLETE_SITES
    title(T6, sprintf('Individual Site SIR Dynamics - Combo %d (T=%.4f, S=%.1f) - Complete Sites Only', ...
                      target_combo, ConnTest(target_combo).thresh, ConnTest(target_combo).shapeParam), ...
          'FontSize', 14, 'FontWeight', 'bold');
else
    title(T6, sprintf('Individual Site SIR Dynamics - Combo %d (T=%.4f, S=%.1f)', ...
                      target_combo, ConnTest(target_combo).thresh, ConnTest(target_combo).shapeParam), ...
          'FontSize', 14, 'FontWeight', 'bold');
end

for idx = 1:n_selected
    site = selected_sites(idx);
    
    nexttile(T6, idx);
    hold on;
    
    % Extract data for this site
    S_site = ConnTest(target_combo).Total_S(:, site);
    I_site = ConnTest(target_combo).Total_I(:, site);
    R_site = ConnTest(target_combo).Total_R(:, site);
    
    % Plot SIR curves
    plot(days_vec, S_site, 'b-', 'LineWidth', 2, 'DisplayName', 'S');
    plot(days_vec, I_site, 'r-', 'LineWidth', 2, 'DisplayName', 'I');
    plot(days_vec, R_site, 'k-', 'LineWidth', 2, 'DisplayName', 'R');
    
    % Add initial cover reference
    yline(CC(site), 'g--', 'LineWidth', 1, 'DisplayName', 'Initial', 'Alpha', 0.5);
    
    % Labels and formatting
    xlabel('Days');
    ylabel('Cover');
    
    % Determine loss category for title
    rank = find(loss_idx == site);
    if rank <= 5
        category = sprintf('HIGH LOSS (Rank %d)', rank);
    else
        category = sprintf('MID LOSS (Rank %d)', rank);
    end
    
    title(sprintf('Site %d - %s\n%.3f loss (%.1f%%)', ...
                  site, category, final_loss(site), pct_loss_per_site(site)), ...
          'FontSize', 9);
    
    legend('Location', 'best', 'FontSize', 7);
    grid on;
    ylim([0, CC(site) * 1.1]);
    hold off;
end

%% NEW: Visualization 7: Spatial location of selected sites

fig7 = figure('Position', [350 350 1000 800]);
hold on;

% Plot all sites in gray
scatter(XY(:,1), XY(:,2), 30, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.3);

% If excluding incomplete sites, mark them differently
if EXCLUDE_INCOMPLETE_SITES
    excluded_idx = ~sites_with_all_classes;
    scatter(XY(excluded_idx,1), XY(excluded_idx,2), 40, [0.5 0.5 0.5], 'x', 'LineWidth', 2);
end

% Plot selected sites in color, sized by loss
colors_selected = [repmat([1 0 0], 5, 1);      % Red for high loss
                   repmat([1 0.5 0], 3, 1)];    % Orange for mid loss
sizes_selected = 200 * (final_loss(selected_sites) / max(final_loss(selected_sites))) + 50;

for idx = 1:n_selected
    site = selected_sites(idx);
    scatter(XY(site,1), XY(site,2), sizes_selected(idx), ...
            colors_selected(idx,:), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    
    % Add site labels
    text(XY(site,1), XY(site,2), sprintf(' %d', site), ...
         'FontSize', 8, 'FontWeight', 'bold', 'Color', 'k');
end

xlabel('Longitude');
ylabel('Latitude');
if EXCLUDE_INCOMPLETE_SITES
    title(sprintf('Location of Selected Sites - Combo %d (T=%.4f, S=%.1f)\nComplete Sites Only (X = Excluded)', ...
                  target_combo, ConnTest(target_combo).thresh, ConnTest(target_combo).shapeParam), ...
          'FontSize', 12, 'FontWeight', 'bold');
    legend('All sites', 'Excluded (incomplete)', 'High loss', 'Mid loss', 'Location', 'best');
else
    title(sprintf('Location of Selected Sites - Combo %d (T=%.4f, S=%.1f)', ...
                  target_combo, ConnTest(target_combo).thresh, ConnTest(target_combo).shapeParam), ...
          'FontSize', 12, 'FontWeight', 'bold');
    legend('All sites', 'High loss', 'Mid loss', 'Location', 'best');
end
axis equal tight;
grid on;
hold off;

%% Print summary for selected sites
fprintf('SELECTED SITE DETAILS\n');
fprintf('=====================\n');
if EXCLUDE_INCOMPLETE_SITES
    fprintf('NOTE: Only showing sites with complete cover (LS, MS, AND HS > 0)\n');
end
fprintf('\n');
fprintf('Site | Rank | Category  | Final Loss | %% Loss | Initial Cover | LS/MS/HS\n');
fprintf('-----|------|-----------|------------|--------|---------------|----------\n');
for idx = 1:n_selected
    site = selected_sites(idx);
    rank = find(loss_idx == site);
    
    if rank <= 5
        category = 'HIGH';
    else
        category = 'MID';
    end
    
    fprintf('%4d | %4d | %9s | %10.4f | %6.1f%% | %13.4f | %.2f/%.2f/%.2f\n', ...
            site, rank, category, final_loss(site), pct_loss_per_site(site), CC(site), ...
            CLS1(site), CMS1(site), CHS1(site));
end
fprintf('\n');

%% Save results
fprintf('Saving connectivity test results...\n');
save(fullfile(tempPath, 'ConnectivityTest_Results.mat'), 'ConnTest', ...
     'thresh_test', 'shapeParam_test', 'c_test', '-v7.3');
fprintf('Results saved to: %s\n', fullfile(tempPath, 'ConnectivityTest_Results.mat'));

fprintf('\n========================================\n');
fprintf('CONNECTIVITY TEST COMPLETE\n');
fprintf('Generated %d figures:\n', 7);
fprintf('  Fig 1: Time series for all combinations\n');
fprintf('  Fig 2: Spatial maps at peak infection\n');
fprintf('  Fig 3: Final impact comparison\n');
fprintf('  Fig 4: Parameter space heatmaps\n');
fprintf('  Fig 5: Spatial maps of %% coral cover removed\n');
fprintf('  Fig 6: Individual site SIR dynamics\n');
fprintf('  Fig 7: Spatial location of selected sites\n');
fprintf('========================================\n\n');

%% Key Insights to Look For:
fprintf('KEY QUESTIONS TO INVESTIGATE:\n');
fprintf('------------------------------\n');
fprintf('1. Threshold effects:\n');
fprintf('   - Lower thresh → more sites can transmit sooner\n');
fprintf('   - Higher thresh → disease spreads only from established outbreaks\n\n');
fprintf('2. ShapeParam effects:\n');
fprintf('   - Negative (concave): Diminishing returns from multiple sources\n');
fprintf('   - Positive (convex): Accelerating effect from multiple sources\n');
fprintf('   - Near-zero (linear): Direct proportional relationship\n\n');
fprintf('3. Combined effects:\n');
fprintf('   - Which combination matches observed spatial patterns?\n');
fprintf('   - Does spatial spread timing align with real data?\n');
fprintf('   - Are "hotspot" sites realistically positioned?\n\n');
fprintf('4. External shutoff mechanism:\n');
fprintf('   - Sites should show rapid local growth after initial seeding\n');
fprintf('   - External influence mainly affects TIMING of arrival\n');
fprintf('   - Intensity driven by local community composition\n\n');
fprintf('5. Site-level dynamics (NEW):\n');
fprintf('   - Do high-loss sites show characteristic outbreak curves?\n');
fprintf('   - When does infection arrive at different sites?\n');
fprintf('   - Is mortality proportional to initial coral cover?\n');
fprintf('========================================\n\n');