%% Load saved outputs

% Setup paths (same as in main script)
projectPath = matlab.project.rootProject().RootFolder;
outputPath = fullfile(projectPath, 'output');
seascapePath = fullfile(outputPath, 'seascape_SIR');

workspace_filename = fullfile(seascapePath, 'seascape_SIR_workspace.mat');

fprintf('\n========================================\n');
fprintf('LOADING SAVED OUTPUTS\n');
fprintf('========================================\n');
fprintf('Loading from:\n  %s\n', workspace_filename);

if exist(workspace_filename, 'file') == 0
    error('Workspace file not found: %s', workspace_filename);
end

tic;
load(workspace_filename, 'outputs');
elapsed_load = toc;

% Get file size
file_info = dir(workspace_filename);
file_size_mb = file_info.bytes / (1024^2);

fprintf('Outputs loaded successfully!\n');
fprintf('  File size: %.1f MB\n', file_size_mb);
fprintf('  Load time: %.1f seconds\n', elapsed_load);
fprintf('  Simulation days: %d\n', outputs.metadata.num_days);
fprintf('  Number of sites: %d\n', outputs.metadata.num_sites);
fprintf('========================================\n\n');

% Extract variables for convenience
unique_IDs = outputs.sites.unique_IDs;
locations = outputs.sites.locations;
N_site = outputs.sites.N_site;
target_reef_IDs = outputs.sites.target_reef_IDs;

S_total_output_days = outputs.totals.S_total;
I_total_output_days = outputs.totals.I_total;
R_total_output_days = outputs.totals.R_total;

num_sites = outputs.metadata.num_sites;
num_days = outputs.metadata.num_days;

fprintf('Variables extracted and ready for plotting!\n\n');






%% Static figure of final state

fprintf('Creating static figure of final state...\n');

% Recreate colormaps (if not already done)
breakpoints = [0 .000000001 .001 .1 1];
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

% Get final timepoint and date
last_valid = num_days;
final_day = outputs.metadata.tspan(1) + last_valid - 1;
final_date = outputs.metadata.simulation_dates(end);

fprintf('Final timepoint: Day %d (%s)\n', final_day, string(final_date, 'dd-MMM-yyyy'));

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

% % Save figure
% figname = sprintf('FinalState_Diagnostic_day%d_ET%s_I0_%s_FS%.1f_FP%.2f', ...
%                   final_day, ...
%                   strrep(num2str(outputs.params.export_thresh, '%.8f'), '.', 'p'), ...
%                   strrep(num2str(outputs.params.I0, '%.10f'), '.', 'p'), ...
%                   outputs.params.flux_scale, ...
%                   outputs.params.flux_shape);
% saveas(fig_final, fullfile(seascapePath, [figname '.png']));
% saveas(fig_final, fullfile(seascapePath, [figname '.fig']));

% fprintf('Static figure saved as: %s\n', figname);

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