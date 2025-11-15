%% Read in observed disease data and create disease hull visualization

% Read the observed data CSV
observedDataFile = fullfile(outputPath, 'combined_coral_data.csv');
obsData = readtable(observedDataFile);

fprintf('========================================\n');
fprintf('OBSERVED DATA LOADED\n');
fprintf('Total observations: %d\n', height(obsData));
fprintf('Date range: %s to %s\n', ...
    string(min(obsData.date), 'dd-MMM-yyyy'), ...
    string(max(obsData.date), 'dd-MMM-yyyy'));
fprintf('Unique locations: %d\n', length(unique(obsData.location)));
fprintf('========================================\n\n');

% Get unique dates for temporal analysis
unique_dates = unique(obsData.date);
unique_dates = sort(unique_dates);

% Filter to 2018-2019
obsData = obsData(year(obsData.date) >= 2018 & year(obsData.date) <= 2019, :);
fprintf('Filtered to 2018-2019: %d observations\n\n', height(obsData));

% Group by month
obsData.year_month = datetime(year(obsData.date), month(obsData.date), 1);
unique_months = unique(obsData.year_month);
unique_months = sort(unique_months);

fprintf('Monthly observation summary (2018-2019):\n');
for i = 1:length(unique_months)
    n_obs = sum(obsData.year_month == unique_months(i));
    n_present = sum(obsData.year_month == unique_months(i) & strcmp(obsData.presence, 'P'));
    fprintf('  %s: %d observations (%d disease-present)\n', ...
            string(unique_months(i), 'MMM-yyyy'), n_obs, n_present);
end
fprintf('\n');

% Create figure showing disease presence over time (monthly)
figure('Name', 'Observed Disease Spread - Monthly 2018-2019', 'Position', [100 100 1400 900]);

% 12 months max = 3x4 grid
n_months = length(unique_months);
n_cols = 4;
n_rows = ceil(n_months / n_cols);

% Track cumulative disease-present sites
cumulative_present_locs = [];

for i = 1:n_months
    fprintf('Processing month %d/%d: %s...\n', i, n_months, string(unique_months(i), 'MMM-yyyy'));
    subplot(n_rows, n_cols, i);
    
    % Get data for this month
    month_mask = obsData.year_month == unique_months(i);
    month_data = obsData(month_mask, :);
    
    % Identify NEW disease-present sites for this month
    present_mask = strcmp(month_data.presence, 'P');
    
    if sum(present_mask) > 0
        % Add new disease-present locations to cumulative set
        new_locs = [month_data.lon(present_mask), month_data.lat(present_mask)];
        cumulative_present_locs = [cumulative_present_locs; new_locs];
        
        % Remove duplicates (same location reported multiple times)
        cumulative_present_locs = unique(cumulative_present_locs, 'rows');
    end
    
    % Plot all cumulative disease-present sites
    if ~isempty(cumulative_present_locs)
        scatter(cumulative_present_locs(:,1), cumulative_present_locs(:,2), 50, 'r', 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1);
        
        % Create disease boundary hull if enough points
        if size(cumulative_present_locs, 1) >= 3
            try
                k = boundary(cumulative_present_locs(:,1), cumulative_present_locs(:,2), 1);
                patch(cumulative_present_locs(k,1), cumulative_present_locs(k,2), [1 0.8 0.8], ...
                      'EdgeColor', 'r', 'FaceAlpha', 0.2, 'LineWidth', 2);
            catch
                % Not enough points for boundary
            end
        end
    end
    
    % Mark seed location (Flat Cay)
    hold on;
    if ~isempty(flat_cay_site_IDs)
        seed_coords = locations(flat_cay_site_IDs(1), :);
        plot(seed_coords(1), seed_coords(2), 'p', 'MarkerSize', 15, ...
             'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    end
    
    title(sprintf('%s\n%d cumulative disease sites', ...
          string(unique_months(i), 'MMM-yyyy'), size(cumulative_present_locs, 1)));
    xlabel('Longitude');
    ylabel('Latitude');
    axis equal tight;
    grid on;
    hold off;
    drawnow;  % Force update to show progress
end

% Add overall legend (only if there's space)
if n_months < (n_rows * n_cols)
    legend_ax = subplot(n_rows, n_cols, n_months + 1);
    axis off;
    hold on;
    scatter(nan, nan, 50, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    plot(nan, nan, 'p', 'MarkerSize', 15, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    legend('Disease Present', 'Seed Site (Flat Cay)', 'Location', 'best');
    hold off;
end

% Save figure
figname = 'ObservedDisease_Monthly_2018-2019';
saveas(gcf, fullfile(seascapePath, [figname '.png']));
saveas(gcf, fullfile(seascapePath, [figname '.fig']));

fprintf('Observed disease figure saved as: %s\n', figname);

% Summary statistics by island
fprintf('\n=== DISEASE OBSERVATIONS BY ISLAND ===\n');
islands = unique(obsData.island);
for i = 1:length(islands)
    island_data = obsData(strcmp(obsData.island, islands{i}), :);
    n_present = sum(strcmp(island_data.presence, 'P'));
    n_absent = sum(strcmp(island_data.presence, 'A'));
    fprintf('%s: %d present, %d absent\n', islands{i}, n_present, n_absent);
end
fprintf('======================================\n\n');