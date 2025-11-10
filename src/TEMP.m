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