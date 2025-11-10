%% Plot infection dynamics for a random site (TWO 4-panel figures)

% Pick a random site that got infected
final_I = Y(end, num_sites*3+1:num_sites*6);
final_P = sum(reshape(final_I, num_sites, 3), 2);
infected_sites = find(final_P > 1e-5);

if isempty(infected_sites)
    fprintf('No sites infected - cannot plot\n');
else
    % Pick random infected site
    random_site = infected_sites(randi(length(infected_sites)));
    
    % Extract time series for this site
    S_LS = Y(:, random_site);
    S_MS = Y(:, num_sites + random_site);
    S_HS = Y(:, 2*num_sites + random_site);
    
    I_LS = Y(:, 3*num_sites + random_site);
    I_MS = Y(:, 4*num_sites + random_site);
    I_HS = Y(:, 5*num_sites + random_site);
    
    R_LS = Y(:, 6*num_sites + random_site);
    R_MS = Y(:, 7*num_sites + random_site);
    R_HS = Y(:, 8*num_sites + random_site);
    
    % ===== FIGURE 1: Full SIR curves =====
    figure('Position', [100 100 1200 800]);
    
    % Panel 1: Total (all groups combined)
    subplot(2,2,1);
    hold on;
    plot(t, S_LS + S_MS + S_HS, 'b-', 'LineWidth', 2);
    plot(t, I_LS + I_MS + I_HS, 'r-', 'LineWidth', 2);
    plot(t, R_LS + R_MS + R_HS, 'k-', 'LineWidth', 2);
    xlabel('Days'); ylabel('Fraction');
    title('Total (All Groups)');
    legend('S', 'I', 'R', 'Location', 'best');
    grid on;
    
    % Panel 2: Low susceptibility
    subplot(2,2,2);
    hold on;
    plot(t, S_LS, 'b-', 'LineWidth', 2);
    plot(t, I_LS, 'r-', 'LineWidth', 2);
    plot(t, R_LS, 'k-', 'LineWidth', 2);
    xlabel('Days'); ylabel('Fraction');
    title(sprintf('Low Susceptibility (N=%.4f)', N_LS(random_site)));
    legend('S', 'I', 'R', 'Location', 'best');
    grid on;
    
    % Panel 3: Moderate susceptibility
    subplot(2,2,3);
    hold on;
    plot(t, S_MS, 'b-', 'LineWidth', 2);
    plot(t, I_MS, 'r-', 'LineWidth', 2);
    plot(t, R_MS, 'k-', 'LineWidth', 2);
    xlabel('Days'); ylabel('Fraction');
    title(sprintf('Moderate Susceptibility (N=%.4f)', N_MS(random_site)));
    legend('S', 'I', 'R', 'Location', 'best');
    grid on;
    
    % Panel 4: High susceptibility
    subplot(2,2,4);
    hold on;
    plot(t, S_HS, 'b-', 'LineWidth', 2);
    plot(t, I_HS, 'r-', 'LineWidth', 2);
    plot(t, R_HS, 'k-', 'LineWidth', 2);
    xlabel('Days'); ylabel('Fraction');
    title(sprintf('High Susceptibility (N=%.4f)', N_HS(random_site)));
    legend('S', 'I', 'R', 'Location', 'best');
    grid on;
    
    % Overall title
    sgtitle(sprintf('FULL SIR DYNAMICS | Site %d (unique_ID: %d) | Total Cover: %.4f', ...
                    random_site, unique_IDs(random_site), N_site(random_site)));
    
    % ===== FIGURE 2: Infected curves only (SAME Y-AXIS) =====
    % Find max infected value across all groups
    max_I = max([max(I_LS + I_MS + I_HS), max(I_LS), max(I_MS), max(I_HS)]);
    
    figure('Position', [150 150 1200 800]);
    
    % Panel 1: Total infected
    subplot(2,2,1);
    plot(t, I_LS + I_MS + I_HS, 'r-', 'LineWidth', 2);
    xlabel('Days'); ylabel('Infected Fraction');
    title('Total Infected (All Groups)');
    ylim([0 max_I]);
    grid on;
    
    % Panel 2: Low susceptibility infected
    subplot(2,2,2);
    plot(t, I_LS, 'r-', 'LineWidth', 2);
    xlabel('Days'); ylabel('Infected Fraction');
    title(sprintf('Low Susceptibility Infected (N=%.4f)', N_LS(random_site)));
    ylim([0 max_I]);
    grid on;
    
    % Panel 3: Moderate susceptibility infected
    subplot(2,2,3);
    plot(t, I_MS, 'r-', 'LineWidth', 2);
    xlabel('Days'); ylabel('Infected Fraction');
    title(sprintf('Moderate Susceptibility Infected (N=%.4f)', N_MS(random_site)));
    ylim([0 max_I]);
    grid on;
    
    % Panel 4: High susceptibility infected
    subplot(2,2,4);
    plot(t, I_HS, 'r-', 'LineWidth', 2);
    xlabel('Days'); ylabel('Infected Fraction');
    title(sprintf('High Susceptibility Infected (N=%.4f)', N_HS(random_site)));
    ylim([0 max_I]);
    grid on;
    
    % Overall title
    sgtitle(sprintf('INFECTED CURVES ONLY | Site %d (unique_ID: %d) | Total Cover: %.4f', ...
                    random_site, unique_IDs(random_site), N_site(random_site)));
    
    fprintf('Plotted site %d (unique_ID: %d)\n', random_site, unique_IDs(random_site));
end