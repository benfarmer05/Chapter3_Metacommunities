% Clear workspace and command window
clear all;
clc;

%% Model Parameters
linewidths = 0.75;
symbsizes = 1.3;
titlesize = 10;
textsize = 9;

%% Run All Scenarios
scenario1 = run_one_scenario(3.64, 15.7, 4.97, true);
scenario1.scenario = repmat({'Nearshore (cover)'}, height(scenario1), 1);

scenario2 = run_one_scenario(3.64, 15.7, 4.97, false);
scenario2.scenario = repmat({'Nearshore (raw)'}, height(scenario2), 1);

scenario3 = run_one_scenario(10, 5, 2, true);
scenario4 = run_one_scenario(5, 2.5, 1, true);
scenario5 = run_one_scenario(2.5, 1.25, 3, true);
scenario6 = run_one_scenario(0.5, 5, 1, true);
scenario7 = run_one_scenario(0.25, 2.5, 0.5, true);
scenario8 = run_one_scenario(0.125, 1.25, 0.25, true);
scenario9 = run_one_scenario(0.1, 1, 0.02, true);

% Combine all scenarios
all_results = [scenario1; scenario2; scenario3; scenario4; scenario5; ...
               scenario6; scenario7; scenario8; scenario9];

%% Plot Full SIR
plot_full_sir(all_results, linewidths, titlesize, textsize);

%% ==================== FUNCTION DEFINITIONS ====================

function dydt = SIR_multi(t, y, p)
    % Extract state variables
    S_LS = y(1); I_LS = y(2); R_LS = y(3);
    S_MS = y(4); I_MS = y(5); R_MS = y(6);
    S_HS = y(7); I_HS = y(8); R_HS = y(9);
    
    % Total infected prevalence
    P = I_LS + I_MS + I_HS;
    
    % Low susceptibility group
    dS_LS_dt = -p.b_LS * S_LS * P / p.N_LS;
    dI_LS_dt = p.b_LS * S_LS * P / p.N_LS - p.g_LS * I_LS;
    dR_LS_dt = p.g_LS * I_LS;
    
    % Moderate susceptibility group
    dS_MS_dt = -p.b_MS * S_MS * P / p.N_MS;
    dI_MS_dt = p.b_MS * S_MS * P / p.N_MS - p.g_MS * I_MS;
    dR_MS_dt = p.g_MS * I_MS;
    
    % High susceptibility group
    dS_HS_dt = -p.b_HS * S_HS * P / p.N_HS;
    dI_HS_dt = p.b_HS * S_HS * P / p.N_HS - p.g_HS * I_HS;
    dR_HS_dt = p.g_HS * I_HS;
    
    dydt = [dS_LS_dt; dI_LS_dt; dR_LS_dt; ...
            dS_MS_dt; dI_MS_dt; dR_MS_dt; ...
            dS_HS_dt; dI_HS_dt; dR_HS_dt];
end

function result = run_one_scenario(cover_LS, cover_MS, cover_HS, use_cover_scale)
    % Prevent division by zero
    epsilon = 1e-6;
    cover_LS = max(cover_LS, epsilon);
    cover_MS = max(cover_MS, epsilon);
    cover_HS = max(cover_HS, epsilon);
    
    % Parameters for THIS scenario only
    beta_LS = 0.03;
    beta_MS = 0.14;
    beta_HS = 2.08;
    gamma_LS = 0.05;
    gamma_MS = 0.55;
    gamma_HS = 3.33;
    
    % Scale to 0-1 if requested
    if use_cover_scale
        N_LS = cover_LS / 100;
        N_MS = cover_MS / 100;
        N_HS = cover_HS / 100;
        initial_infected = 0.01 / 100;
    else
        N_LS = cover_LS;
        N_MS = cover_MS;
        N_HS = cover_HS;
        initial_infected = 0.01;
    end
    
    % Seed infection: HS -> MS -> LS priority
    if cover_HS > epsilon * 10
        I0_LS = 0; I0_MS = 0; I0_HS = initial_infected;
    elseif cover_MS > epsilon * 10
        I0_LS = 0; I0_MS = initial_infected; I0_HS = 0;
    else
        I0_LS = initial_infected; I0_MS = 0; I0_HS = 0;
    end
    
    % Initial state
    state0 = [N_LS - I0_LS; I0_LS; 0; ...
              N_MS - I0_MS; I0_MS; 0; ...
              N_HS - I0_HS; I0_HS; 0];
    
    % Parameters structure
    parms.b_LS = beta_LS; parms.b_MS = beta_MS; parms.b_HS = beta_HS;
    parms.g_LS = gamma_LS; parms.g_MS = gamma_MS; parms.g_HS = gamma_HS;
    parms.N_LS = N_LS; parms.N_MS = N_MS; parms.N_HS = N_HS;
    
    % Run model
    times = 0:0.1:400;
    [t, y] = ode45(@(t,y) SIR_multi(t, y, parms), times, state0);
    
    % Create result table
    result = array2table([t, y], 'VariableNames', ...
        {'time', 'S_LS', 'I_LS', 'R_LS', 'S_MS', 'I_MS', 'R_MS', 'S_HS', 'I_HS', 'R_HS'});
    
    % Add metadata
    total = cover_LS + cover_MS + cover_HS;
    scenario_label = sprintf('L:%.1f M:%.1f H:%.1f (%.1f%%)', ...
                            cover_LS, cover_MS, cover_HS, total);
    result.scenario = repmat({scenario_label}, height(result), 1);
    
    if use_cover_scale
        result.scale_type = repmat({'Cover'}, height(result), 1);
    else
        result.scale_type = repmat({'Raw'}, height(result), 1);
    end
end

function plot_full_sir(all_results, linewidths, titlesize, textsize)
    % Get unique scenarios
    scenarios = unique(all_results.scenario);
    n_scenarios = length(scenarios);
    
    % Create figure
    figure('Position', [100, 100, 1200, 900]);
    
    % Colors and markers
    colors = [30/255, 144/255, 1; ...      % Low - Blue
              1, 215/255, 0; ...           % Moderate - Gold
              1, 20/255, 147/255];         % High - Deep Pink
    
    markers = {'o', '^', 's'};  % Susceptible, Infected, Removed
    marker_sizes = [4, 4, 4];
    
    % Plot each scenario
    for i = 1:n_scenarios
        subplot(3, 3, i);
        hold on;
        
        % Filter data for this scenario
        idx = strcmp(all_results.scenario, scenarios{i});
        data = all_results(idx, :);
        
        % Plot each group and compartment
        groups = {'LS', 'MS', 'HS'};
        compartments = {'S', 'I', 'R'};
        
        for g = 1:3
            for c = 1:3
                var_name = sprintf('%s_%s', compartments{c}, groups{g});
                
                % Plot line with markers (but sparse markers)
                plot(data.time, data.(var_name), '-', ...
                    'Color', colors(g,:), 'LineWidth', linewidths);
                
                % Add sparse markers
                marker_idx = 1:400:length(data.time);
                plot(data.time(marker_idx), data.(var_name)(marker_idx), ...
                    markers{c}, 'MarkerSize', marker_sizes(c), ...
                    'MarkerFaceColor', colors(g,:), ...
                    'MarkerEdgeColor', colors(g,:));
            end
        end
        
        % Formatting
        xlabel('Day of outbreak', 'FontSize', titlesize, 'FontName', 'Georgia');
        ylabel('Percent cover', 'FontSize', titlesize, 'FontName', 'Georgia');
        title(scenarios{i}, 'FontSize', 8, 'FontName', 'Georgia', 'FontWeight', 'normal');
        set(gca, 'FontSize', textsize, 'FontName', 'Georgia');
        box off;
        hold off;
    end
    
    % Add legend to the figure
    subplot(3, 3, 9);
    hold on;
    h = zeros(6, 1);
    labels = cell(6, 1);
    
    % Group colors
    for g = 1:3
        h(g) = plot(nan, nan, '-', 'Color', colors(g,:), 'LineWidth', 2);
    end
    labels(1:3) = {'Low', 'Moderate', 'High'};
    
    % Compartment markers
    for c = 1:3
        h(3+c) = plot(nan, nan, markers{c}, 'MarkerSize', 6, ...
            'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    end
    labels(4:6) = {'Susceptible', 'Infected', 'Removed'};
    
    legend(h, labels, 'Location', 'best', 'FontSize', textsize, 'FontName', 'Georgia');
    axis off;
    
    hold off;
end