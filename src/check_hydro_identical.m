%% NetCDF File Analysis and Comparison Script
% This script extracts, saves, plots, and compares all variables between
% HPC and drive NetCDF files

clear; clc; close all;

%% File paths (update these to your actual file paths)
tempPath = '/Users/benja/Documents/Farmer_Ben_Dissertation/Chapters/Chapter3_Metacommunities/temp/';
hpc_file = fullfile(tempPath, 'HPC_nest_1_20190101000000.nc');
drive_file = fullfile(tempPath, 'drive_nest_1_20190101000000.nc');

% Create output directory for plots and data
output_dir = fullfile(tempPath, 'analysis_output');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Read all variables from both files
fprintf('Reading NetCDF files...\n');

% Get file information
hpc_info = ncinfo(hpc_file);
drive_info = ncinfo(drive_file);

% Extract variable names (excluding dimensions)
var_names = {'zu', 'zv', 'zw', 'zt', 'zs'};
coord_names = {'Longitude', 'Latitude', 'Depth', 'Time'};

%% Read coordinate variables
fprintf('Reading coordinate variables...\n');
coordinates = struct();

for i = 1:length(coord_names)
    coord_name = coord_names{i};
    coordinates.hpc.(coord_name) = ncread(hpc_file, coord_name);
    coordinates.drive.(coord_name) = ncread(drive_file, coord_name);
end

%% Read data variables
fprintf('Reading data variables...\n');
data = struct();

for i = 1:length(var_names)
    var_name = var_names{i};
    fprintf('  Reading %s...\n', var_name);
    data.hpc.(var_name) = ncread(hpc_file, var_name);
    data.drive.(var_name) = ncread(drive_file, var_name);
end

%% Save data to .mat file
fprintf('Saving data to .mat file...\n');
save(fullfile(output_dir, 'netcdf_data_comparison.mat'), 'data', 'coordinates', 'hpc_info', 'drive_info');

%% Compare variables and generate comparison report
fprintf('Comparing variables between HPC and drive sources...\n');
comparison_results = struct();

% Compare coordinates first
fprintf('\n=== COORDINATE COMPARISON ===\n');
for i = 1:length(coord_names)
    coord_name = coord_names{i};
    hpc_coord = coordinates.hpc.(coord_name);
    drive_coord = coordinates.drive.(coord_name);
    
    % Check if identical
    is_identical = isequal(hpc_coord, drive_coord);
    comparison_results.coordinates.(coord_name).identical = is_identical;
    
    if is_identical
        fprintf('%s: IDENTICAL\n', coord_name);
    else
        fprintf('%s: DIFFERENT\n', coord_name);
        max_diff = max(abs(hpc_coord(:) - drive_coord(:)));
        fprintf('  Maximum absolute difference: %.2e\n', max_diff);
        comparison_results.coordinates.(coord_name).max_diff = max_diff;
    end
end

% Compare data variables
fprintf('\n=== DATA VARIABLE COMPARISON ===\n');
for i = 1:length(var_names)
    var_name = var_names{i};
    hpc_var = data.hpc.(var_name);
    drive_var = data.drive.(var_name);
    
    % Check if identical
    is_identical = isequal(hpc_var, drive_var);
    comparison_results.data.(var_name).identical = is_identical;
    
    if is_identical
        fprintf('%s: IDENTICAL\n', var_name);
    else
        fprintf('%s: DIFFERENT\n', var_name);
        max_diff = max(abs(hpc_var(:) - drive_var(:)));
        mean_diff = mean(abs(hpc_var(:) - drive_var(:)));
        fprintf('  Maximum absolute difference: %.2e\n', max_diff);
        fprintf('  Mean absolute difference: %.2e\n', mean_diff);
        comparison_results.data.(var_name).max_diff = max_diff;
        comparison_results.data.(var_name).mean_diff = mean_diff;
        
        % Calculate relative differences for non-zero values
        nonzero_mask = abs(hpc_var) > 1e-10;
        if any(nonzero_mask(:))
            rel_diff = abs((hpc_var(nonzero_mask) - drive_var(nonzero_mask)) ./ hpc_var(nonzero_mask));
            max_rel_diff = max(rel_diff);
            mean_rel_diff = mean(rel_diff);
            fprintf('  Maximum relative difference: %.2e%%\n', max_rel_diff * 100);
            fprintf('  Mean relative difference: %.2e%%\n', mean_rel_diff * 100);
            comparison_results.data.(var_name).max_rel_diff = max_rel_diff;
            comparison_results.data.(var_name).mean_rel_diff = mean_rel_diff;
        end
    end
end

%% Create comprehensive plots
fprintf('\nGenerating plots...\n');

% Set up plotting parameters
depth_indices = [1, 16, 32]; % Surface, mid-depth, bottom
lon = coordinates.hpc.Longitude;
lat = coordinates.hpc.Latitude;
depth = coordinates.hpc.Depth;

% Create longitude and latitude grids for plotting
[LON, LAT] = meshgrid(lon, lat);

%% Plot each variable at different depths
for var_idx = 1:length(var_names)
    var_name = var_names{var_idx};
    
    % Get variable info for titles and units
    var_info_hpc = [];
    for j = 1:length(hpc_info.Variables)
        if strcmp(hpc_info.Variables(j).Name, var_name)
            var_info_hpc = hpc_info.Variables(j);
            break;
        end
    end
    
    if isempty(var_info_hpc)
        continue;
    end
    
    % Extract units and long name
    units = '';
    long_name = var_name;
    for attr_idx = 1:length(var_info_hpc.Attributes)
        if strcmp(var_info_hpc.Attributes(attr_idx).Name, 'units')
            units = var_info_hpc.Attributes(attr_idx).Value;
        elseif strcmp(var_info_hpc.Attributes(attr_idx).Name, 'long_name')
            long_name = var_info_hpc.Attributes(attr_idx).Value;
        end
    end
    
    for depth_idx = 1:length(depth_indices)
        d_idx = depth_indices(depth_idx);
        
        % Extract data at this depth level
        hpc_data_2d = squeeze(data.hpc.(var_name)(:,:,d_idx,1))';
        drive_data_2d = squeeze(data.drive.(var_name)(:,:,d_idx,1))';
        
        % Create figure
        fig = figure('Position', [100, 100, 1400, 1000]);
        
        % Plot HPC data
        subplot(2,2,1);
        pcolor(LON, LAT, hpc_data_2d);
        shading interp;
        colorbar;
        title(sprintf('HPC: %s at %.1fm depth', long_name, depth(d_idx)));
        xlabel('Longitude (°E)');
        ylabel('Latitude (°N)');
        caxis_range = caxis;
        
        % Plot Drive data
        subplot(2,2,2);
        pcolor(LON, LAT, drive_data_2d);
        shading interp;
        colorbar;
        title(sprintf('Drive: %s at %.1fm depth', long_name, depth(d_idx)));
        xlabel('Longitude (°E)');
        ylabel('Latitude (°N)');
        caxis(caxis_range); % Use same color scale
        
        % Plot difference
        subplot(2,2,3);
        diff_data = hpc_data_2d - drive_data_2d;
        pcolor(LON, LAT, diff_data);
        shading interp;
        colorbar;
        title(sprintf('Difference (HPC - Drive): %s', long_name));
        xlabel('Longitude (°E)');
        ylabel('Latitude (°N)');
        
        % Statistics subplot
        subplot(2,2,4);
        axis off;
        
        % Calculate statistics
        max_hpc = max(hpc_data_2d(:));
        min_hpc = min(hpc_data_2d(:));
        mean_hpc = mean(hpc_data_2d(:), 'omitnan');
        
        max_drive = max(drive_data_2d(:));
        min_drive = min(drive_data_2d(:));
        mean_drive = mean(drive_data_2d(:), 'omitnan');
        
        max_diff = max(abs(diff_data(:)));
        mean_diff = mean(abs(diff_data(:)), 'omitnan');
        
        stats_text = {
            sprintf('Variable: %s (%s)', var_name, units);
            sprintf('Depth: %.1f m', depth(d_idx));
            '';
            'HPC Statistics:';
            sprintf('  Min: %.6f', min_hpc);
            sprintf('  Max: %.6f', max_hpc);
            sprintf('  Mean: %.6f', mean_hpc);
            '';
            'Drive Statistics:';
            sprintf('  Min: %.6f', min_drive);
            sprintf('  Max: %.6f', max_drive);
            sprintf('  Mean: %.6f', mean_drive);
            '';
            'Difference Statistics:';
            sprintf('  Max |diff|: %.2e', max_diff);
            sprintf('  Mean |diff|: %.2e', mean_diff);
            '';
            sprintf('Identical: %s', char(string(comparison_results.data.(var_name).identical)));
        };
        
        text(0.05, 0.95, stats_text, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
             'FontSize', 10, 'FontName', 'FixedWidth');
        
        % Save figure
        filename = sprintf('%s_depth_%.0fm_comparison.png', var_name, depth(d_idx));
        saveas(fig, fullfile(output_dir, filename));
        fprintf('  Saved: %s\n', filename);
    end
end

%% Create depth profiles comparison
fprintf('Creating depth profile comparisons...\n');

% Select a representative location (middle of domain)
lon_idx = round(length(lon)/2);
lat_idx = round(length(lat)/2);

for var_idx = 1:length(var_names)
    var_name = var_names{var_idx};
    
    % Extract depth profiles
    hpc_profile = squeeze(data.hpc.(var_name)(lon_idx, lat_idx, :, 1));
    drive_profile = squeeze(data.drive.(var_name)(lon_idx, lat_idx, :, 1));
    
    % Create figure
    fig = figure('Position', [100, 100, 800, 600]);
    
    plot(hpc_profile, -depth, 'b-', 'LineWidth', 2, 'DisplayName', 'HPC');
    hold on;
    plot(drive_profile, -depth, 'r--', 'LineWidth', 2, 'DisplayName', 'Drive');
    
    xlabel(sprintf('%s', var_name));
    ylabel('Depth (m)');
    title(sprintf('Depth Profile Comparison: %s', var_name));
    legend('Location', 'best');
    grid on;
    
    % Save figure
    filename = sprintf('%s_depth_profile_comparison.png', var_name);
    saveas(fig, fullfile(output_dir, filename));
    fprintf('  Saved: %s\n', filename);
end

%% Save comparison results
fprintf('Saving comparison results...\n');
save(fullfile(output_dir, 'comparison_results.mat'), 'comparison_results');

% Write summary report to text file
fid = fopen(fullfile(output_dir, 'comparison_summary.txt'), 'w');
fprintf(fid, 'NetCDF File Comparison Summary\n');
fprintf(fid, '==============================\n\n');

fprintf(fid, 'Files compared:\n');
fprintf(fid, 'HPC: %s\n', hpc_file);
fprintf(fid, 'Drive: %s\n\n', drive_file);

fprintf(fid, 'COORDINATE VARIABLES:\n');
for i = 1:length(coord_names)
    coord_name = coord_names{i};
    if comparison_results.coordinates.(coord_name).identical
        fprintf(fid, '%s: IDENTICAL\n', coord_name);
    else
        fprintf(fid, '%s: DIFFERENT (max diff: %.2e)\n', coord_name, ...
                comparison_results.coordinates.(coord_name).max_diff);
    end
end

fprintf(fid, '\nDATA VARIABLES:\n');
for i = 1:length(var_names)
    var_name = var_names{i};
    if comparison_results.data.(var_name).identical
        fprintf(fid, '%s: IDENTICAL\n', var_name);
    else
        fprintf(fid, '%s: DIFFERENT\n', var_name);
        fprintf(fid, '  Max absolute diff: %.2e\n', comparison_results.data.(var_name).max_diff);
        fprintf(fid, '  Mean absolute diff: %.2e\n', comparison_results.data.(var_name).mean_diff);
        if isfield(comparison_results.data.(var_name), 'max_rel_diff')
            fprintf(fid, '  Max relative diff: %.2e%%\n', comparison_results.data.(var_name).max_rel_diff * 100);
            fprintf(fid, '  Mean relative diff: %.2e%%\n', comparison_results.data.(var_name).mean_rel_diff * 100);
        end
    end
end

fclose(fid);

fprintf('\n=== ANALYSIS COMPLETE ===\n');
fprintf('Results saved in: %s\n', output_dir);
fprintf('- netcdf_data_comparison.mat: All extracted data\n');
fprintf('- comparison_results.mat: Detailed comparison results\n');
fprintf('- comparison_summary.txt: Text summary of differences\n');
fprintf('- Multiple PNG files: Visualizations of each variable\n');

%% Display final summary
fprintf('\n=== FINAL SUMMARY ===\n');
all_identical = true;
for i = 1:length(var_names)
    var_name = var_names{i};
    if ~comparison_results.data.(var_name).identical
        all_identical = false;
        break;
    end
end

if all_identical
    fprintf('ALL VARIABLES ARE IDENTICAL between HPC and drive sources!\n');
else
    fprintf('Some variables show differences between HPC and drive sources.\n');
    fprintf('Check the detailed comparison results and plots for more information.\n');
end