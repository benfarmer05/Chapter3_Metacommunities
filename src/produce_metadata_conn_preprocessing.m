%% CREATE METADATA FILE FOR EXISTING PREPROCESSED EVENTS
% Run this if preprocessing was interrupted and metadata.mat wasn't created
clear; clc;

%% Configuration - MATCH YOUR PREPROCESSING SETTINGS
YEAR = 2019;
QUARTER = 1;

%% Setup paths
projectPath = matlab.project.rootProject().RootFolder;
dataPath = fullfile(projectPath, 'data');
outputPath = fullfile('D:\Dissertation\CMS_traj\output');
quarter_name = sprintf('Q%d_%d', QUARTER, YEAR);
preprocessPath = fullfile(outputPath, 'CMS_traj', quarter_name);

fprintf('=== CREATING METADATA FILE ===\n');
fprintf('Path: %s\n', preprocessPath);

%% Check for unique_events.mat (created early in preprocessing)
unique_events_file = fullfile(preprocessPath, 'unique_events.mat');
if ~exist(unique_events_file, 'file')
    error('unique_events.mat not found. Preprocessing may not have started yet.');
end

fprintf('Loading unique_events.mat...\n');
load(unique_events_file, 'unique_release_dates');
n_events = length(unique_release_dates);
fprintf('Total release events in this quarter: %d\n', n_events);

%% Load reef data for n_locations
centroids = readmatrix(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
n_locations = size(centroids, 1);
fprintf('Reef locations: %d\n', n_locations);

%% Find existing event files
event_files = dir(fullfile(preprocessPath, 'event_*.mat'));
n_processed = length(event_files);
fprintf('Found %d processed event files\n', n_processed);

if n_processed == 0
    error('No event_*.mat files found. Nothing to create metadata for.');
end

%% Create metadata structure
metadata = struct();
metadata.year = YEAR;
metadata.quarter = QUARTER;
metadata.n_events = n_events;
metadata.n_events_processed = n_processed;
metadata.unique_release_dates = unique_release_dates;
metadata.n_locations = n_locations;
metadata.preprocessing_date = datetime('now');
metadata.metadata_created_date = datetime('now');
metadata.variables_saved = {'lon', 'lat', 'depth', 'distance', 'exitcode', 'time'};
metadata.event_files = {event_files.name}';
metadata.preprocessing_complete = (n_processed == n_events);

%% Save metadata
metadata_file = fullfile(preprocessPath, 'metadata.mat');
save(metadata_file, 'metadata', '-v7.3');

fprintf('\n=== SUCCESS ===\n');
fprintf('Created: %s\n', metadata_file);
fprintf('Processed events: %d/%d\n', n_processed, n_events);
if metadata.preprocessing_complete
    fprintf('Status: Preprocessing COMPLETE\n');
else
    fprintf('Status: Preprocessing PARTIAL (%d events remaining)\n', n_events - n_processed);
end
fprintf('\nYou can now run the connectivity analysis script.\n');