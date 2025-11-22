%% Script to create a release file for the CMS, simulating the dispersal and connectivity of SCTLD in the Caribbean
% 31 Oct 2023
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/Release';
clear;clc

%% test version for a 5-day release for Dan's Nov 2023 EPSCoR meeting

% .nc file format: %nest_nestnumber_yyyymmddhhmmss [.nc]
%   here, I'm working with 'nest_1_20190101000000' through
%                          'nest_1_20190101200000'
%   so,                    'nest_1_2019_January_1_00:00' through
%                          'nest_1_2019_January_1_20:00'

% Release file format:
%   Polygon Longitude Latitude Depth Number Year Month Day Second
%   1       277.2     24.66    0     10     2016 1     1   0

repnum = 10; %want to use 20 processors, but split over two release lines
IDs = (1:20)';
lons1 = repelem(295.0456, repnum)'; %Grammanik Bank
lons2 = repelem(295.011888, repnum)'; %Flat Cay
lons = vertcat(lons1, lons2);
lats1 = repelem(18.1951, repnum)'; %Grammanik Bank
lats2 = repelem(18.316191, repnum)';
lats = vertcat(lats1, lats2);

release(:,1) = IDs;
release(:,2) = lons;
release(:,3) = lats;
release(:,4) = 1; %in meters. releasing at surface
release(:,5) = 100; %#/particles. 1000 total per location (2 locations, 100 per release line [10])
release(:,6) = 2019;
release(:,7) = 2; %February
release(:,8) = 4; %4th
release(:,9) = 0; %in seconds. releasing at midnight

writematrix(release, 'ReleaseFile_tester.txt', 'delimiter', '\t');

