%% Script to create a release file for the CMS, simulating the dispersal and connectivity of SCTLD in the Virgin Islands & Puerto Rico
% 21 May 2024

clear;clc

%% test version for couple-day release
% GOALS HERE:
% - Figure out how to code a release from every point, at every hour. Do
% this just for the test files I have created (keep it to just July 1st for
% now).
% - Just release 1 particle per point for now
% - Figure out depth here! That would be awesome and really useful.
% - Need to specify year (2019 for all), month (July for now), day (1 for
% now), and second (this will be 0, 3600, 7200, etc. up to 82800 for the
% 23rd hour of the day)

% .nc file format: %nest_nestnumber_yyyymmddhhmmss [.nc]
%   here, I'm working with 'nest_1_20190101000000' through
%                          'nest_1_20190101200000'
%   so,                    'nest_1_2019_January_1_00:00' through
%                          'nest_1_2019_January_1_20:00'

% Release file format:
%   Polygon Longitude Latitude Depth Number Year Month Day Second
%   1       277.2     24.66    0     10     2016 1     1   0

scriptDrive = '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/Release';

cd(scriptDrive)

%read in release points from GIS output
[relpnts] = readmatrix('points_650_none-on-land.csv');
num_pnts = length(relpnts);

%construct bones of the release structure
repvalue = 365; %releasing daily?
IDs = repmat((1:num_pnts)', repvalue, 1);
lons = repmat(relpnts(:,11), repvalue, 1);
lons = lons+360;
lats = repmat(relpnts(:,12), repvalue, 1);
% days_months = load('Days_Months.txt');
% days = repelem(days_months(:,1),num_pnts,1);
% months = repelem(days_months(:,2),num_pnts,1);

% Generate days and months for a non-leap year (2019)
days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
num_days = sum(days_in_month);  % Total number of days in a year

% Pre-initialize months and days arrays
months = zeros(num_days, 1);  % Pre-initialize months as a column vector of zeros
days = zeros(num_days, 1);    % Pre-initialize days as a column vector of zeros

index = 1;  % Initialize index to keep track of the position in the days/months arrays
for month = 1:12
    days(index:index + days_in_month(month) - 1) = 1:days_in_month(month);  % Generate days for the month
    months(index:index + days_in_month(month) - 1) = month;  % Set the month number
    index = index + days_in_month(month);  % Update index for the next month
end

% Repeat the days and months data to match the number of points
days = repelem(days, num_pnts, 1);
months = repelem(months, num_pnts, 1);


% Fill in release file matrices with pre-allocated columns
release(:,1) = IDs;
release(:,2) = lons;
release(:,3) = lats;
release(:,4) = 1; %in meters. releasing at surface
release(:,5) = 10; %#/particles
release(:,6) = 2019;
release(:,7) = months;
release(:,8) = days;
release(:,9) = 0; %in seconds. releasing at midnight


%write the matrix to the file with the constructed file name
currentDateTime = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
currentDateTimeStr = string(currentDateTime, 'yyyyMMdd_HHmmss');
fileName = "ReleaseFile_tester_" + currentDateTimeStr + ".txt";
writematrix(release, fileName, 'delimiter', '\t');



% %current redacted - very interesting code but may not be necessary with the
% % date range-informed revised code above
% %generate days and months for a non-leap year (2019)
% days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
% num_days = sum(days_in_month); % Total number of days in a year
% 
% %pre-initialize months and days arrays
% months = zeros(num_days, 1);  % Pre-initialize months as a column vector of zeros
% days = zeros(num_days, 1);    % Pre-initialize days as a column vector of zeros
% index = 1; % Initialize index to keep track of the position in the days/months arrays
% 
% %fill in the days and months arrays:
% for month = 1:length(days_in_month)
%     days(index:index + days_in_month(month) - 1) = 1:days_in_month(month); % Generate days for the month
%     months(index:index + days_in_month(month) - 1) = month; % Set the month number
%     index = index + days_in_month(month); % Update index for the next month
% end
% 
% %generate a continuous day count for the entire year
% day_of_year = (1:num_days)';






%% 
[numbers] = xlsread('points.xlsx'); %bring in point data from GIS
numReleasePoints = length(numbers);

% create matrix with enough zeros for all release points and 9 columns
% numReleases = 365*numReleasePoints; %number of total releases (or rows); could be automated to account of #/days in a leap year (e.g. 2020)
numReleases = 12*numReleasePoints; %for once-a-month release. this should be automated
SCTLD_rel = zeros(numReleases,9); %matrix of zeros; days will be 366 in a leap year (e.g. 2020)

% ID# (column 1)
% idOneRep = numbers(:,1);
% id = repmat(idOneRep, 365, 1);
% SCTLD_rel(:,1) = id;

idOneRep = numbers(:,1);
id = repmat(idOneRep, 12, 1); %changed all this stuff to be once a month rather than once a day all year - THIS SHOULD BE AUTOMATED! need to make a variable for #/release events/year (365, 12, etc.)
SCTLD_rel(:,1) = id;

% X/Y (columns 2/3)
% xOneRep = numbers(:,5);
% xCoordinate = repmat(xOneRep, 365, 1);
% SCTLD_rel(:,2) = xCoordinate;

xOneRep = numbers(:,2);
xCoordinate = repmat(xOneRep, 12, 1);
SCTLD_rel(:,2) = xCoordinate;

% yOneRep = numbers(:,6);
% yCoordinate = repmat(yOneRep, 365, 1);
% SCTLD_rel(:,3) = yCoordinate;

yOneRep = numbers(:,3);
yCoordinate = repmat(yOneRep, 12, 1);
SCTLD_rel(:,3) = yCoordinate;

% assign depth (column 4)
% SCTLD_rel(:,4) = interpvalsRelease(:,1); %for depth-adjusted release
SCTLD_rel(:,4) = zeros;

% assign number of particles/release event (I think)
SCTLD_rel(:,5) = 50;

% assign year to each row (column 6)
SCTLD_rel(:,6) = 2008; % should automate integration of multiple years for sure!

% assign month value to each row (column 7 in txtfile)
monthsOneRep = vertcat((ones(31,1)), (2*ones(28,1)), (3*ones(31,1)), (4*ones(30,1)), (5*ones(31,1)), (6*ones(30,1)), (7*ones(31,1)), (8*ones(31,1)), (9*ones(30,1)), (10*ones(31,1)), (11*ones(30,1)), (12*ones(31,1)));
months = repmat(monthsOneRep, numReleasePoints, 1); % "length(numbers)" is #/release points
SCTLD_rel(:,7) = months;

%dumb monthly (rather than daily) release version. again, need to automate
%badly
monthsOneRep = vertcat(1,2,3,4,5,6,7,8,9,10,11,12);
months = repmat(monthsOneRep, numReleasePoints, 1); % "length(numbers)" is #/release points
SCTLD_rel(:,7) = months;

%%

%assign day value to each row (column 8 in txtfile)
Jan = 1:31; Feb = 1:28; Mar = 1:31; Apr = 1:30; May = 1:31; Jun = 1:30; Jul = 1:31; Aug = 1:31; 
Sep = 1:30; Oct = 1:31; Nov = 1:30; Dec = 1:31;

Jan = Jan'; Feb = Feb'; Mar = Mar'; Apr = Apr'; May = May'; Jun = Jun'; Jul = Jul'; Aug = Aug';
Sep = Sep'; Oct = Oct'; Nov = Nov'; Dec = Dec'; %transpose so matrices are vertical

daysOneRep = vertcat(Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec); %creates big column of all the days in each month, sequentially
days = repmat(daysOneRep, numReleasePoints,1); %replicates the big column times #/release points
SCTLD_rel(:,8) = days;

%dumb monthly version
Jan = 1; Feb = 1; Mar = 1; Apr = 1; May = 1; Jun = 1; Jul = 1; Aug = 1; Sep = 1; Oct = 1; Nov = 1; Dec = 1; 

%% careful with this; dlmwrite is apparently just terrible lol
% dlmwrite('release2019.txt', SCTLD_rel, 'delimiter', '\t', 'precision', 7); % bad??
writematrix(SCTLD_rel, 'MAR_ReleaseFile.txt', 'delimiter', '\t'); % /t is tab. also good??

%% make the short release file
% Take a subset of rows in the release file
SCLTD_rel_shorty = SCTLD_rel(1:3050:end,:); %greatly shorten release file
SCLTD_rel_shorty(:,1) = 1:length(SCLTD_rel_shorty); % re-ID the shortened release file
dlmwrite('SCTLD_rel_short.txt', SCLTD_rel_short, 'delimiter', '\t', 'precision', 7);
writematrix(SCLTD_rel_shorty, 'SCTLD_rel_shorty.txt', 'delimiter', '\t');

SCLTD_rel_shorty(:,2) = round(SCLTD_rel_shorty(:,2),8);
SCLTD_rel_shorty(:,3) = round(SCLTD_rel_shorty(:,3),8);
type('SCTLD_rel_short.txt'); %spit out the text file in command window

%% note to self: dear god please go back and fix all of this code i'm dying just looking at it :(
% below was a shotgun approach to running the USVI model for the 2021 RAPID
% meeting. this also nearly destroyed my coding soul. approach with caution
% and keep 6 feet of distance, may cause rapid deterioration of mental
% state
cd '/Users/benja/MATLAB-Drive/USVI_release';

% bring in point data from excel
[numbers, strings, raw] = xlsread('releasepoints_BVI_gonezo.xlsx'); %bring in point data from GIS
numReleasePoints = length(numbers); %number of release points
numReleases = numReleasePoints*12;

IDs = 1:numReleasePoints;
xOneRep = numbers(:,2);
yOneRep = numbers(:,3);
Xs = xOneRep;
Ys = yOneRep;

IDsrep = repmat(IDs',12,1);
Xsrep = repmat(Xs,12,1);
Ysrep = repmat(Ys,12,1);
Col4 = zeros(numReleases,1);
Col5 = 50*ones(numReleases,1);
Col6 = 2010*ones(numReleases,1);
Col8 = ones(numReleases,1);
Col9 = zeros(numReleases,1);

Monthrep = [repmat(1,numReleasePoints,1);
            repmat(2,numReleasePoints,1);
            repmat(3,numReleasePoints,1)
            repmat(4,numReleasePoints,1)
            repmat(5,numReleasePoints,1)
            repmat(6,numReleasePoints,1)
            repmat(7,numReleasePoints,1)
            repmat(8,numReleasePoints,1)
            repmat(9,numReleasePoints,1)
            repmat(10,numReleasePoints,1)
            repmat(11,numReleasePoints,1)
            repmat(12,numReleasePoints,1)];
        
Test = [IDsrep,Xsrep,Ysrep,Col4,Col5,Col6,Monthrep,Col8,Col9];

% % Take a subset of rows in the release file
% SCLTD_rel_shorty = Test(1:numReleasePoints,:); %run full thing, but just for one day
% SCLTD_rel_shorty(:,1) = 1:length(SCLTD_rel_shorty); % re-ID the shortened release file
% SCLTD_rel_shorty(:,2) = round(SCLTD_rel_shorty(:,2),8);
% SCLTD_rel_shorty(:,3) = round(SCLTD_rel_shorty(:,3),8);
% writematrix(SCLTD_rel_shorty, 'SCTLD_rel_short.txt', 'delimiter', '\t');

SCLTD_rel_full = Test; %run full thing, but just for one day
SCLTD_rel_full(:,2) = round(SCLTD_rel_full(:,2),8);
SCLTD_rel_full(:,3) = round(SCLTD_rel_full(:,3),8);
writematrix(SCLTD_rel_full, 'SCTLD_rel_BVI_point_2010.txt', 'delimiter', '\t');


%% Depth extraction

% Note - 28 September 2023
%       This is really confusing. I've tried looking at interp1,
%       interp2, griddedinterpolant, and scatteredinterpolant and none of
%       them really seem to make sense for what I'm doing. a nearest
%       neighbor search may be the most applicable, though - check out this
%       one: https://www.mathworks.com/matlabcentral/answers/446980-how-to-find-closest-lat-lon-on-a-grid-to-an-existing-lat-lon

clear;clc
cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/Release';     
[relpts] = readmatrix('points_650_none-on-land.csv');

cd '/Users/benja/Library/CloudStorage/Box-Box/Big_Projects/Chapter3_Metacommunities/CMS-formatting_Spring2022_to_Summer2023/2023_hydro/Processed/Testing'
ROMS_grid = 'vigrid.nc';
lons = ncread(ROMS_grid, 'lon_rho'); %sticking with rho-grid, as bathymetry ('h') is on there
lats = ncread(ROMS_grid, 'lat_rho');
mask = ncread(ROMS_grid, 'mask_rho');
depths = ncread(ROMS_grid, 'h');
depthcut1 = depths;
depthcut1(depthcut1==0) = NaN;
depthcut1(mask == 0) = NaN;
figure();imagescn(lons', lats', -depthcut1');set(gca, 'ydir', 'normal');hold on;colorbar;axis equal; clim([-50 -1]); hold on

% point = [-65, 18.22];
% plot(point(1), point(2), 'r*')
% [latval, latindex] = min(abs(lats - point(:,2))); %find the lat point that matches up best with my point of interest (smallest diff b/t lat value and my point)
% [lonval, lonindex] = min(abs(lons - point(:,1)));
% chooseme = [lonindex(1), latindex(1)];
% plot(lons(lonindex(1)), lats(latindex(1)), 'g*'); %you can see discrepancy b/t release point (in this case it was random) and hydrodynamics

%for loop covering all release points
% https://stackoverflow.com/questions/13578279/find-closest-point-in-matlab-grid

latsforindexing = lats'; %get the dimensions right so 'min' function works as intended. NOTE - issue here
relptswdepth = zeros(length(relpts), 3);
relptswdepth(relptswdepth == 0) = NaN;
% for i = 1:10 %for testing

% NOTE - for some reason it breaks at release point 5195 currently. wth ???
% 
for i = 1:length(relptswdepth)
    % i = randi([1 length(relptswdepth)]); %for testing

    lats = lats'; %get the dimensions right so 'min' function works as intended. NOTE - issue here


    point = [relpts(i, 11), relpts(i, 12)];
    plot(point(1), point(2), 'r*');hold on
    [lonval, lonindex] = min(abs(lons - point(:,1)));
    [latval, latindex] = min(abs(latsforindexing - point(:,2))); % NOTE - issue here
    closestcoords = [lonindex(1), latindex(1)];
    plot(lons(lonindex(1)), lats(latindex(1)), 'g*');hold on % NOTE - issue here
    relptswdepth(i, 1) = closestcoords(1); %index of closest lon to release point
    relptswdepth(i, 2) = closestcoords(2); %index of closest lon to release point
    relptswdepth(i, 3) = depthcut1(latindex(1),lonindex(1)); %depth value at those coordinate indices. lats up to 536; lons up to 716
end


