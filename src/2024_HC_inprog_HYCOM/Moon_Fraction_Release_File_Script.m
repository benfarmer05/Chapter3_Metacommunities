%% Script to create a release file for the CMS DAMSELFISH
% 14 May 2024

%cd '/Users/hannah/Library/CloudStorage/Box-Box/CMS/Matlab/Relfile_codes/'; % make sure youre in the right directory
clear all;clc % clear the workspace

% .nc file format: nest_nestnumber_yyyymmddhhmmss [.nc]
%   Working with 20 years of hydrodynamics: 2004 - 2023 (these are the
%                                                           years we have the moon fraction for)
%       So          'nest_1_20020101000000' through 
%                   'nest_1_20221231000000'

% Release file format:
%   Polygon Longitude Latitude Depth Number Year Month Day Second
%   1       277.2     24.66    0     10     2006 3     1   0
%% Put all the moon phases into one column
% don't need to do this if you already have them into one column
clear;clc

cd '/Users/hannah/Library/CloudStorage/Box-Box/CMS/Matlab/Release_file_scripts/Moon/years/';

for year = 2001:2023
    file_name = sprintf('%d_moon.txt', year);

    data = dlmread(file_name);

    one_column = data(:);

    output_file = sprintf('%d_one.txt', year);

    fid = fopen(output_file, 'w');
    fprintf(fid, '%f\n', one_column);
    fclose(fid);
end

%% FIRST
    % want to generate the moon fraction to tune the release #'s to that
    % this is for DAMSELFISH and they spawn daily, according to the moon
        % highest spawning at 3rd Q lowest at new moon
% the 0.50 fraction appears on day 22 +/- 2 days, tolerance for 3rd quarter
%           moon set at 0.05
% 0.0 fraction acts as day 1 of the cycle, so every 30 days, there is a new
%           cycle starting with 0

% dan's code: 
%moonfract = repmat([0:.05:1 fliplr([.1:.05:.9])],1,10); % Just an example
% Want peak at .75, trough at .25
%larvfract = [zeros(1,26) moonfract]; % starting point
%plot(moonfract,'-r'); hold on
%plot(larvfract,'-b')
% Then you could multiply the larval fraction by a scalar (or something) to arrive at a number
% of larvae (like round(larvfract*100))

% things i try: 
 
% 3rd quarter moon is at 0.50 visibility
    % so peak at 0.50 and min at 0 (new moon)
% period is 30 days avg (29.53)
    % 0.50 is day ~22
    % 0 is day 0?
clear all;clc

cd '/Users/hannah/Library/CloudStorage/Box-Box/CMS/Matlab/Release_file_scripts/Moon/one_columns/'; % make sure we are in the right directory (this has the correct format of moon fraction)

% pull in the moon fraction files and make them into variables 
for i = 2004:2023 % the years that we have moon data for
    filename = sprintf('%d_one.txt', i); % the file name is based on the text file name
    try 
        eval(['moon_', num2str(i), ' = readmatrix(filename);']); % dynamically generating variable names and assigning values to those variables
        % ['matrix_', num2str(i)] : dynamically generates a string
        %                           representing the variable name  (if i is 2004, it generates the
        %                           string moon_2004
        % readmatrix(filename) : reads the moon data from the file
        %                        specified by 'filename '
    catch
        fprintf('Error reading file: %s\n', filename); % error message if something isn't working
        continue;
    end
end
clear("filename","i");
%%
% use the moon matrices to generate larval # matrices
        % shifting the larval phase by 7 as that is the average time
        %           between moon phases
for i = 2004:2023
    if mod(i, 4) == 0 && (mod(i, 100) ~= 0 || mod(i, 400) == 0); % looking for leap years using math (described more below)
       days_in_year = 366;
   else
       days_in_year = 365;
   end

   x = (1:days_in_year)'; % creates vector x with the # of days in each year

   moon_var = sprintf('moon_%d', i); % constructs a string using sprintf by concatenating the string "moon_" with the current value of i
   moon_year = eval(moon_var); % evaluates the string as a variable name and retrieves its value

   moon_year = [moon_year, x]; % appends the vector x to the end of the variable 'moon_year' (combines data)

   assignin('base', moon_var, moon_year); % assigns the value of moon_year to moon_var (this is where the matrix name comes from)

   larva_var = sprintf('larva_%d', i); % same string variable construction
   larva_year = moon_year;  % copies the content of moon_year to larva_year
   larva_year(:,2) = larva_year(:,2)+7; % 7 is added to every element in the second column of larva_year
   larva_year = larva_year(:,[2 1]); % rearranges the columns so that the day values are x and the larvae #'s are y 

   assignin('base', larva_var, larva_year); % assigns the values to the correct variable
end 

% get rid of the day value column that was added to the moon matrices
for i = 2004:2023
    moon_var = sprintf('moon_%d', i); % string value
    moon_year = evalin('base',moon_var); % retrieves the value of the variable whose name is stored in 'moon_var'
        % evalin evaluates the expression 'base' which indivates the
        %        variable is to be fetched from the base workspace and then
        %        evaluates 'moon_var'

    moon_year(:,2) = []; % deletes the second column in the moon_year matrix
    assignin('base', moon_var, moon_year); % assigns the correct variables to the correct name
end


for i = 2004:2023
    larva_var = sprintf('larva_%d', i);
    larva_matrix = evalin('base', larva_var);


    z_var = sprintf('z_%d', i);
    z_year = larva_matrix(:,1);
    assignin('base', z_var, z_year);
    
    y_var = sprintf('y_%d', i);
    y_year = larva_matrix(:,2);
    assignin('base', y_var, y_year);
end
%%
% FOR PLOTTING
% creating z and y variables to plot the moon and larval fractions
%           throughout a year
for i = 2004:2006
   z_var = sprintf('z_%d', i);
   y_var = sprintf('y_%d', i);
   z = evalin('base', z_var);
   y = evalin('base', y_var);

   figure;
   plot(z,y); hold on
   plot(eval(['moon_', num2str(i)]));
   title(['Plot %d', num2str(i)]);
   legend('Larval Moon Fraction','Day of Year');
end


%%
%cd '/Users/hannah/Library/CloudStorage/Box-Box/CMS/Fish Model/Input Files/Release_Files/'; % make sure we are in the right directory

% so now have the values (relatively) for larval fractions
    % if releasing 50k a year, want to release 4167 a month
        % so release 140 particles a day
%larva_2004 = round(larva_2004 * 300); % that gives around 4566 particles for the month of Feb (is that okay?)
% 4425 for the month of March
% 55722 for the year
% scalar to bring larval #'s up to desired model values (should be changed)
for i = 2004:2023
   y_var = sprintf('y_%d', i);
   y_year = evalin('base', y_var);
   y_year = round(y_year * 300); % multiplies every value of the y_%d matrix by 300

   new_name = sprintf('scaled_larva_%d', i);
   assignin('base', new_name, y_year);

   % filename = sprintf('larvae_values_%d.txt', i); % sprintf: generates a string formatted according to the specified format
    % the %d placeholder in the format is replaced by the value of the year
    % variable
   % writematrix(y_year, filename, 'delimiter', '\t'); 
end
% so now have the larval values for each year in a txt file
% want to pull those values back into matlab and put them into the release
%           file
        % don't think i had to save them, but just in case


% replicates every scaled_larva_%d matrix by 19 to allign with the 19
%           release locations for the FGB
rep_factor = 19;
for i = 2004:2023
    scaled_var = sprintf('scaled_larva_%d', i);
    scaled_year = evalin('base', scaled_var);
    scaled_year = repelem(scaled_year, rep_factor);

    new_name = sprintf('scaled_replicated_%d', i);
    assignin('base', new_name, scaled_year);

    all_data.(new_name) = scaled_year;
end

save('scaled_replicated_all_years.mat','-struct','all_data');



%% RELEASE FILE
cd '/Users/hannah/Library/CloudStorage/Box-Box/CMS/Fish_Model/Input_Files/Damselfish/Damselfish_Release/'; % make sure we are in the right directory


rep_factor = 17; % how many replicates (locations) are there per day in the index file
num_years = 23; % number of years that we are trying to create release files for
depth_values = [14;11;90;12;41;50;41;54;15;41;56;58;54;64;48;68;62];


% loop through each year
for year = 2001:2023 %+ num_years - 1 % calculated range of years that we want to loop over
    leap_year = mod(year, 4) == 0 && (mod(year, 100) ~= 0 || mod(year, 400) == 0); % checks whether or not it is a leap year
    % mod(year, 4) == 0: checks if the  year is divisible by 4 without any
    %                    remainder, if it is then it could be a leap year
    % mod(year, 100) ~= 0: checks if the year is not dividisble by 100
    %                      without any remainder. Years divisible by 100 are not leap years
    %                      unless they are also divisible by 400
    % mod(year, 400) == 0: chekcs if the year is divisible by 400 without
    %                      remainder. Years divisible by 400 are always leap years
    % combined conditions with operators && = combines two locgical
    %                                         expressions and returns true if both are true and false otherwise
    % also || = OR, meaning it returns true if at least one expression is
    %           true and false if both are false

    if leap_year 
        num_day_feb = 29; % if the year is a leap year, give february 29 days
    else
        num_day_feb = 28; % otherwise, february has 28 days
    end

    num_day_month = [31, num_day_feb, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; % assigning the number of days to each month and allowing february change depending on leap year
    
    date_matrix = zeros(sum(num_day_month) * rep_factor, 3); %initializes the date matrix with dimensions determined by the total number of days in the months and the repetition factor
        % sum(num_day_month): calculates the total # of days in all months
        %                     of the year
        % num_day_month: is a vector containing the # of days in each month
        % sum(num)day_month)*rep_factor: multiplies total # days by rep
        %                                factor which indicated how many times each day's data will be
        %                                replicated
        % zeros(sum(num_day_month)*rep_factor,3): creates a matrix of 0's
        %                                         with calculated dimensions. The matrix will have as many rows as
        %                                         the total # of days in all the months multiplied by the rep
        %                                         factor, and 3 columns. Each row represents a day, and the column
        %                                         stores the year, month, and day values. 
    repnum = sum(num_day_month) * rep_factor; % calculates total # of data points (repetitions) based on # of dasy in each month and the rep factor
        % rep_factor: variable represents a factor by which each day's data
        %             is replicated or repeated
        % repnum: stores total # of data points needed based on # of days
        %         in each month and the rep factor
    IDs = (1:repnum)'; % creates a vector of sequential IDs ranging from 1 to repnum. The ' at the end transposes the vector from row to a column vector
% longitudes of each location
%   repnum/19: determines how many times to repeat each location 
%   1: specifies the dimension along which to repeat (along rows)
    lons1 = repelem(265.7017, repnum/17, 1)'; % Elvers-1
    lons2 = repelem(266.1638, repnum/17, 1)'; % Elvers-2
    lons3 = repelem(266.3125, repnum/17, 1)'; % Geyer
    lons4 = repelem(266.3963, repnum/17, 1)'; % Horseshoe
    lons5 = repelem(266.7045, repnum/17, 1)'; % Bright
    lons6 = repelem(266.9393, repnum/17, 1)'; % West Flower Garden
    lons7 = repelem(267.3977, repnum/17, 1)'; % Rankin
    lons8 = repelem(266.4796, repnum/17, 1)'; % McGrail-2
    lons9 = repelem(267.5384, repnum/17, 1)'; % Sidner
    lons10 = repelem(266.5504, repnum/17, 1)'; % Parker-2
    lons11 = repelem(266.5475, repnum/17, 1)'; % Parker-1
    lons12 = repelem(267.0999, repnum/17, 1)'; % McGrail-1
    lons13 = repelem(267.5461, repnum/17, 1)'; % East Flower Garden
    lons14 = repelem(267.6262, repnum/17, 1)'; % Rezak
    lons15 = repelem(267.6239, repnum/17, 1)'; % MacNeil
    lons16 = repelem(267.9719, repnum/17, 1)'; % Bouma
    lons17 = repelem(267.9965, repnum/17, 1)'; % Alderdice
    % lons18 = repelem(265.7021, repnum/19, 1)'; % Stetson
    % lons19 = repelem(267.5404, repnum/19, 1)'; % Sonnier
% concatenate all the lons into one vertical column vector variable
    lons = vertcat(lons1,lons2,lons3,lons4,lons5,lons6,lons7,lons8,lons9,lons10,lons11,lons12,lons13,lons14,lons15,lons16,lons17); % vertically concatenate to bring two locations into one array
    lats1 = repelem(28.1647, repnum/17, 1)'; % Elvers-1
    lats2 = repelem(27.8713, repnum/17, 1)'; % Elvers-2
    lats3 = repelem(27.8331, repnum/17, 1)'; % Geyer
    lats4 = repelem(27.9310, repnum/17, 1)'; % Horseshoe
    lats5 = repelem(27.8915, repnum/17, 1)'; % Bright
    lats6 = repelem(27.8213, repnum/17, 1)'; % West Flower Garden
    lats7 = repelem(27.9655, repnum/17, 1)'; % Rankin
    lats8 = repelem(28.0031, repnum/17, 1)'; % McGrail-2
    lats9 = repelem(28.3377, repnum/17, 1)'; % Sidner
    lats10 = repelem(27.9131, repnum/17, 1)'; % Parker-2
    lats11 = repelem(27.8978, repnum/17, 1)'; % Parker-1
    lats12 = repelem(27.8274, repnum/17, 1)'; % McGrail-1
    lats13 = repelem(28.0581, repnum/17, 1)'; % East Flower Garden
    lats14 = repelem(27.9692, repnum/17, 1)'; % Rezak
    lats15 = repelem(27.9278, repnum/17, 1)'; % MacNeil
    lats16 = repelem(27.9446, repnum/17, 1)'; % Bouma
    lats17 = repelem(28.0838, repnum/17, 1)'; % Alderdice
    % lats18 = repelem(28.166, repnum/19, 1)'; % Stetson
    % lats19 = repelem(28.341, repnum/19, 1)'; % Sonnier
% concatenate all the lats into one vertical column vector variable
    lats = vertcat(lats1,lats2,lats3,lats4,lats5,lats6,lats7,lats8,lats9,lats10,lats11,lats12,lats13,lats14,lats15,lats16,lats17);
    
    % inialize the release matrix with dimensions repnum rows and 9 columns for
    %   the structure of the release file
    release = zeros(repnum, 9);
    
% Repeat the depth values for the total rows
    repeated_depths = repmat(depth_values, repnum / length(depth_values), 1);

    % used to keep track of the current day within the loop, allowing the
    %       code to properly index and assign data to the correct rows in the
    %       release matrix
    day_counter = 1;

% loop through each month
    for month = 1:12 % for months 1 - 12
        num_day_in_month = num_day_month(month); % retrieves the # of days in the current month using the num_day_month vector
        % month variable is used as an idex to access the corresponding
        %       element in the vector which holds the # of days of the month
        for day = 1:num_day_in_month % this loop iterates over the values from 1 to num_day_in_month representing the days in the current month
            start_index = (day_counter - 1) * rep_factor + 1; % calculates the start index for the current day's data in the release matrix
            % based on the current value of day_counter and rep_factor
            end_index = start_index + rep_factor - 1; %calculates end index for the current days data in the release matrix
            % based on start index and rep_factor, ensuring the correct #
            %       of rows is selected for the current day's data
            %larval_data = evalin('base', sprintf('scaled_replicated_%d', year));

%larval_data = evalin('base', sprintf('scaled_replicated_%d', year + num_years - 1));
%fprintf('Year: %d, Elements in larval_data: %d\n', year, numel(larval_data)); % Debugging output
%if numel(larval_data) < repnum
%    error('larval_data for year %d is empty or does not contain enough elements for %d repetitions.', year, repnum);
%end

%            if end_index > numel(larval_data)
%                end_index = numel(larval_data) % adjusts end index to the max value needed
%            end

            %date_matrix = zeros(num_day_in_month * rep_factor, 3); % creates matrix calculated based on # of day sin current month and the rep factor with 3 columns
                % maybe don't need this code twice? not sure

            %date_matrix(start_index:end_index, :) = repmat([year, month, day], rep_factor, 1); % assigns a repeated patten of year, month, day values to a subset of rows in the date_matrix
                % date_matrix(start_index:end_index, :): selects a subset
                %                                        of rows in the matrix starting from the start_index row
                %                                        and ending at the end_index row
                    % : means that all columns are selected
                % repmat([year,month,day],rep_factor,1): replicates the
                %                                        [year, month, day] vector rep factor times along the rows
                %                                        and 1 time along the columns
                    % creates a matrix where each row contains the year,
                    %       month, and day values repeated rep factor times
% prepare data for release matrix
            ids_section = reshape(IDs(start_index:end_index), [], 1); % reshapes IDs from start to end_index into a column vector where the empty [] argument indicates that the # of rows is 
            % automatically calculated based on the # of elements in the selectede subset, and 1 specifies that there is one column 
            lons_sec = reshape(lons(start_index:end_index), [], 1); % same process applied to lons and lats
            lats_sec = reshape(lats(start_index:end_index), [], 1);
            %larva_sec = larval_data((day - 1) * rep_factor + 1 : day * rep_factor);
            %larva_sec = larval_data((day - 1) * rep_factor + 1 : min(day * rep_factor, numel(larval_data)));
            %if numel(larva_sec) < rep_factor
            %    larva_sec = [larva_sec; zeros(rep_factor - numel(larva_sec), 1)];
            %end
% fill the release matrix with location information            
            release(start_index:end_index,1) = ids_section; % assigns values from ids column vector to 1st column of the release matrix, starting from start and ending with end_index row
                % assigns IDs to specified subset of rows in the first
                %       column of release matrix
            release(start_index:end_index,2) = lons_sec;
            release(start_index:end_index,3) = lats_sec;
            release(start_index:end_index,4) = repeated_depths(start_index:end_index); %repeated_depths(start_index:end_index);
            release(start_index:end_index,5) = 43;%larva_sec;
            release(start_index:end_index,6) = year;
            release(start_index:end_index,7) = month;
            release(start_index:end_index,8) = day;
            release(start_index:end_index,9) = 0;
                
            day_counter = day_counter + 1; % this makes sure that day_counter is incremented by 1, ensuring the loop progresses to the next day
        end
    end

    filename = sprintf('damsel_relfile_%d.txt', year); % sprintf: generates a string formatted according to the specified format
    % the %d placeholder in the format is replaced by the value of the year
    % variable
    writematrix(release, filename, 'delimiter', '\t'); % writes the release file in a text format
end 