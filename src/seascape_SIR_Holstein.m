clear; clc

%% Initialize paths
projectPath = matlab.project.rootProject().RootFolder;
tempPath = fullfile(projectPath, 'temp');
dataPath = fullfile(projectPath, 'data');
outputPath = fullfile(projectPath, 'output');

%% Initial Conditions

load reefs_new.mat

%% Create distance matrix - unncessary when connectivity matrices are

% Create less densly connected matrix
% Lat/long (centroids) of habitats
XY = [reefs(:,3) reefs(:,4)];

habs = length(XY); % The number of patches

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % This creates a distance matrix - Not used
% [D2, Wraw] = buildConnectivityKNN(XY(:,1), XY(:,2), 'K', 30, 'KJitter',10,...
%     'LenScaleKm', 25, 'DropProb', 0.05, ...
%     'WeightNoiseCV', 0.25, 'RowNormalize', true);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Load some of Ben's conn data
% load connectivity_2019_Jan01_120000.mat
% D2 = connectivity_results.ConnMatrix_normalized;
% D2(1:size(D2,1)+1:end) = 0; % Remove retention

myFiles = dir("C_*");

% Create connectivity matrices - the below example has months/years - THIS
% NEEDS TO BE UPDATED WITH YOUR MATRICES, BEN!
P = struct();
norm_val = 65*sum(decay_weights); % This will be the numpart*sum(decay)
for i = 1:length(myFiles)
    filenam = myFiles(i).name;
    Con = load(filenam);
    Dat = Con.connectivity_results.calendar_date;
    Mo = month(Dat);
    Dy = day(Dat);
    Yr = year(Dat);

    P(i).DY = Dy;
    P(i).MO = Mo;
    P(i).YR = Yr;
    P(i).Date = Dat;
    conmat = sparse(Con.connectivity_results.ConnMatrix_raw);
    % Still need to remove self retention
    conmat = spdiags(zeros(size(conmat,1),1), 0, conmat);
    conmat = conmat ./ norm_val;

    P(i).full = conmat;
end
    
save P.mat P -v7.3


% P(1).MO = 12;
% P(1).YR = 2018;
% P(1).full = sparse(D2);
% for i = 2:13
%     P(i).MO = i-1;
%     P(i).YR = 2019;
%     P(i).full = sparse(D2);
% end
% for i = 14:25
%     P(i).MO = i-13;
%     P(i).YR = 2020;
%     P(i).full = sparse(D2);
% end

% Define day vector for the ODE - this will change depending on the
% periodicity of connectivity matrices - but the goal is a vector of days
% that indicates where matrices START
% Yrs = [P.YR];
% Mos = [P.MO];
% Dys = [P.DY];
[P.date] = P.Date;
dates = [P.date];
% Pdates = datetime(Yrs,Mos,Dys, ones(size(Yrs)));
ref = datetime(2019,1,1);  
Pdays = days(dates - ref);

%% Initial conditions

% Total coral cover @ each site
CC = reefs(:,5);  

% Initial susceptiple coral cover
CLS1 = reefs(:,6);%*area;
CMS1 = reefs(:,7);%*area;
CHS1 = reefs(:,8);%*area;  

% Initial healthy coral cover at every hab 
SLS1 = CLS1;
SMS1 = CMS1;
SHS1 = CHS1;

% Initial diseased coral at every hab
ILS1 = zeros(habs,1);
IMS1 = zeros(habs,1);
IHS1 = zeros(habs,1);

% Where did the disease start? These indices were provided by BEN.
% Consider CHANGING the initial values of disease. Here I am seeding BOTH
% HS and LS corals. 
% 
% ****NOTE**** It DOES matter if you start the disease in a
% single FLAT grid cell vs 2 or 3 or 4. These cells appear to have
% different connectivity in the first weeks of 2019.

Flat = [find(reefs(:,2)==29088) find(reefs(:,2)==29338) find(reefs(:,2)==29089) find(reefs(:,2)==29339) find(reefs(:,2)==29087)]; % Flay Cay
IHS1(Flat) = .01*CHS1(Flat); % Likely sensitive to how much disease you set initially
SHS1(Flat) = SHS1(Flat)-IHS1(Flat);
IMS1(Flat) = .01*CMS1(Flat); % Likely sensitive to how much disease you set initially
SMS1(Flat) = SMS1(Flat)-IMS1(Flat);
ILS1(Flat) = .01*CLS1(Flat); % Likely sensitive to how much disease you set initially
SLS1(Flat) = SLS1(Flat)-ILS1(Flat);

% Initial dead cover - this used to be "recovered" (hence the R)
RLS1 = zeros(habs,1);
RMS1 = zeros(habs,1);
RHS1 = zeros(habs,1);

%% Parameters
% Most parameters should be unique for each category of coral. This is the
% meat and potatoes of disease dynamics at each site. Note that rates are
% sensitive to time step.
 
% Infection, beta - Ben's paper, multihost, NS
bls = 0.03;
bms = 0.14;
bhs = 2.08;

% Recovery rate, gamma - Ben's paper, multihost, NS
kls = .05; %
kms = .55;% 
khs = 3.33;% 

% Mortality rate of diseased tissue in "recovered" state. This rate is different for each coral
% GROUP. This will be parameterized empirically. THOMAS - 0.1431 - this
% would actually be recovery rate
% mls = 1;
% mms = 1;
% mhs = 1;

% Loss of protection rate - Move from "recovered" back to susceptible %
% Turning this off
% rls = 0;
% rms = 0;
% rhs = 0;

% A (minimum) threshold amount of diseased tissue (cover) a site must have in order
% to transmit disease to other sites. May need to be OPTIMIZED. Or not
% used.
% Some thoughts - if an outgoin or incoming threshold is not used, sites
% can get infected very quickly, even if at very low magnitudes. One
% approach might be to think about % of bottom... another might be to think
% of the stage of a local outbreak...
thresh = 0.0003; % At .001, transmission from Flat is unlikely with current parameterization. At .0005 disease immediatelyy begins to transmit but not super duper fast...

% NO LONGER TRUE -- A (minimum) threshold for incoming disease probability to result in local infection.
% thresh2 = .00000001;

% I've added a reshape that controls the input of T to local DP. 
% c.*(1-exp(-B.*T(:)))/(1-exp(-B));
c = 1; % limits max, ranges 0:1
B = -4; % determines curve. <0 is concave, >0 is convex; set small (.001) for no (linear) reshape. 0 will produce Inf and result in no transmission.

% Save all parameters to a vector for use later...
% pars = [bls; bms; bhs; kls; kms; khs; thresh; sh1; sh2];
pars_simp = [bls; bms; bhs; kls; kms; khs];

%% The ODE Solver - full

% 10/25 - DH - % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% I have updated this ODE - ODEfun_07.m - to use the functions in Ben's
% paper.
%
% I have also updated this code to use a sparse connectivity matrix, which
% has dramtically improved efficiency, particularly when disease
% connectivity (active disease sites) is sparse. It slows quite a lot when
% disease is active at many sites (~>3000).
%
% Note that when incoming disease is not enough to initiate local infection
% ANYWHERE based on thresh2, a warning will display. This is not
% necessarily a bad thing, even if it occurs intermittently. Of course, if
% disease nevery moves in space or time, that's not great.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Initial conditions to supply to the ODE solver. These are arranged in a
% single vertical vector. The order matters! Notice that dead proportion is
% ignored until the end. It is 1-(S+I+R).
Y0 = [SLS1;SMS1;SHS1;ILS1;IMS1;IHS1;RLS1;RMS1;RHS1]; 
habs = length(XY);
save inits_01.mat Y0

% If not using odeWaitbar, turn this off
opts  = odeset('OutputFcn', @odeWaitbar);

% The time span to integrate in days
tspan = [1 365]; 

% The ODE solver. Note that ODEfun_07 is a separate .m file 
clear Y t
tic
[t,Y] = ode45(@(t,Y) ODEfun_08(t,Y,CLS1,CMS1,CHS1,bls,bms,bhs,kls,kms,khs,thresh,c,B,habs,P,Pdays), tspan, Y0, opts);
toc

save RealRun_0001.mat Y Y0 CLS1 CMS1 CHS1 t pars

%% I don't currently have anything to fit against
% Making a loop for at least thresh

threshvec = [.0001:.0002:.0009];
Bvec = [-4:2:4];

count = 0;
Results = struct();
tic
for tv = 1:length(threshvec)
    for bv = 1:length(Bvec)
        count = count+1;

        Y0 = [SLS1;SMS1;SHS1;ILS1;IMS1;IHS1;RLS1;RMS1;RHS1]; 
        habs = length(XY);
        opts  = odeset('OutputFcn', @odeWaitbar);
        tspan = [1 731]; 
        clear Y t
        
        [t,Y] = ode45(@(t,Y) ODEfun_08(t,Y,CLS1,CMS1,CHS1,bls,bms,bhs,kls,kms,khs,threshvec(tv),c,Bvec(bv),habs,P,Pdays), tspan, Y0, opts);
        

        LSHP = Y(:,1:habs);
        MSHP = Y(:,habs+1:habs*2);
        HSHP = Y(:,2*habs+1:habs*3);
        
        LSIP = Y(:,3*habs+1:habs*4);
        MSIP = Y(:,4*habs+1:habs*5);
        HSIP = Y(:,5*habs+1:habs*6);
        
        LSRP = Y(:,6*habs+1:habs*7);
        MSRP = Y(:,7*habs+1:habs*8);
        HSRP = Y(:,8*habs+1:habs*9);
        
        interpLSHP = interp1(t,LSHP,tspan(1):tspan(2));
        interpMSHP = interp1(t,MSHP,tspan(1):tspan(2));
        interpHSHP = interp1(t,HSHP,tspan(1):tspan(2));
        
        interpLSIP = interp1(t,LSIP,tspan(1):tspan(2));
        interpMSIP = interp1(t,MSIP,tspan(1):tspan(2));
        interpHSIP = interp1(t,HSIP,tspan(1):tspan(2));
        
        interpLSRP = interp1(t,LSRP,tspan(1):tspan(2));
        interpMSRP = interp1(t,MSRP,tspan(1):tspan(2));
        interpHSRP = interp1(t,HSRP,tspan(1):tspan(2));
        
        % Just in case
        interpLSHP(interpLSHP<0) = 0;
        interpMSHP(interpMSHP<0) = 0;
        interpHSHP(interpHSHP<0) = 0;
        
        interpLSIP(interpLSIP<0) = 0;
        interpMSIP(interpMSIP<0) = 0;
        interpHSIP(interpHSIP<0) = 0;
        
        interpLSRP(interpLSRP<0) = 0;
        interpMSRP(interpMSRP<0) = 0;
        interpHSRP(interpHSRP<0) = 0;

        Results(count).thresh = threshvec(tv);
        Results(count).B = Bvec(bv);
        Results(count).LSS = interpLSHP;
        Results(count).MSS = interpMSHP;
        Results(count).HSS = interpHSHP;

        Results(count).LSI = interpLSIP;
        Results(count).MSI = interpMSIP;
        Results(count).HSI = interpHSIP;

        Results(count).LSR = interpLSRP;
        Results(count).MSR = interpMSRP;
        Results(count).HSR = interpHSRP;
    end
end
        
toc

T = tiledlayout(5,5,'TileSpacing','compact','Padding','compact');

for i=1:length(Results)
    RTIP = Results(i).LSI + Results(i).MSI + Results(i).HSI;
    nexttile(T,i)
    plot(RTIP)
    title(strcat('thresh = ',num2str(Results(i).thresh),', B = ',num2str(Results(i).B)))
end


%% Need to interpolate results
% The output from the solver is not in days, but is in unequal time steps.
% You need to interpolate back to days..
LSHP = Y(:,1:habs);
MSHP = Y(:,habs+1:habs*2);
HSHP = Y(:,2*habs+1:habs*3);

LSIP = Y(:,3*habs+1:habs*4);
MSIP = Y(:,4*habs+1:habs*5);
HSIP = Y(:,5*habs+1:habs*6);

LSRP = Y(:,6*habs+1:habs*7);
MSRP = Y(:,7*habs+1:habs*8);
HSRP = Y(:,8*habs+1:habs*9);

interpLSHP = interp1(t,LSHP,tspan(1):tspan(2));
interpMSHP = interp1(t,MSHP,tspan(1):tspan(2));
interpHSHP = interp1(t,HSHP,tspan(1):tspan(2));

interpLSIP = interp1(t,LSIP,tspan(1):tspan(2));
interpMSIP = interp1(t,MSIP,tspan(1):tspan(2));
interpHSIP = interp1(t,HSIP,tspan(1):tspan(2));

interpLSRP = interp1(t,LSRP,tspan(1):tspan(2));
interpMSRP = interp1(t,MSRP,tspan(1):tspan(2));
interpHSRP = interp1(t,HSRP,tspan(1):tspan(2));

% Just in case
interpLSHP(interpLSHP<0) = 0;
interpMSHP(interpMSHP<0) = 0;
interpHSHP(interpHSHP<0) = 0;

interpLSIP(interpLSIP<0) = 0;
interpMSIP(interpMSIP<0) = 0;
interpHSIP(interpHSIP<0) = 0;

interpLSRP(interpLSRP<0) = 0;
interpMSRP(interpMSRP<0) = 0;
interpHSRP(interpHSRP<0) = 0;

TIP = interpLSIP + interpMSIP + interpHSIP;
TSP = interpLSHP + interpMSHP + interpHSHP;
TRP = interpLSRP + interpMSRP + interpHSRP;
%% Make results into a movie - 2025
% Just pumping something out for Ben 

% If you have a coastal poly... this code is from a FL application
% Coast = shaperead('JustUS_GEO_02.shp');
% % [xt, yt] = find(interpTIT>thresh, 1, 'last'); % Find the last step where any disease is more than a threshold, like 0.01
% % v = VideoWriter(strcat('DisVid',num2str(now),'.avi'));

% Temporary
TIP = Results(6).LSI + Results(6).MSI + Results(6).HSI;
TSP = Results(6).LSS + Results(6).MSS + Results(6).HSS;
TRP = Results(6).LSR + Results(6).MSR + Results(6).HSR;

f = figure('renderer', 'zbuffer','Position', [10 10 1400 1000]);

set(f,'nextplot','replacechildren'); 

cmap3 = colormap(flipud(autumn(round(max(max(TIP))*1000000)+1)));
cmap3(1,:) = [0 .3 1];
cmap4 = colormap(cool(round(max(CC)*1000)+1));
cmap4(1,:) = [0 .3 1];

edges = [0 .000000001 .001 .1 1];   % mm/hr
Cstart = [
    0.25 .5 0.25;   % light green
    1 1 .2;   % green
    1.00 0.60 0.00;   % yellow
    0.80 0.00 0.80;   % orange   % red
];
Cend = [
    0.00 0.60 0.00;   % green
    1.00 0.60 0.00;   % yellow
    1.00 0.00 0.00;   % orange
    0 0 0;   % red     % magenta
];
cmap_I = stackedColormap(edges, Cstart, Cend, 256, [1 1 1 1]);

edges = [0 .000000001 .001 .1 1];   % mm/hr
Cstart = [
    0.25 .5 0.25;   % light green
    1 1 .2;   % green
    1.00 0.60 0.00;   % yellow
    0.80 0.00 0.80;   % orange   % red
];
Cend = [
    0.00 0.60 0.00;   % green
    1.00 0.60 0.00;   % yellow
    1.00 0.00 0.00;   % orange
    0 0 0;   % red     % magenta
];
cmap_R = stackedColormap(edges, Cstart, Cend, 256, [1 1 1 1]);

T = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
% If you have a coastal polygon, plot first (before loop).
axv1 = nexttile(T,1);
    h1 = scatter(axv1, XY(:,1),XY(:,2),7,TIP(1,:)','filled');
    colormap(axv1,cmap_I);
    clim([0 .01]);
    c1 = colorbar(axv1);
    axis equal
axv2 = nexttile(T,2);
    h2 = scatter(axv2, XY(:,1),XY(:,2),7,TIP(1,:)'./CC(:),'filled');
    colormap(axv2,cmap_I);
    clim([0 .01]);
    c2 = colorbar(axv2);
    axis equal 
    
axv3 = nexttile(T,3);
    h3 = scatter(axv3, XY(:,1),XY(:,2),7,CC(:)-TSP(1,:)','filled');
    c3 = colorbar(axv3);
    clim([0 1]);
    colormap(axv3,cmap_R)
    axis equal
axv4 = nexttile(T,4);
    h4 = scatter(axv4, XY(:,1),XY(:,2),7,(CC(:)-TSP(1,:)')./CC(:),'filled');
    c4 = colorbar(axv4);
    clim([0 1]);
    colormap(axv4,cmap_R)
    axis equal


v = VideoWriter('DisVid_thresh0003B-4');
v.FrameRate = 5;
open(v);

% Threshold for which sites to include when drawing a non-convex hull.
% Logic here is that very low disease incidence would not be easily
% observable... but could set it to 0.
bthresh = 0.001; % .1% of bottom
for k = 1:1:tspan(end)
    Dt=datetime('1-Jan-2019')+k-1;
    
    sick = TRP(k,:) >= bthresh; % arguably it's the DEAD corals that are observed, not the sick corals
    dfXY = XY(sick,:);
    df = boundary(dfXY(:,1),dfXY(:,2),0.7);
    
    % Total infected tissue
    set(h1, 'CData', TIP(k,:).');
    % plot(axv1,dfXY(df,1),dfXY(df,2),'r-');
    p1 = patch(axv1,dfXY(df,1),dfXY(df,2),[.8 .9 1],'EdgeColor','r','FaceAlpha',0.2);
    % h1 = scatter(axv1,XY(:,1),XY(:,2),10,TIP(k,:),'filled');
    title(axv1,'Disease prevalence - proportion of bottom',strcat('t ='," ",datestr(Dt)));
    colormap(axv1,cmap_I);
    clim([0 .01]);
    c1 = colorbar(axv1);
    axis equal

      % Total infected tissue
    set(h2, 'CData', TIP(k,:).'./CC(:));
    % plot(axv2,dfXY(df,1),dfXY(df,2),'r-');
    p2 = patch(axv2,dfXY(df,1),dfXY(df,2),[.8 .9 1],'EdgeColor','r','FaceAlpha',0.2);
    % h1 = scatter(axv1,XY(:,1),XY(:,2),10,TIP(k,:),'filled');
    title(axv2,'Disease prevalence - proportion of living coral',strcat('t ='," ",datestr(Dt)));
    colormap(axv2,cmap_I);
    clim([0 .1]);
    c2 = colorbar(axv2);
    axis equal
     
    % Coral cover loss
    set(h3, 'CData',(CC(:)-TSP(k,:).'));
    % plot(axv3,dfXY(df,1),dfXY(df,2),'r-');
    p3 = patch(axv3,dfXY(df,1),dfXY(df,2),[.8 .9 1],'EdgeColor','r','FaceAlpha',0.2);
    title(axv3,'Total coral cover lost',strcat('t ='," ",datestr(Dt)));
    c3 = colorbar(axv3);
    clim([0 1]);
    colormap(axv3,cmap_R)
    axis equal

     % Coral cover loss
    set(h4, 'CData',(CC(:)-TSP(k,:).')./CC(:));
    % plot(axv4,dfXY(df,1),dfXY(df,2),'r-');
    p4 = patch(axv4,dfXY(df,1),dfXY(df,2),[.8 .9 1],'EdgeColor','r','FaceAlpha',0.2);
    title(axv4,'Proportion coral cover lost',strcat('t ='," ",datestr(Dt)));
    c4 = colorbar(axv4);
    clim([0 1]);
    colormap(axv4,cmap_R)
    axis equal
    
    F = getframe(f);
    writeVideo(v, F);

    delete(p1)
    delete(p2)
    delete(p3)
    delete(p4)
end
close(v); 


%% Just a spacer, below code is older and may or may not be useful
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Just looking at the I curves at random spots. In locations where
% outbreaks don't really happen, it gets pretty interesting.
figure; 
for i = 1:1000 
    plot(TIP(:,randi([1 habs])))
    % axis([1 tspan(2) 0 .0001])
    pause(.75); 
    i = i+1; 
end

figure; 
for i = 1:1000 
    loc = randi([1 habs]);
    plot(interpHSIP(:,loc),'-r'); hold on
    plot(interpLSIP(:,loc),'-g');
    plot(interpMSIP(:,loc),'-b');
    axis([1 tspan(2) 0 .0001])
    hold off
    pause(1); 
    i = i+1; 
end

figure; 
for i = 1:1000 
    loc = randi([1 habs]);
    plot(interpHSHP(:,loc),'-r'); hold on
    plot(interpLSHP(:,loc),'-g');
    plot(interpMSHP(:,loc),'-b');
    plot(interpHSRP(:,loc),'-r');
    plot(interpLSRP(:,loc),'-g');
    plot(interpMSRP(:,loc),'-b');
    axis([1 tspan(2) 0 .05])
    hold off
    pause(1); 
    i = i+1; 
end

figure('Position', [10 10 1000 500]); 
T = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for i = 1:1000 
    loc = randi([1 habs]);
    nexttile(T,1)
        plot(interpHSHP(:,loc),'Color',[1 0 0]); hold on
        plot(interpLSHP(:,loc),'Color',[.1 1 .1]);
        plot(interpMSHP(:,loc),'Color',[0 0 1]);
        plot(interpHSRP(:,loc),'Color',[1 .2 .2]);
        plot(interpLSRP(:,loc),'Color',[.3 .9 .3]);
        plot(interpMSRP(:,loc),'Color',[.2 .2 1]);
        % plot(interpHSIP(:,loc).*10,'Color',[1 0 0]);
        % plot(interpLSIP(:,loc).*10,'Color',[0 1 0]);
        % plot(interpMSIP(:,loc).*10,'Color',[0 0 1]);
        axis([1 tspan(2) -.005 .2])
        hold off

    nexttile(T,2)
        scatter(XY(:,1),XY(:,2),5,'blue','filled'); hold on
        scatter(XY(loc,1),XY(loc,2),25,'red','filled');
        axis equal
        hold off
    pause(2); 
    i = i+1; 
end


%% Some plotting - there are all sorts of plots below - some will be useful, some will not, and some will just confuse you
% Just dummy variables to name outputs
cur = clock;
now = cur(6);

% Save parameters
% csvwrite(strcat('Pars',num2str(now)),pars);



% Total healthy tissue
THT = Y(:,1:habs).*CLS1+Y(:,habs+1:habs*2).*CMS1+Y(:,2*habs+1:habs*3).*CHS1;
% Total infectious tissue
TIT = Y(:,3*habs+1:habs*4).*CLS1+Y(:,4*habs+1:habs*5).*CMS1+Y(:,5*habs+1:habs*6).*CHS1;
% Total dead tissue
TDT = (1-(Y(:,1:habs)+Y(:,3*habs+1:habs*4)+Y(:,6*habs+1:habs*7))).*CLS1+(1-(Y(:,habs+1:habs*2)+Y(:,4*habs+1:habs*5)+Y(:,7*habs+1:habs*8))).*CMS1+(1-(Y(:,2*habs+1:habs*3)+Y(:,5*habs+1:habs*6)+Y(:,8*habs+1:habs*9))).*CHS1;
% Total recovered tissue
TRT = Y(:,6*habs+1:habs*7).*CLS1 + Y(:,7*habs+1:habs*8).*CMS1 + Y(:,8*habs+1:habs*9).*CHS1;

% Total healthy proportion
THP = (Y(:,1:habs)+Y(:,habs+1:habs*2)+Y(:,2*habs+1:habs*3))/3;
% Total infectioys proportion
TIP = (Y(:,3*habs+1:habs*4)+Y(:,4*habs+1:habs*5)+Y(:,5*habs+1:habs*6))/3;
% Total dead proportion
TDP = ((1-(Y(:,1:habs)+Y(:,3*habs+1:habs*4)+Y(:,6*habs+1:habs*7)))+(1-(Y(:,habs+1:habs*2)+Y(:,4*habs+1:habs*5)+Y(:,7*habs+1:habs*8)))+(1-(Y(:,2*habs+1:habs*3)+Y(:,5*habs+1:habs*6)+Y(:,8*habs+1:habs*9))))/3;
% Total recovered proportion
TRP = (Y(:,6*habs+1:habs*7) + Y(:,7*habs+1:habs*8) + Y(:,8*habs+1:habs*9))/3;

% p1 = plot(t,Y(:,3*habs+1:habs*4).*CLS1, 'color',[.3 .6 .1 .2]); hold on;
% % p1.Color(4) = .2;
% p2 = plot(t,Y(:,4*habs+1:habs*5).*CMS1,'color',[.5 .5 0 .2]);
% % p2.Color(4) = .2;
% p3 = plot(t,Y(:,5*habs+1:habs*6).*CHS1,'color', [1 0 0 .2]); hold off;
% % p3.Color(4) = .2; hold off;
 
    % Absolute cover figure
figure('Position', [10 10 700 900]);
    
subplot(8,3,1:3)
    plot(t,Y(:,1:habs).*CLS1+Y(:,habs+1:habs*2).*CMS1+Y(:,2*habs+1:habs*3).*CHS1);
    title('Total healthy tissue')
    ylim([0 .2])
subplot(8,3,4)
    plot(t,Y(:,1:habs).*CLS1);
    title('LS healthy tissue')
    ylim([0 .2])
subplot(8,3,5)
    plot(t,Y(:,habs+1:habs*2).*CMS1);
    title('MS healthy tissue')
    ylim([0 .2])
subplot(8,3,6)
    plot(t,Y(:,2*habs+1:habs*3).*CHS1);
    title('HS healthy tissue')
    ylim([0 .2])

subplot(8,3,7:9)
    plot(t,Y(:,3*habs+1:habs*4).*CLS1+Y(:,4*habs+1:habs*5).*CMS1+Y(:,5*habs+1:habs*6).*CHS1);
    ylim([0 .2])
    title('Total infected tissue')
subplot(8,3,10)
    plot(t,Y(:,3*habs+1:habs*4).*CLS1, 'color',[.3 .6 .1]);
    ylim([0 .2])
    title('LS infected tissue')
subplot(8,3,11)
    plot(t,Y(:,4*habs+1:habs*5).*CMS1,'color',[.6 .2 0]);
    ylim([0 .2])
    title('MS infected tissue')
subplot(8,3,12)
    plot(t,Y(:,5*habs+1:habs*6).*CHS1,'r');
    ylim([0 .2])
    title('HS infected tissue')
    
subplot(8,3,13:15)
    plot(t,(1-(Y(:,1:habs)+Y(:,3*habs+1:habs*4)+Y(:,6*habs+1:habs*7))).*CLS1+(1-(Y(:,habs+1:habs*2)+Y(:,4*habs+1:habs*5)+Y(:,7*habs+1:habs*8))).*CMS1+(1-(Y(:,2*habs+1:habs*3)+Y(:,5*habs+1:habs*6)+Y(:,8*habs+1:habs*9))).*CHS1);
    ylim([0 .2])
    title('Total dead tissue')
subplot(8,3,16)
    plot(t,(1-(Y(:,1:habs)+Y(:,3*habs+1:habs*4)+Y(:,6*habs+1:habs*7))).*CLS1);
    ylim([0 .2])
    title('LS dead tissue')
subplot(8,3,17)
    plot(t,(1-(Y(:,habs+1:habs*2)+Y(:,4*habs+1:habs*5)+Y(:,7*habs+1:habs*8))).*CMS1);
    ylim([0 .2])
    title('MS dead tissue')
subplot(8,3,18)
    plot(t,(1-(Y(:,2*habs+1:habs*3)+Y(:,5*habs+1:habs*6)+Y(:,8*habs+1:habs*9))).*CHS1);
    ylim([0 .2])
    title('HS dead tissue')
      
subplot(8,3,19:21)
    plot(t,Y(:,6*habs+1:habs*7).*CLS1 + Y(:,7*habs+1:habs*8).*CMS1 + Y(:,8*habs+1:habs*9).*CHS1 ); 
    ylim([0 .2])
    title('Total recovered tissue - now protected')
subplot(8,3,22)
    plot(t,Y(:,6*habs+1:habs*7).*CLS1);
    ylim([0 .2])
    title('LS recovered tissue')
subplot(8,3,23)
    plot(t,Y(:,7*habs+1:habs*8).*CMS1);
    ylim([0 .2])
    title('MS recovered tissue')
subplot(8,3,24)
    plot(t,Y(:,8*habs+1:habs*9).*CHS1);
    ylim([0 .2])
    title('HS recovered tissue')
mtit('Absolute tissue','xoff', 0,'yoff', .04)

% print(strcat('AbsTiss',num2str(now),'.pdf'),'-dpdf','-fillpage')

% Proportions
figure('Position', [10 10 700 900]);

subplot(8,3,1:3)
    plot(t,(Y(:,1:habs)+Y(:,habs+1:habs*2)+Y(:,2*habs+1:habs*3))/3);
    title('Total healthy coral (proportion)')
    ylim([0 1])
subplot(8,3,4)
    plot(t,Y(:,1:habs));
    title('LS healthy proportion')
    ylim([0 1])
subplot(8,3,5)
    plot(t,Y(:,habs+1:habs*2));
    title('MS healthy proportion')
    ylim([0 1])
subplot(8,3,6)
    plot(t,Y(:,2*habs+1:habs*3));
    title('HS healthy proportion')
    ylim([0 1])

subplot(8,3,7:9)
    plot(t,(Y(:,3*habs+1:habs*4)+Y(:,4*habs+1:habs*5)+Y(:,5*habs+1:habs*6))/3);
    ylim([0 1])
    title('Total infected proportion')
subplot(8,3,10)
    plot(t,Y(:,3*habs+1:habs*4));
    ylim([0 1])
    title('LS infected proportion')
subplot(8,3,11)
    plot(t,Y(:,4*habs+1:habs*5));
    ylim([0 1])
    title('MS infected proportion')
subplot(8,3,12)
    plot(t,Y(:,5*habs+1:habs*6));
    ylim([0 1])
    title('HS infected proportion')
    
subplot(8,3,13:15)
    plot(t,((1-(Y(:,1:habs)+Y(:,3*habs+1:habs*4)+Y(:,6*habs+1:habs*7)))+(1-(Y(:,habs+1:habs*2)+Y(:,4*habs+1:habs*5)+Y(:,7*habs+1:habs*8)))+(1-(Y(:,2*habs+1:habs*3)+Y(:,5*habs+1:habs*6)+Y(:,8*habs+1:habs*9))))/3);
    ylim([0 1])
    title('Total dead proportion')
subplot(8,3,16)
    plot(t,(1-(Y(:,1:habs)+Y(:,3*habs+1:habs*4)+Y(:,6*habs+1:habs*7))));
    ylim([0 1])
    title('LS dead proportion')
subplot(8,3,17)
    plot(t,(1-(Y(:,habs+1:habs*2)+Y(:,4*habs+1:habs*5)+Y(:,7*habs+1:habs*8))));
    ylim([0 1])
    title('MS dead proportion')
subplot(8,3,18)
    plot(t,(1-(Y(:,2*habs+1:habs*3)+Y(:,5*habs+1:habs*6)+Y(:,8*habs+1:habs*9))));
    ylim([0 1])
    title('HS dead proportion')
      
subplot(8,3,19:21)
    plot(t,(Y(:,6*habs+1:habs*7) + Y(:,7*habs+1:habs*8) + Y(:,8*habs+1:habs*9))/3 ); 
    ylim([0 1])
    title('Total recovered proportion - now protected')
subplot(8,3,22)
    plot(t,Y(:,6*habs+1:habs*7));
    ylim([0 1])
    title('LS recovered proportion')
subplot(8,3,23)
    plot(t,Y(:,7*habs+1:habs*8));
    ylim([0 1])
    title('MS recovered proportion')
subplot(8,3,24)
    plot(t,Y(:,8*habs+1:habs*9));
    ylim([0 1])
    title('HS recovered proportion')
mtit('Proportions','xoff', 0,'yoff', .04)

% print(strcat('Proportions',num2str(now),'.pdf'),'-dpdf','-fillpage')

% Metrics
p = figure();
firstlast = [THT(1,1:habs)' THT(end,1:habs)'];
subplot(1,2,1);
b1 = boxplot(firstlast,'Labels',{'First','Last'});
ylabel('Coral cover');

ts = repmat(t,1,habs);
tdurD = ts.*(TIT>0);
tdurD(tdurD<=0)=NaN;
durD = max(tdurD)-min(tdurD);
durD(isnan(durD)) = 0;

tdurDT = ts.*(TIT>thresh);
tdurDT(tdurDT<=0)=NaN;
durDT = max(tdurDT) - min(tdurDT);
durDT(isnan(durDT)) = 0;

durs = [durD' durDT'];
subplot(1,2,2);
b2 = boxplot(durs,'Labels',{'Any','Above threshold'});
ylabel('Outbreak duration');

mtit('Metrics','xoff', 0,'yoff', .04);

% print(strcat('ODEmetrics',num2str(now),'.pdf'),'-dpdf')

% Metacommunity level SIRDS
figure
plot(t,mean(THP,2),'g'); hold on
plot(t,mean(TIP,2),'r');
plot(t,mean(TDP,2),'k');
plot(t,mean(TRP,2),'b');

%% Some plotting INITIAL PRE RUN

cur = clock;
now = cur(6);

% Save parameters
% csvwrite(strcat('Pars',num2str(now)),pars);

% Absolute cover figure


CLsimp = CLS1(SimpRunInit(:,1))';
CMsimp = CMS1(SimpRunInit(:,1))';
CHsimp = CHS1(SimpRunInit(:,1))';

habs = length(SimpRunInit);
% Total healthy tissue
THT = Y(:,1:habs).*CLsimp + Y(:,habs+1:habs*2).*CMsimp + Y(:,2*habs+1:habs*3).*CHsimp;
% Total infectious tissue
TIT = Y(:,3*habs+1:habs*4).*CLsimp+Y(:,4*habs+1:habs*5).*CMsimp+Y(:,5*habs+1:habs*6).*CHsimp;
% Total dead tissue
TDT = (1-(Y(:,1:habs)+Y(:,3*habs+1:habs*4)+Y(:,6*habs+1:habs*7))).*CLsimp+(1-(Y(:,habs+1:habs*2)+Y(:,4*habs+1:habs*5)+Y(:,7*habs+1:habs*8))).*CMsimp+(1-(Y(:,2*habs+1:habs*3)+Y(:,5*habs+1:habs*6)+Y(:,8*habs+1:habs*9))).*CHsimp;
% Total recovered tissue
TRT = Y(:,6*habs+1:habs*7).*CLsimp + Y(:,7*habs+1:habs*8).*CMsimp + Y(:,8*habs+1:habs*9).*CHsimp;

% Total healthy proportion
THP = (Y(:,1:habs)+Y(:,habs+1:habs*2)+Y(:,2*habs+1:habs*3))/3;
% Total infectioys proportion
TIP = (Y(:,3*habs+1:habs*4)+Y(:,4*habs+1:habs*5)+Y(:,5*habs+1:habs*6))/3;
% Total dead proportion
TDP = ((1-(Y(:,1:habs)+Y(:,3*habs+1:habs*4)+Y(:,6*habs+1:habs*7)))+(1-(Y(:,habs+1:habs*2)+Y(:,4*habs+1:habs*5)+Y(:,7*habs+1:habs*8)))+(1-(Y(:,2*habs+1:habs*3)+Y(:,5*habs+1:habs*6)+Y(:,8*habs+1:habs*9))))/3;
% Total recovered proportion
TRP = (Y(:,6*habs+1:habs*7) + Y(:,7*habs+1:habs*8) + Y(:,8*habs+1:habs*9))/3;

habs+1:habs*7) + Y(:,7*habs+1:habs*8) + Y(:,8*habs+1:habs*9))/3;

interpTHT = interp1(t,THT,1:tspan(2));
interpTIT = interp1(t,TIT,1:tspan(2));
interpTDT = interp1(t,TDT,1:tspan(2));
interpTRT = interp1(t,TRT,1:tspan(2));

interpTHT(interpTHT<0) = 0;
interpTIT(interpTIT<0) = 0;
interpTDT(interpTDT<0) = 0;
interpTRT(interpTRT<0) = 0;

interpTHP = interp1(t,THP,1:tspan(2));
interpTIP = interp1(t,TIP,1:tspan(2));
interpTDP = interp1(t,TDP,1:tspan(2));
interpTRP = interp1(t,TRP,1:tspan(2));

interpTHP(interpTHP<0) = 0;
interpTIP(interpTIP<0) = 0;
interpTDP(interpTDP<0) = 0;
interpTRP(interpTRP<0) = 0;


%% Make results into a movie - 2025
% Just pumping something out for Ben - not sure on T vs P...

% Coast = shaperead('JustUS_GEO_02.shp');
% % [xt, yt] = find(interpTIT>thresh, 1, 'last'); % Find the last step where any disease is more than a threshold, like 0.01
% % v = VideoWriter(strcat('DisVid',num2str(now),'.avi'));

v = VideoWriter('DisVid_Ex');
v.FrameRate = 24;
open(v);

% f = figure('renderer', 'zbuffer');

f = figure('renderer', 'zbuffer','Position', [10 10 700 1000]);

set(f,'nextplot','replacechildren'); 

cmap3 = colormap(flipud(autumn(round(max(max(TIP))*10000)+1)));
cmap3(1,:) = [0 .3 1];
cmap4 = colormap(cool(round(max(CC)*100)+1));
% Make a colormap that is blue at 0 and then immediately, say, pink or gold
% spectrum for infection and prevalence

I_tot = (interpHSIP+interpHSIP+interpLSIP)/3; % This is not right yet
H_tot = (interpHSHP+interpHSHP+interpLSHP)/3; % Neither is this

axv1 = subplot(2,1,1);
    scatter(XY(:,1),XY(:,2),1,'filled')
    c1 = colorbar;
    colormap(axv1,cmap3)
    caxis([0 1])
    axis equal

    
axv2 = subplot(2,1,2);
    scatter(XY(:,1),XY(:,2),1,'filled')
    c2 = colorbar;
    colormap(axv2,cmap4)
    caxis([0 1]);
    axis equal


for k = 1:1:tspan(end)
    Dt=datetime('1-Jan-2019')+k-1;


    h1 = scatter(axv1,XY(:,1),XY(:,2),10,interpTIP(k,:),'filled'); hold on
    axis equal

    % h12 = scatter(axv1,XY(TIP(k,:)>0,1),XY(TIP(k,:)>0,2),10,TIP(k,TIP(k,:)>0),'filled'); hold on

    title(axv1,'Disease prevalence',strcat('t ='," ",datestr(Dt)))
      
    h2 = scatter(axv2, XY(:,1),XY(:,2),10,CC-interpTHP(k,:)','filled'); hold on % This is not right
    axis equal


    title(axv2,'Coral cover lost',strcat('t ='," ",datestr(Dt)))
    
    F(k).frame = getframe(f);
    
    writeVideo(v, F(k).frame);

    delete(h1)
    delete(h2)
    delete(h12)

end

close(v); 
close(f);

clear F