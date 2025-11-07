clear; clc

%% Paths
projectPath = matlab.project.rootProject().RootFolder;
tempPath = fullfile(projectPath, 'temp');
dataPath = fullfile(projectPath, 'data');

%% Load reef data
reefData = readtable(fullfile(dataPath, 'centroids_vertices_FINALFORCMS.csv'));
habs = height(reefData);
CLS1 = reefData.low_coral_cover;
CMS1 = reefData.moderate_coral_cover;
CHS1 = reefData.high_coral_cover;

%% Load cached connectivity
load(fullfile(tempPath, 'P_20190101_to_20190331.mat'), 'P');
dates = [P.Date];
Pdays = days(dates - datetime(2019,1,1));

%% Seed initial infection at Flat Cay sites
Flat = [find(reefData.unique_ID==29088) find(reefData.unique_ID==29338) ...
        find(reefData.unique_ID==29089) find(reefData.unique_ID==29339) ...
        find(reefData.unique_ID==29087)];

ILS1 = zeros(habs,1); IMS1 = zeros(habs,1); IHS1 = zeros(habs,1);
ILS1(Flat) = .01*CLS1(Flat);
IMS1(Flat) = .01*CMS1(Flat);
IHS1(Flat) = .01*CHS1(Flat);

SLS1 = CLS1 - ILS1;
SMS1 = CMS1 - IMS1;
SHS1 = CHS1 - IHS1;

%% Disease parameters
bls = 0.03; bms = 0.14; bhs = 2.08;
kls = .05; kms = .55; khs = 3.33;

%% TEST SCENARIOS - Fixed structure
scenarios = {
    0,      0.001, 0,      1e-10, 1,     'Baseline: No filters';
    0,      0.001, 0.01,   0.0001, 1,     'Only I0 threshold';
    0.0003, 0.001, 0,      1e-10, 1,     'Only connectivity threshold';
    0.0003, -4,    0,      1e-10, 1,     'Conn thresh + nonlinear shape';
    0.0003, 0.001, 0.01,   0.0001, 1,     'Both thresholds';
    0,      0.001, 0,      1e-10, 100,   'No thresh, conn x100';
    0,      0.001, 0,      1e-10, 10000, 'Connectivity x10000';
};

tspan = [1 90];
Y0 = [SLS1;SMS1;SHS1;ILS1;IMS1;IHS1;zeros(habs*3,1)];

Results = struct();

fprintf('\n========== CONNECTIVITY TESTING ==========\n\n');

% Preallocate results
n_scenarios = size(scenarios,1);
Results(n_scenarios).desc = [];

parfor s = 1:n_scenarios
    thresh = scenarios{s,1};
    shapeParam = scenarios{s,2};
    I0_frac = scenarios{s,3};
    tau = scenarios{s,4};
    c_mult = scenarios{s,5};
    desc = scenarios{s,6};
    
    fprintf('Scenario %d: %s\n', s, desc);
    
    tic;
    [t,Y] = ode45(@(t,Y) simple_ODE(t,Y,CLS1,CMS1,CHS1,bls,bms,bhs,kls,kms,khs,...
                                     thresh,shapeParam,I0_frac,tau,c_mult,habs,P,Pdays), ...
                  tspan, Y0);
    runtime = toc;
    
    % Extract results
    TI = Y(:,3*habs+1:4*habs) + Y(:,4*habs+1:5*habs) + Y(:,5*habs+1:6*habs);
    TR = Y(:,6*habs+1:7*habs) + Y(:,7*habs+1:8*habs) + Y(:,8*habs+1:9*habs);
    
    n_infected = sum(max(TI,[],1) > 0);
    total_removed = sum(TR(end,:));
    
    Results(s).desc = desc;
    Results(s).params = scenarios(s,1:5);
    Results(s).n_infected = n_infected;
    Results(s).total_removed = total_removed;
    Results(s).runtime = runtime;
    
    fprintf('  → %d sites infected, %.4f total removed (%.2f sec)\n\n', ...
            n_infected, total_removed, runtime);
end

fprintf('========================================\n\n');

%% Summary table
fprintf('SUMMARY:\n');
fprintf('%-40s %10s %15s\n', 'Scenario', 'Infected', 'Total Removed');
fprintf('%s\n', repmat('-',67,1));
for s = 1:length(Results)
    fprintf('%-40s %10d %15.4f\n', Results(s).desc, ...
            Results(s).n_infected, Results(s).total_removed);
end

%% Simple ODE function
function f = simple_ODE(t,Y,CLS1,CMS1,CHS1,bls,bms,bhs,kls,kms,khs,...
                        thresh,shapeParam,I0_frac,tau,c_mult,habs,P,Pdays)
    persistent last_k Pint
    if isempty(last_k), last_k = -inf; Pint = []; end
    
    % Get connectivity matrix
    k = max(1,find(Pdays <= floor(t),1,'last'));
    if isempty(k), k = 1; end
    if k ~= last_k
        Pint = P(k).full;
        last_k = k;
    end
    
    % Extract states
    SLS = Y(1:habs); SMS = Y(habs+1:2*habs); SHS = Y(2*habs+1:3*habs);
    ILS = Y(3*habs+1:4*habs); IMS = Y(4*habs+1:5*habs); IHS = Y(5*habs+1:6*habs);
    
    % Disease pool
    DP = ILS + IMS + IHS;
    
    % Connectivity influence with threshold
    T = incomingRisk_sparse(Pint, DP, thresh, 'rowsAreSources');
    if abs(shapeParam) < 1e-6
        T = c_mult * T;  % Linear case
    else
        T = c_mult * (1-exp(-shapeParam*T)) / (1-exp(-shapeParam));
    end
    
    % Diagnostic output
    if floor(t) == 1 && t < 1.1
        fprintf('t=%.2f: mean(T)=%.8f, max(T)=%.8f, sum(DP>0)=%d\n', ...
                t, mean(T), max(T), sum(DP>0));
    end
    
    % Sites with disease present
    active = (CLS1>0 | CMS1>0 | CHS1>0) & (DP>0 | T>0);
    
    % Initialize derivatives
    dSLS = zeros(habs,1); dSMS = zeros(habs,1); dSHS = zeros(habs,1);
    dILS = zeros(habs,1); dIMS = zeros(habs,1); dIHS = zeros(habs,1);
    
    if any(active)
        % External transmission
        p_ext_LS = 1 - exp(-bls * T(active));
        p_ext_MS = 1 - exp(-bms * T(active));
        p_ext_HS = 1 - exp(-bhs * T(active));
        
        % Local transmission
        if I0_frac > 0
            site_I = DP(active);
            I0_site = I0_frac * (CLS1(active) + CMS1(active) + CHS1(active));
            beta_mult = 0.5 * (1 + tanh((site_I - I0_site) / tau));
        else
            beta_mult = 1;
        end
        
        fI_LS = DP(active) ./ max(CLS1(active), 1e-10);
        fI_MS = DP(active) ./ max(CMS1(active), 1e-10);
        fI_HS = DP(active) ./ max(CHS1(active), 1e-10);
        
        p_loc_LS = 1 - exp(-beta_mult .* bls .* fI_LS);
        p_loc_MS = 1 - exp(-beta_mult .* bms .* fI_MS);
        p_loc_HS = 1 - exp(-beta_mult .* bhs .* fI_HS);
        
        % Combined transmission
        p_tot_LS = p_loc_LS + p_ext_LS - p_loc_LS.*p_ext_LS;
        p_tot_MS = p_loc_MS + p_ext_MS - p_loc_MS.*p_ext_MS;
        p_tot_HS = p_loc_HS + p_ext_HS - p_loc_HS.*p_ext_HS;
        
        % Apply to active sites
        dSLS(active) = -p_tot_LS .* SLS(active);
        dSMS(active) = -p_tot_MS .* SMS(active);
        dSHS(active) = -p_tot_HS .* SHS(active);
        
        dILS(active) = p_tot_LS .* SLS(active) - kls*ILS(active);
        dIMS(active) = p_tot_MS .* SMS(active) - kms*IMS(active);
        dIHS(active) = p_tot_HS .* SHS(active) - khs*IHS(active);
    end
    
    dRLS = kls*ILS;
    dRMS = kms*IMS;
    dRHS = khs*IHS;
    
    f = [dSLS;dSMS;dSHS;dILS;dIMS;dIHS;dRLS;dRMS;dRHS];
end