function f = ODEfun_SCTLD_seascape(t,Y,CLS1,CMS1,CHS1,bls,bms,bhs,kls,kms,khs,thresh,c,shapeParam,habs,P,Pdays,I0_frac,tau)
% ODEfun_SCTLD_seascape - Multi-host SIR disease dynamics with spatial connectivity
%
% UPDATED: Site-wide activation when ANY class exceeds I0_frac threshold
% - External transmission: always active (no shutoff)
% - Local transmission: activated by smooth sigmoid when ANY class at site exceeds threshold
% - I0_frac: threshold fraction (e.g., 0.01 = 1%) of class population
% - tau: transition width for sigmoid (smaller = sharper activation)
% - Key change: if LS, MS, or HS exceeds their threshold, ALL classes activate locally

persistent last_k Pint_cache
if isempty(last_k)
    last_k = -inf;
    Pint_cache = [];
end

%% Connectivity matrix
tday = floor(t);
k = max(1,find(Pdays <= tday,1,'last'));
if isempty(k); k =1; end

if k ~= last_k
    Pint_cache = P(k).full;
    last_k = k;
end

%% "Initial" conditions for current time step
SLS = Y(1:habs);
SMS = Y(habs+1:habs*2);
SHS = Y(habs*2+1:habs*3);
SLS = SLS(:);
SMS = SMS(:);
SHS = SHS(:);

LS_prez = CLS1>0;
MS_prez = CMS1>0;
HS_prez = CHS1>0;
LS_prez = LS_prez(:);
MS_prez = MS_prez(:);
HS_prez = HS_prez(:);

ILS = Y(habs*3+1:habs*4);
IMS = Y(habs*4+1:habs*5);
IHS = Y(habs*5+1:habs*6);
ILS = ILS(:);
IMS = IMS(:);
IHS = IHS(:);

RLS = Y(habs*6+1:habs*7);
RMS = Y(habs*7+1:habs*8);
RHS = Y(habs*8+1:habs*9);
RLS = RLS(:);
RMS = RMS(:);
RHS = RHS(:);

%% Upstream influence - Calculate incoming disease flux
DP = ILS + IMS + IHS;
DP = DP(:);

% Calculate incoming risk from connectivity
T = incomingRisk_sparse(Pint_cache, DP, thresh, 'rowsAreSources');

% Apply transformation (if shapeParam is near 0.001, this is essentially linear)
T = c.*(1-exp(-shapeParam.*T(:)))/(1-exp(-shapeParam));
T = T(:);

D_prez = DP > 0 | T > 0;
D_prez = D_prez(:);

LS_use = LS_prez & D_prez;
MS_use = MS_prez & D_prez;
HS_use = HS_prez & D_prez;

%% Calculate per-class thresholds for each site
% Each class has its own threshold based on I0_frac
I0_LS = I0_frac * CLS1;  % Threshold for LS class
I0_MS = I0_frac * CMS1;  % Threshold for MS class
I0_HS = I0_frac * CHS1;  % Threshold for HS class

% For each site, determine if ANY class has exceeded its threshold
% This will be used to activate ALL classes at that site
sites_to_check = LS_use | MS_use | HS_use;

% Initialize site activation signal
site_total_I = zeros(habs, 1);
site_min_threshold = zeros(habs, 1);

for idx = 1:habs
    if sites_to_check(idx)
        % Total infected at this site across all classes
        site_total_I(idx) = ILS(idx) + IMS(idx) + IHS(idx);
        
        % Minimum threshold (easiest to cross) among present classes
        thresholds_present = [];
        if LS_prez(idx)
            thresholds_present = [thresholds_present; I0_LS(idx)];
        end
        if MS_prez(idx)
            thresholds_present = [thresholds_present; I0_MS(idx)];
        end
        if HS_prez(idx)
            thresholds_present = [thresholds_present; I0_HS(idx)];
        end
        
        if ~isempty(thresholds_present)
            site_min_threshold(idx) = min(thresholds_present);
        end
    end
end

%% The SIR equations with site-wide threshold activation
% External transmission: always active at base rate beta
% Local transmission: smooth threshold-activated when ANY class exceeds threshold
% All classes at a site share the same activation state

% ===== LS class =====
if any(LS_use)
    DPLS = DP(LS_use);
    TLS = T(LS_use);
    fI_LS = DPLS./CLS1(LS_use);
    
    % External transmission: always active
    p_ext_LS = 1 - exp(-bls * TLS);
    
    % Local transmission: activated by site-wide threshold
    % Use total infected vs minimum threshold at each site
    tau_site = site_min_threshold(LS_use) / 10;  % Adaptive tau based on threshold
    beta_prime_LS = (bls/2) * (1 + tanh((site_total_I(LS_use) - site_min_threshold(LS_use)) ./ tau_site));
    p_loc_LS = 1 - exp(-beta_prime_LS .* fI_LS);
    
    % Combine (no shutoff - both always contribute)
    p_tot_LS = p_loc_LS + p_ext_LS - p_loc_LS.*p_ext_LS;
else
    p_tot_LS = [];
end

% ===== MS class =====
if any(MS_use)
    DPMS = DP(MS_use);
    TMS = T(MS_use);
    fI_MS = DPMS./CMS1(MS_use);
    
    p_ext_MS = 1 - exp(-bms * TMS);
    
    % Use same site-wide activation as LS
    tau_site = site_min_threshold(MS_use) / 10;
    beta_prime_MS = (bms/2) * (1 + tanh((site_total_I(MS_use) - site_min_threshold(MS_use)) ./ tau_site));
    p_loc_MS = 1 - exp(-beta_prime_MS .* fI_MS);
    
    p_tot_MS = p_loc_MS + p_ext_MS - p_loc_MS.*p_ext_MS;
else
    p_tot_MS = [];
end

% ===== HS class =====
if any(HS_use)
    DPHS = DP(HS_use);
    THS = T(HS_use);
    fI_HS = DPHS./CHS1(HS_use);
    
    p_ext_HS = 1 - exp(-bhs * THS);
    
    % Use same site-wide activation as LS and MS
    tau_site = site_min_threshold(HS_use) / 10;
    beta_prime_HS = (bhs/2) * (1 + tanh((site_total_I(HS_use) - site_min_threshold(HS_use)) ./ tau_site));
    p_loc_HS = 1 - exp(-beta_prime_HS .* fI_HS);
    
    p_tot_HS = p_loc_HS + p_ext_HS - p_loc_HS.*p_ext_HS;
else
    p_tot_HS = [];
end

%% Calculate derivatives
dSLS = zeros(length(SLS),1);
if any(LS_use)
    dSLS(LS_use) = -p_tot_LS.*SLS(LS_use);
end

dSMS = zeros(length(SMS),1);
if any(MS_use)
    dSMS(MS_use) = -p_tot_MS.*SMS(MS_use);
end

dSHS = zeros(length(SHS),1);
if any(HS_use)
    dSHS(HS_use) = -p_tot_HS.*SHS(HS_use);
end

dILS = zeros(length(ILS),1);
if any(LS_use)
    dILS(LS_use) = p_tot_LS.*SLS(LS_use) - kls*ILS(LS_use);
end

dIMS = zeros(length(IMS),1);
if any(MS_use)
    dIMS(MS_use) = p_tot_MS.*SMS(MS_use) - kms*IMS(MS_use);
end

dIHS = zeros(length(IHS),1);
if any(HS_use)
    dIHS(HS_use) = p_tot_HS.*SHS(HS_use) - khs*IHS(HS_use);
end

dRLS = kls*ILS;
dRMS = kms*IMS;
dRHS = khs*IHS;

f = [dSLS;dSMS;dSHS;dILS;dIMS;dIHS;dRLS;dRMS;dRHS];
end