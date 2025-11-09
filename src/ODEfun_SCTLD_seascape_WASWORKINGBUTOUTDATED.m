


% %% without '1 -' for transmission
% 
% 
% function f = ODEfun_SCTLD_seascape(t,Y,CLS1,CMS1,CHS1,bls,bms,bhs,kls,kms,khs,thresh,c,shapeParam,habs,P,Pdays)
% % ODEfun_SCTLD_seascape - Multi-host SIR disease dynamics with spatial connectivity
% % 
% % MODIFIED VERSION: Uses rate-based transmission (matching R code) instead of 
% % probability-based transmission (1 - exp(-rate))
% %
% % Each location has 3 variables per susceptibility class: S, I, R
% % Disease prevalence is translated back and forth to disease cover (%).
% 
% persistent last_k Pint_cache
% if isempty(last_k); last_k = -inf; Pint_cache = []; end
% 
% %% Connectivity matrix
% % ODE solvers integrate across non-integer time steps. Thus, will need to
% % interpolate connectivity matrices at time t.
% tday = floor(t);
% k = max(1,find(Pdays <= tday,1,'last'));
% if isempty(k); k =1; end
% 
% if k ~= last_k
%     Pint_cache = P(k).full;
%     last_k = k;
% end
% 
% %% "Initial" conditions for current time step
% % Susceptible
% SLS = Y(1:habs);
% SMS = Y(habs+1:habs*2);
% SHS = Y(habs*2+1:habs*3);
% SLS = SLS(:);
% SMS = SMS(:);
% SHS = SHS(:);
% 
% LS_prez = CLS1>0;
% MS_prez = CMS1>0;
% HS_prez = CHS1>0;
% LS_prez = LS_prez(:);
% MS_prez = MS_prez(:);
% HS_prez = HS_prez(:);
% 
% % Infected
% ILS = Y(habs*3+1:habs*4);
% IMS = Y(habs*4+1:habs*5);
% IHS = Y(habs*5+1:habs*6);
% ILS = ILS(:);
% IMS = IMS(:);
% IHS = IHS(:);
% 
% % Recovering
% RLS = Y(habs*6+1:habs*7);
% RMS = Y(habs*7+1:habs*8);
% RHS = Y(habs*8+1:habs*9);
% RLS = RLS(:);
% RMS = RMS(:);
% RHS = RHS(:);
% 
% %% Upstream influence
% % Infection probability is related to susceptibility, upstream
% % connectivity, and local disease prevalence
% % Calculate transition of disease. Site must have disease over a threshold
% % (defined outside ODE) in order to transmit disease.
% DP = ILS + IMS + IHS;
% DP = DP(:);
% 
% % Estimate disease FLUX, depends on upstream DP, connectivity probability,
% % and thresh (minimum DP(i) to contribute to T(j))
% T = incomingRisk_sparse(Pint_cache, DP, thresh,'rowsAreSources');
% 
% % Reshape T - this modulates how accumulated FLUX translates to local
% % infections. Very sensitive to this.
% T = c.*(1-exp(-shapeParam.*T(:)))/(1-exp(-shapeParam));
% T = T(:);
% 
% D_prez = DP > 0 | T > 0;
% D_prez = D_prez(:);
% 
% LS_use = LS_prez & D_prez;
% MS_use = MS_prez & D_prez;
% HS_use = HS_prez & D_prez;
% 
% %% The SIR equations - RATE-BASED VERSION (matches R code)
% % Changed from probability-based (1 - exp(-rate)) to direct rate formulation
% % Implements: dI/dt = beta * S * (I_local/N + I_external) - gamma * I
% 
% % Initialize all derivatives to zero
% dSLS = zeros(length(SLS),1);
% dSMS = zeros(length(SMS),1);
% dSHS = zeros(length(SHS),1);
% 
% dILS = zeros(length(ILS),1);
% dIMS = zeros(length(IMS),1);
% dIHS = zeros(length(IHS),1);
% 
% % LS group - rate-based transmission
% if any(LS_use)
%     DPLS = DP(LS_use);
%     TLS = T(LS_use);
% 
%     % Local transmission: frequency-dependent (matches R code)
%     % fI = infected fraction within group
%     fI_LS = DPLS ./ CLS1(LS_use);
% 
%     % Rate-based transmission (not probability)
%     % Local: beta * S * (I_total / N)
%     % External: beta * S * T (T is already a rate)
%     rate_loc_LS = bls * fI_LS;
%     rate_ext_LS = bls * TLS;
%     rate_tot_LS = rate_loc_LS + rate_ext_LS;
% 
%     dSLS(LS_use) = -rate_tot_LS .* SLS(LS_use);
%     dILS(LS_use) = rate_tot_LS .* SLS(LS_use) - kls * ILS(LS_use);
% end
% 
% % MS group - rate-based transmission
% if any(MS_use)
%     DPMS = DP(MS_use);
%     TMS = T(MS_use);
% 
%     fI_MS = DPMS ./ CMS1(MS_use);
% 
%     rate_loc_MS = bms * fI_MS;
%     rate_ext_MS = bms * TMS;
%     rate_tot_MS = rate_loc_MS + rate_ext_MS;
% 
%     dSMS(MS_use) = -rate_tot_MS .* SMS(MS_use);
%     dIMS(MS_use) = rate_tot_MS .* SMS(MS_use) - kms * IMS(MS_use);
% end
% 
% % HS group - rate-based transmission
% if any(HS_use)
%     DPHS = DP(HS_use);
%     THS = T(HS_use);
% 
%     fI_HS = DPHS ./ CHS1(HS_use);
% 
%     rate_loc_HS = bhs * fI_HS;
%     rate_ext_HS = bhs * THS;
%     rate_tot_HS = rate_loc_HS + rate_ext_HS;
% 
%     dSHS(HS_use) = -rate_tot_HS .* SHS(HS_use);
%     dIHS(HS_use) = rate_tot_HS .* SHS(HS_use) - khs * IHS(HS_use);
% end
% 
% % Recovery/mortality always applies (even where transmission isn't happening)
% dRLS = kls * ILS;
% dRMS = kms * IMS;
% dRHS = khs * IHS;
% 
% f = [dSLS;dSMS;dSHS;dILS;dIMS;dIHS;dRLS;dRMS;dRHS];
% 
% end






% %% with '1 -' for transmission
% 
% function f = ODEfun_SCTLD_seascape(t,Y,CLS1,CMS1,CHS1,bls,bms,bhs,kls,kms,khs,thresh,c,shapeParam,habs,P,Pdays)
% % ODEfun_SCTLD_seascape - Multi-host SIR disease dynamics with spatial connectivity
% %
% % Each location has 3 variables per susceptibility class: S, I, R
% % Disease prevalence is translated back and forth to disease cover (%).
% %
% % Updated: Renamed from ODEfun_08, changed B parameter to shapeParam
% 
% persistent last_k Pint_cache
% if isempty(last_k); last_k = -inf; Pint_cache = []; end
% 
% %% Connectivity matrix
% % ODE solvers integrate across non-integer time steps. Thus, will need to
% % interpolate connectivity matrices at time t.
% 
% tday = floor(t);
% 
% k = max(1,find(Pdays <= tday,1,'last'));
% if isempty(k); k =1; end
% 
% if k ~= last_k
%     Pint_cache = P(k).full;
%     last_k = k;
% end
% 
% %% "Initial" conditions for current time step
% % Susceptible
% SLS = Y(1:habs);
% SMS = Y(habs+1:habs*2);
% SHS = Y(habs*2+1:habs*3);
% 
% SLS = SLS(:);
% SMS = SMS(:);
% SHS = SHS(:);
% 
% LS_prez = CLS1>0;
% MS_prez = CMS1>0;
% HS_prez = CHS1>0;
% 
% LS_prez = LS_prez(:);
% MS_prez = MS_prez(:);
% HS_prez = HS_prez(:);
% 
% % Infected
% ILS = Y(habs*3+1:habs*4);
% IMS = Y(habs*4+1:habs*5);
% IHS = Y(habs*5+1:habs*6);
% 
% ILS = ILS(:);
% IMS = IMS(:);
% IHS = IHS(:);
% 
% % Recovering
% RLS = Y(habs*6+1:habs*7);
% RMS = Y(habs*7+1:habs*8);
% RHS = Y(habs*8+1:habs*9);
% 
% RLS = RLS(:);
% RMS = RMS(:);
% RHS = RHS(:);
% 
% %% Upstream influence
% % Infection probability is related to susceptibility, upstream
% % connectivity, and local disease prevalence
% % Calculate transition of disease. Site must have disease over a threshold
% % (defined outside ODE) in order to transmit disease.
% DP = ILS + IMS + IHS; 
% DP = DP(:);
% 
% % Estimate disease FLUX, depends on upstream DP, connectivity probability,
% % and thresh (minimum DP(i) to contribute to T(j))
% T = incomingRisk_sparse(Pint_cache, DP, thresh,'rowsAreSources');
% 
% % Reshape T - this modulates how accumulated FLUX translates to local
% % infections. Very sensitive to this.
% T = c.*(1-exp(-shapeParam.*T(:)))/(1-exp(-shapeParam));
% T = T(:);
% 
% D_prez = DP > 0 | T > 0;
% D_prez = D_prez(:);
% 
% LS_use = LS_prez & D_prez;
% MS_use = MS_prez & D_prez;
% HS_use = HS_prez & D_prez;
% 
% %% The SIR equations
% % Adapted from Ben's model
% % Basic SIR with scalable mortality (k)/recovery
% 
% % Currently all categories will be initiated with external infection. To
% % turn off, remove the T_S term from p_tot calcs...
% 
% DPLS = DP(LS_use);
% TLS = T(LS_use);
% 
% fI_LS = DPLS./CLS1(LS_use);
% p_loc_LS = 1 - exp(-bls * fI_LS);
% p_ext_LS = 1 - exp(-bls * TLS);
% 
% p_tot_LS = p_loc_LS+p_ext_LS - p_loc_LS.*p_ext_LS; % This will never go over 1
% 
% DPMS = DP(MS_use);
% TMS = T(MS_use);
% 
% fI_MS = DPMS./CMS1(MS_use);
% p_loc_MS = 1 - exp(-bms * fI_MS);
% p_ext_MS = 1 - exp(-bms * TMS);
% 
% p_tot_MS = p_loc_MS+p_ext_MS - p_loc_MS.*p_ext_MS; % This will never go over 1
% 
% DPHS = DP(HS_use);
% THS = T(HS_use);
% 
% fI_HS = DPHS./CHS1(HS_use);
% p_loc_HS = 1 - exp(-bhs * fI_HS);
% p_ext_HS = 1 - exp(-bhs * THS);
% 
% p_tot_HS = p_loc_HS+p_ext_HS - p_loc_HS.*p_ext_HS; % This will never go over 1
% 
% dSLS = zeros(length(SLS),1);
%     dSLS(LS_use) = -p_tot_LS.*SLS(LS_use);
% dSMS = zeros(length(SMS),1);
%     dSMS(MS_use) = -p_tot_MS.*SMS(MS_use);
% dSHS = zeros(length(SHS),1);
%     dSHS(HS_use) = -p_tot_HS.*SHS(HS_use);
% 
% dILS = zeros(length(ILS),1);
%     dILS(LS_use) = p_tot_LS.*SLS(LS_use) - kls*ILS(LS_use);
% dIMS = zeros(length(IMS),1);
%     dIMS(MS_use) = p_tot_MS.*SMS(MS_use) - kms*IMS(MS_use);
% dIHS = zeros(length(IHS),1);
%     dIHS(HS_use) = p_tot_HS.*SHS(HS_use) - khs*IHS(HS_use);
% 
% dRLS = kls*ILS;
% dRMS = kms*IMS;
% dRHS = khs*IHS;
% 
% f = [dSLS;dSMS;dSHS;dILS;dIMS;dIHS;dRLS;dRMS;dRHS];
% end


%% with shutoff valve for incoming transmission

function f = ODEfun_SCTLD_seascape(t,Y,CLS1,CMS1,CHS1,bls,bms,bhs,kls,kms,khs,thresh,c,shapeParam,habs,P,Pdays)
% ODEfun_SCTLD_seascape - Multi-host SIR disease dynamics with spatial connectivity
%
% Each location has 3 variables per susceptibility class: S, I, R
% Disease prevalence is translated back and forth to disease cover (%).
%
% Updated: Added external transmission shutoff once local infection exceeds threshold
% External transmission helps initiate outbreaks but turns off once sites are self-sustaining

persistent last_k Pint_cache
if isempty(last_k); last_k = -inf; Pint_cache = []; end

%% Connectivity matrix
% ODE solvers integrate across non-integer time steps. Thus, will need to
% interpolate connectivity matrices at time t.
tday = floor(t);
k = max(1,find(Pdays <= tday,1,'last'));
if isempty(k); k =1; end
if k ~= last_k
    Pint_cache = P(k).full;
    last_k = k;
end

%% "Initial" conditions for current time step
% Susceptible
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

% Infected
ILS = Y(habs*3+1:habs*4);
IMS = Y(habs*4+1:habs*5);
IHS = Y(habs*5+1:habs*6);
ILS = ILS(:);
IMS = IMS(:);
IHS = IHS(:);

% Recovering
RLS = Y(habs*6+1:habs*7);
RMS = Y(habs*7+1:habs*8);
RHS = Y(habs*8+1:habs*9);
RLS = RLS(:);
RMS = RMS(:);
RHS = RHS(:);

%% Upstream influence
% Infection probability is related to susceptibility, upstream
% connectivity, and local disease prevalence
% Calculate transition of disease. Site must have disease over a threshold
% (defined outside ODE) in order to transmit disease.
DP = ILS + IMS + IHS;
DP = DP(:);

% Estimate disease FLUX, depends on upstream DP, connectivity probability,
% and thresh (minimum DP(i) to contribute to T(j))
T = incomingRisk_sparse(Pint_cache, DP, thresh,'rowsAreSources');

% Reshape T - this modulates how accumulated FLUX translates to local
% infections. Very sensitive to this.
T = c.*(1-exp(-shapeParam.*T(:)))/(1-exp(-shapeParam));
T = T(:);

D_prez = DP > 0 | T > 0;
D_prez = D_prez(:);
LS_use = LS_prez & D_prez;
MS_use = MS_prez & D_prez;
HS_use = HS_prez & D_prez;

%% The SIR equations with external transmission shutoff
% External transmission turns off once local infection exceeds threshold
% This allows connectivity to initiate outbreaks but not dominate intensity

% Define "kickoff" threshold: 2x typical initial seeding (0.01% of cover)
% Sites above this threshold are considered "self-sustaining"
kickoff_threshold = 2 * 0.0001;  % 0.02% absolute cover

% For LS class
DPLS = DP(LS_use);
TLS = T(LS_use);
fI_LS = DPLS./CLS1(LS_use);
p_loc_LS = 1 - exp(-bls * fI_LS);
p_ext_LS = 1 - exp(-bls * TLS);

% Turn off external if local infection exceeds threshold
local_established_LS = ILS(LS_use) > kickoff_threshold;
p_ext_LS_active = (~local_established_LS) .* p_ext_LS;
p_tot_LS = p_loc_LS + p_ext_LS_active - p_loc_LS.*p_ext_LS_active;

% For MS class
DPMS = DP(MS_use);
TMS = T(MS_use);
fI_MS = DPMS./CMS1(MS_use);
p_loc_MS = 1 - exp(-bms * fI_MS);
p_ext_MS = 1 - exp(-bms * TMS);

local_established_MS = IMS(MS_use) > kickoff_threshold;
p_ext_MS_active = (~local_established_MS) .* p_ext_MS;
p_tot_MS = p_loc_MS + p_ext_MS_active - p_loc_MS.*p_ext_MS_active;

% For HS class
DPHS = DP(HS_use);
THS = T(HS_use);
fI_HS = DPHS./CHS1(HS_use);
p_loc_HS = 1 - exp(-bhs * fI_HS);
p_ext_HS = 1 - exp(-bhs * THS);

local_established_HS = IHS(HS_use) > kickoff_threshold;
p_ext_HS_active = (~local_established_HS) .* p_ext_HS;
p_tot_HS = p_loc_HS + p_ext_HS_active - p_loc_HS.*p_ext_HS_active;

% Calculate derivatives
dSLS = zeros(length(SLS),1);
    dSLS(LS_use) = -p_tot_LS.*SLS(LS_use);
dSMS = zeros(length(SMS),1);
    dSMS(MS_use) = -p_tot_MS.*SMS(MS_use);
dSHS = zeros(length(SHS),1);
    dSHS(HS_use) = -p_tot_HS.*SHS(HS_use);

dILS = zeros(length(ILS),1);
    dILS(LS_use) = p_tot_LS.*SLS(LS_use) - kls*ILS(LS_use);
dIMS = zeros(length(IMS),1);
    dIMS(MS_use) = p_tot_MS.*SMS(MS_use) - kms*IMS(MS_use);
dIHS = zeros(length(IHS),1);
    dIHS(HS_use) = p_tot_HS.*SHS(HS_use) - khs*IHS(HS_use);

dRLS = kls*ILS;
dRMS = kms*IMS;
dRHS = khs*IHS;

f = [dSLS;dSMS;dSHS;dILS;dIMS;dIHS;dRLS;dRMS;dRHS];
end