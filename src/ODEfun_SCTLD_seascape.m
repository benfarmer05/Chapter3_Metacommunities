
function f = ODEfun_SCTLD_seascape(t,Y,N_LS,N_MS,N_HS,b_LS,b_MS,b_HS,g_LS,g_MS,g_HS,export_thresh,flux_scale,flux_shape,num_sites,conn_structs,conn_days,I0,tau)

%% Select correct connectivity matrix at ODE time t
persistent last_conn_index conn_cache
if isempty(last_conn_index); last_conn_index = -inf; conn_cache = []; end

tday = floor(t);
conn_index = max(1, find(conn_days <= tday, 1, 'last'));
if isempty(conn_index); conn_index = 1; end

if conn_index ~= last_conn_index
    conn_cache = conn_structs(conn_index).full;
    last_conn_index = conn_index;
end

%% Extract compartment states at ODE time t
% Susceptible
S_LS = Y(1:num_sites);
S_MS = Y(num_sites+1:num_sites*2);
S_HS = Y(num_sites*2+1:num_sites*3);
S_LS = S_LS(:);
S_MS = S_MS(:);
S_HS = S_HS(:);

% Infected
I_LS = Y(num_sites*3+1:num_sites*4);
I_MS = Y(num_sites*4+1:num_sites*5);
I_HS = Y(num_sites*5+1:num_sites*6);
I_LS = I_LS(:);
I_MS = I_MS(:);
I_HS = I_HS(:);

% Recovered
R_LS = Y(num_sites*6+1:num_sites*7);
R_MS = Y(num_sites*7+1:num_sites*8);
R_HS = Y(num_sites*8+1:num_sites*9);
R_LS = R_LS(:);
R_MS = R_MS(:);
R_HS = R_HS(:);

% Define presence of each compartment
LS_present = N_LS > 0;
MS_present = N_MS > 0;
HS_present = N_HS > 0;
LS_present = LS_present(:);
MS_present = MS_present(:);
HS_present = HS_present(:);

%% Define upstream disease influence
% Local disease pools
P = I_LS + I_MS + I_HS;
P = P(:);

% Disease flux
incoming_flux = incomingRisk_sparse(conn_cache, P, export_thresh, 'rowsAreSources');
incoming_flux = flux_scale .* (1 - exp(-flux_shape .* incoming_flux(:))) / (1 - exp(-flux_shape));
incoming_flux = incoming_flux(:);

disease_present = P > 0 | incoming_flux > 0;
disease_present = disease_present(:);

LS_active = LS_present & disease_present;
MS_active = MS_present & disease_present;
HS_active = HS_present & disease_present;

%% Extract active site values
N_LS_active = N_LS(LS_active);
N_MS_active = N_MS(MS_active);
N_HS_active = N_HS(HS_active);

S_LS_active = S_LS(LS_active);
S_MS_active = S_MS(MS_active);
S_HS_active = S_HS(HS_active);

I_LS_active = I_LS(LS_active);
I_MS_active = I_MS(MS_active);
I_HS_active = I_HS(HS_active);

flux_LS = incoming_flux(LS_active);
flux_MS = incoming_flux(MS_active);
flux_HS = incoming_flux(HS_active);

%% ========== DOBBELAERE-STYLE INTERNAL TRANSMISSION THRESHOLD ==========
% Calculate smooth internal transmission multiplier using tanh sigmoid
% Based on Dobbelaere et al. 2020, Equation 3, page 5:
%   β'(Ij) = (β'₀/2) * (1 + tanh[(Ij - I0)/τ])
%
% I0  = infection threshold (fraction of site that must be infected before
%       exponential within-site spread activates)
% tau = transition steepness (smaller = sharper threshold)
%
% This creates a smooth 0-to-1 multiplier:
%   - Below I0: internal transmission is suppressed (coral resistance)
%   - At I0: half-strength internal transmission
%   - Above I0: full-strength internal transmission (exponential spread)

% Compute total infection fraction at each site (all groups combined)
N_total = N_LS + N_MS + N_HS;
N_total(N_total == 0) = 1; % Avoid division by zero
frac_infected_total = P ./ N_total;

% Smooth sigmoid transition (Dobbelaere Eq. 3)
% tanh ranges from -1 to +1, so 0.5*(1+tanh) ranges from 0 to 1
internal_multiplier = 0.5 * (1 + tanh((frac_infected_total - I0) / tau));
internal_multiplier = internal_multiplier(:);

% Extract multipliers for active sites only
internal_mult_LS = internal_multiplier(LS_active);
internal_mult_MS = internal_multiplier(MS_active);
internal_mult_HS = internal_multiplier(HS_active);

%% ========== SWITCHABLE TRANSMISSION FORMULATION ==========
% Comment/uncomment ONE of the two blocks below

% ===== VERSION 1: RATE-BASED (matches R code) =====
% External transmission (connectivity-driven, always active)
rate_ext_LS = b_LS * flux_LS;
rate_ext_MS = b_MS * flux_MS;
rate_ext_HS = b_HS * flux_HS;

% % ----- ORIGINAL: Internal transmission (always active) -----
% % Local transmission uses TOTAL disease pool P (all groups), not group-specific I
% rate_loc_LS = b_LS * (P(LS_active) ./ N_LS_active);
% rate_loc_MS = b_MS * (P(MS_active) ./ N_MS_active);
% rate_loc_HS = b_HS * (P(HS_active) ./ N_HS_active);

% ----- DOBBELAERE: Internal transmission (threshold-modulated) -----
% Local transmission now multiplied by smooth internal_multiplier
% This suppresses within-site spread until infection exceeds threshold I0
rate_loc_LS = b_LS * (P(LS_active) ./ N_LS_active) .* internal_mult_LS;
rate_loc_MS = b_MS * (P(MS_active) ./ N_MS_active) .* internal_mult_MS;
rate_loc_HS = b_HS * (P(HS_active) ./ N_HS_active) .* internal_mult_HS;

% Total transmission rate (external + internal)
rate_tot_LS = rate_loc_LS + rate_ext_LS;
rate_tot_MS = rate_loc_MS + rate_ext_MS;
rate_tot_HS = rate_loc_HS + rate_ext_HS;

% % ===== VERSION 2: PROBABILITY-BASED (original formulation) =====
% % External transmission probabilities
% p_ext_LS = 1 - exp(-b_LS * flux_LS);
% p_ext_MS = 1 - exp(-b_MS * flux_MS);
% p_ext_HS = 1 - exp(-b_HS * flux_HS);
% 
% % % ----- ORIGINAL: Internal transmission (always active) -----
% % % Local transmission uses TOTAL disease pool P (all groups), not group-specific I
% % fI_LS = (P(LS_active) ./ N_LS_active);
% % fI_MS = (P(MS_active) ./ N_MS_active);
% % fI_HS = (P(HS_active) ./ N_HS_active);
% % 
% % p_loc_LS = 1 - exp(-b_LS * fI_LS);
% % p_loc_MS = 1 - exp(-b_MS * fI_MS);
% % p_loc_HS = 1 - exp(-b_HS * fI_HS);
% 
% % ----- DOBBELAERE: Internal transmission (threshold-modulated) -----
% % Local transmission probabilities now modulated by internal_multiplier
% fI_LS = (P(LS_active) ./ N_LS_active) .* internal_mult_LS;
% fI_MS = (P(MS_active) ./ N_MS_active) .* internal_mult_MS;
% fI_HS = (P(HS_active) ./ N_HS_active) .* internal_mult_HS;
% 
% p_loc_LS = 1 - exp(-b_LS * fI_LS);
% p_loc_MS = 1 - exp(-b_MS * fI_MS);
% p_loc_HS = 1 - exp(-b_HS * fI_HS);
% 
% % Combined probability (accounting for independence)
% rate_tot_LS = p_loc_LS + p_ext_LS - p_loc_LS .* p_ext_LS;
% rate_tot_MS = p_loc_MS + p_ext_MS - p_loc_MS .* p_ext_MS;
% rate_tot_HS = p_loc_HS + p_ext_HS - p_loc_HS .* p_ext_HS;

%% ========================================================

%% Calculate derivatives
dS_LS = zeros(length(S_LS), 1);
    dS_LS(LS_active) = -rate_tot_LS .* S_LS_active;
dS_MS = zeros(length(S_MS), 1);
    dS_MS(MS_active) = -rate_tot_MS .* S_MS_active;
dS_HS = zeros(length(S_HS), 1);
    dS_HS(HS_active) = -rate_tot_HS .* S_HS_active;

% not explicitly taking out dying coral at each site - could be an issue
dI_LS = zeros(length(I_LS), 1);
    dI_LS(LS_active) = rate_tot_LS .* S_LS_active - g_LS * I_LS_active;
dI_MS = zeros(length(I_MS), 1);
    dI_MS(MS_active) = rate_tot_MS .* S_MS_active - g_MS * I_MS_active;
dI_HS = zeros(length(I_HS), 1);
    dI_HS(HS_active) = rate_tot_HS .* S_HS_active - g_HS * I_HS_active;

% % failsafe to remove at every location
% recovery_LS = g_LS * I_LS; %calculate recovery for ALL sites (not just active)
% recovery_MS = g_MS * I_MS; %calculate recovery for ALL sites (not just active)
% recovery_HS = g_HS * I_HS; %calculate recovery for ALL sites (not just active)
% new_infections_LS = zeros(length(I_LS), 1); %calculate new infections only at active sites
% new_infections_LS(LS_active) = rate_tot_LS .* S_LS_active;
% new_infections_MS = zeros(length(I_MS), 1);
% new_infections_MS(MS_active) = rate_tot_MS .* S_MS_active;
% new_infections_HS = zeros(length(I_HS), 1);
% new_infections_HS(HS_active) = rate_tot_HS .* S_HS_active;
% dI_LS = new_infections_LS - recovery_LS; %combine
% dI_MS = new_infections_MS - recovery_MS;
% dI_HS = new_infections_HS - recovery_HS;

dR_LS = g_LS * I_LS;
dR_MS = g_MS * I_MS;
dR_HS = g_HS * I_HS;

f = [dS_LS; dS_MS; dS_HS; dI_LS; dI_MS; dI_HS; dR_LS; dR_MS; dR_HS];

end



% % before attempting to implement Dobbelaere transmission smooth
% function f = ODEfun_SCTLD_seascape(t,Y,N_LS,N_MS,N_HS,b_LS,b_MS,b_HS,g_LS,g_MS,g_HS,export_thresh,flux_scale,flux_shape,num_sites,conn_structs,conn_days)
% 
% %% Select correct connectivity matrix at ODE time t
% persistent last_conn_index conn_cache
% if isempty(last_conn_index); last_conn_index = -inf; conn_cache = []; end
% 
% tday = floor(t);
% conn_index = max(1, find(conn_days <= tday, 1, 'last'));
% if isempty(conn_index); conn_index = 1; end
% 
% if conn_index ~= last_conn_index
%     conn_cache = conn_structs(conn_index).full;
%     last_conn_index = conn_index;
% end
% 
% %% Extract compartment states at ODE time t
% % Susceptible
% S_LS = Y(1:num_sites);
% S_MS = Y(num_sites+1:num_sites*2);
% S_HS = Y(num_sites*2+1:num_sites*3);
% S_LS = S_LS(:);
% S_MS = S_MS(:);
% S_HS = S_HS(:);
% 
% % Infected
% I_LS = Y(num_sites*3+1:num_sites*4);
% I_MS = Y(num_sites*4+1:num_sites*5);
% I_HS = Y(num_sites*5+1:num_sites*6);
% I_LS = I_LS(:);
% I_MS = I_MS(:);
% I_HS = I_HS(:);
% 
% % Recovered
% R_LS = Y(num_sites*6+1:num_sites*7);
% R_MS = Y(num_sites*7+1:num_sites*8);
% R_HS = Y(num_sites*8+1:num_sites*9);
% R_LS = R_LS(:);
% R_MS = R_MS(:);
% R_HS = R_HS(:);
% 
% % Define presence of each compartment
% LS_present = N_LS > 0;
% MS_present = N_MS > 0;
% HS_present = N_HS > 0;
% LS_present = LS_present(:);
% MS_present = MS_present(:);
% HS_present = HS_present(:);
% 
% %% Define upstream disease influence
% % Local disease pools
% P = I_LS + I_MS + I_HS;
% P = P(:);
% 
% % Disease flux
% incoming_flux = incomingRisk_sparse(conn_cache, P, export_thresh, 'rowsAreSources');
% incoming_flux = flux_scale .* (1 - exp(-flux_shape .* incoming_flux(:))) / (1 - exp(-flux_shape));
% incoming_flux = incoming_flux(:);
% 
% disease_present = P > 0 | incoming_flux > 0;
% disease_present = disease_present(:);
% 
% LS_active = LS_present & disease_present;
% MS_active = MS_present & disease_present;
% HS_active = HS_present & disease_present;
% 
% %% Extract active site values
% N_LS_active = N_LS(LS_active);
% N_MS_active = N_MS(MS_active);
% N_HS_active = N_HS(HS_active);
% 
% S_LS_active = S_LS(LS_active);
% S_MS_active = S_MS(MS_active);
% S_HS_active = S_HS(HS_active);
% 
% I_LS_active = I_LS(LS_active);
% I_MS_active = I_MS(MS_active);
% I_HS_active = I_HS(HS_active);
% 
% flux_LS = incoming_flux(LS_active);
% flux_MS = incoming_flux(MS_active);
% flux_HS = incoming_flux(HS_active);
% 
% %% ========== SWITCHABLE TRANSMISSION FORMULATION ==========
% % Comment/uncomment ONE of the two blocks below
% 
% % ===== VERSION 1: RATE-BASED (matches R code) =====
% % Local transmission uses TOTAL disease pool P (all groups), not group-specific I
% rate_loc_LS = b_LS * (P(LS_active) ./ N_LS_active);
% rate_loc_MS = b_MS * (P(MS_active) ./ N_MS_active);
% rate_loc_HS = b_HS * (P(HS_active) ./ N_HS_active);
% 
% rate_ext_LS = b_LS * flux_LS;
% rate_ext_MS = b_MS * flux_MS;
% rate_ext_HS = b_HS * flux_HS;
% 
% rate_tot_LS = rate_loc_LS + rate_ext_LS;
% rate_tot_MS = rate_loc_MS + rate_ext_MS;
% rate_tot_HS = rate_loc_HS + rate_ext_HS;
% 
% % % ===== VERSION 2: PROBABILITY-BASED (original formulation) =====
% % % Local transmission uses TOTAL disease pool P (all groups), not group-specific I
% % fI_LS = (P(LS_active) ./ N_LS_active);
% % fI_MS = (P(MS_active) ./ N_MS_active);
% % fI_HS = (P(HS_active) ./ N_HS_active);
% % 
% % p_loc_LS = 1 - exp(-b_LS * fI_LS);
% % p_loc_MS = 1 - exp(-b_MS * fI_MS);
% % p_loc_HS = 1 - exp(-b_HS * fI_HS);
% % 
% % p_ext_LS = 1 - exp(-b_LS * flux_LS);
% % p_ext_MS = 1 - exp(-b_MS * flux_MS);
% % p_ext_HS = 1 - exp(-b_HS * flux_HS);
% % 
% % rate_tot_LS = p_loc_LS + p_ext_LS - p_loc_LS .* p_ext_LS;
% % rate_tot_MS = p_loc_MS + p_ext_MS - p_loc_MS .* p_ext_MS;
% % rate_tot_HS = p_loc_HS + p_ext_HS - p_loc_HS .* p_ext_HS;
% 
% %% ========================================================
% 
% %% Calculate derivatives
% dS_LS = zeros(length(S_LS), 1);
%     dS_LS(LS_active) = -rate_tot_LS .* S_LS_active;
% dS_MS = zeros(length(S_MS), 1);
%     dS_MS(MS_active) = -rate_tot_MS .* S_MS_active;
% dS_HS = zeros(length(S_HS), 1);
%     dS_HS(HS_active) = -rate_tot_HS .* S_HS_active;
% 
% % not explicitly taking out dying coral at each site - could be an issue
% dI_LS = zeros(length(I_LS), 1);
%     dI_LS(LS_active) = rate_tot_LS .* S_LS_active - g_LS * I_LS_active;
% dI_MS = zeros(length(I_MS), 1);
%     dI_MS(MS_active) = rate_tot_MS .* S_MS_active - g_MS * I_MS_active;
% dI_HS = zeros(length(I_HS), 1);
%     dI_HS(HS_active) = rate_tot_HS .* S_HS_active - g_HS * I_HS_active;
% 
% % % failsafe to remove at every location
% % recovery_LS = g_LS * I_LS; %calculate recovery for ALL sites (not just active)
% % recovery_MS = g_MS * I_MS; %calculate recovery for ALL sites (not just active)
% % recovery_HS = g_HS * I_HS; %calculate recovery for ALL sites (not just active)
% % new_infections_LS = zeros(length(I_LS), 1); %calculate new infections only at active sites
% % new_infections_LS(LS_active) = rate_tot_LS .* S_LS_active;
% % new_infections_MS = zeros(length(I_MS), 1);
% % new_infections_MS(MS_active) = rate_tot_MS .* S_MS_active;
% % new_infections_HS = zeros(length(I_HS), 1);
% % new_infections_HS(HS_active) = rate_tot_HS .* S_HS_active;
% % dI_LS = new_infections_LS - recovery_LS; %combine
% % dI_MS = new_infections_MS - recovery_MS;
% % dI_HS = new_infections_HS - recovery_HS;
% 
% dR_LS = g_LS * I_LS;
% dR_MS = g_MS * I_MS;
% dR_HS = g_HS * I_HS;
% 
% f = [dS_LS; dS_MS; dS_HS; dI_LS; dI_MS; dI_HS; dR_LS; dR_MS; dR_HS];
% 
% end