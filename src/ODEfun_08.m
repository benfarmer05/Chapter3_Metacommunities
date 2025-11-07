function f = ODEfun_08(t,Y,CLS1,CMS1,CHS1,bls,bms,bhs,kls,kms,khs,thresh,c,B,habs,P,Pdays)
% Each location has 3 variables, C (coral cover), S (susceptibility) and D
% (disease prevalence). Disease prevalence is translated back and forth to
% disease cover (%).

% 10/24/25 -DH
% I am playing around with the logic of local vs exogenous infection

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

% % Dead - isnt used.
% DLS = Y(habs*6+1:habs*7);
% DMS = Y(habs*7+1:habs*8);
% DHS = Y(habs*8+1:habs*9);

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
% (definided outside ODE) in order to transmit disease.
% DTP = ILS.*CLS1' + IMS.*CMS1' + IHS.*CHS1'; % mult. prop. by original coral cover
DP = ILS + IMS + IHS; 
DP = DP(:);

% Estimate disease FLUX, depends on upstream DP, connectivity probability,
% and, if using, thresh (minimum DP(i) to contribute to T(j))
T = incomingRisk_sparse(Pint_cache, DP, thresh,'rowsAreSources'); % May need to transpose T, will have to see. 
% maxT = max(T);

% Reshape T - this modulates how accumulated FLUX translates to local
% infections. Very sensitive to this.
% T=shp1.*(1-(1-T(:))).^sh2;
T = c.*(1-exp(-B.*T(:)))/(1-exp(-B));
T=T(:);

% Cumulative probabilities less than thresh2 are zero.
% T(T<thresh2) = 0;

% if sum(T) == 0
%     fprintf(strcat('Warning - Incoming disease does not surpass threshold\n',num2str(maxT)))
% end

D_prez = DP > 0 | T > 0;

D_prez = D_prez(:);

LS_use = LS_prez & D_prez;
MS_use = MS_prez & D_prez;
HS_use = HS_prez & D_prez;

%% The SIR eqns
% I've tried to adapt Ben's model
% Basic SIR with scalable mortality (m)/recovery (k)/loss of protection (r)

% Currently all categories will be initiated with external infection. To
% turn off, remove the T_S term from p_tot calcs...

DPLS = DP(LS_use);
TLS = T(LS_use);

fI_LS = DPLS./CLS1(LS_use);
p_loc_LS = 1 - exp(-bls * fI_LS); % This will never go over 1
p_ext_LS = 1 - exp(-bls * TLS); % Just added this - 

p_tot_LS = p_loc_LS+p_ext_LS - p_loc_LS.*p_ext_LS; % This will never go over 1

% phiLS = ((DPLS(:)+TLS(:))./(1+TLS(:)));

DPMS = DP(MS_use);
TMS = T(MS_use);

fI_MS = DPMS./CMS1(MS_use);
p_loc_MS = 1 - exp(-bms * fI_MS);
p_ext_MS = 1 - exp(-bms * TMS);

p_tot_MS = p_loc_MS+p_ext_MS - p_loc_MS.*p_ext_MS; % This will never go over 1

% phiMS = ((DPMS(:)+TMS(:))./(1+TMS(:)));

DPHS = DP(HS_use);
THS = T(HS_use);

fI_HS = DPHS./CHS1(HS_use);
p_loc_HS = 1 - exp(-bhs * fI_HS);
p_ext_HS = 1 - exp(-bhs * THS);

p_tot_HS = p_loc_HS+p_ext_HS - p_loc_HS.*p_ext_HS; % This will never go over 1


% phiHS = ((DPHS(:)+THS(:))./(1+THS(:)));

dSLS = zeros(length(SLS),1);
    dSLS(LS_use) = -p_tot_LS.*SLS(LS_use); % Pretty sure dividing by N here is causing outbreaks to peter out immediately
dSMS = zeros(length(SMS),1);
    dSMS(MS_use) = -p_tot_MS.*SMS(MS_use);
dSHS = zeros(length(SHS),1);
    dSHS(HS_use) = -p_tot_HS.*SHS(HS_use);

dILS = zeros(length(ILS),1);
    dILS(LS_use) = p_tot_LS.*SLS(LS_use) - kls*ILS(LS_use); %./CLS1(LS_use) - kls*ILS(LS_use);
dIMS = zeros(length(IMS),1);
    dIMS(MS_use) = p_tot_MS.*SMS(MS_use) - kms*IMS(MS_use); %./CMS1(MS_use) - kms*IMS(MS_use);
dIHS = zeros(length(IHS),1);
    dIHS(HS_use) = p_tot_HS.*SHS(HS_use) - khs*IHS(HS_use); %./CHS1(HS_use) - khs*IHS(HS_use);

% dDLS = kls*mls*ILS;
% dDMS = kms*mms*IMS;
% dDHS = khs*mhs*IHS;

dRLS = kls*ILS;
dRMS = kms*IMS;
dRHS = khs*IHS;

% dRLS = (1-mls)*kls*ILS - rls*RLS;
% dRMS = (1-mms)*kms*IMS - rms*RMS;
% dRHS = (1-mhs)*khs*IHS - rhs*RHS;

f = [dSLS;dSMS;dSHS;dILS;dIMS;dIHS;dRLS;dRMS;dRHS]; %dDLS;dDMS;dDHS;

