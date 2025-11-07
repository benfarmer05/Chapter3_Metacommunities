function [P, Wraw] = buildConnectivityKNN(lon, lat, varargin)
%BUILDCONNECTIVITYKNN  Sparse neighbor graph on lon/lat with variable degree & weights
%
% [P, Wraw] = buildConnectivityKNN(lon, lat, 'K', 10, 'KJitter', 3, ...
%                                  'LenScaleKm', 25, 'DropProb', 0.1, ...
%                                  'WeightNoiseCV', 0.25, 'RowNormalize', true)
%
% Inputs:
%   lon, lat       : n×1 (degrees)
%
% Name-Value options (reasonable defaults shown):
%   'K'            : target mean #neighbors per node (default 10)
%   'KJitter'      : degree jitter (std dev for normal, then rounded; default 3)
%   'LenScaleKm'   : distance decay length L in km for weights exp(-d/L) (default 25 km)
%   'DropProb'     : independent random dropout of candidate edges (0–1) (default 0.1)
%   'WeightNoiseCV': coefficient of variation for lognormal weight noise (default 0.25)
%   'RowNormalize' : if true, return row-stochastic P (probabilities) (default true)
%
% Outputs:
%   P    : n×n sparse, row-stochastic if RowNormalize=true (probabilities)
%   Wraw : n×n sparse, raw (unnormalized) weights (before row-normalization)

opts = struct('K',10,'KJitter',3,'LenScaleKm',25,'DropProb',0.10, ...
              'WeightNoiseCV',0.25,'RowNormalize',true,'Seed',[]);
opts = parseOpts(opts, varargin{:});

if ~isempty(opts.Seed), rng(opts.Seed); end

lon = lon(:); lat = lat(:);
n = numel(lon);
assert(numel(lat)==n, 'lon/lat must be same length');

% --- Fast nearest neighbors in 3D unit-sphere space (Euclidean ~ great-circle)
% Convert (lon,lat) to 3D unit vectors
rad = pi/180;
phi = lat*rad; lam = lon*rad;
x = cos(phi).*cos(lam);
y = cos(phi).*sin(lam);
z = sin(phi);
XYZ = [x y z];

% Target neighbor counts per node (>=1). Normal with jitter, clipped.
K = max(1, round(opts.K + opts.KJitter*randn(n,1)));
Kmax = max(K);

% Get Kmax+1 neighbors (includes self)
% Requires Statistics & ML Toolbox for knnsearch. If unavailable, see fallback below.
[idx, dist3D] = knnsearch(XYZ, XYZ, 'K', Kmax+1, 'NSMethod', 'kdtree');
% Remove self (first column should be self index)
idx  = idx(:,2:end);
dist3D = dist3D(:,2:end);

% Convert 3D chord distance to great-circle distance (km)
% chord = 2*sin(theta/2) => theta = 2*asin(chord/2)
theta = 2*asin( min(1, dist3D/2) );      % radians
R_earth_km = 6371;
d_km = R_earth_km * theta;

% Build variable-degree edge lists with optional random dropouts
keepCell_i = cell(n,1);
keepCell_j = cell(n,1);
wCell      = cell(n,1);

% Lognormal noise parameters from desired CV
cv = max(0, opts.WeightNoiseCV);
if cv>0
    sigma2 = log(1 + cv^2);
    mu = -0.5*sigma2;   % so mean(noise) = 1
else
    sigma2 = 0; mu = 0;
end

L = opts.LenScaleKm;
for ii = 1:n
    k_i = K(ii);
    candJ = idx(ii,1:k_i);
    d_i  = d_km(ii,1:k_i);

    % Randomly drop a fraction of edges
    if opts.DropProb>0
        drop = rand(1,k_i) < opts.DropProb;
        candJ = candJ(~drop);
        d_i   = d_i(~drop);
    end
    if isempty(candJ)
        continue
    end

    % Distance-decay base weights
    w = exp(-d_i / L);

    % Optional lognormal multiplicative noise (mean 1)
    if cv>0
        w = w .* exp(mu + sqrt(sigma2)*randn(size(w)));
    end

    % Clip to (0,1) to be probability-friendly
    w = max(0, min(w, 1 - eps));

    keepCell_i{ii} = repmat(ii, numel(candJ), 1);
    keepCell_j{ii} = candJ(:);
    wCell{ii}      = w(:);
end

I = vertcat(keepCell_i{:});
J = vertcat(keepCell_j{:});
Wv = vertcat(wCell{:});

% Raw sparse weights
Wraw = sparse(I, J, Wv, n, n);

% Optional: ensure no self-loops
Wraw = Wraw - spdiags(diag(Wraw), 0, n, n);

% Row-normalize to get probabilities (each row sums to 1 over its kept neighbors)
if opts.RowNormalize
    rowsum = full(sum(Wraw,2));
    scale = rowsum;
    scale(scale==0) = 1;               % avoid divide-by-zero
    P = spdiags(1./scale, 0, n, n) * Wraw;
else
    P = Wraw;
end

end

% ---------- helpers ----------
function opts = parseOpts(opts, varargin)
if mod(numel(varargin),2)~=0
    error('Name-value args must come in pairs');
end
for k = 1:2:numel(varargin)
    name = varargin{k}; val = varargin{k+1};
    if ~isfield(opts, name), error('Unknown option: %s', name); end
    opts.(name) = val;
end
end