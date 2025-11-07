function cmap = stackedColormap(edges, Cstart, Cend, n, gamma)
% Build a stacked-gradient colormap across numeric edges.
% edges: 1x(K+1) sorted vector of data edges (e.g., [-20 -10 0 10 20])
% Cstart: Kx3 RGB start colors for each segment
% Cend:   Kx3 RGB end colors for each segment
% n: total number of colors in the final map (default 256)
% gamma: 1xK easing per segment (default 1's). >1 pushes weight to centers; <1 to ends.

if nargin < 4 || isempty(n), n = 256; end
K = numel(edges) - 1;
assert(K >= 1, 'edges must have at least two values.');
assert(all(size(Cstart) == [K 3]) && all(size(Cend) == [K 3]), ...
    'Cstart and Cend must be Kx3.');

if nargin < 5 || isempty(gamma), gamma = ones(1,K); end
assert(numel(gamma) == K, 'gamma must be length K.');

% Allocate steps per segment proportional to value width (ensure ≥2 each)
widths = diff(edges);
widths(widths < eps) = eps;
steps = max(2, round(n * widths / sum(widths)));

% Adjust to hit exactly n colors
delta = n - sum(steps);
while delta ~= 0
    [~, idx] = max(widths);              % add/remove where widest
    take = sign(delta);
    steps(idx) = steps(idx) + take;
    delta = n - sum(steps);
end

% Build each segment with per-segment easing
chunks = cell(K,1);
for k = 1:K
    m = steps(k);
    t = linspace(0,1,m)';
    g = gamma(k);
    if ~isfinite(g) || g <= 0, g = 1; end
    t = t.^g;  % simple power-easing
    % interpolate
    seg = (1-t).*Cstart(k,:) + t.*Cend(k,:);
    chunks{k} = seg;
    if k < K
        chunks{k} = seg(1:end-1,:);  % drop last to avoid duplicate seams
    end
end
cmap = vertcat(chunks{:});
cmap = max(0, min(1, cmap));       % clamp
end