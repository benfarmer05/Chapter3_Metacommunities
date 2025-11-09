function T = incomingRisk_sparse(Pint, DP, thresh, orientation)

% T = incomingRisk_sparse(Pint, DP, thresh, orientation)
% Computes T_j = 1 - prod_i (1 - DP(i)*Pint(i,j)) over nonzero edges.
% Pint: n×n sparse (rows=sources, cols=destinations by default)
% DP  : n×1 column vector of source “infection intensity”
% thresh: scalar; sources with DP <= thresh are ignored
% orientation (optional): 'rowsAreSources' (default) or 'colsAreSources'
%
% Returns T as n×1 (column). Transpose at call site if you need a row.

    if nargin < 4 || isempty(orientation)
        orientation = 'rowsAreSources';
    end

    % ---- Coerce shapes, validate sizes ----
    if ~issparse(Pint)
        warning('incomingRisk_sparse:ExpectSparse', ...
                'Pint is not sparse; converting (this may be slow).');
        Pint = sparse(Pint);
    end

    [n1, n2] = size(Pint);
    if n1 ~= n2
        error('incomingRisk_sparse:SquareRequired', ...
              'Pint must be square; got %dx%d.', n1, n2);
    end
    n = n1;

    DP = DP(:);                 % force column
    if numel(DP) ~= n
        error('incomingRisk_sparse:SizeMismatch', ...
              'numel(DP)=%d does not match size(Pint,1)=%d.', numel(DP), n);
    end

    % ---- Extract edges ----
    [is, js, vs] = find(Pint);  % all column vectors

    switch lower(orientation)
        case 'rowsaresources'
            % as-is: i = source, j = dest
        case 'colsaresources'
            % flip if your Pint uses columns as sources
            [is, js] = deal(js, is);
        otherwise
            error('incomingRisk_sparse:BadOrientation', ...
                  'orientation must be ''rowsAreSources'' or ''colsAreSources''.');
    end

    % ---- Filter to active sources (DP > thresh) ----
    active = DP(is) > thresh;
    if ~any(active)
        T = zeros(n,1);         % nothing contributing → zero incoming risk
        return
    end
    is = is(active);
    js = js(active);
    vs = vs(active);

    % ---- Edge weights and stability clamps ----
    w = DP(is) .* vs;           % element-wise, both column
    % numerical safety: keep 0 <= w < 1
    w = max(0, min(w, 1 - eps));

    % ---- Accumulate sum_j log(1 - w_ij) ----
    % s(j) = Σ_i log(1 - w_ij)
    s = accumarray(js, log1p(-w), [n, 1], @sum, 0);

    % ---- T = 1 - exp(Σ log(1 - w)) ----
    T = 1 - exp(s);

    % Guard against any rare numeric issues
    bad = ~isfinite(T) | T < 0 | T > 1;
    if any(bad)
        % Replace non-finite with zeros; clamp to [0,1]
        T(~isfinite(T)) = 0;
        T = max(0, min(T, 1));
    end
end