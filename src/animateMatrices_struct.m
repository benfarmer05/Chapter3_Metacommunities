function animateMatrices_struct(S, opts)
% Streaming tween animation for struct array S with fields:
%   S(k).full : sparse/dense matrix (values in [0,1] for logit)
%   S(k).date : datetime (or datenum convertible)
%
% opts fields (all optional):
%   nTween  : # of in-between frames per pair (default 14)
%   fps     : frames per second (default 24)
%   ease    : 'linear'|'cos'|'quad'|'cubic' (default 'cos')
%   transform : 'log10'|'logit'|'asinh'|'none' (default 'log10')
%   eps     : small clamp for log/logit (default 1e-6)
%   clim    : [lo hi] fixed color limits; if empty, auto from pass 1 (default [])
%   export  : filename to write a video (e.g., 'out.mp4'); if empty, no export

% ---- defaults
if nargin < 2, opts = struct; end
def = struct('nTween',14,'fps',24,'ease','cos', ...
             'transform','log10','eps',1e-6,'clim',[], 'export','');
fn = fieldnames(def);
for i=1:numel(fn), if ~isfield(opts,fn{i}), opts.(fn{i}) = def.(fn{i}); end, end

assert(isstruct(S) && isfield(S,'full') && isfield(S,'date'), ...
    'Input must be a struct array with fields .full and .date');

T = numel(S);
assert(T >= 2, 'Need at least two matrices to tween.');
[n1,n2] = size(S(1).full);

% ---- Transform function
switch lower(opts.transform)
    case 'log10'
        tfun = @(X) log10(max(full(X), opts.eps));    % values <= eps -> log10(eps)
    case 'logit'
        % works best when X in [0,1]; clamp to avoid ±Inf
        tfun = @(X) log( min(max(full(X),opts.eps),1-opts.eps) ./ ...
                         max(1 - min(max(full(X),opts.eps),1-opts.eps), opts.eps) );
    case 'asinh'
        tfun = @(X) asinh(full(X));                   % symmetric for small/large
    case 'none'
        tfun = @(X) full(X);
    otherwise
        error('Unknown transform: %s', opts.transform);
end

% ---- Ease function
ease = lower(opts.ease);
switch ease
    case 'linear',  easefun = @(a) a;
    case 'cos',     easefun = @(a) 0.5 - 0.5*cos(pi*a);
    case 'quad',    easefun = @(a) a.^2;
    case 'cubic',   easefun = @(a) a.^3;
    otherwise,      error('Unknown ease: %s', opts.ease);
end

% ---- Pass 1: determine color limits in transform space (streaming)
if isempty(opts.clim)
    lo = inf; hi = -inf;
    for k = 1:T
        Tk = tfun(S(k).full);
        % finite only
        f = isfinite(Tk);
        if any(f,'all')
            lo = min(lo, min(Tk(f),[],'all'));
            hi = max(hi, max(Tk(f),[],'all'));
        end
    end
    if ~isfinite(lo) || ~isfinite(hi)
        lo = 0; hi = 1;  % fallback
    end
    clim = [lo hi];
else
    clim = opts.clim;
end

% ---- Setup figure/video
hfig = figure('Color','w'); ax = axes('Parent',hfig);
axis(ax,'image'); ax.YDir = 'normal';
vid = [];
if ~isempty(opts.export)
    [~,~,ext] = fileparts(opts.export);
    if strcmpi(ext,'.mp4')
        vid = VideoWriter(opts.export, 'MPEG-4');
    else
        vid = VideoWriter(opts.export); % let MATLAB choose
    end
    vid.FrameRate = opts.fps;
    open(vid);
end

% ---- Animate (streaming)
nTween = opts.nTween;
for k = 1:(T-1)
    % endpoints in transform space
    T0 = tfun(S(k).full);
    T1 = tfun(S(k+1).full);

    % ensure consistent size (paranoia)
    if ~isequal(size(T0), [n1 n2]) || ~isequal(size(T1), [n1 n2])
        error('All matrices must be the same size.');
    end

    % include the left endpoint frame
    for j = 0:nTween
        a  = j / (nTween+1);      % 0 .. <1
        aw = easefun(a);
        F  = (1-aw).*T0 + aw.*T1; % interpolate in transform space

        imagesc(ax, F);
        colormap(ax, parula);     %#ok<*MCOL> (use your colormap of choice)
        caxis(ax, clim);
        colorbar(ax);

        % Title: show date range & transform
        t0 = S(k).date;  t1 = S(k+1).date;
        if ~isa(t0,'datetime'), t0 = datetime(t0,'ConvertFrom','datenum'); end
        if ~isa(t1,'datetime'), t1 = datetime(t1,'ConvertFrom','datenum'); end
        title(ax, sprintf('%s → %s  |  %s tween (a=%.2f)', ...
              datestr(t0,'yyyy-mm-dd'), datestr(t1,'yyyy-mm-dd'), ...
              upper(opts.transform), aw));

        drawnow;
        if ~isempty(vid), writeVideo(vid, getframe(hfig)); end
        pause(1/opts.fps);
    end
end

% finally, show the very last slice exactly
Tlast = tfun(S(T).full);
imagesc(ax, Tlast);
colormap(ax, parula); caxis(ax, clim); colorbar(ax);
if ~isempty(vid), writeVideo(vid, getframe(hfig)); close(vid); end
end