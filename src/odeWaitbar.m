function status = odeWaitbar(t, y, flag)
% ODE output function that shows a progress waitbar for any ode* solver.
% Save as: odeWaitbar.m
persistent wb t0 tf

status = 0; % 0 = continue, 1 = stop

switch flag
    case 'init'
        % On init, the solver calls with t = [t0 tf]
        t0 = t(1);
        tf = t(end);
        wb = waitbar(0, 'Integrating...', 'Name', 'ode progress');

    case []
        % Regular step(s): t is a vector; use last time reached
        if ~isempty(t) && ~isempty(tf)
            p = (t(end) - t0) / max(tf - t0, eps);
            p = max(0, min(1, p));

            habs = numel(y(:,end)) / 9; % total sites if 9 compartments
            Istart = 3*habs + 1; Iend = 6*habs;   % indices for infected blocks
            Ycurr = y(:,end);
            active = sum(Ycurr(Istart:Iend) > 1e-6); % active disease sites

            if isvalid(wb)
                waitbar(p, wb, ...
                    sprintf('t = %.1f (%.0f%%)  |  Active compartments: %d', ...
                    t(end), 100*p, active));
            end
        end

    case 'done'
        if ~isempty(wb) && isvalid(wb), close(wb); end
        % Clear persistent state for next run
        wb = []; t0 = []; tf = [];
end
end