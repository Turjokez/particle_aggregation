function [t_days, eps_series, file, used] = load_epsilon_series(matpath, t0, varargin)
% LOAD_EPSILON_SERIES  Load ε(t) from a .mat file (or a folder) and return daily series.
%
% Usage:
%   [t_days, eps_series] = load_epsilon_series('/path/epsilon_daily.mat', 0, ...)
%   [t_days, eps_series] = load_epsilon_series('/path/to/folder', 0, ...)
%
% Key name/value options (all optional):
%   'time_field'   : e.g., 'mtime' or 'S.mtime'
%   'eps_field'    : e.g., 'eps'   or 'S.eps' (Nz×Nt)
%   'depth_field'  : e.g., 'z'     or 'S.z'   (Nz×1) negative up to surface ~0
%   'mld_field'    : e.g., 'mld'   or 'S.mld' (1×Nt) usually negative (depth)
%   'unit'         : 'auto' (default) | 'linear' | 'log10' | 'dB'
%   'agg'          : 'mld' (default) | 'topN' | 'surface' | 'none'
%   'topN'         : 3 (default) when agg='topN' or MLD is NaN
%   'resample_daily' : true (default) — average within day bins
%   'smooth_win'     : 0 or integer (default 0; set 3 or 5 for gentle smoothing)
%
% Returns:
%   t_days      : column vector, days since start + t0
%   eps_series  : column vector, ε in W kg^-1 (linear)
%   file        : file actually used
%   used        : struct of used options/paths (for reproducibility)

% -------------------------- parse options -------------------------------
opt = struct( ...
    'time_field', '', ...
    'eps_field',  '', ...
    'depth_field','', ...
    'mld_field',  '', ...
    'unit', 'auto', ...
    'agg',  'mld', ...
    'topN', 3, ...
    'resample_daily', true, ...
    'smooth_win', 0);
if ~isempty(varargin), opt = parse_opts(opt, varargin{:}); end
if nargin < 2 || isempty(t0), t0 = 0; end

% --------------------- accept folder OR file ----------------------------
if isfolder(matpath)
    cand = [ dir(fullfile(matpath,'*epsilon*.mat')); ...
             dir(fullfile(matpath,'*eps*.mat')); ...
             dir(fullfile(matpath,'*.mat')) ];
    assert(~isempty(cand), 'No .mat files found in folder: %s', matpath);
    file = fullfile(cand(1).folder, cand(1).name);
else
    file = matpath;
end
assert(isfile(file), 'File not found: %s', file);

raw0 = load(file);
raw  = flatten_if_single_struct(raw0);  % allow S.<fields>

% -------------------- fetch fields (plain or dotted) --------------------
time_vec = fetch_field(raw, opt.time_field, {'mtime','time','t','days','day'});
assert(~isempty(time_vec), 'Could not find a time vector. Set ''time_field''.');

% eps is Nz×Nt; depth is Nz×1; mld is 1×Nt
Ez   = fetch_field(raw, opt.eps_field,  {'eps','epsilon','eps_ml','epsilon_ml'});
z    = fetch_field(raw, opt.depth_field,{'z','depth'});
mld  = fetch_field(raw, opt.mld_field,  {'mld'});
assert(~isempty(Ez), 'Could not find epsilon array. Set ''eps_field''.');
assert(~isempty(z),  'Could not find depth vector. Set ''depth_field''.');

% shape/consistency
z = z(:);
if size(Ez,1) ~= numel(z) && size(Ez,2) == numel(z)
    Ez = Ez.'; % transpose if provided Nt×Nz
end
assert(size(Ez,1) == numel(z), 'Size mismatch: size(Ez,1) must equal numel(z).');

% --------------------- convert time to days -----------------------------
t_days_raw = to_days_since_start(time_vec);
t_days_raw = t_days_raw - t_days_raw(1) + t0;  % start at t0

% --------------------- unit conversion to linear ------------------------
Ez = double(Ez);
switch lower(opt.unit)
    case 'linear'
        % nothing
    case 'log10'
        Ez = 10.^Ez;
    case 'db'
        Ez = 10.^(Ez/10);
    case 'auto'
        Ez = auto_unit_convert(Ez);
    otherwise
        error('Unknown unit: %s', opt.unit);
end
Ez(Ez<0) = NaN;  % guardrail

% --------------------- aggregate in vertical ----------------------------
eps_col = aggregate_eps(Ez, z, mld, t_days_raw, opt);

% --------------------- daily resampling & smoothing ---------------------
[t_days, eps_series] = daily_average(t_days_raw, eps_col, opt.resample_daily);
if opt.smooth_win > 0
    eps_series = smoothdata(eps_series, 'movmean', opt.smooth_win, 'omitnan');
end

% final shape
t_days     = t_days(:);
eps_series = eps_series(:);

% record choices
used = opt;
end

% =======================================================================
function t_days = to_days_since_start(time_vec)
    % accepts datetime, datenum-like, or numeric days
    if isdatetime(time_vec)
        t_days = days(time_vec - time_vec(1));
    else
        time_vec = time_vec(:);
        if max(time_vec) > 7e4  % looks like datenum
            dt = datetime(time_vec, 'ConvertFrom','datenum');
            t_days = days(dt - dt(1));
        else
            t_days = time_vec - time_vec(1);
        end
    end
    t_days = double(t_days);
end

function Ez_lin = auto_unit_convert(Ez)
    % Try to infer unit: if mostly negative and >= -300, treat as dB.
    % If typical range [-5, +5], treat as log10.
    % Else assume linear.
    x = Ez(isfinite(Ez));
    if isempty(x)
        Ez_lin = Ez; return;
    end
    q = quantile(x, [0.05 0.5 0.95]);
    if all(x <= 0 | isnan(x)) && q(2) > -300 && q(2) < -10
        Ez_lin = 10.^(Ez/10);   % dB
    elseif q(1) > -10 && q(3) < 10
        Ez_lin = 10.^Ez;        % log10
    else
        Ez_lin = Ez;            % linear
    end
end

function eps_col = aggregate_eps(Ez, z, mld, t, opt)
    Nt = size(Ez,2);
    eps_col = nan(Nt,1);
    have_mld = ~isempty(mld);
    if have_mld
        mld = double(mld(:));   % often negative
    end
    for j = 1:Nt
        col = Ez(:,j);
        if all(~isfinite(col)), continue; end
        switch lower(opt.agg)
            case 'none'
                % take shallowest valid bin
                k = find(isfinite(col), 1, 'last'); % z often negative ascending to 0
                if ~isempty(k), eps_col(j) = col(k); end
            case 'surface'
                % pick z closest to 0
                [~,k] = min(abs(z-0));
                eps_col(j) = col(k);
            case 'topn'
                % mean of the N shallowest finite bins
                [~,ord] = sort(z,'descend'); % shallow first (towards 0)
                ord = ord(isfinite(col(ord)));
                ord = ord(1:min(opt.topN, numel(ord)));
                if ~isempty(ord), eps_col(j) = mean(col(ord),'omitnan'); end
            case 'mld'
                if have_mld && isfinite(mld(j))
                    cap = -abs(mld(j));             % mld negative (depth)
                    idx = z >= cap & z <= 0;        % within mixed layer to surface
                    if any(idx) && any(isfinite(col(idx)))
                        eps_col(j) = mean(col(idx), 'omitnan');
                    else
                        % fallback to topN
                        [~,ord] = sort(z,'descend');
                        ord = ord(isfinite(col(ord)));
                        ord = ord(1:min(opt.topN, numel(ord)));
                        if ~isempty(ord), eps_col(j) = mean(col(ord),'omitnan'); end
                    end
                else
                    % no MLD → topN fallback
                    [~,ord] = sort(z,'descend');
                    ord = ord(isfinite(col(ord)));
                    ord = ord(1:min(opt.topN, numel(ord)));
                    if ~isempty(ord), eps_col(j) = mean(col(ord),'omitnan'); end
                end
            otherwise
                error('Unknown agg option: %s', opt.agg);
        end
    end
end

function [td, y] = daily_average(t, x, do_daily)
    t = t(:); x = x(:);
    ok = isfinite(t) & isfinite(x);
    t = t(ok); x = x(ok);
    if isempty(t)
        td = t; y = x; return;
    end
    if ~do_daily
        td = t; y = x; return;
    end
    d0 = floor(min(t));
    d1 = floor(max(t));
    td = (d0:d1).';
    y  = nan(numel(td),1);
    for k = 1:numel(td)
        sel = (t >= td(k)) & (t < td(k)+1);
        if any(sel), y(k) = mean(x(sel),'omitnan'); end
    end
end

function out = fetch_field(raw, explicit, guesses)
    if ~isempty(explicit)
        out = get_by_path(raw, explicit);
        return;
    end
    out = [];
    if isstruct(raw)
        fns = fieldnames(raw);
        for k = 1:numel(guesses)
            if any(strcmpi(fns, guesses{k}))
                out = raw.(guesses{k});
                return;
            end
        end
    end
end

function v = get_by_path(s, path)
    % supports dotted lookup, e.g. 'S.mtime'
    v = s;
    parts = split(path, '.');
    for i = 1:numel(parts)
        p = strtrim(parts{i});
        if isstruct(v) && isfield(v, p)
            v = v.(p);
        else
            v = [];
            return;
        end
    end
end

function S = flatten_if_single_struct(S)
    % If the .mat contains a single top-level struct (e.g., S),
    % flatten it so you can refer to plain names.
    if isstruct(S)
        f = fieldnames(S);
        if numel(f)==1 && isstruct(S.(f{1}))
            S = S.(f{1});
        end
    end
end

function opt = parse_opts(opt, varargin)
    assert(mod(numel(varargin),2)==0, 'Name/value pairs expected.');
    for i = 1:2:numel(varargin)
        k = varargin{i}; v = varargin{i+1};
        if isfield(opt,k), opt.(k) = v; end
    end
end