classdef MassBalanceBiovolume
% MassBalanceBiovolume
% ------------------------------------------------------------
% Clean "simple budget" in BIOVOLUME space.
%
% M(t)        = column inventory of biovolume      [cm^3/cm^2]
% dM/dt_fd    = finite-difference derivative       [cm^3/cm^2/day]
% EX(t)       = bottom export biovolume flux       [cm^3/cm^2/day]
% PP(t)       = biovolume production from growth   [cm^3/cm^2/day]
%
% Residual:
%   RES_simple = dM/dt_fd - (PP - EX)
%
% IMPORTANT:
% - For biovolume export, prefer:
%     out.output_data.bottom_fluxsect_weighted   (cm^3/m^2/day)
% - The field named bottom_fluxsect_cm3m2d is actually number flux
%   (#/m^2/day) in the current OutputGenerator, so it MUST be converted.
%
% ------------------------------------------------------------

methods(Static)

function mb = compute(sim, out, cfg, opts)
    if nargin < 4, opts = struct(); end
    opts = MassBalanceBiovolume.fillDefaults(opts);

    t = out.time(:);
    Y = out.concentrations;
    nt = numel(t);

    % ---- Section grid ----
    gridS = sim.rhs.getSectionGrid();

    % Ns
    Ns = [];
    try
        if isprop(cfg,'n_sections') && ~isempty(cfg.n_sections)
            Ns = cfg.n_sections;
        end
    catch
    end
    if isempty(Ns)
        r_cm = MassBalanceBiovolume.getRadii_cm(gridS);
        Ns = numel(r_cm);
    end

    % Nz
    is_col = isprop(cfg,'use_column') && cfg.use_column;
    if is_col
        Nz = cfg.getNumLayers();
    else
        Nz = 1;
    end

    dz_cm = cfg.dz * 100; % m -> cm

    % ---- Vbin (cm^3 per particle) ----
    Vbin = MassBalanceBiovolume.getVbin_cm3(gridS, Ns);

    % ---- Inventory M(t): cm^3/cm^2 ----
    M = zeros(nt,1);
    for it = 1:nt
        vflat = Y(it,:).';
        M(it) = MassBalanceBiovolume.integrateColumn(vflat, Vbin, Ns, Nz, dz_cm, "integrated");
    end

    % ---- dM/dt finite difference ----
    dMdt_fd = nan(nt,1);
    dMdt_fd(2:end) = diff(M) ./ diff(t);

    % ---- Export EX(t): cm^3/cm^2/day ----
    [EX, EX_source] = MassBalanceBiovolume.computeExport(sim, out, cfg, Ns, Nz, Vbin, dz_cm, opts);

    % ---- PP(t): cm^3/cm^2/day ----
    [PP, PP_source] = MassBalanceBiovolume.computePPfromGrowth(out, cfg, Ns, Nz, Vbin, dz_cm, opts);

    % ---- Residual ----
    RES_simple = dMdt_fd - (PP - EX);

    % ---- Pack ----
    mb = struct();
    mb.t = t;
    mb.M = M;
    mb.dMdt_fd = dMdt_fd;
    mb.EX = EX;
    mb.PP = PP;
    mb.RES_simple = RES_simple;

    mb.Ns = Ns;
    mb.Nz = Nz;
    mb.dz_cm = dz_cm;
    mb.Vbin = Vbin;

    mb.EX_source = EX_source;
    mb.PP_source = PP_source;

    if opts.do_print
        ok = isfinite(RES_simple) & ((1:nt)' > 1);
        fprintf('MassBalanceBiovolume:\n');
        fprintf('  EX source: %s\n', mb.EX_source);
        fprintf('  PP source: %s\n', mb.PP_source);
        fprintf('  RES_simple median(|.|)=%.3e  max(|.|)=%.3e\n', ...
            median(abs(RES_simple(ok))), max(abs(RES_simple(ok))));
    end
end

% ==========================================================
% Export
% ==========================================================
function [EX, EX_source] = computeExport(sim, out, cfg, Ns, Nz, Vbin, dz_cm, opts)
    t = out.time(:);
    Y = out.concentrations;
    nt = numel(t);

    EX = zeros(nt,1);
    EX_source = "none";

    % Prefer OutputGenerator weighted export (biovolume)
    if opts.use_output_flux
        try
            if isfield(out,'output_data') && ~isempty(out.output_data)
                od = out.output_data;

                % 1) TRUE BIOVOLUME export (cm^3/m^2/day) -> /1e4 => cm^3/cm^2/day
                if isfield(od,'bottom_fluxsect_weighted') && ~isempty(od.bottom_fluxsect_weighted)
                    Fw = od.bottom_fluxsect_weighted;           % [nt x Ns], cm^3/m^2/day
                    EX = sum(max(Fw,0), 2) / 1e4;               % cm^3/cm^2/day
                    EX_source = "out.output_data.bottom_fluxsect_weighted (cm^3/m^2/d -> cm^3/cm^2/d)";
                    return;
                end

                % 2) If only the NUMBER flux exists (really #/m^2/day), convert using Vbin
                if isfield(od,'bottom_fluxsect_cm3m2d') && ~isempty(od.bottom_fluxsect_cm3m2d)
                    Fn = od.bottom_fluxsect_cm3m2d;             % [nt x Ns], actually #/m^2/day
                    EX = (max(Fn,0) * Vbin(:)) / 1e4;           % (cm^3/m^2/day)/1e4
                    EX_source = "out.output_data.bottom_fluxsect_cm3m2d converted using Vbin (treated as #/m^2/d)";
                    return;
                end
            end
        catch
        end
    end

    % Fallback: compute from bottom state + settling velocity if available
    % (NOTE: sim.operators.sink_rate is 1/day loss-rate, not cm/day velocity)
    try
        if isfield(sim.operators,'set_vel_cmday') && ~isempty(sim.operators.set_vel_cmday)
            w_cmday = sim.operators.set_vel_cmday(:); % [Ns x 1]
            for it = 1:nt
                vflat = Y(it,:).';
                N2 = reshape(vflat, [Ns Nz]);
                ybot = N2(:,end); % #/cm^3
                EX(it) = sum( max(ybot,0) .* max(w_cmday,0) .* Vbin ); % cm^3/cm^2/day
            end
            EX_source = "computed from bottom state using set_vel_cmday * Vbin (cm-units)";
            return;
        end
    catch
    end

    % Last fallback: use sink_rate (1/day) as bottom loss * dz to make a flux
    try
        if isfield(sim.operators,'sink_rate') && ~isempty(sim.operators.sink_rate) && Nz > 1
            sink_rate = sim.operators.sink_rate(:); % 1/day  (≈ w/dz)
            for it = 1:nt
                vflat = Y(it,:).';
                N2 = reshape(vflat, [Ns Nz]);
                ybot = N2(:,end); % #/cm^3
                loss = max(ybot,0) .* max(sink_rate,0);         % #/cm^3/day
                EX(it) = sum(loss .* Vbin) * dz_cm;             % cm^3/cm^2/day
            end
            EX_source = "computed from bottom state using sink_rate (1/d) * dz_cm * Vbin";
            return;
        end
    catch
    end

    EX_source = "export unavailable (returned zeros)";
end

% ==========================================================
% PP from growth (Adrian PP)
% ==========================================================
function [PP, PP_source] = computePPfromGrowth(out, cfg, Ns, Nz, Vbin, dz_cm, opts)
    t = out.time(:);
    Y = out.concentrations;
    nt = numel(t);

    mu = MassBalanceBiovolume.getGrowthRatePerDay(cfg); % 1/day
    ib = opts.pp_bin;

    PP = zeros(nt,1);

    if mu == 0
        PP_source = "growth rate is zero (PP=0)";
        return;
    end

    if Nz == 1
        % 0-D box: PP = mu * N(bin) * Vbin
        for it = 1:nt
            vflat = Y(it,:).';
            N2 = reshape(vflat, [Ns Nz]);
            PP(it) = mu * max(N2(ib,1),0) * Vbin(ib);
        end
        PP_source = sprintf("growth(mu=%.3g 1/d) on bin=%d (0-D)", mu, ib);
        return;
    end

    % Column: default is SURFACE layer only (this matches your "surface PP" idea)
    if strcmpi(opts.pp_apply, "surface")
        kz = 1;
        for it = 1:nt
            vflat = Y(it,:).';
            N2 = reshape(vflat, [Ns Nz]);
            PP(it) = mu * max(N2(ib,kz),0) * Vbin(ib) * dz_cm; % cm^3/cm^2/day
        end
        PP_source = sprintf("growth(mu=%.3g 1/d) on bin=%d at surface layer only", mu, ib);
        return;
    end

    % Optional: apply growth in ALL layers (usually not what you want)
    if strcmpi(opts.pp_apply, "all")
        for it = 1:nt
            vflat = Y(it,:).';
            N2 = reshape(vflat, [Ns Nz]);
            PP(it) = mu * sum(max(N2(ib,:),0)) * Vbin(ib) * dz_cm;
        end
        PP_source = sprintf("growth(mu=%.3g 1/d) on bin=%d across ALL layers", mu, ib);
        return;
    end

    % fallback
    PP_source = "PP apply mode not recognized (PP=0)";
end

% ==========================================================
% Helpers
% ==========================================================
function opts = fillDefaults(opts)
    if ~isfield(opts,'use_output_flux'), opts.use_output_flux = true; end
    if ~isfield(opts,'do_print'),        opts.do_print        = true; end
    if ~isfield(opts,'pp_bin'),          opts.pp_bin          = 1;    end
    if ~isfield(opts,'pp_apply'),        opts.pp_apply        = "surface"; end % "surface" or "all"
end

function r_cm = getRadii_cm(gridS)
    r_cm = [];
    try
        if ismethod(gridS,'getConservedRadii')
            r_cm = gridS.getConservedRadii();
        else
            r_cm = gridS.getFractalRadii();
        end
    catch
    end
    r_cm = r_cm(:);
    if isempty(r_cm)
        error('MassBalanceBiovolume: cannot get radii from section grid.');
    end
end

function Vbin = getVbin_cm3(gridS, Ns)
    Vbin = [];
    try
        if ismethod(gridS,'getBinVolumes')
            Vbin = gridS.getBinVolumes();
            Vbin = Vbin(:);
        end
    catch
    end

    if isempty(Vbin) || numel(Vbin) ~= Ns
        r_cm = MassBalanceBiovolume.getRadii_cm(gridS);
        r_cm = r_cm(:);
        if numel(r_cm) ~= Ns
            r_cm = r_cm(1:Ns);
        end
        Vbin = (4/3) * pi * (r_cm.^3);
    end
    Vbin = Vbin(:);
end

function val = integrateColumn(vflat, w, Ns, Nz, dz_cm, inventory_mode)
    if isempty(vflat), val = 0; return; end

    if Nz == 1
        val = sum(vflat(:) .* w(:));
        return;
    end

    N2 = reshape(vflat(:), [Ns Nz]);
    if inventory_mode == "surface"
        val = sum(N2(:,1) .* w(:));
    else
        W2 = repmat(w(:), 1, Nz);
        val = sum(sum(N2 .* W2)) * dz_cm;
    end
end

function mu = getGrowthRatePerDay(cfg)
    mu = 0;
    try
        if isprop(cfg,'growth') && ~isempty(cfg.growth)
            mu = cfg.growth;
            return;
        end
    catch
    end
    try
        if isprop(cfg,'mu') && ~isempty(cfg.mu)
            mu = cfg.mu;
            return;
        end
    catch
    end
end

end % methods(Static)
end % classdef