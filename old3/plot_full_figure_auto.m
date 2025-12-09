function plot_full_figure_auto(outdir, tag)
%PLOT_FULL_FIGURE_AUTO  Pulls res/sim from base and calls plot_full_figure.
% Usage:  plot_full_figure_auto('/path/to/out','baseline')

    if nargin < 1 || isempty(outdir), outdir = pwd; end
    if nargin < 2, tag = 'run'; end
    if ~exist(outdir,'dir'), mkdir(outdir); end

    res = evalin('base','res');
    sim = evalin('base','sim');

    od  = OutputGenerator.spectraAndFluxes(res.time, res.concentrations, sim.grid, sim.config);
    plot_full_figure(res, od, sim.grid, sim.config, outdir, tag);
end