function save_all_figs(outdir, prefix)
%SAVE_ALL_FIGS  Save every open figure in numeric order.
%   save_all_figs(outdir) or save_all_figs(outdir,'tag')

    if nargin < 2, prefix = 'fig'; end
    if ~exist(outdir,'dir'), mkdir(outdir); end

    figs = findall(groot,'Type','figure');
    if isempty(figs), return; end

    nums = arrayfun(@(f) f.Number, figs);
    [~,ord] = sort(nums);                % <-- safe sort for handles
    figs = figs(ord);

    for f = figs'
        % hide axes toolbars when present (no axtoolbar('none') calls)
        ax = findall(f,'Type','axes');
        for a = ax'
            if isprop(a,'Toolbar') && ~isempty(a.Toolbar)
                a.Toolbar.Visible = 'off';
            end
        end
        fn = sprintf('%s_%02d.png', prefix, f.Number);
        exportgraphics(f, fullfile(outdir, fn), 'Resolution', 300);
    end
end