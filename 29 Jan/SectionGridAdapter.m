classdef SectionGridAdapter < handle
    % SECTIONGRIDADAPTER Fallback grid adapter for Disaggregation.applyWithScaling
    % Only used if obj.grid is not a true DerivedGrid.

    properties
        rv
        config
    end

    methods
        function obj = SectionGridAdapter(rv, cfg)
            obj.rv = rv(:);
            obj.config = cfg;
        end

        function r = getConservedRadii(obj)
            r = obj.rv;
        end

        function r = getFractalRadii(obj)
            r = obj.rv;
        end
    end
end
