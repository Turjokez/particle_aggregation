classdef SectionGrid < handle
    %SECTIONGRID Minimal grid wrapper so Disaggregation sees real methods

    properties
        radii_cm
        config
    end

    methods
        function obj = SectionGrid(radii_cm, cfg)
            obj.radii_cm = radii_cm(:);
            obj.config   = cfg;
        end

        function r = getConservedRadii(obj)
            r = obj.radii_cm;
        end

        function r = getFractalRadii(obj)
            r = obj.radii_cm;
        end
    end
end