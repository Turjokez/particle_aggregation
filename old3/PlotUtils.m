classdef PlotUtils
    % Lightweight helpers that keep log-plots and ratios robust.

    methods (Static)
        function y = pos(y, tiny)
            % Return y with non-finite or <= tiny set to NaN (safe for log plots)
            if nargin < 2 || isempty(tiny), tiny = 1e-30; end
            y = double(y);
            y(~isfinite(y) | y <= tiny) = NaN;
        end

        function [x2,y2] = pair_pos(x, y, tiny)
            % Paired version: remove entries where y is not positive/finite
            if nargin < 3, tiny = 1e-30; end
            x = x(:); y = y(:);
            ok = isfinite(x) & isfinite(y) & (y > tiny);
            if any(ok)
                x2 = x(ok); y2 = y(ok);
            else
                x2 = NaN; y2 = NaN;
            end
        end

        function r = ratio_safe(num, den, tiny)
            % Safe ratio num./den; returns NaN where den<=tiny or either is non-finite
            if nargin < 3, tiny = 1e-30; end
            num = double(num); den = double(den);
            r = nan(size(num));
            ok = isfinite(num) & isfinite(den) & (den > tiny);
            r(ok) = num(ok) ./ den(ok);
            r(~isfinite(r) | r <= 0) = NaN;   % keep log-friendly
        end

        function Zp = nonneg(Z)
            % For linear 3-D surfaces: clamp negatives to 0 and drop non-finites
            Zp = double(Z);
            Zp(~isfinite(Zp)) = NaN;
            Zp(Zp < 0) = 0;
        end
    end
end