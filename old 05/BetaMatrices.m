classdef BetaMatrices < handle
    % BETAMATRICES Container for coagulation kernel beta matrices
    %
    % b25 is a convenience matrix used by the RHS:
    %     b25 = b2 + b5   (total loss rate from section i)
    %
    % Make sure BetaAssembler always constructs b25 with this convention.

    properties
        b1;   % Gain in section i by collisions of (i-1) & j < i
        b2;   % Gain in section i by collisions of i & j < i
        b3;   % Loss from section i by collisions of i & j < i
        b4;   % Loss from section i by collisions with itself
        b5;   % Loss from section i by collisions of i & j > i
        b25;  % Convenience sum: b2 + b5 (used in CoagulationRHS)
    end

    methods
        function obj = BetaMatrices()
            % Constructor — initialize empty matrices
            obj.b1  = [];
            obj.b2  = [];
            obj.b3  = [];
            obj.b4  = [];
            obj.b5  = [];
            obj.b25 = [];
        end

        function validate(obj)
            % Validate that required matrices exist and share dimensions
            req = {'b1','b2','b3','b4','b5'};
            for k = 1:numel(req)
                if isempty(obj.(req{k}))
                    error('BetaMatrices:%sEmpty','%s is empty', req{k}); %#ok<SPERR>
                end
            end
            sz = size(obj.b1);
            for k = 2:numel(req)
                if ~isequal(size(obj.(req{k})), sz)
                    error('BetaMatrices:SizeMismatch','%s size mismatch', req{k});
                end
            end
            if ~isempty(obj.b25) && ~isequal(size(obj.b25), sz)
                error('BetaMatrices:SizeMismatch','b25 size mismatch');
            end
        end

        function n = getNumSections(obj)
            % Number of sections
            if isempty(obj.b1)
                n = 0;
            else
                n = size(obj.b1,1);
            end
        end

        function tf = isSymmetric(obj)
            % Basic symmetry check (self-collision block)
            tf = ~isempty(obj.b4) && issymmetric(obj.b4);
        end

        function displaySummary(obj)
            % Print a compact summary
            n = obj.getNumSections();
            fprintf('BetaMatrices: %d sections | fields: b1,b2,b3,b4,b5', n);
            if ~isempty(obj.b25), fprintf(',b25'); end
            fprintf('\n');
            if ~isempty(obj.b1)
                fprintf('  dims: %dx%d\n', size(obj.b1,1), size(obj.b1,2));
                fprintf('  b1 range:  [%.2e, %.2e]\n', min(obj.b1(:)), max(obj.b1(:)));
                fprintf('  b2 range:  [%.2e, %.2e]\n', min(obj.b2(:)), max(obj.b2(:)));
                fprintf('  b3 range:  [%.2e, %.2e]\n', min(obj.b3(:)), max(obj.b3(:)));
                fprintf('  b4 range:  [%.2e, %.2e]\n', min(obj.b4(:)), max(obj.b4(:)));
                fprintf('  b5 range:  [%.2e, %.2e]\n', min(obj.b5(:)), max(obj.b5(:)));
                if ~isempty(obj.b25)
                    fprintf('  b25 range: [%.2e, %.2e]\n', min(obj.b25(:)), max(obj.b25(:)));
                end
                fprintf('  symmetric b4: %s\n', string(obj.isSymmetric()));
            end
        end
    end
end