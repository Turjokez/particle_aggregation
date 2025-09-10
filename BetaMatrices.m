classdef BetaMatrices < handle
    %BETAMATRICES Container for coagulation kernel beta matrices
    
    properties
        b1;     % Gain in section i by collisions of (i-1) & j < i
        b2;     % Gain in section i by collisions of i & j < i
        b3;     % Loss from section i by collisions of i & j < i
        b4;     % Loss from section i by collisions with itself
        b5;     % Loss from section i by collisions of i & j > i
        b25;    % Net coagulation effect: b2 - b3 - b4 - b5
    end
    
    methods
        function obj = BetaMatrices()
            %BETAMATRICES Constructor - initialize empty matrices
            obj.b1 = [];
            obj.b2 = [];
            obj.b3 = [];
            obj.b4 = [];
            obj.b5 = [];
            obj.b25 = [];
        end
        
        function validate(obj)
            %VALIDATE Validate that all matrices have consistent dimensions
            matrices = {obj.b1, obj.b2, obj.b3, obj.b4, obj.b5, obj.b25};
            
            % Check that all matrices exist and have same size
            if isempty(obj.b1)
                error('Beta matrices not initialized');
            end
            
            [rows, cols] = size(obj.b1);
            for i = 2:length(matrices)
                if ~isequal(size(matrices{i}), [rows, cols])
                    error('Beta matrix %d has inconsistent dimensions', i);
                end
            end
        end
        
        function n = getNumSections(obj)
            %GETNUMSECTIONS Get number of sections from matrix dimensions
            if isempty(obj.b1)
                n = 0;
            else
                n = size(obj.b1, 1);
            end
        end
        
        function result = isSymmetric(obj)
            %ISSYMMETRIC Check if beta matrices represent symmetric coagulation
            % This is a basic check - more sophisticated validation may be needed
            if isempty(obj.b1)
                result = false;
                return;
            end
            
            % Check if b4 is symmetric (self-coagulation)
            result = issymmetric(obj.b4);
        end
        
        function displaySummary(obj)
            %DISPLAYSUMMARY Display summary of beta matrices
            fprintf('Beta Matrices Summary:\n');
            fprintf('  Number of sections: %d\n', obj.getNumSections());
            fprintf('  Matrix dimensions: %dx%d\n', size(obj.b1));
            fprintf('  Symmetric: %s\n', mat2str(obj.isSymmetric()));
            
            % Display some statistics
            if ~isempty(obj.b1)
                fprintf('  b1 range: [%.2e, %.2e]\n', min(obj.b1(:)), max(obj.b1(:)));
                fprintf('  b2 range: [%.2e, %.2e]\n', min(obj.b2(:)), max(obj.b2(:)));
                fprintf('  b3 range: [%.2e, %.2e]\n', min(obj.b3(:)), max(obj.b3(:)));
                fprintf('  b4 range: [%.2e, %.2e]\n', min(obj.b4(:)), max(obj.b4(:)));
                fprintf('  b5 range: [%.2e, %.2e]\n', min(obj.b5(:)), max(obj.b5(:)));
                fprintf('  b25 range: [%.2e, %.2e]\n', min(obj.b25(:)), max(obj.b25(:)));
            end
        end
    end
end
