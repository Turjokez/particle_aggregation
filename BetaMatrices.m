classdef BetaMatrices < handle
    % BETAMATRICES Container for coagulation kernel beta matrices
    
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
            % BETAMATRICES Constructor - initialize empty matrices
            obj.b1 = [];  
            obj.b2 = [];  
            obj.b3 = [];  
            obj.b4 = [];  
            obj.b5 = [];  
            obj.b25 = []; 
        end
        
        function validate(obj)
            % VALIDATE Validate that all matrices have consistent dimensions
            matrices = {obj.b1, obj.b2, obj.b3, obj.b4, obj.b5, obj.b25};  % List all matrices
            
            % Check that all matrices exist and have same size
            if isempty(obj.b1)  % If b1 is empty, throw an error
                error('Beta matrices not initialized');
            end
            
            [rows, cols] = size(obj.b1);  % Get the size of b1
            for i = 2:length(matrices)  % Loop through all matrices
                if ~isequal(size(matrices{i}), [rows, cols])  % Check if matrix dimensions match
                    error('Beta matrix %d has inconsistent dimensions', i);  % Throw error if dimensions mismatch
                end
            end
        end
        
        function n = getNumSections(obj)
            % GETNUMSECTIONS Get number of sections from matrix dimensions
            if isempty(obj.b1)  % If b1 is empty, return 0
                n = 0;
            else
                n = size(obj.b1, 1);  % Return the number of rows (sections) in b1
            end
        end
        
        function result = isSymmetric(obj)
            % ISSYMMETRIC Check if beta matrices represent symmetric coagulation
            % This is a basic check - more sophisticated validation may be needed
            if isempty(obj.b1)  % If b1 is empty, return false
                result = false;
                return;
            end
            
            % Check if b4 is symmetric (self-coagulation)
            result = issymmetric(obj.b4);  % Check symmetry of b4 matrix
        end
        
        function displaySummary(obj)
            % DISPLAYSUMMARY Display summary of beta matrices
            fprintf('Beta Matrices Summary:\n');
            fprintf('  Number of sections: %d\n', obj.getNumSections());  % Print number of sections
            fprintf('  Matrix dimensions: %dx%d\n', size(obj.b1));  % Print matrix dimensions
            fprintf('  Symmetric: %s\n', mat2str(obj.isSymmetric()));  % Print if matrices are symmetric
            
            % Display some statistics
            if ~isempty(obj.b1)  % If b1 is not empty, display range of values
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
