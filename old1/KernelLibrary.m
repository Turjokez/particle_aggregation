classdef KernelLibrary < handle
    %KERNELLIBRARY Collection of coagulation kernel implementations
    
    methods (Static)
        function b = brownian(r, rcons, param)
            %BROWNIAN Brownian motion kernel
            % r = 2xN vector of particle radii [cm]
            % rcons and param are dummy variables for efficiency
            
            b = (2 + r(1,:)./r(2,:) + r(2,:)./r(1,:));
        end
        
        function b = curvilinearDS(r, rcons, param)
            %CURVILINEARDS Curvilinear differential sedimentation kernel
            % r = column vector of particle radii [cm]
            % rcons = column vector of conserved volume radii [cm]
            % param = parameters including r_to_rg and setcon
            
            r_small = min(r) * param.r_to_rg;
            rcons_3 = rcons .* rcons .* rcons;
            
            b = 0.5 * pi * abs(rcons_3(1,:)./r(1,:) - rcons_3(2,:)./r(2,:)) .* r_small .* r_small;
        end
        
        function b = curvilinearShear(r, rcons, param)
            %CURVILINEARSHEAR Curvilinear shear kernel
            % r = column vector of particle radii [cm]
            % param = parameters including r_to_rg
            
            rg = (r(1,:) + r(2,:)) * param.r_to_rg;
            
            p = min(r) ./ max(r);
            p1 = 1.0 + p;
            p5 = p1 .* p1 .* p1 .* p1 .* p1;
            
            efficiency = 1.0 - (1.0 + 5.0*p + 2.5*p.*p) ./ p5;
            
            b = sqrt(8.0*pi/15.0) * efficiency .* rg .* rg .* rg;
        end
        
        function b = fractalDS(r, rcons, param)
            %FRACTALDS Fractal differential sedimentation kernel
            % r = column vector of particle radii [cm]
            % rcons = column vector of conserved volume radii [cm]
            % param = parameters including r_to_rg
            
            c1 = 0.984;  % Constant from Li and Logan
            
            rg = (r(1,:) + r(2,:)) * param.r_to_rg;
            r_ratio = min(r(1,:)./r(2,:), r(2,:)./r(1,:));
            rcons_3 = rcons .* rcons .* rcons;
            
            b = pi * abs(rcons_3(1,:)./r(1,:) - rcons_3(2,:)./r(2,:)) .* rg .* rg;
            b = b .* r_ratio.^c1;
        end
        
        function b = fractalShear(r, rcons, param)
            %FRACTALSHEAR Fractal shear kernel
            % r = column vector of particle radii [cm]
            % param = parameters including r_to_rg
            
            c1 = 0.785;  % Constant from Li and Logan
            r_ratio = min(r(1,:)./r(2,:), r(2,:)./r(1,:));
            rg = (r(1,:) + r(2,:)) * param.r_to_rg;
            
            b = 1.3 * rg .* rg .* rg;
            b = b .* r_ratio.^c1;
        end
        
        function b = rectilinearDS(r, rcons, param)
            %RECTILINEARDS Rectilinear differential sedimentation kernel
            % r = column vector of particle radii [cm]
            % rcons = column vector of conserved volume radii [cm]
            % param = parameters including r_to_rg
            
            rg = (r(1,:) + r(2,:)) * param.r_to_rg;
            rcons_3 = rcons .* rcons .* rcons;
            
            b = pi * abs(rcons_3(1,:)./r(1,:) - rcons_3(2,:)./r(2,:)) .* rg .* rg;
        end
        
        function b = rectilinearShear(r, rcons, param)
            %RECTILINEARSHEAR Rectilinear shear kernel
            % r = column vector of particle radii [cm]
            % param = parameters including r_to_rg
            
            rg = (r(1,:) + r(2,:)) * param.r_to_rg;
            b = 1.3 * rg .* rg .* rg;
        end
        
        function kernel = getKernel(kernelName)
            %GETKERNEL Get kernel function handle by name
            switch kernelName
                case 'KernelBrown'
                    kernel = @KernelLibrary.brownian;
                case 'KernelCurDS'
                    kernel = @KernelLibrary.curvilinearDS;
                case 'KernelCurSh'
                    kernel = @KernelLibrary.curvilinearShear;
                case 'KernelFracDS'
                    kernel = @KernelLibrary.fractalDS;
                case 'KernelFracSh'
                    kernel = @KernelLibrary.fractalShear;
                case 'KernelRectDS'
                    kernel = @KernelLibrary.rectilinearDS;
                case 'KernelRectSh'
                    kernel = @KernelLibrary.rectilinearShear;
                otherwise
                    error('Unknown kernel: %s', kernelName);
            end
        end
    end
end
