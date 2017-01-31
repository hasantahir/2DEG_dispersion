function plotDielectricCylinderTotalFieldUnderPlaneWave(varargin)
% Plot the total field as a plane wave is scattered by a dielectric
% cylinder. The plane wave propagates from left to right, i.e. in the
% positive x direction.
%
% 
%                     ^ y
%                     |
%                     |____>  x
%
% Optional:
%
% 'radius         radius of the dielectric cylinder
%                 (relative to the background
%                 wavelength). Default to 1.
%
% 'cylinder'      object of DielectricMaterial. Default to
%                 DielectricMaterial(4.0, 0.0).
%
% 'background'    object of DielectricMaterial. Default to
%                 DielectricMaterial(1.0, 0.0).
%
% 'frequency'     scalar to denote the frequency of the plane
%                 wave. Default to 1 Hz.
%
% 'domain_size'   4x1 vector to specify the size of the domain in
%                 the format  [x_min; x_max; y_min; y_max]
%                 (relative to the background wavelength). Default
%                 to [-3; 3; -2; 2].
%
% 'resolution'    2x1 vector to specify the resolution of the domain in
%                 the format [x_res; y_res] (relative to the
%                 background wavelength). Default to [0.05; 0.05].
%
% 'figure_name'   figure name. Default to [].
%
% 'workspace_name'  workspace name to store the
%                   data. Default to [].
%
    
% $Author:: kzhu                                          $
% $Rev:: 1617                                             $
% $Date:: 2011-03-27 13:55:53 -0400 (Sun, 27 Mar 2011)    $

    p = inputParser;
    p.addParamValue('radius', 1.0, @isnumeric);
    p.addParamValue('cylinder', DielectricMaterial(4.0, 0.0), @isobject);
    p.addParamValue('background', DielectricMaterial(1.0, 0.0), @isobject);
    p.addParamValue('frequency', 1.0, @isnumeric);
    p.addParamValue('domain_size', [-3; +3; -2; +2], @isvector);
    p.addParamValue('resolution', [0.05; 0.05], @isvector);
    p.addParamValue('dB', 0, @(x)x==0 || x==1);
    p.addParamValue('figure_name', '', @ischar);
    p.addParamValue('workspace_name', '', @ischar);
    p.addParamValue('verbose', 0, @(x)x==0||x==1);
    p.parse(varargin{:});
    if (p.Results.verbose)
        disp(p.Results);
    end
    
    wavelength = getElectromagneticWavelength(p.Results.background, ...
                                              p.Results.frequency);
    x          = p.Results.domain_size(1):p.Results.resolution(1):p.Results.domain_size(2);
    y          = p.Results.domain_size(3):p.Results.resolution(2):p.Results.domain_size(4);
    [x, y]     = ndgrid(x, y);
    
    Ez         = zeros(size(x));
    H_rho      = zeros(size(x));
    H_phi      = zeros(size(x));
    for iX = 1:size(x, 1);
        for iY = 1:size(y, 2);
            [Ez(iX,iY), H_rho(iX,iY), H_phi(iX,iY)] = getDielectricCylinderFieldUnderPlaneWave(p.Results.radius*wavelength, ...
                                                              p.Results.cylinder, ...
                                                              p.Results.background, ...
                                                              [x(iX,iY); y(iX, iY)]*wavelength, ...
                                                              p.Results.frequency);
        end
    end

    phi = cart2pol(x, y);
    Hx  = cos(phi).*H_rho - sin(phi).*H_phi;
    Hy  = sin(phi).*H_rho + cos(phi).*H_phi;

    save;
    for iX = 1:size(x,1)
        for iY = 1:size(y,2)
            sensor_location = [0; -y(iX, iY); x(iX,iY)];
            if (norm(sensor_location) > p.Results.radius) % outside the cylinder
                [Ez_inc, Hy_inc] = getPlaneWaveUsingCylindricalExpansion(p.Results.background, ...
                                                                  sensor_location*wavelength, ... 
                                                                  p.Results.frequency);
                % getPlaneWaveUsingCylindricalExpansion assumes +z
                % progating wave, but we have +x propagating wave here.
                Hy_inc           = -Hy_inc; 

                Ez_tot(iX,iY)    = Ez(iX,iY)+Ez_inc;
                Hy_tot(iX,iY)    = Hy(iX,iY)+Hy_inc;
            else % inside
                Ez_tot(iX,iY)    = Ez(iX,iY);
                Hy_tot(iX,iY)    = Hy(iX,iY);
            end
            
            % plane wave has not Hx component.
            Hx_tot(iX,iY)    = Hx(iX,iY);
        end%iY
    end%iX
    
    E        = zeros([size(Ez_tot) 3]);
    E(:,:,3) = Ez_tot;
    H        = zeros(size(E));
    H(:,:,1) = Hx_tot;
    H(:,:,2) = Hy_tot;
    
    S        = cross(E, conj(H)); % poynting vector
    clear Ez H_rho H_phi E H; % save some space

    save;
    
    set(figure, 'color', 'white');
    if (p.Results.dB)
        subplot(221); imagesc(x(:,1), y(1,:), 20*log10(abs(Ez_tot))'); 
        plotCylinder(p.Results.radius); title('|E_z| (V/m)');
        quiver(x(1:4:end,1:4:end), y(1:4:end,1:4:end), ...
               real(S(1:4:end,1:4:end,1)), real(S(1:4:end,1:4:end,2)), ...
               'k-','linewidth', 1); 
        
        
        subplot(223); imagesc(x(:,1), y(1,:), 20*log10(abs(Hx_tot))'); 
        plotCylinder(p.Results.radius); title('|H_x| (A/m)');
        subplot(224); imagesc(x(:,1), y(1,:), 20*log10(abs(Hy_tot))'); 
        plotCylinder(p.Results.radius); title('|H_y| (A/m)');
    else
        subplot(221); imagesc(x(:,1), y(1,:), abs(Ez_tot)'); 
        plotCylinder(p.Results.radius); title('|E_z| (V/m)');
        quiver(x(1:4:end,1:4:end), y(1:4:end,1:4:end), ...
               real(S(1:4:end,1:4:end,1)), real(S(1:4:end,1:4:end,2)), ...
               'k-','linewidth', 1); 
        
        subplot(223); imagesc(x(:,1), y(1,:), abs(Hx_tot)'); 
        plotCylinder(p.Results.radius); title('|H_x| (A/m)');
        subplot(224); imagesc(x(:,1), y(1,:), abs(Hy_tot)'); 
        plotCylinder(p.Results.radius); title('|H_y| (A/m)');
    end
    subplot(222); imagesc(x(:,1), y(1,:), angle(Ez_tot)'); 
    plotCylinder(p.Results.radius); title('Phase (rad)')

    if (~isempty(p.Results.figure_name))
        saveas(gcf, p.Results.figure_name, 'fig');
    end

    if (~isempty(p.Results.workspace_name))
        save(p.Results.workspace_name, 'x', 'y', 'Ez_tot')
    end

end


