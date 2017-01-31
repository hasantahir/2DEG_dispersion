function plotConductingCylinderTotalFieldUnderLineSource(varargin)
% Plot the total field, as an infinitely-long line
% source radiates a dielectric cylinder.
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
% 'background'    object of DielectricMaterial. Default to
%                 DielectricMaterial(1.0, 0.0).
%
% 'source_location'  2x1 vector to specify the location of the
%                    current source in the format [x; y] (relative to the
%                    wavelength). Default to [-2.5; -1.5].
%
% 'frequency'     scalar to denote the frequency of the radiated field.
%                 Default to 1 Hz.
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
% $Rev:: 1426                                             $
% $Date:: 2011-01-23 21:57:19 -0500 (Sun, 23 Jan 2011)    $

    p = inputParser;
    p.addParamValue('radius', 1.0, @isnumeric);
    p.addParamValue('background', DielectricMaterial(1.0, 0.0), @isobject);
    p.addParamValue('source_location', [-2.5; -2.5], @isvector);
    p.addParamValue('frequency', 1.0, @isnumeric);
    p.addParamValue('domain_size', [-5; +3; -5; +2], @isvector);
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
            [Ez(iX,iY), H_rho(iX,iY), H_phi(iX,iY)] = getConductingCylinderFieldUnderLineSource(p.Results.radius*wavelength, ...
                                                              p.Results.background, ...
                                                              p.Results.source_location*wavelength, ...
                                                              [x(iX,iY); y(iX, iY)]*wavelength, ...
                                                              p.Results.frequency);
        end
    end
    
    for iX = 1:size(x,1)
        for iY = 1:size(y,2)
            sensor_location = [x(iX,iY); y(iX, iY)];
            if (norm(sensor_location) > p.Results.radius) % outside the conducting cylinder
                [Ez_inc, H_rho_inc, H_phi_inc] = getCylindricalWaveUsingCylindricalExpansion(p.Results.background, ...
                                                                  p.Results.source_location*wavelength, ...
                                                                  sensor_location*wavelength, ... 
                                                                  p.Results.frequency);
                Ez(iX,iY)    = Ez(iX,iY)   +Ez_inc;
                H_rho(iX,iY) = H_rho(iX,iY)+H_rho_inc;
                H_phi(iX,iY) = H_phi(iX,iY)+H_phi_inc;
            end
        end%iY
    end%iX
    
    phi     = cart2pol(x, y);
    Hx  = cos(phi).*H_rho - sin(phi).*H_phi;
    Hy  = sin(phi).*H_rho + cos(phi).*H_phi;
    
    E        = zeros([size(Ez) 3]);
    E(:,:,3) = Ez;
    H        = zeros(size(E));
    H(:,:,1) = Hx;
    H(:,:,2) = Hy;
    
    S        = cross(E, conj(H)); % poynting vector
    clear E H; % save some space

    set(figure, 'color', 'white');
    if (p.Results.dB)
        subplot(221); imagesc(x(:,1), y(1,:), 20*log10(abs(Ez))'); 
        plotCylinder(p.Results.radius); title('|E_z| (V/m)');
        quiver(x(1:8:end,1:8:end), y(1:8:end,1:8:end), ...
               real(S(1:8:end,1:8:end,1)), real(S(1:8:end,1:8:end,2)), ...
               'k-','linewidth', 1); 

        subplot(223); imagesc(x(:,1), y(1,:), 20*log10(abs(Hx))'); 
        plotCylinder(p.Results.radius); title('|H_x| (A/m)');
        subplot(224); imagesc(x(:,1), y(1,:), 20*log10(abs(Hy))'); 
        plotCylinder(p.Results.radius); title('|H_y| (A/m)');
    else
        subplot(221); imagesc(x(:,1), y(1,:), abs(Ez)'); 
        plotCylinder(p.Results.radius); title('|E_z| (V/m)');
        quiver(x(1:8:end,1:8:end), y(1:8:end,1:8:end), ...
               real(S(1:8:end,1:8:end,1)), real(S(1:8:end,1:8:end,2)), ...
               'k-','linewidth', 1); 
        
        subplot(223); imagesc(x(:,1), y(1,:), abs(Hx)'); 
        plotCylinder(p.Results.radius); title('|H_x| (A/m)');
        subplot(224); imagesc(x(:,1), y(1,:), abs(Hy)'); 
        plotCylinder(p.Results.radius); title('|H_y| (A/m)');
    end
    subplot(222); imagesc(x(:,1), y(1,:), angle(Ez)'); 
    plotCylinder(p.Results.radius); title('Phase (rad)')

    if (~isempty(p.Results.figure_name))
        saveas(gcf, p.Results.figure_name, 'fig');
    end

    if (~isempty(p.Results.workspace_name))
        save(p.Results.workspace_name, 'x', 'y', 'Ez', 'H_rho', 'H_phi');
    end

end


