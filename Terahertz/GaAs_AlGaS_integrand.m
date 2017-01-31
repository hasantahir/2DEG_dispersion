clear all
close all;tic

lambda = 1;%2.99792458e-4; % 1 THz
% (refractiveindex.info)
eps1 =  2.4 - 1i*0.1; % Upper layer - GaAs
eps2 = 2.3 - 1i*0.01 ; % Lower Layer - AlGaAs
% eps_silver = -265.06 - 1i*29.436;  % @ 2500 nm
load em_constants.mat % Contains varepsilon, mu and c
eps_0 = epsilon_0;
c = 1/sqrt(mu_0*eps_0);
omega = 2*pi*c/lambda; % angular frequency
eta_0 = sqrt(mu_0/eps_0);
k1 = omega*sqrt(mu_0*eps_0); % propagation constant of air
k2 = omega * sqrt(mu_0*eps_0*eps1); % propagation constant of silver
%%
len = 5e2; % Vector Length
%%
kxx = horzcat(linspace(0*k1,.5*k1,len/4),...
    linspace(.5*k1 + 1e-6, k1,len/4),...
    linspace(k1 + 1e-6, 1.5*k1,len/4),...
    linspace(1.5*k1 + 1e-6, 1e3*k1,len/4)); % 1e-6 used to avoid duplicates in the array

kxy = horzcat(linspace(-1e2*k1,-.5*k1,len/4),...
    linspace(-.5*k1 + 1e-6, -.1*k1,len/4),...
    linspace(-.1*k1 + 1e-6, 0,len/4),...
    linspace(0 + 1e-6, 1e1*k1,len/4));

% Find the branch point location on the kx x-axis
dif = abs(kxx - k1);
bp1_loc = find(dif == min(dif)); % Index in kxx with the nearest value of k_air

% x = horzcat(linspace(1e-2*lambda,1e0*lambda,len/2),...
%     linspace(1e0*lambda,1e4*lambda,len/2));
x = 1e0;
H_c1 = zeros(length(x),1); % Initialize the magnetic field vector
H_c2 = zeros(length(x),1); % Initialize the magnetic field vector
%% Define Contour

%               | k1 |
%               |    |
%             C1|    |c2
%               |    |
%       ------> |    V <--------
%       bottom  |    | top sheet
%               ^    |
% Define real and imaginary start points for c1
%
% Path location in terms of re_kx in kx plane
c1_start_real = kxx(bp1_loc);
c1_start_imag = -kxx(bp1_loc)*1e1;
diff = abs(kxy - c1_start_imag);
c1_sty_loc = find(diff == min(diff)); % Starting point on the Im_kx axis

%
% Define real and imaginary end points for c1
%
c1_end_real = kxx(bp1_loc);
diff = abs(kxy - 0);
c1_eny_loc = find(diff == min(diff)); % Ending point on the Im_kx axis
% c1_end_imag = kxy(c1_eny_loc) - 1e-6;
c1_end_imag = 0;
%
% Make C1 Contour
c1_real = c1_start_real*ones(len*1e3,1);
c1_imag = linspace(c1_start_imag, c1_end_imag,len*1e3);
%
c1 = horzcat(c1_real, c1_imag'); % concatenate real and imaginary parts
% kx_c1 = c1_real(:,1) + 1i*c1_imag(1,:)';
kx_c1 = transpose(-.6:1e-3:0);
kx_c1 = k1 + 1i*kx_c1;
% Define real and imaginary start points for c2
%
c2_start_real = kxx(bp1_loc);
c2_start_imag = 0;
%
% Define real and imaginary end points for c2
%
c2_end_real = kxx(bp1_loc);
c2_end_imag = -kxx(bp1_loc)*1e1;
%
% Make c2 Contour
c2_real = c2_start_real*ones(len*1e3,1);

c2_imag = linspace(c2_start_imag, c2_end_imag,len*1e3);
c2 = horzcat(c2_real, c2_imag');
% kx_c2 = c2_real(:,1) + 1i*c2_imag(1,:)';
% kx_c2 = transpose(1:1e-4:10);
kx_c2 = transpose(0:-1e-3:-.6);
kx_c2 = k1 + 1i*kx_c2;
%% Define Green's function


kz_1 = @(kx) sqrt(k1^2 - kx.^2);
kz_2 = @(kx) sqrt(k2^2 - kx.^2);
D = @(kz_1, kz_2) kz_2/eps1 + kz_1/eps2;
% G = @(kz_1, kz_2) eps_0./D;

% On Contour C1
kz1_c1 = kz_1(kx_c1);
kz2_c1 = kz_2(kx_c1);

% Define D and G
% D_c1 = D(kz1_c1, kz2_c1);
% G_c1 = 1./D_c1;

% On Contour c2
kz1_c2 = kz_1(kx_c2);
kz2_c2 = kz_2(kx_c2);

% Define D and G
% D_c2 = D(kz1_c2, kz2_c2);
% G_c2 = 1./D_c2;

%%
% Branch Cut Curve
hyp_silver = imag(k2^2)./(2*kxx); % Hyperbolic cruve for silver
hyp_silver_re = +kxx.^2 - real(k2^2); % Hyperbolic cruve for silver
% Intersection of Branch cut with vertical cut
y_int1 = imag(k1^2)/(2*k2);
y_int2 = imag(k2^2)/(2*k1);


%% Integrate
% C1 lies totally on the bottom sheet of k1
% C1 partially lies on the bottom sheet of k2 ( until y_int)
%    and partially on the top sheet

% Integrate on left edge

for j = 1 : length(kx_c1)
    
    % Enfore bottom sheet on k2
    if abs(imag(kx_c1(j))) > abs(y_int1)
        
        
        % Enfore bottom sheet on k1
        % Im(kz1) > 0, Re(kz1) > 0
        if real(kz1_c1(j)) < 0
            kz1_c1(j) = -conj(kz1_c1(j));
        end
        
        if imag(kz1_c1(j)) < 0
            kz1_c1(j) = conj(kz1_c1(j));
        end
    else
        % Top sheet of kz1
        % Im(kz1) < 0, Re(kz1) > 0
        
        if real(kz1_c1(j)) < 0
            kz1_c1(j) = -conj(kz1_c1(j));
        end
        
        if imag(kz1_c1(j)) > 0
            kz1_c1(j) = conj(kz1_c1(j));
        end
        
    end
    
    if abs(imag(kx_c1(j))) > abs(y_int1)
        
        
        % Enfore bottom sheet on k2
        % Im(kz2) > 0, Re(kz2) > 0
        if real(kz2_c1(j)) < 0
            kz2_c1(j) = -conj(kz2_c1(j));
        end
        
        if imag(kz2_c1(j)) < 0
            kz2_c1(j) = conj(kz2_c1(j));
        end
    else
        % Top sheet of kz2
        % Im(kz2) < 0, Re(kz2) > 0
        
        if real(kz2_c1(j)) < 0
            kz2_c1(j) = -conj(kz2_c1(j));
        end
        
        if imag(kz2_c1(j)) > 0
            kz2_c1(j) = conj(kz2_c1(j));
        end
        
    end
end

D_c1 = D(kz1_c1, kz2_c1);
G_c1 = 1./D_c1;

% Right Edge
for j = 1 : length(kx_c2)
    
    % Enfore top sheet on k1
    if abs(imag(kx_c2(j))) < abs(y_int1)
        
        % Enfore top sheet on k1
        % Im(kz1) < 0, Re(kz1) > 0
        if real(kz1_c2(j)) < 0
            kz1_c2(j) = -conj(kz1_c2(j));
        end
        
        if imag(kz1_c2(j)) > 0
            kz1_c2(j) = conj(kz1_c2(j));
        end
    else
        % Bottom sheet of kz1
        % Im(kz1) > 0, Re(kz1) > 0
        
        if real(kz1_c2(j)) < 0
            kz1_c2(j) = -conj(kz1_c2(j));
        end
        
        if imag(kz1_c2(j)) < 0
            kz1_c2(j) = conj(kz1_c2(j));
        end
        
    end
    
    if abs(imag(kx_c2(j))) < abs(y_int2)
        
        % Enfore top sheet on k2
        % Im(kz2) < 0, Re(kz2) > 0
        if real(kz2_c2(j)) < 0
            kz2_c2(j) = -conj(kz2_c2(j));
        end
        
        if imag(kz2_c2(j)) > 0
            kz2_c2(j) = conj(kz2_c2(j));
        end
    else
        % Bottom sheet of kz2
        % Im(kz2) > 0, Re(kz2) > 0
        
        if real(kz2_c2(j)) < 0
            kz2_c2(j) = -conj(kz2_c2(j));
        end
        
        if imag(kz2_c2(j)) < 0
            kz2_c2(j) = conj(kz2_c2(j));
        end
        
    end
end

D_c2 = D(kz1_c2, kz2_c2);
G_c2 = 1./D_c2;

%% Plot the integrands
%
figure(1)
integrand1 = G_c1.*exp(-1i*kx_c1*1e1);
plot(imag(kx_c1), real(integrand1),'LineWidth',1.4)
hold on
plot(imag(kx_c1), imag(integrand1),'LineWidth',1.4)
set(gcf,'Color','white');

% Create ylabel
ylabel('Integrand along $C_1$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create xlabel
xlabel('$\Im k_x$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create Legend
legend({'Real Part', 'Imaginary Part'},...
    'FontSize',10,...
    'Interpreter','latex',...
    'Location','SouthWest');

% Create Title
title('Integrand along the left edge of vertical integration path',...
    'FontSize',12,...
    'Interpreter','latex');

% Create Annotation box
str = '$\frac{x}{\lambda_0} = 1\times 10^1$';
dim = [.2 .2 .2 .2];
annotation('textbox',dim,'String',str,...
    'FitBoxToText','on',...
    'FontSize',10,...
    'Interpreter','latex');

grid on
% Save as a tikz file
% cleanfigure();
% matlab2tikz('filename',sprintf('figures/source/integrand_vertical_left.tex'),'showInfo', false)

figure(2)
integrand2 = G_c2.*exp(-1i*kx_c2*1e1);
plot(imag(kx_c2), real(integrand2),'LineWidth',1.4)
hold on
plot(imag(kx_c2), imag(integrand2),'LineWidth',1.4)
set(gcf,'Color','white');

% Create ylabel
ylabel('Integrand along $C_2$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create xlabel
xlabel('$\Im k_x$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create Legend
legend({'Real Part', 'Imaginary Part'},...
    'FontSize',10,...
    'Interpreter','latex',...
    'Location','SouthWest');

% Create Title
title('Integrand along the left edge of vertical integration path',...
    'FontSize',12,...
    'Interpreter','latex');

% Create Annotation box
str = '$\frac{x}{\lambda_0} = 1\times 10^1$';
dim = [.2 .2 .2 .2];
annotation('textbox',dim,'String',str,...
    'FitBoxToText','on',...
    'FontSize',10,...
    'Interpreter','latex');

grid on
% Save as a tikz file
% cleanfigure();
% matlab2tikz('filename',sprintf('figures/source/integrand_vertical_right.tex'),'showInfo', false)

% figure(3)
% hold on
% plot(kxx,hyp_silver,'LineWidth',1.4)
% plot(kxx,hyp_silver_re,'LineWidth',1.4)
% set(gcf,'Color','white');
%
% ylabel('Integrand along C2',...
%     'HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'FontSize',12,...
%     'Interpreter','latex');
%
% % Create xlabel
% xlabel('$\Im k_x$',...
%     'HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'FontSize',12,...
%     'Interpreter','latex');
% grid on; box on
% legend({'Real Part', 'Imaginary Part'},...
%     'FontSize',10,...
%     'Interpreter','latex',...
%     'Location','SouthWest');
% xlim([0 10])
% ylim([-50 50])
%%
toc