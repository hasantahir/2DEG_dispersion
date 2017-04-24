% plot images
load microscopic_free_space.mat
% load microscopic_plasmonic.mat
% load microscopic_plasmonic_hi.mat

%% Surface plot
figure(4)
axis tight
h = surf(2*x,2*y,abs(ifz),'LineStyle', 'none',...
    'FaceColor', 'interp');  
axis equal
% h = imagesc(x,y,z);  
%

% plt_min = min(min(abs(plt)));
% plt_max = max(max(abs(plt)));
% light('Position',[-10 10 -10],'Style','infinite')
% Set z-axis to logscale
set(gca,'Xscale','Lin','Yscale','Lin','Zscale','Lin')

% Use Brewermap color schemes
colormap(brewermap([],'Spectral')) %RdYlGn 'RdYlBu' Spectral
% colormap(brewermap([],'RdYlGn')) %RdYlGn 'RdYlBu' Spectral
% colormap viridis
% Create some lighting for nice figures
% lighting gouraud
material metal
% shading interp
set(gca,...
    'box','on',...
    'FontName','times new roman',...
    'FontSize',15);
grid on
set(gcf, 'color','white');
view([0 90])
axis tight
axis([0 500 0 500])
% ylim([1e-2 1])
% Labeling
xlabel('$x (nm)$','interpreter','latex')
ylabel('$y (nm)$','interpreter','latex')
% caxis([abs((plt_min)) .04*abs(log(plt_max))])


% if plt == abs(f_8)
%     title('$J_0(z)$', 'interpreter', 'latex');
%     matlab2tikz('filename',sprintf('free_space_sim.tex'),...
%         'showInfo', false,'floatFormat','%.3f');
    print(gcf,'-dpng','free_space_sim','-r1200')

%     matlab2tikz('filename',sprintf('plasmonic_sim.tex'),...
%         'showInfo', false,'floatFormat','%.3f');
%     print(gcf,'-dpng','plasmonic_sim','-r1200')