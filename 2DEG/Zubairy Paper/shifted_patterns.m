clear all;
load shifted_plots;




close all;clf;figure(1)
% 
N = 3;
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = plot(xaxis, plot1./max(plot1), 'linewidth',1.4);
hold on
% h1 = plot(abs(imag(root)), f*1e-12, 'linewidth',1.4);
h2 = plot(xaxis, plot4./max(plot4), 'linewidth',1.4, 'linestyle','-');
h3 = plot(xaxis, plot7./max(plot7), 'linewidth',1.4, 'linestyle','-');
set(gcf,'Color','white');
set(gca,'FontName','times new roman','FontSize',15);
set(gca,'FontName','times new roman','FontSize',15,'YScale', 'lin','XScale', 'lin') % Set axes fonts to Times New Roman
%
xlabel('$x~(\mathrm{\u m})$','interpreter','latex')
ylabel('$|E/E_{max}| $','interpreter','latex')
legend([h1 h2 h3], {'$\chi = 0^\circ$', '$\chi = 72^\circ$',  '$\chi = 126^\circ$' },...
    'location','northeast','Orientation','horizontal',...
    'interpreter','latex','FontSize',15);
set(gca,'FontName','times new roman','FontSize',15) % Set axes fonts to Times New Roman
% grid on
box on
xlim([ 0 1.94])
ylim([ 0 1.3])
cleanfigure();
matlab2tikz('filename',sprintf('shifted.tex'));
