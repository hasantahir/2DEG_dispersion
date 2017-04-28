clear all;
load lo_hi_patterns;




close all;clf;figure(1)
% 
N = 2;
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = plot(xaxis, low./max(low), 'linewidth',1.4);
hold on
% h1 = plot(abs(imag(root)), f*1e-12, 'linewidth',1.4);
% h2 = plot(xaxis, hi./max(hi), 'linewidth',1.4, 'linestyle','-');
h2 = plot(xaxis, highest./max(highest), 'linewidth',1.2, 'linestyle','-');
set(gcf,'Color','white');
set(gca,'FontName','times new roman','FontSize',15);
set(gca,'FontName','times new roman','FontSize',15,'YScale', 'lin','XScale', 'lin') % Set axes fonts to Times New Roman
%
xlim([ 0 1.94])
ylim([ 0 1.3])
% grid on
xlabel('$x~(\mathrm{\u m})$','interpreter','latex')
ylabel('$|E/E_{max}| $','interpreter','latex')
legend([h1 h2], {'$\chi = 0^\circ$', '$\chi = 72^\circ$'},...
    'location','northeast','Orientation','horizontal',...
    'interpreter','latex','FontSize',15);
set(gca,'FontName','times new roman','FontSize',15) % Set axes fonts to Times New Roman
% grid on
box on
cleanfigure();
matlab2tikz('filename',sprintf('lo_hi.tex'));
