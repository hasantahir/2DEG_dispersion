load standing_Wave.mat x y
close all;clf;figure(1)
% 
N = 2;
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = plot(x,y, 'linewidth',1.4);
set(gcf,'Color','white');
set(gca,'FontName','times new roman','FontSize',15);
set(gca,'FontName','times new roman','FontSize',15,'YScale', 'lin','XScale', 'lin') % Set axes fonts to Times New Roman
% xlim([1 10])
% ylim([-35 55])
%
xlabel('$x (\mathrm{\mu m})$','interpreter','latex')
ylabel('$\vert E \vert $','interpreter','latex')
set(gca,'FontName','times new roman','FontSize',15) % Set axes fonts to Times New Roman
% grid on
box on
% % cleanfigure();
matlab2tikz('filename',sprintf('standing.tex'));