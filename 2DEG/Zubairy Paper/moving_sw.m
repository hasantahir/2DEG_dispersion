% load comsol_data_standing_wave_25THZ.mat
close all;clf;figure(1)
% % 
% N = 5;
% axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
% h1 = plot(x,y1, 'linewidth',1.4);
% hold on
% % h2 = plot(x,y2, 'linewidth',1.4);
% h3 = plot(x,y3, 'linewidth',1.4);
% % h4 = plot(x,y4, 'linewidth',1.4);
% h5 = plot(x,y5, 'linewidth',1.4);
% 
% set(gcf,'Color','white');
% set(gca,'FontName','times new roman','FontSize',15);
% set(gca,'FontName','times new roman','FontSize',15,'YScale', 'lin','XScale', 'lin') % Set axes fonts to Times New Roman
% % legend([h1 h2 h3 h4 h5], {'$20 ^\circ$', '$40^\circ$','$60^\circ$','$80^\circ$','$120^\circ$'},...
% %     'location','northeast','interpreter','latex','FontSize',15);
% legend([h1 h3 h5], {'$20 ^\circ$', '$60^\circ$','$120^\circ$'},...
%     'location','northeast','interpreter','latex','FontSize',15);
% % legend([h1 h2 ], {'bla', 'blu'},...
% %     'location','northeast','interpreter','latex','FontSize',15);
% % xlim([1 10])
% % ylim([-35 55])
% %
% xlabel('$x (\mathrm{\mu m})$','interpreter','latex')
% ylabel('$\vert E \vert $','interpreter','latex')
% set(gca,'FontName','times new roman','FontSize',15) % Set axes fonts to Times New Roman
% % grid on
% box on
% % % cleanfigure();
% matlab2tikz('filename',sprintf('standing_shifted.tex'));


%%
% for no port excitation
% load no_port.mat dd
x = dd(:,1);
y = dd(:,2);
N = 1;
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = plot(x,y, 'linewidth',1.4);


set(gcf,'Color','white');
set(gca,'FontName','times new roman','FontSize',15);
set(gca,'FontName','times new roman','FontSize',15,'YScale', 'lin','XScale', 'lin') % Set axes fonts to Times New Roman

% legend([h1 h3 h5], {'$20 ^\circ$', '$60^\circ$','$120^\circ$'},...
%     'location','northeast','interpreter','latex','FontSize',15);
xlabel('$x (\mathrm{\mu m})$','interpreter','latex')
ylabel('$\vert E \vert $','interpreter','latex')
set(gca,'FontName','times new roman','FontSize',15) % Set axes fonts to Times New Roman
% grid on
box on
% % cleanfigure();
matlab2tikz('filename',sprintf('no_port.tex'));