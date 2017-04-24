close all
clf
% Classical Interband Drude-Sommerfeld Model of Gold
c = 3e8;
w =  2*pi*linspace(3e14,1.5e15 ,1e4);
l = 3e8./w;
wp = 45e14;
gamma = 8.35e-16;
ll = 450e-9;
w0 = 2*pi*c/ll;
e_d = 1;
% e_g = 1 + wp^2./((w0^2 - w.^2) - 1i*gamma*w); 

% Drude Critical points model
% Vial and Laroche
% 
ep_inf = 1.1431;
wd = 1.3202e16;
gamma = 1.0805e14;
A1 = .26698;
phi1 = -1.2371;
W1 = 3.8711e15;
T1 = 4.4642e14;

A2 = 3.0834;
phi2 = -1.0968;
W2 = 4.1684e15;
T2 = 2.3555e15;

A = [A1 A2];
phi = [phi1 phi2];
W = [W1 W2];
T = [T1 T2];

e_g = ep_inf - wd^2./(w.^2 + 1i*gamma*w);
for i = 1 : length(A)
    e_g = e_g + A(i)*W(i)*(exp(1i*phi(i))./(W(i) - w - 1i*T(i)) + ...
        exp(-1i*phi(i))./(W(i) + w + 1i*T(i)));
end

figure(1)
plot(w,real(e_g));
hold on
plot(w,imag(e_g));
k0 = 2*pi./l;
k_sp = k0.*sqrt(e_g*e_d./(e_g+e_d));
k_x_real = k0.*sqrt(real(e_g)*e_d./(real(e_g)+e_d));
k_x_imag = k0.*imag(e_g)./(2*real(e_g).^2).*sqrt(real(e_g)*e_d./(real(e_g)+e_d));
beta = real(k_sp);
alpha = imag(k_sp);

figure(2)
% h = plot(beta,w/w0);
h = plot((k_x_real),w/w0);
hold on
plot((k_x_imag),w/w0);
% plot(alpha,w/w0);
plot(k0,w/wp,':k',...
                'LineWidth',1) 
set(gcf,'Color','white');
set(h,'Color','black','LineWidth',1.4)
% set(gca,'FontSize',10, 'FontName' , 'times new roman','YScale', 'log','XScale', 'log')
xlabel('Wavenumber, $k$ [arb. units]',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',10,...
    'Interpreter','latex');

ylabel('$\omega/\omega_p$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',10,...
    'Interpreter','latex');


%     title('Dispersion Diagram of Air-Gold Interface','FontSize',11)
    
%     axis([.250e-3 .5e-3 0 8])

% xlim([0 1]*1e9)
% ylim([0 5])

legend({'SPP','Light Line'}...
    ,'Location','northeast','FontWeight','bold',...
    'FontSize',10,...
    'Interpreter','latex')            
            
set(gca,'box','on','ticklength',[0.01 0.01])
% line([4.6 4.6],[1e-10 1e10],'--k','LineWidth','1')
    hold off
% export_fig dispersion_gold '-pdf' '-png' -nocrop -r300 -native -painters -transparent -q101