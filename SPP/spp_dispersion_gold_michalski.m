close all;clear all
clf
% Classical Interband Drude-Sommerfeld Model of Gold
c = 3e8;
w =  2*pi*linspace(1e14,1.5e15 ,1e4);
% w =  2*pi*linspace(5e14,6.0e14 ,1e4);
f = w/(2*pi);
l = 3e8./f;


% wp = 45e14;
% gamma = 8.35e-16;
% ll = 450e-9;
% w0 = 2*pi*c/ll;
% e_d = 4;
% e_g = 1 + wp^2./((w0^2 - w.^2) - 1i*gamma*w); 

% Drude Critical points model
% Vial and Laroche
% Au 1.1431 1.3202E16 1.0805E14 
% 0.26698 ?1.2371 3.8711E15 
% 4.4642E14 3.0834 ?1.0968 4.1684E15 2.3555E15
% 
% ep_inf = 5.1431; % Check This value
e = 1.0;
d = 9.6320e2;
s = 1i*w;

% Data Palik
load gold_JC.mat
N = n + 1i*k;
epsilon = N.^2;

p1 = -8.0129e-2;
p2 = -4.7259e-1 + 1i*2.7018;
p3 = -1.8003 +1i*2.7991;
p4 = -7.4424e-1 + 1i*4.4474;

c1 = -9.6287e2;
c2 = 7.2435e-1 - 1i*1.4189;
c3 = 9.7416 - 1i*3.4040e-2;
c4 = -7.753e-1 - 1i*4.2546e-1;


p = [p1 p2 p3 p4];
c = [c1 c2 c3 c4];


e_g = e + d./s;
for i = 1 : length(p)
    e_g = e_g + c(i)./(s - p(i));
end

figure(1)
plot(l,real(e_g));
hold on
plot(l,imag(e_g));
plot(lambda*1e-6, real(epsilon),'o');
plot(lambda*1e-6, imag(epsilon),'o')
k0 = 2*pi./l;
k_sp = k0.*sqrt(e_g*e_d./(e_g+e_d));
% k_x_real = w/c.*sqrt(real(e_g)*e_d./(real(e_g)+e_d));
% k_x_imag = w/c.*imag(e_g)./(2*(real(e_g)).^2).*(real(e_g)*e_d./(real(e_g)+e_d)).^(3/2);
beta = real(k_sp);
alpha = imag(k_sp);

figure(2)
h = plot(beta,f);
hold on
% h = plot(real(k_x_real),f);

% plot(imag(k_x_imag),f);
plot(alpha,f);
plot(sqrt(e_d)*k0,f,':k',...
                'LineWidth',1) 
set(gcf,'Color','white');
xlim([0 1]*1e8)
% ylim([0 5])
% set(h,'Color','black','LineWidth',1.4)
set(gca,'FontSize',10, 'FontName' , 'times new roman','YScale', 'log','XScale', 'log')
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



legend({'SPP','Light Line'}...
    ,'Location','northeast','FontWeight','bold',...
    'FontSize',10,...
    'Interpreter','latex')            
            
set(gca,'box','on','ticklength',[0.01 0.01])
% line([4.6 4.6],[1e-10 1e10],'--k','LineWidth','1')
    hold off
% export_fig dispersion_gold '-pdf' '-png' -nocrop -r300 -native -painters -transparent -q101