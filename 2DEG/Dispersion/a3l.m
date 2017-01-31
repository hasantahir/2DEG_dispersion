% File name: a3L.m
% Analysis for TE modes
% First, we conduct test plots to find ranges of beta where possible
% solutions exists. This is done in several steps:
% 1. Plots are done treating kappa_f as an independent variable
% 2. Ranges of kappa_f are determined where there are zeros of functions
% 3. Corresponding ranges of beta are determined
% 4. Searches are performed to find propagation constants
%
clear all
% Definition of structure
n_f = 1.50; % ref. index of film layer
n_s = 1.45; % ref. index of substrate
n_c = 1.40; % ref. index of cladding
lambda = 1.0; % wavelength in microns
h = 5.0; % thickness of film layer in microns
a_c = 2*h; % thickness of cladding region
a_s = 2*h; % thickness of substrate region
%
k = 2*pi/lambda; % wave number
kappa_f = 0:0.01:3.0; % establish range of kappa_f
beta_temp = sqrt((n_f*k)^2 - kappa_f.^2);
beta_min = min(beta_temp);
beta_max = max(beta_temp);
% Before searches, we plot search function versus beta
beta = beta_min:0.001:beta_max; % establish range of beta
N = beta./k;
ff = func_asym(beta,n_c,n_s,n_f,k,h);
plot(beta,ff)
xlabel('\beta','FontSize',22);
ylabel('Search function','FontSize',22);
ylim([-10.0 10.0])
grid on
pause
close all
% From the above plot, one must choose proper search range for each mode.
% Search numbers provided below are only for the waveguide defined above.
% For different waveguide, one must choose different ranges
% for searches.
beta0 = fzero(@(beta) func_asym(beta,n_c,n_s,n_f,k,h),[9.40 9.41])
beta1 = fzero(@(beta) func_asym(beta,n_c,n_s,n_f,k,h),[9.35 9.37])
beta2 = fzero(@(beta) func_asym(beta,n_c,n_s,n_f,k,h),[9.27 9.29])
beta3 = fzero(@(beta) func_asym(beta,n_c,n_s,n_f,k,h),[9.17 9.18])
% Plot of field profiles
A_s = 1.0;
thickness = h + a_c + a_s;
beta_field = beta0; % Select appropriate propagation constant for plotting
gamma_s = sqrt(beta_field^2 - (n_s*k)^2);
gamma_c = sqrt(beta_field^2 - (n_c*k)^2);
kappa_f = sqrt((n_f*k)^2 - beta_field^2);
%
% In the formulas below for electric field E_y we have shifted
% x-coordinate by a_s
% We also 'reversed' direction of plot in the substrate region
NN = 100;
delta = thickness/NN;
x = 0.0:delta:thickness; % coordinates of plot points
x_t = 0;
for i=1:NN+1
x_t(i+1)= x_t(i) + delta;
if (x_t(i)<=a_s);
E_y(i) = A_s*exp(gamma_s*(x_t(i)-a_s));
elseif (a_s<=x_t(i)) && (x_t(i)<=a_s+h);
E_y(i) = A_s*(cos(kappa_f*(x_t(i)-a_s))+...
gamma_s*sin(kappa_f*(x_t(i)-a_s))/kappa_f);
else (a_s+h<=x_t(i)) & (x_t(i)<=thickness);
E_y(i) = A_s*(cos(kappa_f*h)+gamma_s*sin(kappa_f*h)/kappa_f)...
*exp(-gamma_c*(x_t(i)-h-a_s));
end
end
%
h=plot(x,E_y);
% add text on x-axix and y-axis and size of x and y labels
xlabel('x (microns)','FontSize',22);
ylabel('TE electric field','FontSize',22);
set(h,'LineWidth',1.5); % new thickness of plotting lines
set(gca,'FontSize',22); % new size of tick marks on both axes
grid on
pause
close all