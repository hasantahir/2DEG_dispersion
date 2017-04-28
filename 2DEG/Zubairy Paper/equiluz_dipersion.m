%% Equiluz dispersion relation
clear all
e0 = 8.85e-12; % Free-space permittivity
c = 3e8; % speed of light
e = 1.60218e-19; % electron charge
me = 9.109e-31; % Electron mass

w = linspace(.001e12, 2e12, 1e1); % Frequency in Hz
% w = .3e12;
% Structure properties
d = 100e-6;
N = 2e12*1e4;
es = 12;
eo = 3.7;
ms = .2*me;


root = [];
for i = 1 : length(w)
%     r = [];
    % Nakayama conductivity
    sigma = 1i*e^2*N./(ms*w(i));
    
    % Eguiluz polarizability
    X = sigma./(-1i*w(i));
    
    Bs = @(q) sqrt(q.^2 - (w(i)/c).^2*es);
    Bo = @(q) sqrt(q.^2 - (w(i)/c).^2*eo);
    
    f = @(q) 4*pi*X + es./Bs(q) + e0.*coth(Bo(q)*d)./Bo(q);
    
    r = newtzero(f,1e2)
    fval = f(r);
    s = min((fval));
    D = find(f(r) == s);
    if length(D) > 1
        
        if real(r(D(1)) > 0)
            root = vertcat(root,r(D(1)));
        else
            root = vertcat(root,r(D(2)));
        end
    else
        root = vertcat(root,r(D(1)));
    end

end


% close all
plot(real(root)/c,w);
hold on
plot(w/c,w);
% set(gca,'Xscale','Log','Yscale','Log');