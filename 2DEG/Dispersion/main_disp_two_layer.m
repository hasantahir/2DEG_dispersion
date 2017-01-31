close all;clear all
% Define frequency sweep
w = 2*pi*linspace(1e10,1e15, 1e4);
% w = 2*pi*linspace(1e7,10e12, 1e5);


c = 3e8;
ep0_0 = 12.9;
ep0_inf = 11.0;
lambda_l = 1/(292.77e2);
w_l = 2*pi*c/lambda_l;
lambda_t = 1/(268.7e2);
w_t = 2*pi*c/lambda_t;
lambda_gt = 1/(2.4e2);
g_t = 2*pi*c/lambda_gt;
k_t = w_t/c;

% InAs
% c = 3e8;
% ep0_0 = 12.9;
% ep0_inf = 11.8;
% w_t = 4.12e13;
% w_l = 4.58e13;
% k_t = w_t/c;
% w_p = 2.95e13;
% W_p = 2.61e12;
% m = .021*m_e;
% Constant Parameters
h = 6.626e-34/(2*pi);
c = 3e8;
e2 = (2.3068e-28);
% e = 1.60218e-19; 
m_e = 9.109e-31; % Electron mass

% Calculate Plasma Frequency of GaAs
% N = .9e15*1e4;
Ns = 5e12*1e4;
m = .063*m_e;

% W_p = (4*pi*Ns*e2/m/c); % Square of plasma frequency
% w_p = sqrt(W_p)


% ep0 = ep0_inf*(w.^2 - w_l^2)./(w.^2 - w_t^2) - w_p^2./w.^2;
% Get Lorentz-based materials
ep2 = Lorentz(ep0_0, ep0_inf, w_t, g_t,w);
sigma = Ns*e2./(1i*m*w);
W_p = (4*pi*Ns*e2/m); % Square of plasma frequency
w_p = sqrt(W_p);

% w = 1e11;
i = 500;
% for i = 1 : length(w)
    ep1 = 1;
    k1 = (w(i)/c).^2;
    k2 = (w(i)/c).^2.*ep2(i);
    
    kz1 = @(kx) sqrt(k1.^2 - kx.^2);
    kz2 = @(kx) sqrt(k2.^2 - kx.^2);
    func = @(kx) ep2(i)./kz1(kx) + 1./kz2(kx) - 1i*Ns*e2./(m.*w(i).^2);
    
    lxlim = k1;
    uxlim = k2;
    num = 100;
    p = linspace(lxlim,uxlim,num);
    root = [];
    for i = 1 : length(p)
        r = newtzero(func,p(i));
        root = vertcat(root,r);
    end
    
    % Sort the array
    root = sort(root);
    roots = root((real(root))>1e3);% & abs(imag(root))<tol);
    if ~isempty(roots)
    cnt = 1;  % Counter for while loop.
    
    while ~isempty(roots)
        vct = abs(roots - roots(1)) < 1e-4; % Minimum spacing between roots.
        C = roots(vct);  % C has roots grouped close together.
        [idx,idx] = min(abs(func(C)));  % Pick the best root per group.
        rt(cnt) = C(idx); %  Most root vectors are small.
        roots(vct) = []; % Deplete the pool of roots.
        cnt = cnt + 1;  % Increment the counter.
    end
    roots = sort(rt).';  % return a nice, sorted column vector
end
% end
