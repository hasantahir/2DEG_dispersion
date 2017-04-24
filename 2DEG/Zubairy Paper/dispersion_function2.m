clear all
w = 2*pi*.5e12;

e_0 = 8.85e-12;
h = 6.626e-34/(2*pi);
c = 3e8;
e = 1.60218e-19; 
m_e = 9.109e-31; % Electron mass

% Material properties
Ns = 1e11 *1e4;
e_d = 9.6;
ms = .067*m_e; % Effective mass
mu = 1e6*1e-4;
tau = 1.14e-10;

% layer thicknesses
d = 200e-9;

k0 = w/c;
kz = linspace(1e6,1e7,1e3);

root = [];
for i = 1 : length(kz)
    r = newtzero(@dispersion,kz(i));
    root = vertcat(root,r);
end
% Sort the array
root = sort(root);
% Clean up roots by weeding out too close values
if ~isempty(root)
    cnt = 1;  % Counter for while loop.
    
    while ~isempty(root)
        vct = abs(root - root(1)) < 1e-1; % Minimum spacing between roots.
        C = root(vct);  % C has roots grouped close together.
        [idx,idx] = min(abs(dispersion(C)));  % Pick the best root per group.
        rt(cnt) = C(idx); %  Most root vectors are small.
        root(vct) = []; % Deplete the pool of roots.
        cnt = cnt + 1;  % Increment the counter.
    end
    root = sort(rt).';  % return a nice, sorted column vector
end