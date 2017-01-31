%% Calculate 2DEG condunctivity

%% \sigma = N e \mu
% Sigma = Electrical Conductivity
% e = Electron Charge
% mu = Electron Mobility

% From [1] 
% [1]  Computational model of 2DEG mobility in AlGaN/GaN heterostructures
% Karine Abgaryan*, Ilya Mutigullin, and Dmitry Reviznikov
% Phys. Status Solidi C 12, No. 4–5, 460–465 (2015) / DOI 10.1002/pssc.201400200

mu = 1e6; % cm^2/(Vs)
e = 1.602e-19; % C 
N = 1e13; % cm^-3
sigma = N * e * mu; % C/(cm * Vs)
