%% Evaluate complex integrals
% This program invokes MUPAD functions by using the command
% feval(symengine,'MuPAD_Function',x1,...,xn);
% where symengine is a symbolic engine
% MUPAD_function is the exact format of the mupad function
% x1,...,xn are the variable involved
syms r theta
p = .9;
% s = symengine
% fun = @(r) (r.^(p-1)./(r+1));
% fun = @(r) sqrt(1-r.^2)./(1+r.^2);
% fun = @(r) log(r)./(1+r.^2).^2;
% fun = @(r) sqrt(r)./(r.^2 + 6*r + 8);
% fun = @(r) r.^(3/4).*(3-r).^(1/4)./(5-r);
% fun = @(r) 1./sqrt(4*r.^2 + 4*r + 3);
fun = @(r) r.^(-p)./(r+1);
%%  Poisson Kernel (Lookup wikipedia) r = .4;;
a = 0;
b = inf;
% % solve(f)
% % MaxCalls = infinity;
% d = evalin(s, 'numeric::quadrature',f, 'r = 0 .. infinity');
% d = vpa(d)
% DIGITS:= 50:
% numeric::quadrature((r^(.9-1)/(r+1)), r = 0 .. infinity, MaxCalls = infinity)
q = quadgk(fun,a,b,...
    'MaxIntervalCount',20e5,...
    'RelTol',1e-6,...
    'AbsTol',1e-6)
Q = integral(fun,a,b,...
    'RelTol',1e-12,...
    'AbsTol',1e-12)
% I = pi/(sin(pi*p));
% I = pi*(sqrt(2)-1);
% I = -pi/4
% I = pi*(1-1/sqrt(2))
% I = pi/(2*sqrt(2))*(17-40^(3/4))
% I = -1i*pi
I = pi/sin(p*pi)
I - q
ezplot3(r,r,fun)