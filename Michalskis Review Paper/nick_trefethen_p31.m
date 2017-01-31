x = chebfun('x');

f = sin(12*x).*exp(-x);         % A function on [-1, 1]
g = max(f, 1./(x+2));           % The max of f and 1./(x+2)
plot(g)                         % A function with discontinuous derivative
sum(g)                          % The integral of g
plot(diff(g))                   % The derivative of g
h = g + x - .8;                 % A function with several roots in [-1, 1]
rr = roots(h);                  % Compute the roots of h
plot(h, 'k', rr, h(rr), 'ro')   % Plot h and its roots

figure(2)
sym y
t = @(y) sin(12*y).*exp(-y);
h =  max(t,1./(y+2));
ezplot(h)
hold on
symsum(h)
plot(diff(h))
i = h + y - .8;
rrr = roots(i);
ezplot(i,'k', rrr, i(rrr), 'ro')