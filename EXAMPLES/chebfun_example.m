close all;
x = chebfun('x');
f = @(x) x.^(2).*sin(2*exp(2*sin(2*exp(2*x.^2))));
fc = chebfun(f,1e4);
LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize';
figure, plot(fc,LW,1.2)
title('Function f',FS,16)
format long
Ichebfun = sum(fc)
Npts = length(fc);
[s,w] = chebpts(Npts);
Iclenshawcurtis = w*f(s)
[s,w] = legpts(Npts);
Igauss = w*f(s)
figure, tic, err = [];
NN = 10:10:500;
Igauss_konrod = quadgk(f,-1,1,'RelTol',1e-10,'AbsTol',1e-15,'MaxIntervalCount',1e15)
for Npts = NN
    [s,w] = legpts(Npts);
    Igauss = w*f(s);
    err = [err abs(Igauss-Ichebfun)];
end
semilogy(NN,err,'.-',LW,1,MS,16), grid on
ylim([1e-18 1])
xlabel('Npts',FS,12), ylabel('Error',FS,12)
title('Gauss quadrature convergence',FS,16), toc
hold on; tic, err = [];
for Npts = NN
    [s,w] = chebpts(Npts);
    Iclenshawcurtis = w*f(s);
    err = [err abs(Iclenshawcurtis-Ichebfun)];
end
semilogy(NN,err,'.-r',LW,1,MS,16)
title('Gauss and Clenshaw-Curtis',FS,16)
legend('Gauss','Clenshaw-Curtis','location','southwest'), toc
