w = linspace(.1e12, 2e12, 1e2); % Frequency in Hz
load root 
load root1
load root2
load root3
load root4

plot(real(root),w)
hold on
plot(real(root1),w)
plot(real(root2),w)
plot(real(root3),w)
plot(real(root4),w)


