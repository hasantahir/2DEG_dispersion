function plotCylinder(radius)
    shading('interp'); 
    colorbar;
    hold('on')
    theta=linspace(0,2*pi,50);
    plot(radius*cos(theta), radius*sin(theta),'--w','linewidth', 1.0)
    axis('equal', 'tight');
    set(gca,'fontname', 'arial', 'fontsize',8)
    xlabel('x (\lambda_{0})');
    ylabel('y (\lambda_{0})');
end

