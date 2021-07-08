function visualize_coil_(coil)

    % Visualizes the coil structure and the ports
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    figure
    
    pos = zeros(length(coil.port),3);
    s = zeros(length(coil.port),1);
    
    for i = 1:length(coil.port)
        pos(i,:) = coil.node(:,coil.edge(1,coil.index == coil.port(i).t(1)));
        s(i) = i;
    end
    
    s = string(s);
    
    TR = triangulation(coil.elem(1:3,:)',coil.node');
    
    trisurf(TR,'FaceColor',[0.72 0.45 0.2],'LineStyle','none','FaceAlpha',.8); hold on;
    xlabel('x','interpreter','latex');
    ylabel('y','interpreter','latex');
    zlabel('z','interpreter','latex');
    text(pos(:,1),pos(:,2),pos(:,3),s,'HorizontalAlignment','left','Color','b','FontSize',10);
    axis equal;
    grid on;
    
end