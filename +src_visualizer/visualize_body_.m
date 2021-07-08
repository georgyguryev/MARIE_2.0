function visualize_body_(xd,yd,zd,id)

    % Visualizes the body model structure
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    figure();
    
    P = [xd(id) yd(id) zd(id)];
    k = boundary(P,1);
    
    trisurf(k,P(:,1),P(:,2),P(:,3),'LineStyle','none');
    colormap([linspace(0/256,30/256,256)', linspace(144/256,191/256,256)', ones(256,1)]);
    camlight headlight;
    lighting gouraud;
    xlabel('x','interpreter','latex');
    ylabel('y','interpreter','latex');
    zlabel('z','interpreter','latex');
    axis equal;
    grid on;
    
end