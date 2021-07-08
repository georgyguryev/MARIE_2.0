function visualize_normals_(coil)

    % Visualizes the normals of the triangular mesh of the coil structure
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    figure
    
    TR = triangulation(coil.elem(1:3,:)',coil.node');
    F = faceNormal(TR);
    P = incenter(TR);
    
    trisurf(TR,'FaceColor',[0.72 0.45 0.2],'LineStyle','none','FaceAlpha',.8); hold on;
    quiver3(P(:,1),P(:,2),P(:,3),F(:,1),F(:,2),F(:,3),1,'color','r');
    xlabel('x','interpreter','latex');
    ylabel('y','interpreter','latex');
    zlabel('z','interpreter','latex');
    axis equal;
    grid on;
end