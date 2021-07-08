function[coil_out] = COIL_Displacement_(Cnt,rot,coil_out)
% function COIL_Displacement_(Cnt, rot, coil_out)
% Displaces the coil structure
% Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    coil_out.node(1,:) = coil_out.node(1,:)+Cnt(1);
    coil_out.node(2,:) = coil_out.node(2,:)+Cnt(2);
    coil_out.node(3,:) = coil_out.node(3,:)+Cnt(3);
    
    mesh_center = sum(coil_out.node,2)/size(coil_out.node,2);
    rot_x = cosd(rot)*(coil_out.node(1,:)-mesh_center(1))  - sind(rot) *(coil_out.node(2,:)-mesh_center(2)) + mesh_center(1);
    rot_y = sind(rot) *(coil_out.node(1,:)-mesh_center(1)) + cosd(rot)*(coil_out.node(2,:)-mesh_center(2)) + mesh_center(2);
    coil_out.node(1,:) = rot_x;
    coil_out.node(2,:) = rot_y;
    
    [Ct,Ln,Pn] = Mesh_CLP(coil_out.node,coil_out.elem);
    
    coil_out.Ct = Ct;
    coil_out.Ln = Ln;
    coil_out.Pn = Pn;
    
end