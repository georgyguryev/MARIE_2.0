function [handle] = visualize_currents_(obj,COIL,E,port_cut,freq_cut)

    % Visualizes the currents from one port on the coil structure
    % Author: Ioannis P. Georgakis, Moscow, 2019
    
    handle = figure();
    
    % crop port and frequency
    Jc =  squeeze(E(:,port_cut,freq_cut));
    
    first_node  = [3 1 2];
    second_node = [2 3 1];
    
    j_cart = zeros( size(COIL.elem,2),3);
    
    for jj = 1: size(COIL.elem,2)
        
        rp = COIL.Ct(:,jj);
        r1 = COIL.node(:,COIL.elem(1,jj));
        r2 = COIL.node(:,COIL.elem(2,jj));
        r3 = COIL.node(:,COIL.elem(3,jj));
        
        Ae = obj.triangle_area(r1,r2,r3);
        Ae1 = obj.triangle_area(rp,r2,r3);
        Ae2 = obj.triangle_area(r1,rp,r3);
        Ae3 = obj.triangle_area(r1,r2,rp);
        
        zs = [Ae1 Ae2 Ae3]./Ae;
        
        ed = abs(COIL.etod(:,jj));
        
        lo = [r2(1)-r3(1) r3(1)-r1(1) r1(1)-r2(1); r2(2)-r3(2) r3(2)-r1(2) r1(2)-r2(2); r2(3)-r3(3) r3(3)-r1(3) r1(3)-r2(3)];
        
        for id = 1:3
            
            if COIL.index(ed(id))~=0
                
                i_1 = second_node(id);
                
                i_2 = first_node(id);
                
                L_i = sqrt(lo(1,id)^2 + lo(2,id)^2 + lo(3,id)^2);
                
                j_cart(jj,:) = j_cart(jj,:) + (L_i/(2*Ae)) * (zs(i_2)*lo(:,i_1)' - zs(i_1)*lo(:,i_2)') * sign(COIL.etod(id,jj)) * Jc(COIL.index(ed(id)),1);
                
            end
        end
    end
    
    j_cart = real(j_cart);
    
    quiver3(COIL.Ct(1,:)',COIL.Ct(2,:)',COIL.Ct(3,:)',j_cart(:,1),j_cart(:,2),j_cart(:,3),2,'k');
    xlabel('x','interpreter','latex');
    ylabel('y','interpreter','latex');
    zlabel('z','interpreter','latex');
    axis equal;
    grid on;    
    
    
end