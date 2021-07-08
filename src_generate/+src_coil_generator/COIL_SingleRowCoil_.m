function[mshfile] = COIL_SingleRowCoil_(rad1,rad2,dist,len,z,t,meshing,shield)

    % Generates an elliptical single row decoupled coil with 8 ports
    % Author: Ilias I. Giannakopoulos, Atlanta, GA, 2019
    
    %% Write in .geo file
    
    if shield.flag == 1
        geofile  = strcat('./data/coils/geo_files/SingleRow_decoupled_shield_semiaxis_x_',num2str(rad1*100),'cm_semiaxis_y_',num2str(rad2*100),'cm_distance_',num2str(dist*100),'cm_len_',num2str(len*100),'cm_displace_',num2str(z*100),'cm_mesh_',num2str(meshing),'.geo');
        mshfile = strcat('./data/coils/msh_files/SingleRow_decoupled_shield_semiaxis_x_',num2str(rad1*100),'cm_semiaxis_y_',num2str(rad2*100),'cm_distance_',num2str(dist*100),'cm_len_',num2str(len*100),'cm_displace_',num2str(z*100),'cm_mesh_',num2str(meshing),'.msh');
    else
        geofile  = strcat('./data/coils/geo_files/SingleRow_decoupled_semiaxis_x_',num2str(rad1*100),'cm_semiaxis_y_',num2str(rad2*100),'cm_distance_',num2str(dist*100),'cm_len_',num2str(len*100),'cm_displace_',num2str(z*100),'cm_mesh_',num2str(meshing),'.geo');
        mshfile = strcat('./data/coils/msh_files/SingleRow_decoupled_semiaxis_x_',num2str(rad1*100),'cm_semiaxis_y_',num2str(rad2*100),'cm_distance_',num2str(dist*100),'cm_len_',num2str(len*100),'cm_displace_',num2str(z*100),'cm_mesh_',num2str(meshing),'.msh');
    end
    
    fileID = fopen(geofile,'w');
    fprintf(fileID,'lc = %f; \n\n',meshing);   
    
    %% Create Coil   
    %% Points
    % Centers of ellipse
    fprintf(fileID,'Point( 1000 )   = {%f, %f, %f, lc}; \n',[0; 0; len/2+z]);
    fprintf(fileID,'Point( 2000 )   = {%f, %f, %f, lc}; \n',[0; 0; len/2-t+z]);
    fprintf(fileID,'Point( 3000 )   = {%f, %f, %f, lc}; \n',[0; 0; t-len/2+z]);
    fprintf(fileID,'Point( 4000 )   = {%f, %f, %f, lc}; \n',[0; 0; -len/2+z]);
    if rad1>rad2
        focal = sqrt(rad1^2-rad2^2);
        fprintf(fileID,'Point( 1001 )   = {%f, %f, %f, lc}; \n',[focal; 0; len/2+z]);
        fprintf(fileID,'Point( 2001 )   = {%f, %f, %f, lc}; \n',[focal; 0; len/2-t+z]);
        fprintf(fileID,'Point( 3001 )   = {%f, %f, %f, lc}; \n',[focal; 0; t-len/2+z]);
        fprintf(fileID,'Point( 4001 )   = {%f, %f, %f, lc}; \n',[focal; 0; -len/2+z]);
    else 
        focal = sqrt(rad2^2-rad1^2);
        fprintf(fileID,'Point( 1001 )   = {%f, %f, %f, lc}; \n',[0; focal; len/2+z]);
        fprintf(fileID,'Point( 2001 )   = {%f, %f, %f, lc}; \n',[0; focal; len/2-t+z]);
        fprintf(fileID,'Point( 3001 )   = {%f, %f, %f, lc}; \n',[0; focal; t-len/2+z]);
        fprintf(fileID,'Point( 4001 )   = {%f, %f, %f, lc}; \n',[0; focal; -len/2+z]);
    end
    % James Ivory and Bessel formula
    h = (rad1-rad2)^2/(rad1+rad2)^2;
    % sum_ell = 0;
    % for n = 2:100
    %     sum_ell = sum_ell + (fact2(2*n-3)/(2^n*factorial(n)))^2*h^n; 
    % end
    % C = pi*(rad1+rad2)*(1+h/4+sum_ell);
    % Srinivasa Ramanujan formulas
    C = pi*(rad1+rad2)*(1+3*h/(10+sqrt(4-3*h)));
    
    skip = round(dist/C*1000);
    keep = round(1/8*1000);
    t_step = round(t/C*1000);
    dim_4 = round(keep/4);
    dim_2 = round(keep/2);
    dim_34 = round(3*keep/4);
    dim_t_step = keep-t_step;
    
    [x_all,y_all] = get_points_ellipse(rad1,rad2,(skip+keep)*8);
%     plot(x_all,y_all,'r-');
    x = zeros(8,keep);
    y = zeros(8,keep);
    pos1 = zeros(8,2);
    pos2 = zeros(8,2);
    pos3 = zeros(8,2);
    pos4 = zeros(8,2);
    pos5 = zeros(8,2);
    pos6 = zeros(8,2);
    pos7 = zeros(8,2);
     for i = 1:8
        x(i,:) = x_all(skip*(i-1)+keep*(i-1)+1:skip*(i-1)+keep*(i-1)+keep);
        y(i,:) = y_all(skip*(i-1)+keep*(i-1)+1:skip*(i-1)+keep*(i-1)+keep);
     end
     for i = 1:8
         pos1(i,1) = x(i,1);
         pos1(i,2) = y(i,1);
         pos2(i,1) = x(i,dim_4);
         pos2(i,2) = y(i,dim_4);
         pos3(i,1) = x(i,dim_2);
         pos3(i,2) = y(i,dim_2);
         pos4(i,1) = x(i,dim_34);
         pos4(i,2) = y(i,dim_34);
         pos5(i,1) = x(i,keep);
         pos5(i,2) = y(i,keep);
         pos6(i,1) = x(i,t_step);
         pos6(i,2) = y(i,t_step);
         pos7(i,1) = x(i,dim_t_step);
         pos7(i,2) = y(i,dim_t_step);
     end
    
    for i = 1:8
        % Outer
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+1     ; pos1(i,1); pos1(i,2); len/2+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+8     ; pos2(i,1); pos2(i,2); len/2+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+2     ; pos3(i,1); pos3(i,2); len/2+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+12   ; pos4(i,1); pos4(i,2); len/2+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+3     ; pos5(i,1); pos5(i,2); len/2+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+4     ; pos5(i,1); pos5(i,2); len/4+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+5     ; pos5(i,1); pos5(i,2); 0+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+6     ; pos5(i,1); pos5(i,2); -len/4+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+7     ; pos5(i,1); pos5(i,2); -len/2+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+9     ; pos4(i,1); pos4(i,2); -len/2+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+10   ; pos3(i,1); pos3(i,2); -len/2+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+11   ; pos2(i,1); pos2(i,2); -len/2+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+13   ; pos1(i,1); pos1(i,2); -len/2+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+14   ; pos1(i,1); pos1(i,2); -len/4+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+15   ; pos1(i,1); pos1(i,2); 0+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+16   ; pos1(i,1); pos1(i,2); len/4+z]);
        % Inner
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+17   ; pos6(i,1); pos6(i,2); len/2-t+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+24   ; pos2(i,1); pos2(i,2); len/2-t+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+18   ; pos3(i,1); pos3(i,2); len/2-t+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+28   ; pos4(i,1); pos4(i,2); len/2-t+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+19   ; pos7(i,1); pos7(i,2); len/2-t+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+20   ; pos7(i,1); pos7(i,2); len/4+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+21   ; pos7(i,1); pos7(i,2); 0+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+22   ; pos7(i,1); pos7(i,2); -len/4+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+23   ; pos7(i,1); pos7(i,2); -len/2+t+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+25   ; pos4(i,1); pos4(i,2); -len/2+t+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+26   ; pos3(i,1); pos3(i,2); -len/2+t+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+27   ; pos2(i,1); pos2(i,2); -len/2+t+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+29   ; pos6(i,1); pos6(i,2); -len/2+t+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+30   ; pos6(i,1); pos6(i,2); -len/4+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+31   ; pos6(i,1); pos6(i,2); 0+z]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+32   ; pos6(i,1); pos6(i,2); len/4+z]);
    end
    %% Ellipses
    for i = 1:8
        % Outer
        fprintf(fileID,'Ellipse( %d ) = {%d,%d,%d,%d}; \n',[48*(i-1)+1  ;32*(i-1)+1; 1000 ;1001;32*(i-1)+8]);
        fprintf(fileID,'Ellipse( %d ) = {%d,%d,%d,%d}; \n',[48*(i-1)+8  ;32*(i-1)+12; 1000 ;1001;32*(i-1)+2]);
        fprintf(fileID,'Ellipse( %d ) = {%d,%d,%d,%d}; \n',[48*(i-1)+2  ;32*(i-1)+2; 1000 ;1001;32*(i-1)+8]);
        fprintf(fileID,'Ellipse( %d ) = {%d,%d,%d,%d}; \n',[48*(i-1)+24  ;32*(i-1)+12; 1000 ;1001;32*(i-1)+3]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+3  ;32*(i-1)+3 ;32*(i-1)+4]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+4  ;32*(i-1)+4; 32*(i-1)+5]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+5  ;32*(i-1)+5  ;32*(i-1)+6]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+6  ;32*(i-1)+6  ;32*(i-1)+7]);
        fprintf(fileID,'Ellipse( %d ) = {%d,%d,%d,%d}; \n',[48*(i-1)+7  ;32*(i-1)+7; 4000  ;4001;32*(i-1)+9]);
        fprintf(fileID,'Ellipse( %d ) = {%d,%d,%d,%d}; \n',[48*(i-1)+9  ;32*(i-1)+9; 4000 ;4001;32*(i-1)+10]);
        fprintf(fileID,'Ellipse( %d ) = {%d,%d,%d,%d}; \n',[48*(i-1)+10  ;32*(i-1)+10; 4000 ;4001;32*(i-1)+11]);
        fprintf(fileID,'Ellipse( %d ) = {%d,%d,%d,%d}; \n',[48*(i-1)+11  ;32*(i-1)+11; 4000 ;4001;32*(i-1)+13]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+13  ;32*(i-1)+13;  32*(i-1)+14]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+14  ;32*(i-1)+14;  32*(i-1)+15]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+15  ;32*(i-1)+15;  32*(i-1)+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+16  ;32*(i-1)+16;  32*(i-1)+1]);
        % Inner
        fprintf(fileID,'Ellipse( %d ) = {%d,%d,%d,%d}; \n',[48*(i-1)+17  ;32*(i-1)+1+16 ; 2000 ;2001;32*(i-1)+24]);
        fprintf(fileID,'Ellipse( %d ) = {%d,%d,%d,%d}; \n',[48*(i-1)+12  ;32*(i-1)+28; 2000 ;2001;32*(i-1)+2+16]);
        fprintf(fileID,'Ellipse( %d ) = {%d,%d,%d,%d}; \n',[48*(i-1)+18  ;32*(i-1)+2+16; 2000 ;2001;32*(i-1)+24]);
        fprintf(fileID,'Ellipse( %d ) = {%d,%d,%d,%d}; \n',[48*(i-1)+28  ;32*(i-1)+28; 2000 ;2001;32*(i-1)+3+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+19  ;32*(i-1)+3+16 ;32*(i-1)+4+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+20  ;32*(i-1)+4+16; 32*(i-1)+5+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+21  ;32*(i-1)+5+16  ;32*(i-1)+6+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+22  ;32*(i-1)+6+16  ;32*(i-1)+7+16]);
        fprintf(fileID,'Ellipse( %d ) = {%d,%d,%d,%d}; \n',[48*(i-1)+23  ;32*(i-1)+7+16; 3000  ;3001;32*(i-1)+9+16]);
        fprintf(fileID,'Ellipse( %d ) = {%d,%d,%d,%d}; \n',[48*(i-1)+25  ;32*(i-1)+9+16; 3000 ;3001;32*(i-1)+10+16]);
        fprintf(fileID,'Ellipse( %d ) = {%d,%d,%d,%d}; \n',[48*(i-1)+26  ;32*(i-1)+10+16; 3000 ;3001;32*(i-1)+11+16]);
        fprintf(fileID,'Ellipse( %d ) = {%d,%d,%d,%d}; \n',[48*(i-1)+27  ;32*(i-1)+11+16; 3000 ;3001;32*(i-1)+13+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+29  ;32*(i-1)+13+16;  32*(i-1)+14+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+30  ;32*(i-1)+14+16;  32*(i-1)+15+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+31  ;32*(i-1)+15+16;  32*(i-1)+16+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+32  ;32*(i-1)+16+16;  32*(i-1)+1+16]);
        % Angles
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+33  ;32*(i-1)+1  ;32*(i-1)+1+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+34  ;32*(i-1)+3  ;32*(i-1)+3+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+35  ;32*(i-1)+7  ;32*(i-1)+7+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+36  ;32*(i-1)+13  ;32*(i-1)+13+16]);
        % Feed ports and lumped elements
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+37  ;32*(i-1)+2  ;32*(i-1)+2+16]); % Port
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+38  ;32*(i-1)+4  ;32*(i-1)+4+16]); % Capacitor right top
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+39  ;32*(i-1)+5  ;32*(i-1)+5+16]); % Capacitor right middle
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+40  ;32*(i-1)+6  ;32*(i-1)+6+16]); % Capacitor right bottom
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+41  ;32*(i-1)+12  ;32*(i-1)+28]); % Capacitor top left
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+42  ;32*(i-1)+9  ;32*(i-1)+9+16]); % Capacitor 4
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+43  ;32*(i-1)+10  ;32*(i-1)+10+16]); % Capacitor 3
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+44  ;32*(i-1)+11  ;32*(i-1)+11+16]); % Capacitor 2
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+45  ;32*(i-1)+8  ;32*(i-1)+24]); % Capacitor top right
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+46  ;32*(i-1)+14  ;32*(i-1)+14+16]); % Capacitor left top
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+47  ;32*(i-1)+15  ;32*(i-1)+15+16]); % Capacitor left middle
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+48  ;32*(i-1)+16  ;32*(i-1)+16+16]); % Capacitor left bottom
        
    end
    
    %% Line Loops
    for i = 1:8
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+1   ;   -(48*(i-1)+1)   ;   -(48*(i-1)+45)    ; 48*(i-1)+17      ; 48*(i-1)+33]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+2   ;   48*(i-1)+2      ;      48*(i-1)+45    ; -(48*(i-1)+18)   ; -(48*(i-1)+37)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+3   ;   -(48*(i-1)+3)   ;   -(48*(i-1)+38)    ; 48*(i-1)+3+16  ; 48*(i-1)+34]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+4   ;   -(48*(i-1)+4)   ;   -(48*(i-1)+39)    ; 48*(i-1)+4+16  ; 48*(i-1)+38]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+5   ;   -(48*(i-1)+5)   ;   -(48*(i-1)+40)    ; 48*(i-1)+5+16  ; 48*(i-1)+39]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+6   ;   -(48*(i-1)+6)   ;   -(48*(i-1)+35)    ; 48*(i-1)+6+16  ; 48*(i-1)+40]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+7   ;   48*(i-1)+8      ;      48*(i-1)+37    ; -(48*(i-1)+12)   ; -(48*(i-1)+41)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+8   ;   -(48*(i-1)+7 )  ;   -(48*(i-1)+42)    ; 48*(i-1)+8+15  ; 48*(i-1)+35]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+9   ;   -(48*(i-1)+9 )  ;   -(48*(i-1)+43)    ; 48*(i-1)+9+16  ; 48*(i-1)+42]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+10 ;   -(48*(i-1)+10) ;   -(48*(i-1)+44)    ; 48*(i-1)+10+16 ; 48*(i-1)+43]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+11 ;   -(48*(i-1)+11) ;   -(48*(i-1)+36)    ; 48*(i-1)+11+16 ; 48*(i-1)+44]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+12 ;   -(48*(i-1)+24) ;   -(48*(i-1)+34)    ; 48*(i-1)+12+16 ; 48*(i-1)+41]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+13 ;   -(48*(i-1)+13) ;   -(48*(i-1)+46)    ; 48*(i-1)+13+16 ; 48*(i-1)+36]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+14 ;   -(48*(i-1)+14) ;   -(48*(i-1)+47)    ; 48*(i-1)+14+16 ; 48*(i-1)+46]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+15 ;   -(48*(i-1)+15) ;   -(48*(i-1)+48)    ; 48*(i-1)+15+16 ; 48*(i-1)+47]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+16 ;   -(48*(i-1)+16) ;   -(48*(i-1)+33)    ; 48*(i-1)+16+16 ; 48*(i-1)+48]);
    end
    
    
    %% Surfaces
    for i=1:8
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*(i-1)+1  ;16*(i-1)+1  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*(i-1)+2  ;16*(i-1)+2  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*(i-1)+3  ;16*(i-1)+3  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*(i-1)+4  ;16*(i-1)+4  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*(i-1)+5  ;16*(i-1)+5  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*(i-1)+6  ;16*(i-1)+6  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*(i-1)+7  ;16*(i-1)+7  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*(i-1)+8  ;16*(i-1)+8  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*(i-1)+9  ;16*(i-1)+9  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*(i-1)+10  ;16*(i-1)+10  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*(i-1)+11  ;16*(i-1)+11  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*(i-1)+12  ;16*(i-1)+12  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*(i-1)+13  ;16*(i-1)+13  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*(i-1)+14  ;16*(i-1)+14  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*(i-1)+15  ;16*(i-1)+15  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*(i-1)+16  ;16*(i-1)+16  ]);
    end
    
    % Physical Lines
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[i      ;48*(i-1)+37]); % Port
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[8+i  ;48*(i-1)+45]); % Capacitor right top
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[16+i;48*(i-1)+41]); % Capacitor left top
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[24+i;48*(i-1)+42]); % Capacitor right bottom
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[32+i;48*(i-1)+44]); % Capacitor left bottom
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[40+i;48*(i-1)+38]); % Capacitor right up
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[48+i;48*(i-1)+40]); % Capacitor right down
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[56+i;48*(i-1)+48]); % Capacitor left up
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[64+i;48*(i-1)+46]); % Capacitor left down
    end
    
    
    if shield.flag == 1
        fprintf(fileID,'lc_s = %f; \n\n',shield.meshing);   
        
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[32*8+1  ; 0    ; 0    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[32*8+2  ; 0    ; 0    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[32*8+3  ; shield.radius    ; 0    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[32*8+4  ; 0    ; shield.radius    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[32*8+5  ; -shield.radius    ; 0    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[32*8+6  ; 0    ; -shield.radius    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[32*8+7  ; shield.radius    ; 0    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[32*8+8  ; 0    ; shield.radius    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[32*8+9  ; -shield.radius    ; 0    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[32*8+10  ; 0    ; -shield.radius    ; -shield.length/2]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*8+1; 32*8+3;32*8+7 ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*8+2; 32*8+4;32*8+8 ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*8+3; 32*8+5;32*8+9 ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*8+4; 32*8+6;32*8+10 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*8+5; 32*8+3;32*8+1;32*8+4 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*8+6; 32*8+4;32*8+1;32*8+5 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*8+7; 32*8+5;32*8+1;32*8+6 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*8+8; 32*8+6;32*8+1;32*8+3 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*8+9; 32*8+7;32*8+2;32*8+8 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*8+10; 32*8+8;32*8+2;32*8+9 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*8+11; 32*8+9;32*8+2;32*8+10 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*8+12; 32*8+10;32*8+2;32*8+7 ]);
        
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*8+1;   48*8+1 ;   48*8+9    ; -(48*8+2); -(48*8+5)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*8+2;   48*8+2 ;   48*8+10    ; -(48*8+3); -(48*8+6)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*8+3;   48*8+3 ;   48*8+11    ; -(48*8+4); -(48*8+7)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*8+4;   48*8+4 ;   48*8+12    ; -(48*8+1); -(48*8+8)]);
        
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*8+1  ;16*8+1  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*8+2  ;16*8+2  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*8+3  ;16*8+3  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[16*8+4  ;16*8+4  ]);
        
        fprintf(fileID,'Physical Surface(2)  = {16*8+1,16*8+2,16*8+3,16*8+4};\n\n');
        
    end
    
    %% Physical Surfaces
    fprintf(fileID,'Physical Surface(1)  = {');
    for i=1:8
        fprintf(fileID,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d',[...
            16*(i-1)+1,16*(i-1)+2,16*(i-1)+3,16*(i-1)+4,16*(i-1)+5,16*(i-1)+6,16*(i-1)+7,16*(i-1)+8,...
            16*(i-1)+9,16*(i-1)+10,16*(i-1)+11,16*(i-1)+12,16*(i-1)+13,16*(i-1)+14,16*(i-1)+15,16*(i-1)+16]);
        if i~=8
            fprintf(fileID,',');
        end
    end
    fprintf(fileID,'}; \n\n');
    
    command = sprintf('./data/coils/GMSH/bin/gmsh %s -1 -2 -o %s', geofile,mshfile);
    [status,cmdout] = system(command); 
    
    %% Add comments and close
    fprintf(fileID,'Coherence Mesh; \n\n');
    fprintf(fileID,'View "comments" { \n');
    fprintf(fileID,'T2(10, -10, 0){ "Copyright (C) Ilias Giannakopoulos" }; \n');
    fprintf(fileID,'T2(10, 15, 0){ StrCat("File created on ", Today) }; }; \n\n');
    fclose(fileID);

end

function [x,y]=get_points_ellipse(rad1,rad2,points_ellipse)

     t = linspace(0,2*pi,points_ellipse);
    if rad1>rad2
        x1 = sqrt(rad1^2-rad2^2);
        y1 = 0;
        x2 = -sqrt(rad1^2-rad2^2);
        y2 = 0;
         X = rad1*cos(t);
         Y = rad2*sin(t);
    else
        x1 = 0;
        y1 = sqrt(rad2^2-rad1^2);
        x2 = 0;
        y2 = -sqrt(rad2^2-rad1^2);
         X = rad2*cos(t);
         Y = rad1*sin(t);
    end
    
     w = atan2(y2-y1,x2-x1);
     x = (x1+x2)/2 + X*cos(w) - Y*sin(w);
     y = (y1+y2)/2 + X*sin(w) + Y*cos(w);
     
     if rad1 == rad2
        t = linspace(0,2*pi,points_ellipse);
        x = rad1*cos(t);
        y = rad2*sin(t);
     end
        
end