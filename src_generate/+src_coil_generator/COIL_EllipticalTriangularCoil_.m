function[mshfile] = COIL_EllipticalTriangularCoil_(a,b,len,t,meshing,shield)

% Generates an octagonal elliptical triangular coil with 8 ports
% Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    %% Write in .geo file
    
    if shield.flag == 1
        geofile  = strcat('./data/coils/geo_files/EllipticalTriagular_shield_decoupled_rad_x_',num2str(a*100),'cm_rad_y_',num2str(b*100),'cm_len_',num2str(len*100),'cm_mesh_',num2str(meshing),'.geo');
        mshfile = strcat('./data/coils/msh_files/EllipticalTriagular_shield_decoupled_rad_x_',num2str(a*100),'cm_rad_y_',num2str(b*100),'cm_len_',num2str(len*100),'cm_mesh_',num2str(meshing),'.msh');
    else
        geofile  = strcat('./data/coils/geo_files/EllipticalTriagular_decoupled_rad_x_',num2str(a*100),'cm_rad_y_',num2str(b*100),'cm_len_',num2str(len*100),'cm_mesh_',num2str(meshing),'.geo');
        mshfile = strcat('./data/coils/msh_files/EllipticalTriagular_decoupled_rad_x_',num2str(a*100),'cm_rad_y_',num2str(b*100),'cm_len_',num2str(len*100),'cm_mesh_',num2str(meshing),'.msh');
    end
    
    fileID = fopen(geofile,'w');
    fprintf(fileID,'lc = %f; \n\n',meshing);   
    
    %% Create Coil   
    
    %% Points
    pos=get_points_ellipse(a,b);
    
    for i = 1:8
        j = mod(i,8)+1;
        mid_x = (pos(i,1)+pos(j,1))/2;
        mid_y = (pos(i,2)+pos(j,2))/2;
        [x_off1,x_off2,y_off1,y_off2] = get_points_line(pos(i,1),pos(i,2),pos(j,1),pos(j,2),t);
        [x_off3,x_off4,y_off3,y_off4] = get_points_line(pos(i,1),pos(i,2),pos(j,1),pos(j,2),t/2);
        d = sqrt((pos(i,1)-x_off2)^2+(pos(i,2)-y_off2)^2);
        if mod(i,2) == 1
            [mid_x_1,~,mid_y_1,~] = get_points_line(pos(i,1),pos(i,2),x_off2,y_off2,1/4*d);
            [mid_x_1_s,~,mid_y_1_s,~] = get_points_line(pos(i,1),pos(i,2),x_off2,y_off2,t/2);
            [mid_x_2,~,mid_y_2,~] = get_points_line(pos(i,1),pos(i,2),x_off2,y_off2,3/4*d);
%             [mid_x_2_s,~,mid_y_2_s,~] = get_points_line(pos(i,1),pos(i,2),x_off2,y_off2,d-t/2);
            [mid_x_3,~,mid_y_3,~] = get_points_line(x_off1,y_off1,pos(j,1),pos(j,2),1/4*d);
%             [mid_x_3_s,~,mid_y_3_s,~] = get_points_line(x_off1,y_off1,pos(j,1),pos(j,2),t/2);
            [mid_x_4,~,mid_y_4,~] = get_points_line(x_off1,y_off1,pos(j,1),pos(j,2),3/4*d);
            [mid_x_4_s,~,mid_y_4_s,~] = get_points_line(x_off1,y_off1,pos(j,1),pos(j,2),d-t/2);
        else
            [~,mid_x_1,~,mid_y_1] = get_points_line(x_off1,y_off1,pos(j,1),pos(j,2),1/4*d);
            [~,mid_x_1_s,~,mid_y_1_s] = get_points_line(x_off1,y_off1,pos(j,1),pos(j,2),t/2);
            [~,mid_x_2,~,mid_y_2] = get_points_line(x_off1,y_off1,pos(j,1),pos(j,2),3/4*d);
%             [~,mid_x_2_s,~,mid_y_2_s] = get_points_line(x_off1,y_off1,pos(j,1),pos(j,2),d-t/2);
            [~,mid_x_3,~,mid_y_3] = get_points_line(pos(i,1),pos(i,2),x_off2,y_off2,1/4*d);
%             [~,mid_x_3_s,~,mid_y_3_s] = get_points_line(pos(i,1),pos(i,2),x_off2,y_off2,t/2);
            [~,mid_x_4,~,mid_y_4] = get_points_line(pos(i,1),pos(i,2),x_off2,y_off2,3/4*d);
            [~,mid_x_4_s,~,mid_y_4_s] = get_points_line(pos(i,1),pos(i,2),x_off2,y_off2,d-t/2);
        end
        mid_z_1 = (t-len/2)/2;
        mid_z_2 = (len/2-t)/2;
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+1    ; pos(i,1)  ; pos(i,2) ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+2    ; pos(i,1)  ; pos(i,2) ; len/2   ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+3    ; pos(i,1)  ; pos(i,2) ; t-len/2 ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+4    ; pos(i,1)  ; pos(i,2) ; len/2-t ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+5    ; mid_x    ; mid_y    ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+6    ; mid_x    ; mid_y    ; len/2   ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+7    ; mid_x    ; mid_y    ; t-len/2 ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+8    ; mid_x    ; mid_y    ; len/2-t ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+9    ; x_off2    ; y_off2    ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+10  ; x_off2    ; y_off2    ; len/2   ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+11  ; x_off2    ; y_off2    ; t-len/2 ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+12  ; x_off2    ; y_off2    ; len/2-t ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+13  ; x_off1    ; y_off1    ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+14  ; x_off1    ; y_off1    ; len/2   ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+15  ; x_off1    ; y_off1    ; t-len/2 ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+16  ; x_off1    ; y_off1    ; len/2-t ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+17  ; mid_x_1    ; mid_y_1    ; mid_z_1 ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+18  ; mid_x_2    ; mid_y_2    ; mid_z_2 ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+19  ; mid_x_3    ; mid_y_3    ; mid_z_1 ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+20  ; mid_x_4    ; mid_y_4    ; mid_z_2 ]);
        
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+21  ; x_off4    ; y_off4    ; len/2   ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+22  ; x_off4    ; y_off4    ; len/2-t ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+23  ; x_off3    ; y_off3    ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+24  ; x_off3    ; y_off3    ; t-len/2 ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+25  ; x_off4    ; y_off4    ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+26  ; x_off4    ; y_off4    ; t-len/2 ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+27  ; x_off3    ; y_off3    ; len/2   ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+28  ; x_off3    ; y_off3    ; len/2-t ]);
        
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+29    ; pos(i,1)  ; pos(i,2) ; 2*t-len/2 ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+30    ; pos(i,1)  ; pos(i,2) ; len/2-2*t ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+31  ; mid_x_1_s    ; mid_y_1_s    ; 2*t-len/2 ]);
%         fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+32  ; mid_x_3_s    ; mid_y_3_s    ; 2*t-len/2 ]);
%         fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+33  ; mid_x_2_s    ; mid_y_2_s    ; len/2-2*t ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[34*(i-1)+34  ; mid_x_4_s    ; mid_y_4_s    ; len/2-2*t ]);
        
    end
    %% Lines
    for i = 1:8
        j = mod(i,8)+1;
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+2 ; 34*(i-1)+13;34*(i-1)+5]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+3 ; 34*(i-1)+5;34*(i-1)+9]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+5 ; 34*(i-1)+3;34*(i-1)+24]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+32 ; 34*(i-1)+24;34*(i-1)+15]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+6 ; 34*(i-1)+15;34*(i-1)+7]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+7 ; 34*(i-1)+7;34*(i-1)+11]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+8 ; 34*(i-1)+11;34*(i-1)+26]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+33 ; 34*(i-1)+26;34*(j-1)+3]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+10 ; 34*(i-1)+14;34*(i-1)+6]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+11 ; 34*(i-1)+6;34*(i-1)+10]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+13 ; 34*(i-1)+4;34*(i-1)+28]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+38 ; 34*(i-1)+28;34*(i-1)+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+14 ; 34*(i-1)+16;34*(i-1)+8]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+15 ; 34*(i-1)+8;34*(i-1)+12]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+16 ; 34*(i-1)+12;34*(i-1)+22]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+34 ; 34*(i-1)+22;34*(j-1)+4]);
        
        if mod(i,2) == 1
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+12 ; 34*(i-1)+10;34*(i-1)+21]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+9 ; 34*(i-1)+2;34*(i-1)+14]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+1 ; 34*(i-1)+23;34*(i-1)+13]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+4 ; 34*(i-1)+9;34*(j-1)+1]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+17 ; 34*(i-1)+31;34*(i-1)+17]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+18 ; 34*(i-1)+17;34*(i-1)+18]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+19 ; 34*(i-1)+18;34*(i-1)+12]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+20 ; 34*(i-1)+34;34*(i-1)+20]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+21 ; 34*(i-1)+20;34*(i-1)+19]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+22 ; 34*(i-1)+19;34*(i-1)+15]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+30 ; 34*(i-1)+23;34*(i-1)+24]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+24 ; 34*(i-1)+21;34*(i-1)+22]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+26 ; 34*(i-1)+2;34*(i-1)+4]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+35 ; 34*(i-1)+3;34*(i-1)+29]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+36 ; 34*(i-1)+29;34*(i-1)+31]);
            if i == 1
                fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+37 ; 34*(i-1)+29;34*(7)+31]);
            else
                fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+37 ; 34*(i-1)+29;34*(i-2)+31]);
            end
        else
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+12 ; 34*(i-1)+10;34*(j-1)+2]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+9 ; 34*(i-1)+27;34*(i-1)+14]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+1 ; 34*(i-1)+1;34*(i-1)+13]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+4 ; 34*(i-1)+9;34*(i-1)+25]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+17 ; 34*(i-1)+31;34*(i-1)+17]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+18 ; 34*(i-1)+17;34*(i-1)+18]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+19 ; 34*(i-1)+18;34*(i-1)+16]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+20 ; 34*(i-1)+20;34*(i-1)+34]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+21 ; 34*(i-1)+19;34*(i-1)+20]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+22 ; 34*(i-1)+11;34*(i-1)+19]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+30 ; 34*(i-1)+1;34*(i-1)+3]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+31 ; 34*(i-1)+25;34*(i-1)+26]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+25 ; 34*(i-1)+27;34*(i-1)+28]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+35 ; 34*(i-1)+4;34*(i-1)+30]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+36 ; 34*(i-1)+30;34*(i-1)+34]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+37 ; 34*(i-1)+30;34*(i-2)+34]);
        end
        
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+23 ; 34*(i-1)+7;34*(i-1)+5]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+27 ; 34*(i-1)+6;34*(i-1)+8]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+28 ; 34*(i-1)+17;34*(i-1)+19]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*(i-1)+29 ; 34*(i-1)+18;34*(i-1)+20]);
    end
   %% Line Loops
    for i = 1:8
        j = mod(i,8)+1;
        if mod(i,2) == 1
            fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d,%d,%d}; \n',[7*(i-1)+1;38*(i-1)+1;38*(i-1)+2;-(38*(i-1)+23);-(38*(i-1)+6);-(38*(i-1)+32);-(38*(i-1)+30)]);
            fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d,%d,%d,%d}; \n',[7*(i-1)+2;38*(i-1)+23;38*(i-1)+3;38*(i-1)+4;38*(j-1)+30;-(38*(i-1)+33);-(38*(i-1)+8);-(38*(i-1)+7)]);
            fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d,%d,%d,%d}; \n',[7*(i-1)+3;-(38*(i-1)+35);-(38*(i-1)+36);-(38*(i-1)+17);-(38*(i-1)+28);-(38*(i-1)+22);38*(i-1)+32;38*(i-1)+5]);
            fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[7*(i-1)+4;38*(i-1)+28;-(38*(i-1)+21);-(38*(i-1)+29);-(38*(i-1)+18)]);
            fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d,%d,%d,%d}; \n',[7*(i-1)+5;38*(i-1)+29;-(38*(i-1)+20);-(38*(j-1)+37);-(38*(j-1)+35);-(38*(i-1)+34);-(38*(i-1)+16);-(38*(i-1)+19)]);
            fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d,%d,%d,%d}; \n',[7*(i-1)+6;38*(i-1)+26;38*(i-1)+13;38*(i-1)+38;38*(i-1)+14;-(38*(i-1)+27);-(38*(i-1)+10);-(38*(i-1)+9)]);
            fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d,%d,%d}; \n',[7*(i-1)+7;38*(i-1)+27;38*(i-1)+15;38*(i-1)+16;-(38*(i-1)+24);-(38*(i-1)+12);-(38*(i-1)+11);]);
        else
            fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d,%d,%d,%d}; \n',[7*(i-1)+1;38*(i-1)+1;38*(i-1)+2;-(38*(i-1)+23);-(38*(i-1)+6);-(38*(i-1)+32);-(38*(i-1)+5);-(38*(i-1)+30)]);
            fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d,%d,%d}; \n',[7*(i-1)+2;38*(i-1)+23;38*(i-1)+3;38*(i-1)+4;38*(i-1)+31;-(38*(i-1)+8);-(38*(i-1)+7)]);
            fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d,%d,%d,%d}; \n',[7*(i-1)+3;38*(i-1)+17;38*(i-1)+28;-(38*(i-1)+22);38*(i-1)+8;38*(i-1)+33;38*(j-1)+35;38*(j-1)+37]);
            fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[7*(i-1)+4;-(38*(i-1)+28);38*(i-1)+18;38*(i-1)+29;-(38*(i-1)+21)]);
            fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d,%d,%d,%d}; \n',[7*(i-1)+5;-(38*(i-1)+29);-(38*(i-1)+20);38*(i-1)+36;38*(i-1)+35;-(38*(i-1)+13);-(38*(i-1)+38);38*(i-1)+19]);
            fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d,%d,%d}; \n',[7*(i-1)+6;38*(i-1)+25;38*(i-1)+38;38*(i-1)+14;-(38*(i-1)+27);-(38*(i-1)+10);-(38*(i-1)+9)]);
            fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d,%d,%d,%d}; \n',[7*(i-1)+7;38*(i-1)+27;38*(i-1)+15;38*(i-1)+16;38*(i-1)+34;-(38*(j-1)+26);-(38*(i-1)+12);-(38*(i-1)+11)]);
        end
    end
     
    %% Surfaces
    for i=1:8
        fprintf(fileID,'Plane Surface( %d ) = {%d}; \n',[7*(i-1)+1  ;7*(i-1)+1  ]);
        fprintf(fileID,'Plane Surface( %d ) = {%d}; \n',[7*(i-1)+2  ;7*(i-1)+2  ]);
        fprintf(fileID,'Plane Surface( %d ) = {%d}; \n',[7*(i-1)+3  ;7*(i-1)+3  ]);
        fprintf(fileID,'Plane Surface( %d ) = {%d}; \n',[7*(i-1)+4  ;7*(i-1)+4  ]);
        fprintf(fileID,'Plane Surface( %d ) = {%d}; \n',[7*(i-1)+5  ;7*(i-1)+5  ]);
        fprintf(fileID,'Plane Surface( %d ) = {%d}; \n',[7*(i-1)+6  ;7*(i-1)+6  ]);
        fprintf(fileID,'Plane Surface( %d ) = {%d}; \n',[7*(i-1)+7  ;7*(i-1)+7  ]);
    end
    
    %% Physical Lines
    fprintf(fileID,'Physical Line( 1 ) = {%d}; \n',26+76*0);
    fprintf(fileID,'Physical Line( 2 ) = {%d}; \n',68+76*0);
    fprintf(fileID,'Physical Line( 3 ) = {%d}; \n',26+76*1);
    fprintf(fileID,'Physical Line( 4 ) = {%d}; \n',68+76*1);
    fprintf(fileID,'Physical Line( 5 ) = {%d}; \n',26+76*2);
    fprintf(fileID,'Physical Line( 6 ) = {%d}; \n',68+76*2);
    fprintf(fileID,'Physical Line( 7 ) = {%d}; \n',26+76*3);
    fprintf(fileID,'Physical Line( 8 ) = {%d}; \n',68+76*3);
    fprintf(fileID,'Physical Line( 9 ) = {%d}; \n',27+76*0);
    fprintf(fileID,'Physical Line( 10 ) = {%d}; \n',27+76*1);
    fprintf(fileID,'Physical Line( 11 ) = {%d}; \n',27+76*2);
    fprintf(fileID,'Physical Line( 12 ) = {%d}; \n',27+76*3);
    fprintf(fileID,'Physical Line( 13 ) = {%d}; \n',65+76*0);
    fprintf(fileID,'Physical Line( 14 ) = {%d}; \n',65+76*1);
    fprintf(fileID,'Physical Line( 15 ) = {%d}; \n',65+76*2);
    fprintf(fileID,'Physical Line( 16 ) = {%d}; \n',65+76*3);
    fprintf(fileID,'Physical Line( 17 ) = {%d}; \n',61+76*0);
    fprintf(fileID,'Physical Line( 18 ) = {%d}; \n',61+76*1);
    fprintf(fileID,'Physical Line( 19 ) = {%d}; \n',61+76*2);
    fprintf(fileID,'Physical Line( 20 ) = {%d}; \n',61+76*3);
    fprintf(fileID,'Physical Line( 21 ) = {%d}; \n',23+76*0);
    fprintf(fileID,'Physical Line( 22 ) = {%d}; \n',23+76*1);
    fprintf(fileID,'Physical Line( 23 ) = {%d}; \n',23+76*2);
    fprintf(fileID,'Physical Line( 24 ) = {%d}; \n',23+76*3);
    fprintf(fileID,'Physical Line( 25 ) = {%d}; \n',28+76*0);
    fprintf(fileID,'Physical Line( 26 ) = {%d}; \n',28+76*1);
    fprintf(fileID,'Physical Line( 27 ) = {%d}; \n',28+76*2);
    fprintf(fileID,'Physical Line( 28 ) = {%d}; \n',28+76*3);
    fprintf(fileID,'Physical Line( 29 ) = {%d}; \n',29+76*0);
    fprintf(fileID,'Physical Line( 30 ) = {%d}; \n',29+76*1);
    fprintf(fileID,'Physical Line( 31 ) = {%d}; \n',29+76*2);
    fprintf(fileID,'Physical Line( 32 ) = {%d}; \n',29+76*3);
    fprintf(fileID,'Physical Line( 33 ) = {%d}; \n',66+76*0);
    fprintf(fileID,'Physical Line( 34 ) = {%d}; \n',66+76*1);
    fprintf(fileID,'Physical Line( 35 ) = {%d}; \n',66+76*2);
    fprintf(fileID,'Physical Line( 36 ) = {%d}; \n',66+76*3);
    fprintf(fileID,'Physical Line( 37 ) = {%d}; \n',67+76*0);
    fprintf(fileID,'Physical Line( 38 ) = {%d}; \n',67+76*1);
    fprintf(fileID,'Physical Line( 39 ) = {%d}; \n',67+76*2);
    fprintf(fileID,'Physical Line( 40 ) = {%d}; \n',67+76*3);
    fprintf(fileID,'Physical Line( 41 ) = {%d}; \n',46+76*3);
    fprintf(fileID,'Physical Line( 42 ) = {%d}; \n',32+76*0);
    fprintf(fileID,'Physical Line( 43 ) = {%d}; \n',46+76*0);
    fprintf(fileID,'Physical Line( 44 ) = {%d}; \n',32+76*1);
    fprintf(fileID,'Physical Line( 45 ) = {%d}; \n',46+76*1);
    fprintf(fileID,'Physical Line( 46 ) = {%d}; \n',32+76*2);
    fprintf(fileID,'Physical Line( 47 ) = {%d}; \n',46+76*2);
    fprintf(fileID,'Physical Line( 48 ) = {%d}; \n',32+76*3);
    fprintf(fileID,'Physical Line( 49 ) = {%d}; \n',16+76*0);
    fprintf(fileID,'Physical Line( 50 ) = {%d}; \n',76+76*0);
    fprintf(fileID,'Physical Line( 51 ) = {%d}; \n',16+76*1);
    fprintf(fileID,'Physical Line( 52 ) = {%d}; \n',76+76*1);
    fprintf(fileID,'Physical Line( 53 ) = {%d}; \n',16+76*2);
    fprintf(fileID,'Physical Line( 54 ) = {%d}; \n',76+76*2);
    fprintf(fileID,'Physical Line( 55 ) = {%d}; \n',16+76*3);
    fprintf(fileID,'Physical Line( 56 ) = {%d}; \n',76+76*3);
    
    if shield.flag == 1
        fprintf(fileID,'lc_s = %f; \n\n',shield.meshing);   
        
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[34*8+1  ; 0    ; 0    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[34*8+2  ; 0    ; 0    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[34*8+3  ; shield.radius    ; 0    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[34*8+4  ; 0    ; shield.radius    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[34*8+5  ; -shield.radius    ; 0    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[34*8+6  ; 0    ; -shield.radius    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[34*8+7  ; shield.radius    ; 0    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[34*8+8  ; 0    ; shield.radius    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[34*8+9  ; -shield.radius    ; 0    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[34*8+10  ; 0    ; -shield.radius    ; -shield.length/2]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*8+1; 34*8+3;34*8+7 ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*8+2; 34*8+4;34*8+8 ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*8+3; 34*8+5;34*8+9 ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[38*8+4; 34*8+6;34*8+10 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[38*8+5; 34*8+3;34*8+1;34*8+4 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[38*8+6; 34*8+4;34*8+1;34*8+5 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[38*8+7; 34*8+5;34*8+1;34*8+6 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[38*8+8; 34*8+6;34*8+1;34*8+3 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[38*8+9; 34*8+7;34*8+2;34*8+8 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[38*8+10; 34*8+8;34*8+2;34*8+9 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[38*8+11; 34*8+9;34*8+2;34*8+10 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[38*8+12; 34*8+10;34*8+2;34*8+7 ]);
        
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[7*8+1;   38*8+1 ;   38*8+9    ; -(38*8+2); -(38*8+5)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[7*8+2;   38*8+2 ;   38*8+10    ; -(38*8+3); -(38*8+6)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[7*8+3;   38*8+3 ;   38*8+11    ; -(38*8+4); -(38*8+7)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[7*8+4;   38*8+4 ;   38*8+12    ; -(38*8+1); -(38*8+8)]);
        
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[7*8+1  ;7*8+1  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[7*8+2  ;7*8+2  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[7*8+3  ;7*8+3  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[7*8+4  ;7*8+4  ]);
        
        fprintf(fileID,'Physical Surface(2)  = {7*8+1,7*8+2,7*8+3,7*8+4};\n\n');
        
    end
    
    %% Physical Surfaces
    fprintf(fileID,'Physical Surface(1)  = {');
    for i=1:8
        fprintf(fileID,'%d,%d,%d,%d,%d,%d,%d',[...
            7*(i-1)+1,7*(i-1)+2,7*(i-1)+3,7*(i-1)+4,7*(i-1)+5,7*(i-1)+6,7*(i-1)+7]);
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

function [xout_1,xout_2,yout_1,yout_2] = get_points_line(x1,y1,x2,y2,t)

    d = sqrt((x2-x1)^2+(y2-y1)^2);
    t = t/d;
    xout_1 = (1-t)*x1+t*x2;
    xout_2 = (1-t)*x2+t*x1;

    yout_1 = (1-t)*y1+t*y2;
    yout_2 = (1-t)*y2+t*y1;
        
end

function pos=get_points_ellipse(a,b)
    theta_=linspace(0,pi/2,1001);
    theta=[theta_(1:end-1) theta_(1:end-1)+pi/2 theta_(1:end-1)+pi theta_(1:end-1)+3/2*pi];
    x=a*cos(theta);
    y=b*sin(theta);
    arc_length = norm(sqrt(diff(x).^2+diff(y).^2),1);
    arc_lengths_vec=arc_length*1/8*linspace(0,7,8);
    k=2; 
    d=0;
    pos(1,:)=[x(1) y(1)];
    alpha = zeros(1,length(arc_lengths_vec));
    alpha(1)=-90;
    for i=2:length(arc_lengths_vec)
        while d<arc_lengths_vec(i)
            d = norm(sqrt(diff(x(1:k)).^2+diff(y(1:k)).^2),1);
            k=k+1;
        end
        pos(i,:)=[x(k-1) y(k-1)]; %position of coils
        alpha(i)=atand((y(k-2)-y(k))/(x(k-2)-x(k)));% tangent, for the angle
        if y(k-1)>0 && x(k-1)>0
            alpha(i)=alpha(i);
        elseif y(k-1)>0 && x(k-1)<0
            alpha(i)=alpha(i);
        elseif y(k-1)<0 && x(k-1)<0
            alpha(i)=alpha(i)+180;
        elseif y(k-1)<0 && x(k-1)>0
            alpha(i)=alpha(i)+180;
        end
    end
end