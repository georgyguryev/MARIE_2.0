function[mshfile] = COIL_SelfCoil_(rad,dist,len,t,meshing,shield)
   
    % Generates a self decoupled coil with 8 ports
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    %% Write in .geo file
    
    if shield.flag == 1
        geofile  = strcat('./data/coils/geo_files/Self_decoupled_shield_rad_',num2str(rad*100),'cm_distance_',num2str(dist*100),'cm_len_',num2str(len*100),'cm_mesh_',num2str(meshing),'.geo');
        mshfile = strcat('./data/coils/msh_files/Self_decoupled_shield_rad_',num2str(rad*100),'cm_distance_',num2str(dist*100),'cm_len_',num2str(len*100),'cm_mesh_',num2str(meshing),'.msh');
    else
        geofile  = strcat('./data/coils/geo_files/Self_decoupled_rad_',num2str(rad*100),'cm_distance_',num2str(dist*100),'cm_len_',num2str(len*100),'cm_mesh_',num2str(meshing),'.geo');
        mshfile = strcat('./data/coils/msh_files/Self_decoupled_rad_',num2str(rad*100),'cm_distance_',num2str(dist*100),'cm_len_',num2str(len*100),'cm_mesh_',num2str(meshing),'.msh');
    end
    
    fileID = fopen(geofile,'w');
    fprintf(fileID,'lc = %f; \n\n',meshing);   
    
    %% Create Coil   
    %% Points
    dim = (2*pi*rad-8*dist)/8;
    offset = dim/2;
    x = [0;dim+dist;2*(dim+dist);3*(dim+dist);4*(dim+dist);5*(dim+dist);6*(dim+dist);7*(dim+dist)];
    x = x + offset;
    
    for i = 1:8
        % Outer
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+1     ; rad*cos((x(i))/rad)  ; rad*sin((x(i))/rad) ; len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+2     ; rad*cos((x(i)+dim/2)/rad)  ; rad*sin((x(i)+dim/2)/rad) ; len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+3     ; rad*cos((x(i)+dim)/rad)  ; rad*sin((x(i)+dim)/rad) ; len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+4     ; rad*cos((x(i)+dim)/rad)  ; rad*sin((x(i)+dim)/rad) ; len/4  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+5     ; rad*cos((x(i)+dim)/rad)  ; rad*sin((x(i)+dim)/rad) ; 0  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+6     ; rad*cos((x(i)+dim)/rad)  ; rad*sin((x(i)+dim)/rad) ; -len/4  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+7     ; rad*cos((x(i)+dim)/rad)  ; rad*sin((x(i)+dim)/rad) ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+8     ; rad*cos((x(i)+5*dim/6)/rad)  ; rad*sin((x(i)+5*dim/6)/rad) ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+9     ; rad*cos((x(i)+2*dim/3)/rad)  ; rad*sin((x(i)+2*dim/3)/rad) ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+10     ; rad*cos((x(i)+dim/2)/rad)  ; rad*sin((x(i)+dim/2)/rad) ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+11     ; rad*cos((x(i)+dim/3)/rad)  ; rad*sin((x(i)+dim/3)/rad) ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+12     ; rad*cos((x(i)+dim/6)/rad)  ; rad*sin((x(i)+dim/6)/rad) ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+13     ; rad*cos((x(i))/rad)  ; rad*sin((x(i))/rad) ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+14     ; rad*cos((x(i))/rad)  ; rad*sin((x(i))/rad) ; -len/4  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+15     ; rad*cos((x(i))/rad)  ; rad*sin((x(i))/rad) ; 0  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+16     ; rad*cos((x(i))/rad)  ; rad*sin((x(i))/rad) ; len/4  ]);
        % Inner
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+17     ; rad*cos((x(i)+t)/rad)  ; rad*sin((x(i)+t)/rad) ; len/2-t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+18     ; rad*cos((x(i)+dim/2)/rad)  ; rad*sin((x(i)+dim/2)/rad) ; len/2-t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+19     ; rad*cos((x(i)-t+dim)/rad)  ; rad*sin((x(i)-t+dim)/rad) ; len/2-t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+20     ; rad*cos((x(i)-t+dim)/rad)  ; rad*sin((x(i)-t+dim)/rad) ; len/4  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+21     ; rad*cos((x(i)-t+dim)/rad)  ; rad*sin((x(i)-t+dim)/rad) ; 0  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+22     ; rad*cos((x(i)-t+dim)/rad)  ; rad*sin((x(i)-t+dim)/rad) ; -len/4  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+23     ; rad*cos((x(i)-t+dim)/rad)  ; rad*sin((x(i)-t+dim)/rad) ; -len/2+t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+24     ; rad*cos((x(i)+5*dim/6)/rad)  ; rad*sin((x(i)+5*dim/6)/rad) ; -len/2+t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+25     ; rad*cos((x(i)+2*dim/3)/rad)  ; rad*sin((x(i)+2*dim/3)/rad) ; -len/2+t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+26     ; rad*cos((x(i)+dim/2)/rad)  ; rad*sin((x(i)+dim/2)/rad) ; -len/2+t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+27     ; rad*cos((x(i)+dim/3)/rad)  ; rad*sin((x(i)+dim/3)/rad) ; -len/2+t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+28     ; rad*cos((x(i)+dim/6)/rad)  ; rad*sin((x(i)+dim/6)/rad) ; -len/2+t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+29     ; rad*cos((x(i)+t)/rad)  ; rad*sin((x(i)+t)/rad) ; -len/2+t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+30     ; rad*cos((x(i)+t)/rad)  ; rad*sin((x(i)+t)/rad) ; -len/4  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+31     ; rad*cos((x(i)+t)/rad)  ; rad*sin((x(i)+t)/rad) ; 0  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[32*(i-1)+32     ; rad*cos((x(i)+t)/rad)  ; rad*sin((x(i)+t)/rad) ; len/4  ]);
    end
    fprintf(fileID,'Point( 1000 )   = {%f, %f, %f, lc}; \n',[0; 0; len/2]);
    fprintf(fileID,'Point( 2000 )   = {%f, %f, %f, lc}; \n',[0; 0; len/2-t]);
    fprintf(fileID,'Point( 3000 )   = {%f, %f, %f, lc}; \n',[0; 0; t-len/2]);
    fprintf(fileID,'Point( 4000 )   = {%f, %f, %f, lc}; \n',[0; 0; -len/2]);
    %% Circles
    for i = 1:8
        % Outer
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*(i-1)+1  ;32*(i-1)+1     ; 1000 ;32*(i-1)+2]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*(i-1)+2  ;32*(i-1)+2; 1000 ;32*(i-1)+3]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+3  ;32*(i-1)+3 ;32*(i-1)+4]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+4  ;32*(i-1)+4; 32*(i-1)+5]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+5  ;32*(i-1)+5  ;32*(i-1)+6]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+6  ;32*(i-1)+6  ;32*(i-1)+7]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*(i-1)+7  ;32*(i-1)+7; 4000  ;32*(i-1)+8]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*(i-1)+8  ;32*(i-1)+8; 4000 ;32*(i-1)+9]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*(i-1)+9  ;32*(i-1)+9; 4000 ;32*(i-1)+10]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*(i-1)+10  ;32*(i-1)+10; 4000 ;32*(i-1)+11]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*(i-1)+11  ;32*(i-1)+11; 4000 ;32*(i-1)+12]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*(i-1)+12  ;32*(i-1)+12; 4000 ;32*(i-1)+13]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+13  ;32*(i-1)+13;  32*(i-1)+14]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+14  ;32*(i-1)+14;  32*(i-1)+15]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+15  ;32*(i-1)+15;  32*(i-1)+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+16  ;32*(i-1)+16;  32*(i-1)+1]);
        % Inner
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*(i-1)+17  ;32*(i-1)+1+16     ; 1000 ;32*(i-1)+2+16]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*(i-1)+18  ;32*(i-1)+2+16; 1000 ;32*(i-1)+3+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+19  ;32*(i-1)+3+16 ;32*(i-1)+4+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+20  ;32*(i-1)+4+16; 32*(i-1)+5+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+21  ;32*(i-1)+5+16  ;32*(i-1)+6+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+22  ;32*(i-1)+6+16  ;32*(i-1)+7+16]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*(i-1)+23  ;32*(i-1)+7+16; 4000  ;32*(i-1)+8+16]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*(i-1)+24  ;32*(i-1)+8+16; 4000 ;32*(i-1)+9+16]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*(i-1)+25  ;32*(i-1)+9+16; 4000 ;32*(i-1)+10+16]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*(i-1)+26  ;32*(i-1)+10+16; 4000 ;32*(i-1)+11+16]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*(i-1)+27  ;32*(i-1)+11+16; 4000 ;32*(i-1)+12+16]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[48*(i-1)+28  ;32*(i-1)+12+16; 4000 ;32*(i-1)+13+16]);
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
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+38  ;32*(i-1)+4  ;32*(i-1)+4+16]); % Inductor right top
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+39  ;32*(i-1)+5  ;32*(i-1)+5+16]); % Inductor right middle
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+40  ;32*(i-1)+6  ;32*(i-1)+6+16]); % Inductor right bottom
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+41  ;32*(i-1)+8  ;32*(i-1)+8+16]); % Capacitor 5
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+42  ;32*(i-1)+9  ;32*(i-1)+9+16]); % Capacitor 4
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+43  ;32*(i-1)+10  ;32*(i-1)+10+16]); % Capacitor 3
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+44  ;32*(i-1)+11  ;32*(i-1)+11+16]); % Capacitor 2
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+45  ;32*(i-1)+12  ;32*(i-1)+12+16]); % Capacitor 1
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+46  ;32*(i-1)+14  ;32*(i-1)+14+16]); % Inductor left top
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+47  ;32*(i-1)+15  ;32*(i-1)+15+16]); % Inductor left middle
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[48*(i-1)+48  ;32*(i-1)+16  ;32*(i-1)+16+16]); % Inductor left bottom
        
    end
    
    %% Line Loops
    for i = 1:8
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+1   ;   -(48*(i-1)+1)   ;   -(48*(i-1)+37)    ; 48*(i-1)+1+16  ; 48*(i-1)+33]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+2   ;   -(48*(i-1)+2)   ;   -(48*(i-1)+34)    ; 48*(i-1)+2+16  ; 48*(i-1)+37]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+3   ;   -(48*(i-1)+3)   ;   -(48*(i-1)+38)    ; 48*(i-1)+3+16  ; 48*(i-1)+34]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+4   ;   -(48*(i-1)+4)   ;   -(48*(i-1)+39)    ; 48*(i-1)+4+16  ; 48*(i-1)+38]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+5   ;   -(48*(i-1)+5)   ;   -(48*(i-1)+40)    ; 48*(i-1)+5+16  ; 48*(i-1)+39]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+6   ;   -(48*(i-1)+6)   ;   -(48*(i-1)+35)    ; 48*(i-1)+6+16  ; 48*(i-1)+40]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+7   ;   -(48*(i-1)+7)   ;   -(48*(i-1)+41)    ; 48*(i-1)+7+16  ; 48*(i-1)+35]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+8   ;   -(48*(i-1)+8)   ;   -(48*(i-1)+42)    ; 48*(i-1)+8+16  ; 48*(i-1)+41]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+9   ;   -(48*(i-1)+9)   ;   -(48*(i-1)+43)    ; 48*(i-1)+9+16  ; 48*(i-1)+42]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+10 ;   -(48*(i-1)+10) ;   -(48*(i-1)+44)    ; 48*(i-1)+10+16  ; 48*(i-1)+43]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+11 ;   -(48*(i-1)+11) ;   -(48*(i-1)+45)    ; 48*(i-1)+11+16  ; 48*(i-1)+44]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+12 ;   -(48*(i-1)+12) ;   -(48*(i-1)+36)    ; 48*(i-1)+12+16  ; 48*(i-1)+45]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+13 ;   -(48*(i-1)+13) ;   -(48*(i-1)+46)    ; 48*(i-1)+13+16  ; 48*(i-1)+36]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+14 ;   -(48*(i-1)+14) ;   -(48*(i-1)+47)    ; 48*(i-1)+14+16  ; 48*(i-1)+46]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+15 ;   -(48*(i-1)+15) ;   -(48*(i-1)+48)    ; 48*(i-1)+15+16  ; 48*(i-1)+47]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[16*(i-1)+16 ;   -(48*(i-1)+16) ;   -(48*(i-1)+33)    ; 48*(i-1)+16+16  ; 48*(i-1)+48]);
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
    
    %% Physical Lines
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[i      ;48*(i-1)+37]); % Port
    end
%     for i=1:8
%         fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[8+i  ;48*(i-1)+38]); % Inductor right top
%     end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[8+i;48*(i-1)+39]); % Inductor right middle
    end
%     for i=1:8
%         fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[24+i;48*(i-1)+40]); % Inductor right bottom
%     end
%     for i=1:8
%         fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[32+i;48*(i-1)+46]); % Inductor left top
%     end
%     for i=1:8
%         fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[40+i;48*(i-1)+47]); % Inductor left middle
%     end
%     for i=1:8
%         fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[48+i;48*(i-1)+48]); % Inductor left bottom
%     end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[16+i;48*(i-1)+41]); % Capacitor 5
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[24+i;48*(i-1)+42]); % Capacitor 4
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[32+i;48*(i-1)+43]); % Capacitor 3
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[40+i;48*(i-1)+44]); % Capacitor 2
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[48+i;48*(i-1)+45]); % Capacitor 1
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