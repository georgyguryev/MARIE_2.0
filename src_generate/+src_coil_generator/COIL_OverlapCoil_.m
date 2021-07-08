function[mshfile] = COIL_OverlapCoil_(rad1,rad2,over,len,t,meshing,shield)
   
    % Generates an overlapping coil with 8 ports
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    %% Write in .geo file
    
    if shield.flag == 1
        geofile  = strcat('./data/coils/geo_files/Overlap_decoupled_shield_rad_inner_',num2str(rad1*100),'cm_rad_outter_',num2str(rad2*100),'cm_overlap_',num2str(over*100),'cm_len_',num2str(len*100),'cm_mesh_',num2str(meshing),'.geo');
        mshfile = strcat('./data/coils/msh_files/Overlap_decoupled_shield_rad_inner_',num2str(rad1*100),'cm_rad_outter_',num2str(rad2*100),'cm_overlap_',num2str(over*100),'cm_len_',num2str(len*100),'cm_mesh_',num2str(meshing),'.msh');
    else
        geofile  = strcat('./data/coils/geo_files/Overlap_decoupled_rad_inner_',num2str(rad1*100),'cm_rad_outter_',num2str(rad2*100),'cm_overlap_',num2str(over*100),'cm_len_',num2str(len*100),'cm_mesh_',num2str(meshing),'.geo');
        mshfile = strcat('./data/coils/msh_files/Overlap_decoupled_rad_inner_',num2str(rad1*100),'cm_rad_outter_',num2str(rad2*100),'cm_overlap_',num2str(over*100),'cm_len_',num2str(len*100),'cm_mesh_',num2str(meshing),'.msh');
    end
    
    fileID = fopen(geofile,'w');
    fprintf(fileID,'lc = %f; \n\n',meshing);   
    
    %% Create Coil   
    
    %% Points
    dim = [(2*pi*rad1+8*over)/8;(2*pi*rad2+8*over)/8;(2*pi*rad1+8*over)/8;(2*pi*rad2+8*over)/8;(2*pi*rad1+8*over)/8;(2*pi*rad2+8*over)/8;(2*pi*rad1+8*over)/8;(2*pi*rad2+8*over)/8];
    rad = [rad1;rad2;rad1;rad2;rad1;rad2;rad1;rad2];
    offset = dim(1)/2;
    x = [0;dim(2)-over;2*(dim(1)-over);3*(dim(2)-over);4*(dim(1)-over);5*(dim(2)-over);6*(dim(1)-over);7*(dim(2)-over)];
    x = x + offset*0.8;
    
    for i = 1:8
        % Outer
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+1     ; rad(i)*cos((x(i))/rad(i))  ; rad(i)*sin((x(i))/rad(i)) ; len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+2     ; rad(i)*cos((x(i)+dim(i)/3)/rad(i))  ; rad(i)*sin((x(i)+dim(i)/3)/rad(i)) ; len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+3     ; rad(i)*cos((x(i)+2*dim(i)/3)/rad(i))  ; rad(i)*sin((x(i)+2*dim(i)/3)/rad(i)) ; len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+4     ; rad(i)*cos((x(i)+dim(i))/rad(i))  ; rad(i)*sin((x(i)+dim(i))/rad(i)) ; len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+5     ; rad(i)*cos((x(i)+dim(i))/rad(i))  ; rad(i)*sin((x(i)+dim(i))/rad(i)) ; len/2-len/3  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+6     ; rad(i)*cos((x(i)+dim(i))/rad(i))  ; rad(i)*sin((x(i)+dim(i))/rad(i)) ; len/2-2*len/3  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+7     ; rad(i)*cos((x(i)+dim(i))/rad(i))  ; rad(i)*sin((x(i)+dim(i))/rad(i)) ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+8     ; rad(i)*cos((x(i)+3*dim(i)/4)/rad(i))  ; rad(i)*sin((x(i)+3*dim(i)/4)/rad(i)) ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+9     ; rad(i)*cos((x(i)+dim(i)/2)/rad(i))  ; rad(i)*sin((x(i)+dim(i)/2)/rad(i)) ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+10   ; rad(i)*cos((x(i)+dim(i)/4)/rad(i))  ; rad(i)*sin((x(i)+dim(i)/4)/rad(i)) ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+11   ; rad(i)*cos((x(i))/rad(i))  ; rad(i)*sin((x(i))/rad(i)) ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+12   ; rad(i)*cos((x(i))/rad(i))  ; rad(i)*sin((x(i))/rad(i)) ; len/2-2*len/3  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+13   ; rad(i)*cos((x(i))/rad(i))  ; rad(i)*sin((x(i))/rad(i)) ; len/2-len/3  ]);
        % Inner
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+14    ; rad(i)*cos((x(i)+t)/rad(i))  ; rad(i)*sin((x(i)+t)/rad(i)) ; len/2-t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+15    ; rad(i)*cos((x(i)+dim(i)/3)/rad(i))  ; rad(i)*sin((x(i)+dim(i)/3)/rad(i)) ; len/2-t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+16    ; rad(i)*cos((x(i)+2*dim(i)/3)/rad(i))  ; rad(i)*sin((x(i)+2*dim(i)/3)/rad(i)) ; len/2-t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+17    ; rad(i)*cos((x(i)+dim(i)-t)/rad(i))  ; rad(i)*sin((x(i)+dim(i)-t)/rad(i)) ; len/2-t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+18    ; rad(i)*cos((x(i)+dim(i)-t)/rad(i))  ; rad(i)*sin((x(i)+dim(i)-t)/rad(i)) ; len/2-len/3  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+19    ; rad(i)*cos((x(i)+dim(i)-t)/rad(i))  ; rad(i)*sin((x(i)+dim(i)-t)/rad(i)) ; len/2-2*len/3  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+20    ; rad(i)*cos((x(i)+dim(i)-t)/rad(i))  ; rad(i)*sin((x(i)+dim(i)-t)/rad(i)) ; -len/2+t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+21    ; rad(i)*cos((x(i)+3*dim(i)/4)/rad(i))  ; rad(i)*sin((x(i)+3*dim(i)/4)/rad(i)) ; -len/2+t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+22    ; rad(i)*cos((x(i)+dim(i)/2)/rad(i))  ; rad(i)*sin((x(i)+dim(i)/2)/rad(i)) ; -len/2+t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+23    ; rad(i)*cos((x(i)+dim(i)/4)/rad(i))  ; rad(i)*sin((x(i)+dim(i)/4)/rad(i)) ; -len/2+t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+24    ; rad(i)*cos((x(i)+t)/rad(i))  ; rad(i)*sin((x(i)+t)/rad(i)) ; -len/2+t  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+25    ; rad(i)*cos((x(i)+t)/rad(i))  ; rad(i)*sin((x(i)+t)/rad(i)) ; len/2-2*len/3  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+26    ; rad(i)*cos((x(i)+t)/rad(i))  ; rad(i)*sin((x(i)+t)/rad(i)) ; len/2-len/3  ]);
    end
    fprintf(fileID,'Point( 1000 )   = {%f, %f, %f, lc}; \n',[0; 0; len/2]);
    fprintf(fileID,'Point( 2000 )   = {%f, %f, %f, lc}; \n',[0; 0; len/2-t]);
    fprintf(fileID,'Point( 3000 )   = {%f, %f, %f, lc}; \n',[0; 0; t-len/2]);
    fprintf(fileID,'Point( 4000 )   = {%f, %f, %f, lc}; \n',[0; 0; -len/2]);
    %% Circles
    for i = 1:8
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*(i-1)+1  ;26*(i-1)+1; 1000 ;26*(i-1)+2]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*(i-1)+2  ;26*(i-1)+2; 1000 ;26*(i-1)+3]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*(i-1)+3  ;26*(i-1)+3; 1000 ;26*(i-1)+4]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+4  ;26*(i-1)+4; 26*(i-1)+5]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+5  ;26*(i-1)+5  ;26*(i-1)+6]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+6  ;26*(i-1)+6  ;26*(i-1)+7]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*(i-1)+7  ;26*(i-1)+7; 4000  ;26*(i-1)+8]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*(i-1)+8  ;26*(i-1)+8; 4000 ;26*(i-1)+9]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*(i-1)+9  ;26*(i-1)+9; 4000 ;26*(i-1)+10]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*(i-1)+10  ;26*(i-1)+10; 4000 ;26*(i-1)+11]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+11  ;26*(i-1)+11;  26*(i-1)+12]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+12  ;26*(i-1)+12;  26*(i-1)+13]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+13  ;26*(i-1)+13;  26*(i-1)+1]);
        
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*(i-1)+14  ;26*(i-1)+14;  2000 ;26*(i-1)+15]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*(i-1)+15  ;26*(i-1)+15;  2000 ;26*(i-1)+16]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*(i-1)+16  ;26*(i-1)+16;  2000 ;26*(i-1)+17]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+17  ;26*(i-1)+17;  26*(i-1)+18]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+18  ;26*(i-1)+18;  26*(i-1)+19]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+19  ;26*(i-1)+19;  26*(i-1)+20]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*(i-1)+20  ;26*(i-1)+20;  3000 ;26*(i-1)+21]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*(i-1)+21  ;26*(i-1)+21;  3000 ;26*(i-1)+22]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*(i-1)+22  ;26*(i-1)+22;  3000 ;26*(i-1)+23]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*(i-1)+23  ;26*(i-1)+23;  3000 ;26*(i-1)+24]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+24  ;26*(i-1)+24;  26*(i-1)+25]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+25  ;26*(i-1)+25;  26*(i-1)+26]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+26  ;26*(i-1)+26;  26*(i-1)+14]);
        
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+27  ;26*(i-1)+2  ;26*(i-1)+15]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+28  ;26*(i-1)+3  ;26*(i-1)+16]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+29  ;26*(i-1)+5  ;26*(i-1)+18]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+30  ;26*(i-1)+6  ;26*(i-1)+19]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+31  ;26*(i-1)+8  ;26*(i-1)+21]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+32  ;26*(i-1)+9  ;26*(i-1)+22]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+33  ;26*(i-1)+10  ;26*(i-1)+23]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+34  ;26*(i-1)+12  ;26*(i-1)+25]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+35  ;26*(i-1)+13  ;26*(i-1)+26]);
        
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+36  ;26*(i-1)+1  ;26*(i-1)+14]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+37  ;26*(i-1)+4  ;26*(i-1)+17]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+38  ;26*(i-1)+7  ;26*(i-1)+20]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*(i-1)+39  ;26*(i-1)+11  ;26*(i-1)+24]);
        
    end
    
    %% Line Loops
    for i = 1:8
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[13*(i-1)+1   ;   -(39*(i-1)+1)   ;   -(39*(i-1)+27)    ; 39*(i-1)+14  ; 39*(i-1)+36]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[13*(i-1)+2   ;   -(39*(i-1)+2)   ;   -(39*(i-1)+28)    ; 39*(i-1)+15  ; 39*(i-1)+27]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[13*(i-1)+3   ;   -(39*(i-1)+3)   ;   -(39*(i-1)+37)    ; 39*(i-1)+16  ; 39*(i-1)+28]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[13*(i-1)+4   ;   -(39*(i-1)+4)   ;   -(39*(i-1)+29)    ; 39*(i-1)+17  ; 39*(i-1)+37]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[13*(i-1)+5   ;   -(39*(i-1)+5)   ;   -(39*(i-1)+30)    ; 39*(i-1)+18  ; 39*(i-1)+29]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[13*(i-1)+6   ;   -(39*(i-1)+6)   ;   -(39*(i-1)+38)    ; 39*(i-1)+19  ; 39*(i-1)+30]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[13*(i-1)+7   ;   -(39*(i-1)+7)   ;   -(39*(i-1)+31)    ; 39*(i-1)+20  ; 39*(i-1)+38]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[13*(i-1)+8   ;   -(39*(i-1)+8)   ;   -(39*(i-1)+32)    ; 39*(i-1)+21  ; 39*(i-1)+31]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[13*(i-1)+9   ;   -(39*(i-1)+9)   ;   -(39*(i-1)+33)    ; 39*(i-1)+22  ; 39*(i-1)+32]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[13*(i-1)+10 ;   -(39*(i-1)+10) ;   -(39*(i-1)+39)    ; 39*(i-1)+23  ; 39*(i-1)+33]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[13*(i-1)+11 ;   -(39*(i-1)+11) ;   -(39*(i-1)+34)    ; 39*(i-1)+24  ; 39*(i-1)+39]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[13*(i-1)+12 ;   -(39*(i-1)+12) ;   -(39*(i-1)+35)    ; 39*(i-1)+25  ; 39*(i-1)+34]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[13*(i-1)+13 ;   -(39*(i-1)+13) ;   -(39*(i-1)+36)    ; 39*(i-1)+26  ; 39*(i-1)+35]);
    end
    
    
    %% Surfaces
    for i=1:8
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[13*(i-1)+1  ;13*(i-1)+1  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[13*(i-1)+2  ;13*(i-1)+2  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[13*(i-1)+3  ;13*(i-1)+3  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[13*(i-1)+4  ;13*(i-1)+4  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[13*(i-1)+5  ;13*(i-1)+5  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[13*(i-1)+6  ;13*(i-1)+6  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[13*(i-1)+7  ;13*(i-1)+7  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[13*(i-1)+8  ;13*(i-1)+8  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[13*(i-1)+9  ;13*(i-1)+9  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[13*(i-1)+10  ;13*(i-1)+10  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[13*(i-1)+11  ;13*(i-1)+11  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[13*(i-1)+12  ;13*(i-1)+12  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[13*(i-1)+13  ;13*(i-1)+13  ]);
    end
    
    %% Physical Lines
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[i      ;39*(i-1)+32]); % Port
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[8+i  ;39*(i-1)+33]); % Capacitor left of port
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[16+i;39*(i-1)+31]); % Capacitor right of port
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[24+i;39*(i-1)+27]); % Capacitor top left
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[32+i;39*(i-1)+28]); % Capacitor top right
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[40+i;39*(i-1)+35]); % Capacitor left top
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[48+i;39*(i-1)+34]); % Capacitor left bottom
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[56+i;39*(i-1)+29]); % Capacitor right top
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[64+i;39*(i-1)+30]); % Capacitor right bottom
    end
    
    
    if shield.flag == 1
        fprintf(fileID,'lc_s = %f; \n\n',shield.meshing);   
        
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[26*8+1  ; 0    ; 0    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[26*8+2  ; 0    ; 0    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[26*8+3  ; shield.radius    ; 0    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[26*8+4  ; 0    ; shield.radius    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[26*8+5  ; -shield.radius    ; 0    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[26*8+6  ; 0    ; -shield.radius    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[26*8+7  ; shield.radius    ; 0    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[26*8+8  ; 0    ; shield.radius    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[26*8+9  ; -shield.radius    ; 0    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[26*8+10  ; 0    ; -shield.radius    ; -shield.length/2]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*8+1; 26*8+3;26*8+7 ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*8+2; 26*8+4;26*8+8 ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*8+3; 26*8+5;26*8+9 ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[39*8+4; 26*8+6;26*8+10 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*8+5; 26*8+3;26*8+1;26*8+4 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*8+6; 26*8+4;26*8+1;26*8+5 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*8+7; 26*8+5;26*8+1;26*8+6 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*8+8; 26*8+6;26*8+1;26*8+3 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*8+9; 26*8+7;26*8+2;26*8+8 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*8+10; 26*8+8;26*8+2;26*8+9 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*8+11; 26*8+9;26*8+2;26*8+10 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[39*8+12; 26*8+10;26*8+2;26*8+7 ]);
        
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[13*8+1;   39*8+1 ;   39*8+9    ; -(39*8+2); -(39*8+5)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[13*8+2;   39*8+2 ;   39*8+10    ; -(39*8+3); -(39*8+6)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[13*8+3;   39*8+3 ;   39*8+11    ; -(39*8+4); -(39*8+7)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[13*8+4;   39*8+4 ;   39*8+12    ; -(39*8+1); -(39*8+8)]);
        
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[13*8+1  ;13*8+1  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[13*8+2  ;13*8+2  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[13*8+3  ;13*8+3  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[13*8+4  ;13*8+4  ]);
        
        fprintf(fileID,'Physical Surface(2)  = {13*8+1,13*8+2,13*8+3,13*8+4};\n\n');
        
    end
    
    %% Physical Surfaces
    fprintf(fileID,'Physical Surface(1)  = {');
    for i=1:8
        fprintf(fileID,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d',[...
            13*(i-1)+1,13*(i-1)+2,13*(i-1)+3,13*(i-1)+4,13*(i-1)+5,13*(i-1)+6,13*(i-1)+7,13*(i-1)+8,...
            13*(i-1)+9,13*(i-1)+10,13*(i-1)+11,13*(i-1)+12,13*(i-1)+13]);
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