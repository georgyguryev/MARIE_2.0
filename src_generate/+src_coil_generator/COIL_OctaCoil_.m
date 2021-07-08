function[mshfile] = COIL_OctaCoil_(a,b,len,t,meshing,shield)
   
    % Generates an octagonal elliptical birdcage coil with 8 ports
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    %% Write in .geo file
    
    if shield.flag == 1
        geofile  = strcat('./data/coils/geo_files/Octa_shield_decoupled_rad_x_',num2str(a*100),'cm_rad_y_',num2str(b*100),'cm_len_',num2str(len*100),'cm_mesh_',num2str(meshing),'.geo');
        mshfile = strcat('./data/coils/msh_files/Octa_shield_decoupled_rad_x_',num2str(a*100),'cm_rad_y_',num2str(b*100),'cm_len_',num2str(len*100),'cm_mesh_',num2str(meshing),'.msh');
    else
        geofile  = strcat('./data/coils/geo_files/Octa_decoupled_rad_x_',num2str(a*100),'cm_rad_y_',num2str(b*100),'cm_len_',num2str(len*100),'cm_mesh_',num2str(meshing),'.geo');
        mshfile = strcat('./data/coils/msh_files/Octa_decoupled_rad_x_',num2str(a*100),'cm_rad_y_',num2str(b*100),'cm_len_',num2str(len*100),'cm_mesh_',num2str(meshing),'.msh');
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
        [x_off1,x_off2,y_off1,y_off2] = get_points_line(mid_x,mid_y,pos(j,1),pos(j,2),t/2);
        if x_off1<mid_x && x_off1>pos(i,1)
            mid_x1 = (pos(i,1)+x_off1)/2;
            mid_y1 = (pos(i,2)+y_off1)/2;
            mid_x2 = (pos(j,1)+x_off2)/2;
            mid_y2 = (pos(j,2)+y_off2)/2;
        else
            mid_x2 = (pos(i,1)+x_off2)/2;
            mid_y2 = (pos(i,2)+y_off2)/2;
            mid_x1 = (pos(j,1)+x_off1)/2;
            mid_y1 = (pos(j,2)+y_off1)/2;
        end
        mid_z1 = (len/2-t)/2;
        mid_z2 = (t-len/2)/2;
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+1    ; pos(i,1)  ; pos(i,2) ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+2    ; pos(i,1)  ; pos(i,2) ; len/2   ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+3    ; pos(i,1)  ; pos(i,2) ; t-len/2 ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+4    ; pos(i,1)  ; pos(i,2) ; len/2-t ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+5    ; x_off1    ; y_off1    ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+6    ; x_off1    ; y_off1    ; len/2   ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+7    ; x_off1    ; y_off1    ; t-len/2 ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+8    ; x_off1    ; y_off1    ; len/2-t ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+9    ; x_off2    ; y_off2    ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+10  ; x_off2    ; y_off2    ; len/2   ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+11  ; x_off2    ; y_off2    ; t-len/2 ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+12  ; x_off2    ; y_off2    ; len/2-t ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+13  ; mid_x1  ; mid_y1  ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+14  ; mid_x1  ; mid_y1  ; len/2   ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+15  ; mid_x1  ; mid_y1  ; t-len/2 ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+16  ; mid_x1  ; mid_y1  ; len/2-t ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+17  ; mid_x2  ; mid_y2  ; -len/2  ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+18  ; mid_x2  ; mid_y2  ; len/2   ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+19  ; mid_x2  ; mid_y2  ; t-len/2 ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+20  ; mid_x2  ; mid_y2  ; len/2-t ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+21  ; x_off1    ; y_off1    ; 0        ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+22  ; x_off2    ; y_off2    ; 0        ]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+23  ; x_off1    ; y_off1    ; mid_z1]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+24  ; x_off2    ; y_off2    ; mid_z1]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+25  ; x_off1    ; y_off1    ; mid_z2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc}; \n',[26*(i-1)+26  ; x_off2    ; y_off2    ; mid_z2]);
    end
    %% Lines
    for i = 1:8
        j = mod(i,8)+1;
        mid_x = (pos(i,1)+pos(j,1))/2;
        mid_y = (pos(i,2)+pos(j,2))/2;
        [x_off1,~,~,~] = get_points_line(mid_x,mid_y,pos(j,1),pos(j,2),t/2);
        if x_off1<mid_x && x_off1>pos(i,1)
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+1  ;26*(i-1)+1  ;26*(i-1)+13]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+2  ;26*(i-1)+2  ;26*(i-1)+14]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+3  ;26*(i-1)+3  ;26*(i-1)+15]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+4  ;26*(i-1)+4  ;26*(i-1)+16]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+5  ;26*(i-1)+13;26*(i-1)+5  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+6  ;26*(i-1)+14;26*(i-1)+6  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+7  ;26*(i-1)+15;26*(i-1)+7  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+8  ;26*(i-1)+16;26*(i-1)+8  ]);   
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+9  ;26*(i-1)+5  ;26*(i-1)+9  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+10;26*(i-1)+6  ;26*(i-1)+10]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+11;26*(i-1)+7  ;26*(i-1)+11]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+12;26*(i-1)+8  ;26*(i-1)+12]); 
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+13;26*(i-1)+17;26*(i-1)+9  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+14;26*(i-1)+18;26*(i-1)+10]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+15;26*(i-1)+19;26*(i-1)+11]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+16;26*(i-1)+20;26*(i-1)+12]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+17;26*(i-1)+17;26*(j-1)+1  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+18;26*(i-1)+18;26*(j-1)+2  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+19;26*(i-1)+19;26*(j-1)+3  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+20;26*(i-1)+20;26*(j-1)+4  ]);  
        else
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+1  ;26*(i-1)+1  ;26*(i-1)+17]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+2  ;26*(i-1)+2  ;26*(i-1)+18]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+3  ;26*(i-1)+3  ;26*(i-1)+19]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+4  ;26*(i-1)+4  ;26*(i-1)+20]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+5  ;26*(i-1)+17;26*(i-1)+9  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+6  ;26*(i-1)+18;26*(i-1)+10]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+7  ;26*(i-1)+19;26*(i-1)+11]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+8  ;26*(i-1)+20;26*(i-1)+12]);   
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+9  ;26*(i-1)+9  ;26*(i-1)+5  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+10;26*(i-1)+10;26*(i-1)+6  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+11;26*(i-1)+11;26*(i-1)+7  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+12;26*(i-1)+12;26*(i-1)+8  ]); 
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+13;26*(i-1)+13;26*(i-1)+5  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+14;26*(i-1)+14;26*(i-1)+6  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+15;26*(i-1)+15;26*(i-1)+7  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+16;26*(i-1)+16;26*(i-1)+8  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+17;26*(i-1)+13;26*(j-1)+1  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+18;26*(i-1)+14;26*(j-1)+2  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+19;26*(i-1)+15;26*(j-1)+3  ]);
            fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+20;26*(i-1)+16;26*(j-1)+4  ]); 
        end
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+21;26*(i-1)+1  ;26*(i-1)+3  ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+22;26*(i-1)+2  ;26*(i-1)+4  ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+23;26*(i-1)+13;26*(i-1)+15]); % Physical Line: Capacitor bottom
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+24;26*(i-1)+14;26*(i-1)+16]); % Physical Line: Capacitor top
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+25;26*(i-1)+17;26*(i-1)+19]); % Physical Line: Capacitor bottom
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+26;26*(i-1)+18;26*(i-1)+20]); % Physical Line: Capacitor top
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+27;26*(i-1)+8  ;26*(i-1)+23]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+28;26*(i-1)+23;26*(i-1)+21]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+29;26*(i-1)+21;26*(i-1)+25]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+30;26*(i-1)+25;26*(i-1)+7  ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+31;26*(i-1)+12;26*(i-1)+24]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+32;26*(i-1)+24;26*(i-1)+22]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+33;26*(i-1)+22;26*(i-1)+26]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+34;26*(i-1)+26;26*(i-1)+11]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+35;26*(i-1)+21;26*(i-1)+22]); % Physical Line: Port
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+36;26*(i-1)+23;26*(i-1)+24]); % Physical Line: Capacitor next to port top
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+37;26*(i-1)+25;26*(i-1)+26]); % Physical Line: Capacitor next to port bottom
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+38;26*(i-1)+5  ;26*(i-1)+7  ]); 
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+39;26*(i-1)+6  ;26*(i-1)+8  ]); 
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+40;26*(i-1)+9  ;26*(i-1)+11]); 
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*(i-1)+41;26*(i-1)+10;26*(i-1)+12]); 
    end
    
    %% Line Loops
    for i = 1:8
        j = mod(i,8)+1;
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[14*(i-1)+1  ;     41*(i-1)+1      ;    41*(i-1)+25  ; -(41*(i-1)+3)   ; -(41*(i-1)+21)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[14*(i-1)+2  ;   -(41*(i-1)+2)     ;  -(41*(i-1)+26) ;   41*(i-1)+4    ;   41*(i-1)+22 ]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[14*(i-1)+3  ;   -(41*(i-1)+25)   ;  -(41*(i-1)+7)   ;   41*(i-1)+40  ;   41*(i-1)+5   ]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[14*(i-1)+4  ;     41*(i-1)+26    ;    41*(i-1)+8    ; -(41*(i-1)+41) ; -(41*(i-1)+6)  ]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[14*(i-1)+5  ;     41*(i-1)+38    ;    41*(i-1)+9    ; -(41*(i-1)+40) ; -(41*(i-1)+11) ]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[14*(i-1)+6  ;   -(41*(i-1)+39)   ;  -(41*(i-1)+10) ;   41*(i-1)+41  ;   41*(i-1)+12  ]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[14*(i-1)+7  ;   -(41*(i-1)+38)   ;  -(41*(i-1)+13) ;   41*(i-1)+23  ;   41*(i-1)+15  ]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[14*(i-1)+8  ;     41*(i-1)+39    ;  -(41*(i-1)+16) ; -(41*(i-1)+24) ;   41*(i-1)+14  ]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[14*(i-1)+9  ;   -(41*(i-1)+23)   ;  -(41*(i-1)+19) ;   41*(j-1)+21  ;   41*(i-1)+17  ]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[14*(i-1)+10;     41*(i-1)+24    ;    41*(i-1)+20  ; -(41*(j-1)+22) ; -(41*(i-1)+18) ]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[14*(i-1)+11;     41*(i-1)+11    ;    41*(i-1)+34  ;   41*(i-1)+37  ; -(41*(i-1)+30) ]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[14*(i-1)+12;   -(41*(i-1)+12)   ;   41*(i-1)+31   ; -(41*(i-1)+36) ; -(41*(i-1)+27) ]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[14*(i-1)+13;   -(41*(i-1)+37)   ;  -(41*(i-1)+29) ;   41*(i-1)+35  ;   41*(i-1)+33  ]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[14*(i-1)+14;     41*(i-1)+36    ;    41*(i-1)+32  ; -(41*(i-1)+35) ; -(41*(i-1)+28) ]);
    end
    
    %% Surfaces
    for i=1:8
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[14*(i-1)+1  ;14*(i-1)+1  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[14*(i-1)+2  ;14*(i-1)+2  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[14*(i-1)+3  ;14*(i-1)+3  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[14*(i-1)+4  ;14*(i-1)+4  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[14*(i-1)+5  ;14*(i-1)+5  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[14*(i-1)+6  ;14*(i-1)+6  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[14*(i-1)+7  ;14*(i-1)+7  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[14*(i-1)+8  ;14*(i-1)+8  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[14*(i-1)+9  ;14*(i-1)+9  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[14*(i-1)+10;14*(i-1)+10]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[14*(i-1)+11;14*(i-1)+11]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[14*(i-1)+12;14*(i-1)+12]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[14*(i-1)+13;14*(i-1)+13]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[14*(i-1)+14;14*(i-1)+14]);
    end
    
    %% Physical Lines
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[i      ;41*(i-1)+35]); % Port
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[8+i  ;41*(i-1)+36]); % Top capacitor of port
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[16+i;41*(i-1)+37]); % Bottom capacitor of port
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[24+i;41*(i-1)+24]); % Capacitor top
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[32+i;41*(i-1)+26]); % Capacitor top
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[40+i;41*(i-1)+23]); % Capacitor bottom
    end
    for i=1:8
        fprintf(fileID,'Physical Line( %d ) = {%d}; \n',[48+i;41*(i-1)+25]); % Capacitor bottom
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
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*8+1; 26*8+3;26*8+7 ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*8+2; 26*8+4;26*8+8 ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*8+3; 26*8+5;26*8+9 ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[41*8+4; 26*8+6;26*8+10 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[41*8+5; 26*8+3;26*8+1;26*8+4 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[41*8+6; 26*8+4;26*8+1;26*8+5 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[41*8+7; 26*8+5;26*8+1;26*8+6 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[41*8+8; 26*8+6;26*8+1;26*8+3 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[41*8+9; 26*8+7;26*8+2;26*8+8 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[41*8+10; 26*8+8;26*8+2;26*8+9 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[41*8+11; 26*8+9;26*8+2;26*8+10 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[41*8+12; 26*8+10;26*8+2;26*8+7 ]);
        
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[14*8+1;   41*8+1 ;   41*8+9    ; -(41*8+2); -(41*8+5)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[14*8+2;   41*8+2 ;   41*8+10    ; -(41*8+3); -(41*8+6)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[14*8+3;   41*8+3 ;   41*8+11    ; -(41*8+4); -(41*8+7)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[14*8+4;   41*8+4 ;   41*8+12    ; -(41*8+1); -(41*8+8)]);
        
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[14*8+1  ;14*8+1  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[14*8+2  ;14*8+2  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[14*8+3  ;14*8+3  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[14*8+4  ;14*8+4  ]);
        
        fprintf(fileID,'Physical Surface(2)  = {14*8+1,14*8+2,14*8+3,14*8+4};\n\n');
        
    end
    
    %% Physical Surfaces
    fprintf(fileID,'Physical Surface(1)  = {');
    for i=1:8
        fprintf(fileID,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d',[...
            14*(i-1)+1,14*(i-1)+2,14*(i-1)+3,14*(i-1)+4,14*(i-1)+5,14*(i-1)+6,14*(i-1)+7,14*(i-1)+8,...
            14*(i-1)+9,14*(i-1)+10,14*(i-1)+11,14*(i-1)+12,14*(i-1)+13,14*(i-1)+14]);
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
    
    m = (y2-y1)/(x2-x1);
    temp = sqrt(t^2/(m^2+1));
    if temp+x1<x2 && temp+x1>x1
        xout_1 = temp+x1;
        xout_2 = -temp+x1;
        yout_1 = m*(xout_1-x1)+y1;
        yout_2 = m*(xout_2-x1)+y1;
    else
        xout_1 = -temp+x1;
        xout_2 = temp+x1;
        yout_1 = m*(xout_1-x1)+y1;
        yout_2 = m*(xout_2-x1)+y1;
    end
    
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