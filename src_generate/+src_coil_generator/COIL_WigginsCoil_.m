function[mshfile] = COIL_WigginsCoil_(r,w,lex,t,meshing,shield)

    % Generates a Wiggins stadium coil with 8 ports
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    if shield.flag ==1
        geofile  = strcat('./data/coils/geo_files/Wiggins_decoupled_shield_rad_',num2str(r*100),'cm_wide_',num2str(w*100),'cm_len_',num2str(lex*100),'cm_width_',num2str(t*100),'cm_mesh_',num2str(meshing),'.geo');
        mshfile = strcat('./data/coils/msh_files/Wiggins_decoupled_shield_rad_',num2str(r*100),'cm_wide_',num2str(w*100),'cm_len_',num2str(lex*100),'cm_width_',num2str(t*100),'cm_mesh_',num2str(meshing),'.msh');
    else
        geofile  = strcat('./data/coils/geo_files/Wiggins_decoupled_rad_',num2str(r*100),'cm_wide_',num2str(w*100),'cm_len_',num2str(lex*100),'cm_width_',num2str(t*100),'cm_mesh_',num2str(meshing),'.geo');
        mshfile = strcat('./data/coils/msh_files/Wiggins_decoupled_rad_',num2str(r*100),'cm_wide_',num2str(w*100),'cm_len_',num2str(lex*100),'cm_width_',num2str(t*100),'cm_mesh_',num2str(meshing),'.msh');
    end
    
    fileID = fopen(geofile,'w');
    fprintf(fileID,'lc = %f; \n\n',meshing);   
    
    L = 20;
    L1 = ceil(2*L*(pi*r)/(pi*r+w));
    L2 = 2*L - L1;
    L1 = L - L2;
    lin = lex - t;
   
    %% Points
    % Centers
    fprintf(fileID,'Point( 1 ) = {%f, %f, %f, lc}; \n',[w/2 ; 0 ; lex/2]);
    fprintf(fileID,'Point( 2 ) = {%f, %f, %f, lc}; \n',[w/2 ; 0 ; lin/2]);
    fprintf(fileID,'Point( 3 ) = {%f, %f, %f, lc}; \n',[w/2 ; 0 ; -lin/2]);
    fprintf(fileID,'Point( 4 ) = {%f, %f, %f, lc}; \n',[w/2 ; 0 ; -lex/2]);
    fprintf(fileID,'Point( 5 ) = {%f, %f, %f, lc}; \n',[-w/2 ; 0 ; lex/2]);
    fprintf(fileID,'Point( 6 ) = {%f, %f, %f, lc}; \n',[-w/2 ; 0 ; lin/2]);
    fprintf(fileID,'Point( 7 ) = {%f, %f, %f, lc}; \n',[-w/2 ; 0 ; -lin/2]);
    fprintf(fileID,'Point( 8 ) = {%f, %f, %f, lc}; \n',[-w/2 ; 0 ; -lex/2]);
    
    % Bottom Ports
    fprintf(fileID,'Point( 9 ) = {%f, %f, %f, lc}; \n',[r+w/2 ; 0 ; -lex/2]);
    fprintf(fileID,'Point( 10 ) = {%f, %f, %f, lc}; \n',[r+w/2 ; 0 ; -lin/2]);
    fprintf(fileID,'Point( 11 ) = {%f, %f, %f, lc}; \n',[0 ; r ; -lex/2]);
    fprintf(fileID,'Point( 12 ) = {%f, %f, %f, lc}; \n',[0 ; r ; -lin/2]);
    fprintf(fileID,'Point( 13 ) = {%f, %f, %f, lc}; \n',[-r-w/2 ; 0 ; -lex/2]);
    fprintf(fileID,'Point( 14 ) = {%f, %f, %f, lc}; \n',[-r-w/2 ; 0 ; -lin/2]);
    fprintf(fileID,'Point( 15 ) = {%f, %f, %f, lc}; \n',[0 ; -r ; -lex/2]);
    fprintf(fileID,'Point( 16 ) = {%f, %f, %f, lc}; \n',[0 ; -r ; -lin/2]);
    
    theta = 5/18*pi;
    % Top Ports
    fprintf(fileID,'Point( 17 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; r*sin(theta)  ; lex/2]);
    fprintf(fileID,'Point( 18 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; r*sin(theta)  ; lin/2]);
    fprintf(fileID,'Point( 19 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; r*sin(theta); lex/2]);
    fprintf(fileID,'Point( 20 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; r*sin(theta); lin/2]);
    fprintf(fileID,'Point( 21 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; -r*sin(theta)  ; lex/2]);
    fprintf(fileID,'Point( 22 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; -r*sin(theta)  ; lin/2]);
    fprintf(fileID,'Point( 23 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; -r*sin(theta); lex/2]);
    fprintf(fileID,'Point( 24 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; -r*sin(theta)  ; lin/2]);
    
    tweak = 0.01315;
    % Points for mutual inductors top
    [x1,y1,x2,y2] = Tweak(r,0,r,tweak);    
    fprintf(fileID,'Point( 25 ) = {%f, %f, %f, lc}; \n',[x1+w/2 ; y1; lex/2]);
    fprintf(fileID,'Point( 26 ) = {%f, %f, %f, lc}; \n',[x1+w/2 ; y1 ; lin/2]);
    fprintf(fileID,'Point( 27 ) = {%f, %f, %f, lc}; \n',[x2+w/2 ; y2 ; lex/2]);
    fprintf(fileID,'Point( 28 ) = {%f, %f, %f, lc}; \n',[x2+w/2 ; y2 ; lin/2]);
    [x1,y1,x2,y2] = Tweak(0,r,r,tweak);    
    fprintf(fileID,'Point( 29 ) = {%f, %f, %f, lc}; \n',[x1 ; y1; lex/2]);
    fprintf(fileID,'Point( 30 ) = {%f, %f, %f, lc}; \n',[x1 ; y1 ; lin/2]);
    fprintf(fileID,'Point( 31 ) = {%f, %f, %f, lc}; \n',[x2 ; y2 ; lex/2]);
    fprintf(fileID,'Point( 32 ) = {%f, %f, %f, lc}; \n',[x2 ; y2 ; lin/2]);
    [x1,y1,x2,y2] = Tweak(-r,0,r,tweak);    
    fprintf(fileID,'Point( 33 ) = {%f, %f, %f, lc}; \n',[-x1-w/2 ; y1; lex/2]);
    fprintf(fileID,'Point( 34 ) = {%f, %f, %f, lc}; \n',[-x1-w/2 ; y1 ; lin/2]);
    fprintf(fileID,'Point( 35 ) = {%f, %f, %f, lc}; \n',[-x2-w/2 ; y2 ; lex/2]);
    fprintf(fileID,'Point( 36 ) = {%f, %f, %f, lc}; \n',[-x2-w/2 ; y2 ; lin/2]);
    [x1,y1,x2,y2] = Tweak(0,-r,r,tweak);    
    fprintf(fileID,'Point( 37 ) = {%f, %f, %f, lc}; \n',[x1 ; y1; lex/2]);
    fprintf(fileID,'Point( 38 ) = {%f, %f, %f, lc}; \n',[x1 ; y1 ; lin/2]);
    fprintf(fileID,'Point( 39 ) = {%f, %f, %f, lc}; \n',[x2 ; y2 ; lex/2]);
    fprintf(fileID,'Point( 40 ) = {%f, %f, %f, lc}; \n',[x2 ; y2 ; lin/2]);
    
    theta = 5/18*pi;
    % Points for mutual inductors bottom
    [x1,y1,x2,y2] = Tweak(r*cos(theta),r*sin(theta),r,tweak);    
    fprintf(fileID,'Point( 41 ) = {%f, %f, %f, lc}; \n',[x1+w/2 ; y1; -lex/2]);
    fprintf(fileID,'Point( 42 ) = {%f, %f, %f, lc}; \n',[x1+w/2 ; y1 ; -lin/2]);
    fprintf(fileID,'Point( 43 ) = {%f, %f, %f, lc}; \n',[x2+w/2 ; y2 ; -lex/2]);
    fprintf(fileID,'Point( 44 ) = {%f, %f, %f, lc}; \n',[x2+w/2 ; y2 ; -lin/2]);
    [x1,y1,x2,y2] = Tweak(-r*cos(theta),r*sin(theta),r,tweak); 
    fprintf(fileID,'Point( 45 ) = {%f, %f, %f, lc}; \n',[-x1-w/2  ; y1; -lex/2]);
    fprintf(fileID,'Point( 46 ) = {%f, %f, %f, lc}; \n',[-x1-w/2  ; y1 ; -lin/2]);
    fprintf(fileID,'Point( 47 ) = {%f, %f, %f, lc}; \n',[-x2-w/2  ; y2 ; -lex/2]);
    fprintf(fileID,'Point( 48 ) = {%f, %f, %f, lc}; \n',[-x2-w/2  ; y2 ; -lin/2]);
    [x1,y1,x2,y2] = Tweak(-r*cos(theta),-r*sin(theta),r,tweak); 
    fprintf(fileID,'Point( 49 ) = {%f, %f, %f, lc}; \n',[-x1-w/2 ; y1; -lex/2]);
    fprintf(fileID,'Point( 50 ) = {%f, %f, %f, lc}; \n',[-x1-w/2 ; y1 ; -lin/2]);
    fprintf(fileID,'Point( 51 ) = {%f, %f, %f, lc}; \n',[-x2-w/2 ; y2 ; -lex/2]);
    fprintf(fileID,'Point( 52 ) = {%f, %f, %f, lc}; \n',[-x2-w/2 ; y2 ; -lin/2]);
    [x1,y1,x2,y2] = Tweak(r*cos(theta),-r*sin(theta),r,tweak); 
    fprintf(fileID,'Point( 53 ) = {%f, %f, %f, lc}; \n',[x1+w/2 ; y1; -lex/2]);
    fprintf(fileID,'Point( 54 ) = {%f, %f, %f, lc}; \n',[x1+w/2 ; y1 ; -lin/2]);
    fprintf(fileID,'Point( 55 ) = {%f, %f, %f, lc}; \n',[x2+w/2 ; y2 ; -lex/2]);
    fprintf(fileID,'Point( 56 ) = {%f, %f, %f, lc}; \n',[x2+w/2 ; y2 ; -lin/2]);
    
    theta = pi/6;
    % Capacitors on top ring
    fprintf(fileID,'Point( 57 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; r*sin(theta)  ; lex/2]);
    fprintf(fileID,'Point( 58 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; r*sin(theta)  ; lin/2]);
    fprintf(fileID,'Point( 59 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; r*sin(theta); lex/2]);
    fprintf(fileID,'Point( 60 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; r*sin(theta); lin/2]);
    fprintf(fileID,'Point( 61 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; -r*sin(theta)  ; lex/2]);
    fprintf(fileID,'Point( 62 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; -r*sin(theta)  ; lin/2]);
    fprintf(fileID,'Point( 63 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; -r*sin(theta)  ; lex/2]);
    fprintf(fileID,'Point( 64 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; -r*sin(theta)  ; lin/2]);
    theta = 5*pi/12;
    fprintf(fileID,'Point( 65 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; r*sin(theta)  ; lex/2]);
    fprintf(fileID,'Point( 66 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; r*sin(theta)  ; lin/2]);
    fprintf(fileID,'Point( 67 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; r*sin(theta); lex/2]);
    fprintf(fileID,'Point( 68 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; r*sin(theta); lin/2]);
    fprintf(fileID,'Point( 69 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; -r*sin(theta)  ; lex/2]);
    fprintf(fileID,'Point( 70 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; -r*sin(theta)  ; lin/2]);
    fprintf(fileID,'Point( 71 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; -r*sin(theta)  ; lex/2]);
    fprintf(fileID,'Point( 72 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; -r*sin(theta)  ; lin/2]);
    
    theta = pi/6;
    % Capacitors on bottom ring
    fprintf(fileID,'Point( 73 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; r*sin(theta)  ; -lex/2]);
    fprintf(fileID,'Point( 74 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; r*sin(theta)  ; -lin/2]);
    fprintf(fileID,'Point( 75 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; r*sin(theta); -lex/2]);
    fprintf(fileID,'Point( 76 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; r*sin(theta); -lin/2]);
    fprintf(fileID,'Point( 77 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; -r*sin(theta)  ; -lex/2]);
    fprintf(fileID,'Point( 78 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; -r*sin(theta)  ; -lin/2]);
    fprintf(fileID,'Point( 79 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; -r*sin(theta)   ; -lex/2]);
    fprintf(fileID,'Point( 80 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; -r*sin(theta)  ; -lin/2]);
    theta = 5*pi/12;
    fprintf(fileID,'Point( 81 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; r*sin(theta)  ; -lex/2]);
    fprintf(fileID,'Point( 82 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; r*sin(theta)  ; -lin/2]);
    fprintf(fileID,'Point( 83 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; r*sin(theta); -lex/2]);
    fprintf(fileID,'Point( 84 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; r*sin(theta); -lin/2]);
    fprintf(fileID,'Point( 85 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; -r*sin(theta)  ; -lex/2]);
    fprintf(fileID,'Point( 86 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; -r*sin(theta)  ; -lin/2]);
    fprintf(fileID,'Point( 87 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; -r*sin(theta)   ; -lex/2]);
    fprintf(fileID,'Point( 88 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; -r*sin(theta)  ; -lin/2]);
    
    % Break points from Circle to Lines
    fprintf(fileID,'Point( 89 ) = {%f, %f, %f, lc}; \n',[w/2; r ; -lex/2]);
    fprintf(fileID,'Point( 90 ) = {%f, %f, %f, lc}; \n',[w/2; r ; -lin/2]);
    fprintf(fileID,'Point( 91 ) = {%f, %f, %f, lc}; \n',[w/2; r ; lex/2]);
    fprintf(fileID,'Point( 92 ) = {%f, %f, %f, lc}; \n',[w/2; r ; lin/2]);
    fprintf(fileID,'Point( 93 ) = {%f, %f, %f, lc}; \n',[-w/2; r ; -lex/2]);
    fprintf(fileID,'Point( 94 ) = {%f, %f, %f, lc}; \n',[-w/2; r ; -lin/2]);
    fprintf(fileID,'Point( 95 ) = {%f, %f, %f, lc}; \n',[-w/2; r ; lex/2]);
    fprintf(fileID,'Point( 96 ) = {%f, %f, %f, lc}; \n',[-w/2; r ; lin/2]);
    fprintf(fileID,'Point( 97 ) = {%f, %f, %f, lc}; \n',[-w/2; -r ; -lex/2]);
    fprintf(fileID,'Point( 98 ) = {%f, %f, %f, lc}; \n',[-w/2; -r ; -lin/2]);
    fprintf(fileID,'Point( 99 ) = {%f, %f, %f, lc}; \n',[-w/2; -r ; lex/2]);
    fprintf(fileID,'Point( 100 ) = {%f, %f, %f, lc}; \n',[-w/2; -r ; lin/2]);
    fprintf(fileID,'Point( 101 ) = {%f, %f, %f, lc}; \n',[w/2; -r ; -lex/2]);
    fprintf(fileID,'Point( 102 ) = {%f, %f, %f, lc}; \n',[w/2; -r ; -lin/2]);
    fprintf(fileID,'Point( 103 ) = {%f, %f, %f, lc}; \n',[w/2; -r ; lex/2]);
    fprintf(fileID,'Point( 104 ) = {%f, %f, %f, lc}; \n',[w/2; -r ; lin/2]);
    
    tweak = 2*0.01315;
    % Assist points for mutual inductors top
    [x1,y1,x2,y2] = Tweak(r,0,r,tweak);    
    fprintf(fileID,'Point( 105 ) = {%f, %f, %f, lc}; \n',[x1+w/2 ; y1; lex/2]);
    fprintf(fileID,'Point( 106 ) = {%f, %f, %f, lc}; \n',[x1+w/2 ; y1 ; lin/2]);
    fprintf(fileID,'Point( 107 ) = {%f, %f, %f, lc}; \n',[x2+w/2 ; y2 ; lex/2]);
    fprintf(fileID,'Point( 108 ) = {%f, %f, %f, lc}; \n',[x2+w/2 ; y2 ; lin/2]);
    [x1,y1,x2,y2] = Tweak(0,r,r,tweak);    
    fprintf(fileID,'Point( 109 ) = {%f, %f, %f, lc}; \n',[x1 ; y1; lex/2]);
    fprintf(fileID,'Point( 110 ) = {%f, %f, %f, lc}; \n',[x1 ; y1 ; lin/2]);
    fprintf(fileID,'Point( 111 ) = {%f, %f, %f, lc}; \n',[x2 ; y2 ; lex/2]);
    fprintf(fileID,'Point( 112 ) = {%f, %f, %f, lc}; \n',[x2 ; y2 ; lin/2]);
    [x1,y1,x2,y2] = Tweak(-r,0,r,tweak);    
    fprintf(fileID,'Point( 113 ) = {%f, %f, %f, lc}; \n',[-x1-w/2 ; y1; lex/2]);
    fprintf(fileID,'Point( 114 ) = {%f, %f, %f, lc}; \n',[-x1-w/2 ; y1 ; lin/2]);
    fprintf(fileID,'Point( 115 ) = {%f, %f, %f, lc}; \n',[-x2-w/2 ; y2 ; lex/2]);
    fprintf(fileID,'Point( 116 ) = {%f, %f, %f, lc}; \n',[-x2-w/2 ; y2 ; lin/2]);
    [x1,y1,x2,y2] = Tweak(0,-r,r,tweak);    
    fprintf(fileID,'Point( 117 ) = {%f, %f, %f, lc}; \n',[x1 ; y1; lex/2]);
    fprintf(fileID,'Point( 118 ) = {%f, %f, %f, lc}; \n',[x1 ; y1 ; lin/2]);
    fprintf(fileID,'Point( 119 ) = {%f, %f, %f, lc}; \n',[x2 ; y2 ; lex/2]);
    fprintf(fileID,'Point( 120 ) = {%f, %f, %f, lc}; \n',[x2 ; y2 ; lin/2]);
    
    theta = 5/18*pi;
    % Assist points for mutual inductors bottom
    [x1,y1,x2,y2] = Tweak(r*cos(theta),r*sin(theta),r,tweak);    
    fprintf(fileID,'Point( 121 ) = {%f, %f, %f, lc}; \n',[x1+w/2 ; y1; -lex/2]);
    fprintf(fileID,'Point( 122 ) = {%f, %f, %f, lc}; \n',[x1+w/2 ; y1 ; -lin/2]);
    fprintf(fileID,'Point( 123 ) = {%f, %f, %f, lc}; \n',[x2+w/2 ; y2 ; -lex/2]);
    fprintf(fileID,'Point( 124 ) = {%f, %f, %f, lc}; \n',[x2+w/2 ; y2 ; -lin/2]);
    [x1,y1,x2,y2] = Tweak(-r*cos(theta),r*sin(theta),r,tweak); 
    fprintf(fileID,'Point( 125 ) = {%f, %f, %f, lc}; \n',[-x1-w/2  ; y1; -lex/2]);
    fprintf(fileID,'Point( 126 ) = {%f, %f, %f, lc}; \n',[-x1-w/2  ; y1 ; -lin/2]);
    fprintf(fileID,'Point( 127 ) = {%f, %f, %f, lc}; \n',[-x2-w/2  ; y2 ; -lex/2]);
    fprintf(fileID,'Point( 128 ) = {%f, %f, %f, lc}; \n',[-x2-w/2  ; y2 ; -lin/2]);
    [x1,y1,x2,y2] = Tweak(-r*cos(theta),-r*sin(theta),r,tweak); 
    fprintf(fileID,'Point( 129 ) = {%f, %f, %f, lc}; \n',[-x1-w/2 ; y1; -lex/2]);
    fprintf(fileID,'Point( 130 ) = {%f, %f, %f, lc}; \n',[-x1-w/2 ; y1 ; -lin/2]);
    fprintf(fileID,'Point( 131 ) = {%f, %f, %f, lc}; \n',[-x2-w/2 ; y2 ; -lex/2]);
    fprintf(fileID,'Point( 132 ) = {%f, %f, %f, lc}; \n',[-x2-w/2 ; y2 ; -lin/2]);
    [x1,y1,x2,y2] = Tweak(r*cos(theta),-r*sin(theta),r,tweak); 
    fprintf(fileID,'Point( 133 ) = {%f, %f, %f, lc}; \n',[x1+w/2 ; y1; -lex/2]);
    fprintf(fileID,'Point( 134 ) = {%f, %f, %f, lc}; \n',[x1+w/2 ; y1 ; -lin/2]);
    fprintf(fileID,'Point( 135 ) = {%f, %f, %f, lc}; \n',[x2+w/2 ; y2 ; -lex/2]);
    fprintf(fileID,'Point( 136 ) = {%f, %f, %f, lc}; \n',[x2+w/2 ; y2 ; -lin/2]);
    
    % Helices
    theta = 5/18*pi;
    [x1,y1,~,~] = Tweak(r,0,r,tweak);    
    [x2,y2,~,~] = Tweak(r*cos(theta),r*sin(theta),r,tweak);
    theta_start = atan(y2/x2);
    theta_end = atan(y1/x1);
    theta_steps = linspace(theta_start,theta_end,L);
    x = r*cos(theta_steps);
    y = r*sin(theta_steps);
    z = linspace(-lin/2,lin/2,L);
    for i =2:L-1
        fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i-1+136; x(i)+w/2 ; y(i) ; z(i)]); 
    end    
    theta = 5/18*pi;
    [~,~,x1,y1] = Tweak(r,0,r,tweak);    
    [~,~,x2,y2] = Tweak(r*cos(theta),r*sin(theta),r,tweak);
    theta_start = atan(y2/x2);
    theta_end = atan(y1/x1);
    theta_steps = linspace(theta_start,theta_end,L);
    x = r*cos(theta_steps);
    y = r*sin(theta_steps);
    z = linspace(-lin/2,lin/2,L);
    for i =2:L-1
        fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i-1+154; x(i)+w/2 ; y(i) ; z(i)]); 
    end
     theta = 5/18*pi;
    [~,~,x1,y1] = Tweak(r,0,r,tweak);    
    [~,~,x2,y2] = Tweak(r*cos(theta),-r*sin(theta),r,tweak); 
    theta_start = atan(y2/x2);
    theta_end = atan(y1/x1);
    theta_steps = linspace(theta_start,theta_end,L);
    x = r*cos(theta_steps);
    y = r*sin(theta_steps);
    z = linspace(-lin/2,lin/2,L);
    for i =2:L-1
        fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i-1+172; x(i)+w/2 ; y(i) ; z(i)]); 
    end
    theta = 5/18*pi;
    [x1,y1,~,~] = Tweak(r,0,r,tweak);    
    [x2,y2,~,~] = Tweak(r*cos(theta),-r*sin(theta),r,tweak); 
    theta_start = atan(y2/x2);
    theta_end = atan(y1/x1);
    theta_steps = linspace(theta_start,theta_end,L);
    x = r*cos(theta_steps);
    y = r*sin(theta_steps);
    z = linspace(-lin/2,lin/2,L);
    for i =2:L-1
        fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i-1+190; x(i)+w/2 ; y(i) ; z(i)]); 
    end
    
    theta = 5/18*pi;
    [x1,y1,~,~] = Tweak(-r,0,r,tweak);    
    [x2,y2,~,~] = Tweak(-r*cos(theta),r*sin(theta),r,tweak);
    theta_start = atan(y2/x2);
    theta_end = atan(y1/x1);
    theta_steps = linspace(theta_start,theta_end,L);
    x = -r*cos(theta_steps);
    y = r*sin(theta_steps);
    z = linspace(-lin/2,lin/2,L);
    for i =2:L-1
        fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i-1+208; x(i)-w/2 ; y(i) ; z(i)]); 
    end
    theta = 5/18*pi;
    [~,~,x1,y1] = Tweak(-r,0,r,tweak);    
    [~,~,x2,y2] = Tweak(-r*cos(theta),r*sin(theta),r,tweak);
    theta_start = atan(y2/x2);
    theta_end = atan(y1/x1);
    theta_steps = linspace(theta_start,theta_end,L);
    x = -r*cos(theta_steps);
    y = r*sin(theta_steps);
    z = linspace(-lin/2,lin/2,L);
    for i =2:L-1
        fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i-1+226; x(i)-w/2 ; y(i) ; z(i)]); 
    end
    theta = 5/18*pi;
    [~,~,x1,y1] = Tweak(-r,0,r,tweak);    
    [~,~,x2,y2] = Tweak(-r*cos(theta),-r*sin(theta),r,tweak); 
    theta_start = atan(y2/x2);
    theta_end = atan(y1/x1);
    theta_steps = linspace(theta_start,theta_end,L);
    x = -r*cos(theta_steps);
    y = r*sin(theta_steps);
    z = linspace(-lin/2,lin/2,L);
    for i =2:L-1
        fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i-1+244; x(i)-w/2 ; y(i) ; z(i)]); 
    end
    theta = 5/18*pi;
    [x1,y1,~,~] = Tweak(-r,0,r,tweak);    
    [x2,y2,~,~] = Tweak(-r*cos(theta),-r*sin(theta),r,tweak); 
    theta_start = atan(y2/x2);
    theta_end = atan(y1/x1);
    theta_steps = linspace(theta_start,theta_end,L);
    x = -r*cos(theta_steps);
    y = r*sin(theta_steps);
    z = linspace(-lin/2,lin/2,L);
    for i =2:L-1
        fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i-1+262; x(i)-w/2 ; y(i) ; z(i)]); 
    end
        theta = 5/18*pi;
    % Assist Points for helices
    fprintf(fileID,'Point( 281 ) = {%f, %f, %f, lc}; \n',[-r-w/2 ; 0 ; lin/2-0.004]); 
    fprintf(fileID,'Point( 282 ) = {%f, %f, %f, lc}; \n',[r+w/2 ; 0 ; lin/2-0.004]); 
    fprintf(fileID,'Point( 283 ) = {%f, %f, %f, lc}; \n',[0 ; r ; lin/2-0.004]); 
    fprintf(fileID,'Point( 284 ) = {%f, %f, %f, lc}; \n',[0 ; -r ; lin/2-0.004]); 
    fprintf(fileID,'Point( 285 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; r*sin(theta); -lin/2+0.004]); 
    fprintf(fileID,'Point( 286 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; r*sin(theta) ; -lin/2+0.004]); 
    fprintf(fileID,'Point( 287 ) = {%f, %f, %f, lc}; \n',[-r*cos(theta)-w/2  ; -r*sin(theta) ; -lin/2+0.004]); 
    fprintf(fileID,'Point( 288 ) = {%f, %f, %f, lc}; \n',[r*cos(theta)+w/2 ; -r*sin(theta); -lin/2+0.004]); 

    theta = 5/18*pi;
    [x1,y1,~,~] = Tweak(-r*cos(theta),r*sin(theta),r,tweak); 
    theta_start = atan(-y1/x1);
    theta_end = pi/2;
    theta_steps = linspace(theta_start,theta_end,L1);
    x = r*cos(theta_steps);
    y = r*sin(theta_steps);
    z = linspace(-lin/2,lin/2,L);
    for i =2:L1
        fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i-1+288; -x(i)-w/2 ; -y(i) ; z(i)]); 
    end
    for i =3:L1
        fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i-2+303; -x(i)-w/2 ; -y(i) ; z(i-1)]); 
    end
    
    theta = 5/18*pi;
    [~,~,x1,y1] = Tweak(r*cos(theta),r*sin(theta),r,tweak);   
    theta_start = atan(y1/x1);
    theta_end = pi/2;
    theta_steps = linspace(theta_start,theta_end,L1);
    x = r*cos(theta_steps);
    y = r*sin(theta_steps);
    z = linspace(-lin/2,lin/2,L);
    for i =2:L1
        fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i-1+317; x(i)+w/2 ; y(i) ; z(i)]); 
    end
    for i =3:L1
        fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i-1+332; x(i)+w/2 ; y(i) ; z(i-1)]); 
    end  
    
    theta = 5/18*pi;
    [~,~,x1,y1] = Tweak(-r*cos(theta),-r*sin(theta),r,tweak); 
    theta_start = atan(y1/x1);
    theta_end = pi/2;
    theta_steps = linspace(theta_start,theta_end,L1);
    x = r*cos(theta_steps);
    y = r*sin(theta_steps);
    z = linspace(-lin/2,lin/2,L);
    for i =2:L1
        fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i-1+347; -x(i)-w/2 ; y(i) ; z(i)]); 
    end
    for i =3:L1
        fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i-2+362; -x(i)-w/2 ; y(i) ; z(i-1)]); 
    end  
    
    theta = 5/18*pi;
    [x1,y1,~,~] = Tweak(r*cos(theta),-r*sin(theta),r,tweak); 
    theta_start = atan(y1/x1);
    theta_end = -pi/2;
    theta_steps = linspace(theta_start,theta_end,L1);
    x = r*cos(theta_steps);
    y = r*sin(theta_steps);
    z = linspace(-lin/2,lin/2,L);
    for i =2:L1
        fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i-1+376; x(i)+w/2 ; y(i) ; z(i)]); 
    end
    for i =3:L1
        fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i-1+391; x(i)+w/2 ; y(i) ; z(i-1)]); 
    end  
    
    %% Lines
    fprintf(fileID,'Line(1) = {91,111}; \n');    
    fprintf(fileID,'Line(2) = {111,31}; \n');   
    fprintf(fileID,'Line(3) = {29,109}; \n');   
    fprintf(fileID,'Line(4) = {109,95}; \n');   
    fprintf(fileID,'Line(5) = {99,119}; \n');    
    fprintf(fileID,'Line(6) = {119,39}; \n');   
    fprintf(fileID,'Line(7) = {37,117}; \n');   
    fprintf(fileID,'Line(8) = {117,103}; \n');   
    fprintf(fileID,'Line(9) = {96,110}; \n');    
    fprintf(fileID,'Line(10) = {110,30}; \n');   
    fprintf(fileID,'Line(11) = {32,112}; \n');   
    fprintf(fileID,'Line(12) = {112,92}; \n');   
    fprintf(fileID,'Line(13) = {104,118}; \n');    
    fprintf(fileID,'Line(14) = {118,38}; \n');   
    fprintf(fileID,'Line(15) = {40,120}; \n');   
    fprintf(fileID,'Line(16) = {120,100}; \n');   
    fprintf(fileID,'Circle(17) = {107,1,27}; \n');    
    fprintf(fileID,'Circle(18) = {25,1,105}; \n');   
    fprintf(fileID,'Circle(19) = {106,2,26}; \n');   
    fprintf(fileID,'Circle(20) = {28,2,108}; \n');   
    fprintf(fileID,'Circle(21) = {113,5,33}; \n');    
    fprintf(fileID,'Circle(22) = {35,5,115}; \n');   
    fprintf(fileID,'Circle(23) = {116,6,36}; \n');   
    fprintf(fileID,'Circle(24) = {34,6,114}; \n');   
    fprintf(fileID,'Circle(25) = {126,7,46}; \n');    
    fprintf(fileID,'Circle(26) = {48,7,128}; \n');   
    fprintf(fileID,'Circle(27) = {127,8,47}; \n');   
    fprintf(fileID,'Circle(28) = {45,8,125}; \n');   
    fprintf(fileID,'Line(29) = {98,16}; \n');    
    fprintf(fileID,'Line(30) = {16,102}; \n');   
    fprintf(fileID,'Line(31) = {101,15}; \n');   
    fprintf(fileID,'Line(32) = {15,97}; \n');      
    fprintf(fileID,'Circle(33) = {136,3,56}; \n');      
    fprintf(fileID,'Circle(34) = {54,3,134}; \n');      
    fprintf(fileID,'Circle(35) = {133,4,53}; \n');      
    fprintf(fileID,'Circle(36) = {55,4,135}; \n');      
    fprintf(fileID,'Circle(37) = {124,3,44}; \n');      
    fprintf(fileID,'Circle(38) = {42,3,122}; \n');      
    fprintf(fileID,'Circle(39) = {121,4,41}; \n');      
    fprintf(fileID,'Circle(40) = {43,4,123}; \n');      
    fprintf(fileID,'Line(41) = {90,12}; \n');      
    fprintf(fileID,'Line(42) = {12,94}; \n');      
    fprintf(fileID,'Line(43) = {93,11}; \n');   
    fprintf(fileID,'Line(44) = {11,89}; \n');      
    fprintf(fileID,'Circle(45) = {130,7,50}; \n');      
    fprintf(fileID,'Circle(46) = {52,7,132}; \n');      
    fprintf(fileID,'Circle(47) = {131,8,51}; \n');   
    fprintf(fileID,'Circle(48) = {49,8,129}; \n');   
    fprintf(fileID,'Circle( 49 ) = {105, 1, 57}; \n');
    fprintf(fileID,'Circle( 50 ) = {57, 1, 17}; \n');       
    fprintf(fileID,'Circle( 51 ) = {17, 1, 65}; \n');    
    fprintf(fileID,'Circle( 52 ) = {65, 1, 91}; \n');    
    fprintf(fileID,'Circle( 53 ) = {95, 5, 67}; \n');
    fprintf(fileID,'Circle( 54 ) = {67, 5, 19}; \n');       
    fprintf(fileID,'Circle( 55 ) = {19, 5, 59}; \n');    
    fprintf(fileID,'Circle( 56 ) = {59, 5, 113}; \n');
    fprintf(fileID,'Circle( 57 ) = {115, 5, 61}; \n');
    fprintf(fileID,'Circle( 58 ) = {61, 5, 21}; \n');       
    fprintf(fileID,'Circle( 59 ) = {21, 5, 69}; \n');    
    fprintf(fileID,'Circle( 60 ) = {69, 5, 99}; \n');    
    fprintf(fileID,'Circle( 61 ) = {103, 1, 71}; \n');
    fprintf(fileID,'Circle( 62 ) = {71, 1, 23}; \n');       
    fprintf(fileID,'Circle( 63 ) = {23, 1, 63}; \n');    
    fprintf(fileID,'Circle( 64 ) = {63, 1, 107}; \n');    
    fprintf(fileID,'Circle( 65 ) = {58,2,106 }; \n');
    fprintf(fileID,'Circle( 66 ) = {18,2,58}; \n');       
    fprintf(fileID,'Circle( 67 ) = {66,2,18}; \n');    
    fprintf(fileID,'Circle( 68 ) = {92,2,66}; \n');    
    fprintf(fileID,'Circle( 69 ) = {68,6,96}; \n');
    fprintf(fileID,'Circle( 70 ) = {20,6,68}; \n');       
    fprintf(fileID,'Circle( 71 ) = {60,6,20}; \n');    
    fprintf(fileID,'Circle( 72 ) = {114,6,60}; \n');
    fprintf(fileID,'Circle( 73 ) = {62,6,116}; \n');
    fprintf(fileID,'Circle( 74 ) = {22,6,62}; \n');       
    fprintf(fileID,'Circle( 75 ) = {70,6,22}; \n');    
    fprintf(fileID,'Circle( 76 ) = {100,6,70}; \n');    
    fprintf(fileID,'Circle( 77 ) = {72,2,104}; \n');
    fprintf(fileID,'Circle( 78 ) = {24,2,72}; \n');       
    fprintf(fileID,'Circle( 79 ) = {64,2,24}; \n');    
    fprintf(fileID,'Circle( 80 ) = {108,2,64}; \n'); 
    fprintf(fileID,'Circle( 81 ) = {94,7,84}; \n'); 
    fprintf(fileID,'Circle( 82 ) = {84,7,130}; \n'); 
    fprintf(fileID,'Circle( 83 ) = {132,7,76}; \n'); 
    fprintf(fileID,'Circle( 84 ) = {76,7,14}; \n'); 
    fprintf(fileID,'Circle( 85 ) = {14,7,78}; \n'); 
    fprintf(fileID,'Circle( 86 ) = {78,7,126}; \n'); 
    fprintf(fileID,'Circle( 87 ) = {128,7,86}; \n'); 
    fprintf(fileID,'Circle( 88 ) = {86,7,98}; \n'); 
    fprintf(fileID,'Circle( 89 ) = {102,3,88}; \n'); 
    fprintf(fileID,'Circle( 90 ) = {88,3,136}; \n'); 
    fprintf(fileID,'Circle( 91 ) = {134,3,80}; \n'); 
    fprintf(fileID,'Circle( 92 ) = {80,3,10}; \n'); 
    fprintf(fileID,'Circle( 93 ) = {10,3,74}; \n'); 
    fprintf(fileID,'Circle( 94 ) = {74,3,124}; \n'); 
    fprintf(fileID,'Circle( 95 ) = {122,3,82}; \n'); 
    fprintf(fileID,'Circle( 96 ) = {82,3,90}; \n');
    fprintf(fileID,'Circle( 97 ) = {83,8,93}; \n'); 
    fprintf(fileID,'Circle( 98 ) = {129,8,83}; \n'); 
    fprintf(fileID,'Circle( 99 ) = {75,8,131}; \n'); 
    fprintf(fileID,'Circle( 100 ) = {13,8,75}; \n'); 
    fprintf(fileID,'Circle( 101 ) = {77,8,13}; \n'); 
    fprintf(fileID,'Circle( 102 ) = {125,8,77}; \n'); 
    fprintf(fileID,'Circle( 103 ) = {85,8,127}; \n'); 
    fprintf(fileID,'Circle( 104 ) = {97,8,85}; \n'); 
    fprintf(fileID,'Circle( 105 ) = {87,4,101}; \n'); 
    fprintf(fileID,'Circle( 106 ) = {135,4,87}; \n'); 
    fprintf(fileID,'Circle( 107 ) = {79,4,133}; \n'); 
    fprintf(fileID,'Circle( 108 ) = {9,4,79}; \n'); 
    fprintf(fileID,'Circle( 109 ) = {73,4,9}; \n'); 
    fprintf(fileID,'Circle( 110 ) = {123,4,73}; \n'); 
    fprintf(fileID,'Circle( 111 ) = {81,4,121}; \n'); 
    fprintf(fileID,'Circle( 112 ) = {89,4,81}; \n');
    fprintf(fileID,'Line(113) = {35,36}; \n');
    fprintf(fileID,'Line(114) = {115,116}; \n');    
    fprintf(fileID,'Line(115) = {61,62}; \n');      
    fprintf(fileID,'Line(116) = {21,22}; \n');
    fprintf(fileID,'Line(117) = {69,70}; \n');    
    fprintf(fileID,'Line(118) = {99,100}; \n');    
    fprintf(fileID,'Line(119) = {119,120}; \n');
    fprintf(fileID,'Line(120) = {39,40}; \n');    
    fprintf(fileID,'Line(121) = {37,38}; \n');    
    fprintf(fileID,'Line(122) = {117,118}; \n');
    fprintf(fileID,'Line(123) = {103,104}; \n');    
    fprintf(fileID,'Line(124) = {71,72}; \n');    
    fprintf(fileID,'Line(125) = {23,24}; \n');
    fprintf(fileID,'Line(126) = {63,64}; \n');    
    fprintf(fileID,'Line(127) = {107,108}; \n');    
    fprintf(fileID,'Line(128) = {27,28}; \n');
    fprintf(fileID,'Line(129) = {25,26}; \n');    
    fprintf(fileID,'Line(130) = {105,106}; \n');    
    fprintf(fileID,'Line(131) = {57,58}; \n');    
    fprintf(fileID,'Line(132) = {17,18}; \n');    
    fprintf(fileID,'Line(133) = {65,66}; \n');
    fprintf(fileID,'Line(134) = {91,92}; \n');    
    fprintf(fileID,'Line(135) = {111,112}; \n');    
    fprintf(fileID,'Line(136) = {31,32}; \n');
    fprintf(fileID,'Line(137) = {29,30}; \n');    
    fprintf(fileID,'Line(138) = {109,110}; \n');    
    fprintf(fileID,'Line(139) = {95,96}; \n');    
    fprintf(fileID,'Line(140) = {67,68}; \n');    
    fprintf(fileID,'Line(141) = {19,20}; \n');
    fprintf(fileID,'Line(142) = {59,60}; \n');
    fprintf(fileID,'Line(143) = {113,114}; \n');
    fprintf(fileID,'Line(144) = {33,34}; \n');
    fprintf(fileID,'Line(145) = {42,41}; \n');
    fprintf(fileID,'Line(146) = {122,121}; \n');    
    fprintf(fileID,'Line(147) = {82,81}; \n');      
    fprintf(fileID,'Line(148) = {90,89}; \n');
    fprintf(fileID,'Line(149) = {12,11}; \n');    
    fprintf(fileID,'Line(150) = {94,93}; \n');    
    fprintf(fileID,'Line(151) = {84,83}; \n');
    fprintf(fileID,'Line(152) = {130,129}; \n');    
    fprintf(fileID,'Line(153) = {50,49}; \n');    
    fprintf(fileID,'Line(154) = {52,51}; \n');
    fprintf(fileID,'Line(155) = {132,131}; \n');    
    fprintf(fileID,'Line(156) = {76,75}; \n');    
    fprintf(fileID,'Line(157) = {14,13}; \n');
    fprintf(fileID,'Line(158) = {78,77}; \n');    
    fprintf(fileID,'Line(159) = {126,125}; \n');    
    fprintf(fileID,'Line(160) = {46,45}; \n');
    fprintf(fileID,'Line(161) = {48,47}; \n');    
    fprintf(fileID,'Line(162) = {128,127}; \n');    
    fprintf(fileID,'Line(163) = {86,85}; \n');
    fprintf(fileID,'Line(164) = {98,97}; \n');    
    fprintf(fileID,'Line(165) = {16,15}; \n');    
    fprintf(fileID,'Line(166) = {102,101}; \n');
    fprintf(fileID,'Line(167) = {88,87}; \n');    
    fprintf(fileID,'Line(168) = {136,135}; \n');    
    fprintf(fileID,'Line(169) = {56,55}; \n');
    fprintf(fileID,'Line(170) = {54,53}; \n');    
    fprintf(fileID,'Line(171) = {134,133}; \n');    
    fprintf(fileID,'Line(172) = {80,79}; \n');
    fprintf(fileID,'Line(173) = {10,9}; \n');    
    fprintf(fileID,'Line(174) = {74,73}; \n');    
    fprintf(fileID,'Line(175) = {124,123}; \n');
    fprintf(fileID,'Line(176) = {44,43}; \n');
    fprintf(fileID,'Circle( 177 ) = {28,2,26}; \n'); 
    fprintf(fileID,'Line( 178 ) = {30,32}; \n'); 
    fprintf(fileID,'Circle( 179 ) = {36,2,34}; \n'); 
    fprintf(fileID,'Line( 180 ) = {38,40}; \n'); 
    fprintf(fileID,'Circle( 181 ) = {42,3,44}; \n'); 
    fprintf(fileID,'Circle( 182 ) = {52,7,50}; \n'); 
    fprintf(fileID,'Circle( 183 ) = {48,7,46}; \n'); 
    fprintf(fileID,'Circle( 184 ) = {54,3,56}; \n'); 
    fprintf(fileID,'Line( 185 ) = {268,250}; \n'); 
    fprintf(fileID,'Line( 186 ) = {276,258}; \n'); 
    fprintf(fileID,'Line( 187 ) = {222,240}; \n'); 
    fprintf(fileID,'Line( 188 ) = {214,232}; \n');  
    fprintf(fileID,'Line( 189 ) = {168,150}; \n'); 
    fprintf(fileID,'Line( 190 ) = {160,142}; \n');  
    fprintf(fileID,'Line( 191 ) = {186,204}; \n'); 
    fprintf(fileID,'Line( 192 ) = {178,196}; \n'); 
    fprintf(fileID,'Spline( 193 ) = {286,263,264,265,266,267,268}; \n');
    fprintf(fileID,'Spline( 194 ) = {268,269,270,271,272,273,274,275,276}; \n');
    fprintf(fileID,'Spline( 195 ) = {276,277,278,279,280,114}; \n'); 
    fprintf(fileID,'Spline( 196 ) = {281,262,261,260,259,258}; \n');
    fprintf(fileID,'Spline( 197 ) = {258,257,256,255,254,253,252,251,250}; \n');
    fprintf(fileID,'Spline( 198 ) = {250,249,248,247,246,245,132}; \n'); 
    fprintf(fileID,'Spline( 199 ) = {116,244,243,242,241,240}; \n');
    fprintf(fileID,'Spline( 200 ) = {240,239,238,237,236,235,234,233,232}; \n');
    fprintf(fileID,'Spline( 201 ) = {232,231,230,229,228,227,287}; \n'); 
    fprintf(fileID,'Spline( 202 ) = {126,209,210,211,212,213,214}; \n');
    fprintf(fileID,'Spline( 203 ) = {214,215,216,217,218,219,220,221,222}; \n');
    fprintf(fileID,'Spline( 204 ) = {222,223,224,225,226,281}; \n'); 
    fprintf(fileID,'Spline( 205 ) = {108,190,189,187,186}; \n');
    fprintf(fileID,'Spline( 206 ) = {186,185,184,183,182,181,180,179,178}; \n');
    fprintf(fileID,'Spline( 207 ) = {178,177,176,175,174,173,288}; \n'); 
    fprintf(fileID,'Spline( 208 ) = {134,191,192,193,194,195,196}; \n');
    fprintf(fileID,'Spline( 209 ) = {196,197,198,199,200,201,202,203,204}; \n');
    fprintf(fileID,'Spline( 210 ) = {204,205,206,207,208,282}; \n'); 
    fprintf(fileID,'Spline( 211 ) = {282,172,171,170,169,168}; \n');
    fprintf(fileID,'Spline( 212 ) = {168,167,166,165,164,163,162,161,160}; \n');
    fprintf(fileID,'Spline( 213 ) = {160,159,158,157,156,155,124}; \n'); 
    fprintf(fileID,'Spline( 214 ) = {285,137,138,139,140,141,142}; \n');
    fprintf(fileID,'Spline( 215 ) = {142,143,144,145,146,147,148,149,150}; \n');
    fprintf(fileID,'Spline( 216 ) = {150,151,152,153,154,106}; \n'); 
    fprintf(fileID,'Line( 217 ) = {347,331}; \n'); 
    fprintf(fileID,'Line( 218 ) = {339,323}; \n'); 
    fprintf(fileID,'Line( 219 ) = {376,361}; \n'); 
    fprintf(fileID,'Line( 220 ) = {368,353}; \n'); 
    fprintf(fileID,'Line( 221 ) = {302,317}; \n'); 
    fprintf(fileID,'Line( 222 ) = {294,309}; \n'); 
    fprintf(fileID,'Line( 223 ) = {390,406}; \n'); 
    fprintf(fileID,'Line( 224 ) = {382,398}; \n'); 
    fprintf(fileID,'Spline( 225 ) = {287,289,290,291,292,293,294}; \n'); 
    fprintf(fileID,'Spline( 226 ) = {294,295,296,297,298,299,300,301,302}; \n'); 
    fprintf(fileID,'Spline( 227 ) = {302,303,120}; \n'); 
    fprintf(fileID,'Spline( 228 ) = {284,317}; \n'); 
    fprintf(fileID,'Spline( 229 ) = {317,316,315,314,313,312,311,310,309}; \n'); 
    fprintf(fileID,'Spline( 230 ) = {309,308,307,306,305,304,128}; \n'); 
    fprintf(fileID,'Spline( 231 ) = {118,391,390}; \n'); 
    fprintf(fileID,'Spline( 232 ) = {390,389,388,387,386,385,384,383,382}; \n'); 
    fprintf(fileID,'Spline( 233 ) = {382,381,380,379,378,377,288}; \n'); 
    fprintf(fileID,'Spline( 234 ) = {136,393,394,395,396,397,398}; \n'); 
    fprintf(fileID,'Spline( 235 ) = {398,399,400,401,402,403,404,405,406}; \n'); 
    fprintf(fileID,'Spline( 236 ) = {406,284}; \n'); 
    fprintf(fileID,'Spline( 237 ) = {285,318,319,320,321,322,323}; \n'); 
    fprintf(fileID,'Spline( 238 ) = {323,324,325,326,327,328,329,330,331}; \n'); 
    fprintf(fileID,'Spline( 239 ) = {331,332,112}; \n'); 
    fprintf(fileID,'Spline( 240 ) = {283,347}; \n'); 
    fprintf(fileID,'Spline( 241 ) = {347,346,345,344,343,342,341,340,339}; \n'); 
    fprintf(fileID,'Spline( 242 ) = {339,338,337,336,335,334,122}; \n'); 
    fprintf(fileID,'Spline( 243 ) = {286,348,349,350,351,352,353}; \n'); 
    fprintf(fileID,'Spline( 244 ) = {353,354,355,356,357,358,359,360,361}; \n'); 
    fprintf(fileID,'Spline( 245 ) = {361,362,110}; \n'); 
    fprintf(fileID,'Spline( 246 ) = {283,376}; \n'); 
    fprintf(fileID,'Spline( 247 ) = {376,375,374,373,372,371,370,369,368}; \n'); 
    fprintf(fileID,'Spline( 248 ) = {368,367,366,365,364,363,130}; \n'); 
    fprintf(fileID,'Line( 349 ) = {114,281}; \n'); 
    fprintf(fileID,'Line( 350 ) = {281,116}; \n'); 
    fprintf(fileID,'Line( 351 ) = {108,282}; \n'); 
    fprintf(fileID,'Line( 352 ) = {282,106}; \n'); 
    fprintf(fileID,'Line( 353 ) = {112,283}; \n'); 
    fprintf(fileID,'Line( 354 ) = {283,110}; \n'); 
    fprintf(fileID,'Line( 355 ) = {120,284}; \n'); 
    fprintf(fileID,'Line( 356 ) = {284,118}; \n'); 
    fprintf(fileID,'Line( 357 ) = {124,285}; \n'); 
    fprintf(fileID,'Line( 358 ) = {285,122}; \n'); 
    fprintf(fileID,'Line( 359 ) = {130,286}; \n'); 
    fprintf(fileID,'Line( 360 ) = {286,132}; \n'); 
    fprintf(fileID,'Line( 361 ) = {126,287}; \n'); 
    fprintf(fileID,'Line( 362 ) = {287,128}; \n'); 
    fprintf(fileID,'Line( 363 ) = {136,288}; \n'); 
    fprintf(fileID,'Line( 364 ) = {288,134}; \n'); 
    
    %% Line Loops
    fprintf(fileID,'Line Loop( 1 ) = {-18,-130,-19,129}; \n');
    fprintf(fileID,'Line Loop( 2 ) = {-49,-131,-65,130}; \n');
    fprintf(fileID,'Line Loop( 3 ) = {-50,-132,-66,131}; \n');
    fprintf(fileID,'Line Loop( 4 ) = {-51,-133,-67,132}; \n');
    fprintf(fileID,'Line Loop( 5 ) = {-52,-134,-68,133}; \n');
    fprintf(fileID,'Line Loop( 6 ) = {-1,-135,-12,134}; \n');
    fprintf(fileID,'Line Loop( 7 ) = {-2,-136,-11,135}; \n');
    fprintf(fileID,'Line Loop( 8 ) = {-3,-138,-10,137}; \n');
    fprintf(fileID,'Line Loop( 9 ) = {-4,-139,-9,138}; \n');
    fprintf(fileID,'Line Loop( 10 ) = {-53,-140,-69,139}; \n');
    fprintf(fileID,'Line Loop( 11 ) = {-54,-141,-70,140}; \n');
    fprintf(fileID,'Line Loop( 12 ) = {-55,-142,-71,141}; \n');
    fprintf(fileID,'Line Loop( 13 ) = {-56,-143,-72,142}; \n');
    fprintf(fileID,'Line Loop( 14 ) = {-21,-144,-24,143}; \n');
    fprintf(fileID,'Line Loop( 15 ) = {-22,-114,-23,113}; \n');
    fprintf(fileID,'Line Loop( 16 ) = {-57,-115,-73,114}; \n');
    fprintf(fileID,'Line Loop( 17 ) = {-58,-116,-74,115}; \n');
    fprintf(fileID,'Line Loop( 18 ) = {-59,-117,-75,116}; \n');
    fprintf(fileID,'Line Loop( 19 ) = {-60,-118,-76,117}; \n');
    fprintf(fileID,'Line Loop( 20 ) = {-5,-119,-16,118}; \n');
    fprintf(fileID,'Line Loop( 21 ) = {-6,-120,-15,119}; \n');
    fprintf(fileID,'Line Loop( 22 ) = {-7,-122,-14,121}; \n');
    fprintf(fileID,'Line Loop( 23 ) = {-8,-123,-13,122}; \n');
    fprintf(fileID,'Line Loop( 24 ) = {-61,-124,-77,123}; \n');
    fprintf(fileID,'Line Loop( 25 ) = {-62,-125,-78,124}; \n');
    fprintf(fileID,'Line Loop( 26 ) = {-63,-126,-79,125}; \n');
    fprintf(fileID,'Line Loop( 27 ) = {-64,-127,-80,126}; \n');
    fprintf(fileID,'Line Loop( 28 ) = {-17,-128,-20,127}; \n');
    fprintf(fileID,'Line Loop( 29 ) = {-38,-146,-39,145}; \n');
    fprintf(fileID,'Line Loop( 30 ) = {-95,-147,-111,146}; \n');
    fprintf(fileID,'Line Loop( 31 ) = {-96,-148,-112,147}; \n');
    fprintf(fileID,'Line Loop( 32 ) = {-41,-149,-44,148}; \n');
    fprintf(fileID,'Line Loop( 33 ) = {-42,-150,-43,149}; \n');
    fprintf(fileID,'Line Loop( 34 ) = {-81,-151,-97,150}; \n');
    fprintf(fileID,'Line Loop( 35 ) = {-82,-152,-98,151}; \n');
    fprintf(fileID,'Line Loop( 36 ) = {-45,-153,-48,152}; \n');
    fprintf(fileID,'Line Loop( 37 ) = {-46,-155,-47,154}; \n');
    fprintf(fileID,'Line Loop( 38 ) = {-83,-156,-99,155}; \n');
    fprintf(fileID,'Line Loop( 39 ) = {-84,-157,-100,156}; \n');
    fprintf(fileID,'Line Loop( 40 ) = {-85,-158,-101,157}; \n');
    fprintf(fileID,'Line Loop( 41 ) = {-86,-159,-102,158}; \n');
    fprintf(fileID,'Line Loop( 42 ) = {-25,-160,-28,159}; \n');
    fprintf(fileID,'Line Loop( 43 ) = {-26,-162,-27,161}; \n');
    fprintf(fileID,'Line Loop( 44 ) = {-87,-163,-103,162}; \n');
    fprintf(fileID,'Line Loop( 45 ) = {-88,-164,-104,163}; \n');
    fprintf(fileID,'Line Loop( 46 ) = {-29,-165,-32,164}; \n');
    fprintf(fileID,'Line Loop( 47 ) = {-30,-166,-31,165}; \n');
    fprintf(fileID,'Line Loop( 48 ) = {-89,-167,-105,166}; \n');
    fprintf(fileID,'Line Loop( 49 ) = {-90,-168,-106,167}; \n');
    fprintf(fileID,'Line Loop( 50 ) = {-33,-169,-36,168}; \n');
    fprintf(fileID,'Line Loop( 51 ) = {-34,-171,-35,170}; \n');
    fprintf(fileID,'Line Loop( 52 ) = {-91,-172,-107,171}; \n');
    fprintf(fileID,'Line Loop( 53 ) = {-92,-173,-108,172}; \n');
    fprintf(fileID,'Line Loop( 54 ) = {-93,-174,-109,173}; \n');
    fprintf(fileID,'Line Loop( 55 ) = {-94,-175,-110,174}; \n');
    fprintf(fileID,'Line Loop( 56 ) = {-37,-176,-40,175}; \n');
    fprintf(fileID,'Line Loop( 57 ) = {351,352,20,-177,19}; \n');
    fprintf(fileID,'Line Loop( 58 ) = {353,354,11,178,10}; \n');
    fprintf(fileID,'Line Loop( 59 ) = {349,350,24,179,23}; \n');
    fprintf(fileID,'Line Loop( 60 ) = {355,356,15,180,14}; \n');
    fprintf(fileID,'Line Loop( 61 ) = {-357,-358,38,-181,37}; \n');
    fprintf(fileID,'Line Loop( 62 ) = {-359,-360,46,-182,45}; \n');
    fprintf(fileID,'Line Loop( 63 ) = {-361,-362,26,-183,25}; \n');
    fprintf(fileID,'Line Loop( 64 ) = {-363,-364,34,-184,33}; \n');
    fprintf(fileID,'Line Loop( 65 ) = {362,-230,-222,-225}; \n');
    fprintf(fileID,'Line Loop( 66 ) = {222,-226,-221,-229}; \n');
    fprintf(fileID,'Line Loop( 67 ) = {221,-227,-355,-228}; \n');
    fprintf(fileID,'Line Loop( 68 ) = {-356,-231,-223,-236}; \n');
    fprintf(fileID,'Line Loop( 69 ) = {223,-232,-224,-235}; \n');
    fprintf(fileID,'Line Loop( 70 ) = {224,-234,363,-233}; \n');
    fprintf(fileID,'Line Loop( 71 ) = {364,208,-192,207}; \n');
    fprintf(fileID,'Line Loop( 72 ) = {192,209,-191,206}; \n');
    fprintf(fileID,'Line Loop( 73 ) = {191,210,-351,205}; \n');
    fprintf(fileID,'Line Loop( 74 ) = {-352,216,189,211}; \n');
    fprintf(fileID,'Line Loop( 75 ) = {-189,215,190,212}; \n');
    fprintf(fileID,'Line Loop( 76 ) = {-190,214,357,213}; \n');
    fprintf(fileID,'Line Loop( 77 ) = {358,-242,218,-237}; \n');
    fprintf(fileID,'Line Loop( 78 ) = {-218,-241,217,-238}; \n');
    fprintf(fileID,'Line Loop( 79 ) = {-217,-240,-353,-239}; \n');
    fprintf(fileID,'Line Loop( 80 ) = {-354,245,219,246}; \n');
    fprintf(fileID,'Line Loop( 81 ) = {-219,244,220,247}; \n');
    fprintf(fileID,'Line Loop( 82 ) = {-220,243,359,248}; \n');
    fprintf(fileID,'Line Loop( 83 ) = {360,-198,-185,-193}; \n');
    fprintf(fileID,'Line Loop( 84 ) = {185,-197,-186,-194}; \n');
    fprintf(fileID,'Line Loop( 85 ) = {186,-196,-349,-195}; \n');
    fprintf(fileID,'Line Loop( 86 ) = {-350,-199,187,-204}; \n');
    fprintf(fileID,'Line Loop( 87 ) = {-187,-200,188,-203}; \n');
    fprintf(fileID,'Line Loop( 88 ) = {-188,-201,361,-202}; \n');
        
    if shield.flag == 1
        fprintf(fileID,'lc_s = %f; \n\n',shield.meshing);   
        
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[500+1  ; 0    ; 0    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[500+2  ; 0    ; 0    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[500+3  ; shield.radius    ; 0    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[500+4  ; 0    ; shield.radius    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[500+5  ; -shield.radius    ; 0    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[500+6  ; 0    ; -shield.radius    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[500+7  ; shield.radius    ; 0    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[500+8  ; 0    ; shield.radius    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[500+9  ; -shield.radius    ; 0    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[500+10  ; 0    ; -shield.radius    ; -shield.length/2]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[500+1; 500+3;500+7 ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[500+2; 500+4;500+8 ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[500+3; 500+5;500+9 ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[500+4; 500+6;500+10 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[500+5; 500+3;500+1;500+4 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[500+6; 500+4;500+1;500+5 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[500+7; 500+5;500+1;500+6 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[500+8; 500+6;500+1;500+3 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[500+9; 500+7;500+2;500+8 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[500+10; 500+8;500+2;500+9 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[500+11; 500+9;500+2;500+10 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[500+12; 500+10;500+2;500+7 ]);
        
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[500+1;   500+1 ;   500+9    ; -(500+2); -(500+5)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[500+2;   500+2 ;   500+10    ; -(500+3); -(500+6)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[500+3;   500+3 ;   500+11    ; -(500+4); -(500+7)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[500+4;   500+4 ;   500+12    ; -(500+1); -(500+8)]);
        
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[500+1  ;500+1  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[500+2  ;500+2  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[500+3  ;500+3  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[500+4  ;500+4  ]);
        
        fprintf(fileID,'Physical Surface(2)  = {500+1,500+2,500+3,500+4};\n\n');
        
    end
    
    %% Ruled Surfaces
    for i = 1:56
        fprintf(fileID,'Ruled Surface( %i ) = {%i}; \n',i,i);
    end
    for i = 57:64
        fprintf(fileID,'Plane Surface( %i ) = {%i}; \n',i,i);
    end
    for i = 65:88
        fprintf(fileID,'Ruled Surface( %i ) = {%i}; \n',i,i);
    end
    
    %% Physical Lines
    % Sources
    fprintf(fileID,'Physical Line (1) = {173}; \n');    
    fprintf(fileID,'Physical Line (2) = {132}; \n');      
    fprintf(fileID,'Physical Line (3) = {149}; \n');       
    fprintf(fileID,'Physical Line (4) = {141}; \n');      
    fprintf(fileID,'Physical Line (5) = {157}; \n');    
    fprintf(fileID,'Physical Line (6) = {116}; \n');       
    fprintf(fileID,'Physical Line (7) = {165}; \n');        
    fprintf(fileID,'Physical Line (8) = {125}; \n');
    % Capacitors
    fprintf(fileID,'Physical Line (9) = {174}; \n');    
    fprintf(fileID,'Physical Line (10) = {147}; \n');      
    fprintf(fileID,'Physical Line (11) = {151}; \n');       
    fprintf(fileID,'Physical Line (12) = {156}; \n');      
    fprintf(fileID,'Physical Line (13) = {158}; \n');    
    fprintf(fileID,'Physical Line (14) = {163}; \n');       
    fprintf(fileID,'Physical Line (15) = {167}; \n');        
    fprintf(fileID,'Physical Line (16) = {172}; \n');
    fprintf(fileID,'Physical Line (17) = {131}; \n');    
    fprintf(fileID,'Physical Line (18) = {133}; \n');      
    fprintf(fileID,'Physical Line (19) = {140}; \n');       
    fprintf(fileID,'Physical Line (20) = {142}; \n');      
    fprintf(fileID,'Physical Line (21) = {115}; \n');    
    fprintf(fileID,'Physical Line (22) = {117}; \n');       
    fprintf(fileID,'Physical Line (23) = {124}; \n');        
    fprintf(fileID,'Physical Line (24) = {126}; \n');
    fprintf(fileID,'Physical Line (25) = {189}; \n');    
    fprintf(fileID,'Physical Line (26) = {190}; \n');      
    fprintf(fileID,'Physical Line (27) = {218}; \n');       
    fprintf(fileID,'Physical Line (28) = {217}; \n');      
    fprintf(fileID,'Physical Line (29) = {219}; \n');    
    fprintf(fileID,'Physical Line (30) = {220}; \n');       
    fprintf(fileID,'Physical Line (31) = {185}; \n');        
    fprintf(fileID,'Physical Line (32) = {186}; \n');
    fprintf(fileID,'Physical Line (33) = {187}; \n');    
    fprintf(fileID,'Physical Line (34) = {188}; \n');      
    fprintf(fileID,'Physical Line (35) = {222}; \n');       
    fprintf(fileID,'Physical Line (36) = {221}; \n');      
    fprintf(fileID,'Physical Line (37) = {223}; \n');    
    fprintf(fileID,'Physical Line (38) = {224}; \n');       
    fprintf(fileID,'Physical Line (39) = {192}; \n');        
    fprintf(fileID,'Physical Line (40) = {191}; \n');
    % Mutual Inductors
    fprintf(fileID,'Physical Line (41) = {37}; \n');
    fprintf(fileID,'Physical Line (42) = {38}; \n');
    fprintf(fileID,'Physical Line (43) = {11}; \n');
    fprintf(fileID,'Physical Line (44) = {10}; \n');
    fprintf(fileID,'Physical Line (45) = {45}; \n');
    fprintf(fileID,'Physical Line (46) = {46}; \n');
    fprintf(fileID,'Physical Line (47) = {24}; \n');
    fprintf(fileID,'Physical Line (48) = {23}; \n');
    fprintf(fileID,'Physical Line (49) = {25}; \n');
    fprintf(fileID,'Physical Line (50) = {26}; \n');
    fprintf(fileID,'Physical Line (51) = {15}; \n');
    fprintf(fileID,'Physical Line (52) = {14}; \n');
    fprintf(fileID,'Physical Line (53) = {33}; \n');
    fprintf(fileID,'Physical Line (54) = {34}; \n');
    fprintf(fileID,'Physical Line (55) = {20}; \n');
    fprintf(fileID,'Physical Line (56) = {19}; \n');    
    
    
    %% Physical Surfaces
    fprintf(fileID,'Physical Surface( 1 )  = ');
    fprintf(fileID,'{1,2,3,4,5,6,7,8,9,10,');
    fprintf(fileID, '11,12,13,14,15,16,17,18,19,20,');
    fprintf(fileID, '21,22,23,24,25,26,27,28,29,30,');
    fprintf(fileID, '31,32,33,34,35,36,37,38,39,40,');
    fprintf(fileID, '41,42,43,44,45,46,47,48,49,50,');
    fprintf(fileID, '51,52,53,54,55,56,57,58,59,60,');
    fprintf(fileID, '61,62,63,64,65,66,67,68,69,70,');
    fprintf(fileID, '71,72,73,74,75,76,77,78,79,80,');
    fprintf(fileID, '81,82,83,84,85,86,87,88}; \n\n');
    
    fprintf(fileID,'allSurfaces[] = Surface "*"; \n');
    fprintf(fileID,'Rotate {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{ allSurfaces[] } ; } \n');
    
    fprintf(fileID,'Coherence Mesh; \n\n');
    
    %% Add comments and Close
    fprintf(fileID,'View "comments" { \n');
    fprintf(fileID,'T2(10, -10, 0){ "Copyright (C) Ilias Giannakopoulos" }; \n');
    fprintf(fileID,'T2(10, 15, 0){ StrCat("File created on ", Today) }; }; \n\n');
    fclose(fileID);  
    
    command = sprintf('./data/coils/GMSH/bin/gmsh %s -1 -2 -o %s', geofile,mshfile);
    [status,cmdout] = system(command); 
    
end

function[x1,y1,x2,y2] = Tweak(x0,y0,a,tweak)    

    pos = atan(y0/x0);
    x1 = a*cos(pos+tweak);
    y1 = a*sin(pos+tweak);
    x2 = a*cos(pos-tweak);
    y2 = a*sin(pos-tweak);  
    
end
