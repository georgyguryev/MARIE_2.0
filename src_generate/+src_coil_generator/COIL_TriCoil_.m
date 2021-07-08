function[mshfile] = COIL_TriCoil_(a,le,meshing,shield)

    % Generates a Wiggins triangular coil with 8 ports
    % Author: Ilias I. Giannakopoulos, Moscow, 2018
    
    %% Parameters
    
    % default
    % a = 0.152/2;
    % le = 0.146;
    % meshing = 0.01;
    
    t = 0.01;
    yoff = 0;
    yoff2 = t;
    xoff = t/2;
    li = le-2*t;
    h = li+t;
    l = 2*pi*a;
    len =20;
    xt = sqrt(2)/2*a;
    yt = sqrt(2)/2*a;
    
    if shield.flag ==1
        geofile = strcat('./data/coils/geo_files/Triangular_decoupled_shield_rad_',num2str(a*100),'cm_len_',num2str(le*100),'cm_mesh_',num2str(meshing),'.geo');
        mshfile = strcat('./data/coils/msh_files/Triangular_decoupled_shield_rad_',num2str(a*100),'cm_len_',num2str(le*100),'cm_mesh_',num2str(meshing),'.msh');
    else
        geofile = strcat('./data/coils/geo_files/Triangular_decoupled_rad_',num2str(a*100),'cm_len_',num2str(le*100),'cm_mesh_',num2str(meshing),'.geo');
        mshfile = strcat('./data/coils/msh_files/Triangular_decoupled_rad_',num2str(a*100),'cm_len_',num2str(le*100),'cm_mesh_',num2str(meshing),'.msh');
    end
    
    
    %% Write in .geo file
    fileID = fopen(geofile,'w');
    fprintf(fileID,'lc = %f; \n\n',meshing);   
    
    %% Points
    fprintf(fileID,'Point( 1001 ) = {%f, %f, %f, lc}; \n',[0 ; 0 ; 0-le/2]);
    fprintf(fileID,'Point( 1006 ) = {%f, %f, %f, lc}; \n',[0 ; 0 ; t-le/2]);
    fprintf(fileID,'Point( 1011 ) = {%f, %f, %f, lc}; \n',[0 ; 0 ; h-le/2]);
    fprintf(fileID,'Point( 1016 ) = {%f, %f, %f, lc}; \n',[0 ; 0 ; le-le/2]);
    [x1,y1,x2,y2] = Points_on_circle(sqrt(2)/2*a,sqrt(2)/2*a,a);       
    fprintf(fileID,'Point( 1002 ) = {%f, %f, %f, lc}; \n',[x1; y1; 0-le/2]);
    fprintf(fileID,'Point( 1007 ) = {%f, %f, %f, lc}; \n',[x1; y1; t-le/2]);
    fprintf(fileID,'Point( 1021 ) = {%f, %f, %f, lc}; \n',[x2; y2; 0-le/2]);
    fprintf(fileID,'Point( 1025 ) = {%f, %f, %f, lc}; \n',[x2; y2; t-le/2]);
    [x1,y1,x2,y2] = Points_on_circle(-sqrt(2)/2*a,sqrt(2)/2*a,a);      
    fprintf(fileID,'Point( 1003 ) = {%f, %f, %f, lc}; \n',[x1; y1; 0-le/2]);
    fprintf(fileID,'Point( 1008 ) = {%f, %f, %f, lc}; \n',[x1; y1; t-le/2]);
    fprintf(fileID,'Point( 1022 ) = {%f, %f, %f, lc}; \n',[x2; y2; 0-le/2]);
    fprintf(fileID,'Point( 1026 ) = {%f, %f, %f, lc}; \n',[x2; y2; t-le/2]);
    [x1,y1,x2,y2] = Points_on_circle(-sqrt(2)/2*a,-sqrt(2)/2*a,a);      
    fprintf(fileID,'Point( 1004 ) = {%f, %f, %f, lc}; \n',[x1; y1; 0-le/2]);
    fprintf(fileID,'Point( 1009 ) = {%f, %f, %f, lc}; \n',[x1; y1; t-le/2]);
    fprintf(fileID,'Point( 1023 ) = {%f, %f, %f, lc}; \n',[x2; y2; 0-le/2]);
    fprintf(fileID,'Point( 1027 ) = {%f, %f, %f, lc}; \n',[x2; y2; t-le/2]);
    [x1,y1,x2,y2] = Points_on_circle(sqrt(2)/2*a,-sqrt(2)/2*a,a);        
    fprintf(fileID,'Point( 1005 ) = {%f, %f, %f, lc}; \n',[x1; y1; 0-le/2]);
    fprintf(fileID,'Point( 1010 ) = {%f, %f, %f, lc}; \n',[x1; y1; t-le/2]);
    fprintf(fileID,'Point( 1024 ) = {%f, %f, %f, lc}; \n',[x2; y2; 0-le/2]);
    fprintf(fileID,'Point( 1028 ) = {%f, %f, %f, lc}; \n',[x2; y2; t-le/2]);
    [x1,y1,x2,y2] = Points_on_circle(a,0,a);     
    fprintf(fileID,'Point( 1012 ) = {%f, %f, %f, lc}; \n',[x1; y1; h-le/2]);
    fprintf(fileID,'Point( 1017 ) = {%f, %f, %f, lc}; \n',[x1; y1; le-le/2]);
    fprintf(fileID,'Point( 1029 ) = {%f, %f, %f, lc}; \n',[x2; y2; h-le/2]);
    fprintf(fileID,'Point( 1033 ) = {%f, %f, %f, lc}; \n',[x2; y2; le-le/2]);
    [x1,y1,x2,y2] = Points_on_circle(0,a,a);       
    fprintf(fileID,'Point( 1013 ) = {%f, %f, %f, lc}; \n',[x1 ; y1; h-le/2]);
    fprintf(fileID,'Point( 1018 ) = {%f, %f, %f, lc}; \n',[x1 ; y1; le-le/2]);
    fprintf(fileID,'Point( 1030 ) = {%f, %f, %f, lc}; \n',[x2 ; y2; h-le/2]);
    fprintf(fileID,'Point( 1034 ) = {%f, %f, %f, lc}; \n',[x2 ; y2; le-le/2]);
    [x1,y1,x2,y2] = Points_on_circle(-a,0,a);  
    fprintf(fileID,'Point( 1014 ) = {%f, %f, %f, lc}; \n',[x1; y1; h-le/2]);
    fprintf(fileID,'Point( 1019 ) = {%f, %f, %f, lc}; \n',[x1; y1; le-le/2]);
    fprintf(fileID,'Point( 1031 ) = {%f, %f, %f, lc}; \n',[-x2; -y2; h-le/2]);
    fprintf(fileID,'Point( 1035 ) = {%f, %f, %f, lc}; \n',[-x2; -y2; le-le/2]);
    [x1,y1,x2,y2] = Points_on_circle(0,-a,a);   
    fprintf(fileID,'Point( 1015 ) = {%f, %f, %f, lc}; \n',[x1; y1; h-le/2]);
    fprintf(fileID,'Point( 1020 ) = {%f, %f, %f, lc}; \n',[x1; y1; le-le/2]);    
    fprintf(fileID,'Point( 1032 ) = {%f, %f, %f, lc}; \n',[-x2; -y2; h-le/2]);
    fprintf(fileID,'Point( 1036 ) = {%f, %f, %f, lc}; \n',[-x2; -y2; le-le/2]);    
    
    fprintf(fileID,'Point( 2001 ) = {%f, %f, %f, lc}; \n',[a; 0; 0-le/2]);
    fprintf(fileID,'Point( 3001 ) = {%f, %f, %f, lc}; \n',[a; 0; t-le/2]);
    fprintf(fileID,'Point( 2002 ) = {%f, %f, %f, lc}; \n',[0; a; 0-le/2]);
    fprintf(fileID,'Point( 3002 ) = {%f, %f, %f, lc}; \n',[0; a; t-le/2]);
    fprintf(fileID,'Point( 2003 ) = {%f, %f, %f, lc}; \n',[-a; 0; 0-le/2]);
    fprintf(fileID,'Point( 3003 ) = {%f, %f, %f, lc}; \n',[-a; 0; t-le/2]);
    fprintf(fileID,'Point( 2004 ) = {%f, %f, %f, lc}; \n',[0; -a; 0-le/2]);
    fprintf(fileID,'Point( 3004 ) = {%f, %f, %f, lc}; \n',[0; -a; t-le/2]);
    fprintf(fileID,'Point( 2005 ) = {%f, %f, %f, lc}; \n',[xt; yt; h-le/2]);
    fprintf(fileID,'Point( 3005 ) = {%f, %f, %f, lc}; \n',[xt; yt; le-le/2]);
    fprintf(fileID,'Point( 2006 ) = {%f, %f, %f, lc}; \n',[-xt; yt; h-le/2]);
    fprintf(fileID,'Point( 3006 ) = {%f, %f, %f, lc}; \n',[-xt; yt; le-le/2]);
    fprintf(fileID,'Point( 2007 ) = {%f, %f, %f, lc}; \n',[-xt; -yt; h-le/2]);
    fprintf(fileID,'Point( 3007 ) = {%f, %f, %f, lc}; \n',[-xt; -yt; le-le/2]);
    fprintf(fileID,'Point( 2008 ) = {%f, %f, %f, lc}; \n',[xt; -yt; h-le/2]);
    fprintf(fileID,'Point( 3008 ) = {%f, %f, %f, lc}; \n',[xt; -yt; le-le/2]);
    [x1,y1,x2,y2] = Position_of_capacitors(sqrt(2)/2*a,sqrt(2)/2*a,a);       
    fprintf(fileID,'Point( 4001 ) = {%f, %f, %f, lc}; \n',[x1; y1; 0-le/2]);
    fprintf(fileID,'Point( 5001 ) = {%f, %f, %f, lc}; \n',[x1; y1; t-le/2]);
    fprintf(fileID,'Point( 6001 ) = {%f, %f, %f, lc}; \n',[x2; y2; 0-le/2]);
    fprintf(fileID,'Point( 7001 ) = {%f, %f, %f, lc}; \n',[x2; y2; t-le/2]);
    [x1,y1,x2,y2] = Position_of_capacitors(-sqrt(2)/2*a,sqrt(2)/2*a,a);      
    fprintf(fileID,'Point( 4002 ) = {%f, %f, %f, lc}; \n',[x1; y1; 0-le/2]);
    fprintf(fileID,'Point( 5002 ) = {%f, %f, %f, lc}; \n',[x1; y1; t-le/2]);
    fprintf(fileID,'Point( 6002 ) = {%f, %f, %f, lc}; \n',[x2; y2; 0-le/2]);
    fprintf(fileID,'Point( 7002 ) = {%f, %f, %f, lc}; \n',[x2; y2; t-le/2]);
    [x1,y1,x2,y2] = Position_of_capacitors(-sqrt(2)/2*a,-sqrt(2)/2*a,a);      
    fprintf(fileID,'Point( 4003 ) = {%f, %f, %f, lc}; \n',[x1; y1; 0-le/2]);
    fprintf(fileID,'Point( 5003 ) = {%f, %f, %f, lc}; \n',[x1; y1; t-le/2]);
    fprintf(fileID,'Point( 6003 ) = {%f, %f, %f, lc}; \n',[x2; y2; 0-le/2]);
    fprintf(fileID,'Point( 7003 ) = {%f, %f, %f, lc}; \n',[x2; y2; t-le/2]);
    [x1,y1,x2,y2] = Position_of_capacitors(sqrt(2)/2*a,-sqrt(2)/2*a,a);        
    fprintf(fileID,'Point( 4004 ) = {%f, %f, %f, lc}; \n',[x1; y1; 0-le/2]);
    fprintf(fileID,'Point( 5004 ) = {%f, %f, %f, lc}; \n',[x1; y1; t-le/2]);
    fprintf(fileID,'Point( 6004 ) = {%f, %f, %f, lc}; \n',[x2; y2; 0-le/2]);
    fprintf(fileID,'Point( 7004 ) = {%f, %f, %f, lc}; \n',[x2; y2; t-le/2]);
    [x1,y1,x2,y2] = Position_of_capacitors(a,0,a);     
    fprintf(fileID,'Point( 4005 ) = {%f, %f, %f, lc}; \n',[x1; y1; h-le/2]);
    fprintf(fileID,'Point( 5005 ) = {%f, %f, %f, lc}; \n',[x1; y1; le-le/2]);
    fprintf(fileID,'Point( 6005 ) = {%f, %f, %f, lc}; \n',[x2; y2; h-le/2]);
    fprintf(fileID,'Point( 7005 ) = {%f, %f, %f, lc}; \n',[x2; y2; le-le/2]);
    [x1,y1,x2,y2] = Position_of_capacitors(0,a,a);       
    fprintf(fileID,'Point( 4006 ) = {%f, %f, %f, lc}; \n',[x1 ; y1; h-le/2]);
    fprintf(fileID,'Point( 5006 ) = {%f, %f, %f, lc}; \n',[x1 ; y1; le-le/2]);
    fprintf(fileID,'Point( 6006 ) = {%f, %f, %f, lc}; \n',[x2 ; y2; h-le/2]);
    fprintf(fileID,'Point( 7006 ) = {%f, %f, %f, lc}; \n',[x2 ; y2; le-le/2]);
    [x1,y1,x2,y2] = Position_of_capacitors(-a,0,a);  
    fprintf(fileID,'Point( 4007 ) = {%f, %f, %f, lc}; \n',[x1; y1; h-le/2]);
    fprintf(fileID,'Point( 5007 ) = {%f, %f, %f, lc}; \n',[x1; y1; le-le/2]);
    fprintf(fileID,'Point( 6007 ) = {%f, %f, %f, lc}; \n',[-x2; -y2; h-le/2]);
    fprintf(fileID,'Point( 7007 ) = {%f, %f, %f, lc}; \n',[-x2; -y2; le-le/2]);
    [x1,y1,x2,y2] = Position_of_capacitors(0,-a,a);   
    fprintf(fileID,'Point( 4008 ) = {%f, %f, %f, lc}; \n',[x1; y1; h-le/2]);
    fprintf(fileID,'Point( 5008 ) = {%f, %f, %f, lc}; \n',[x1; y1; le-le/2]);    
    fprintf(fileID,'Point( 6008 ) = {%f, %f, %f, lc}; \n',[-x2; -y2; h-le/2]);
    fprintf(fileID,'Point( 7008 ) = {%f, %f, %f, lc}; \n',[-x2; -y2; le-le/2]);      
    % First Helix top
    %Points
    x1 = 0;
    y1 = h;
    x2 = l/8;
    y2 = t;
    % Transform
    x1 = x1+ xoff;
    y1 = y1 - yoff;
    y2 = y2 + yoff + yoff2;
    y1 = y1 - le/2;
    y2 = y2 - le/2;
    m = (y2-y1)/(x2-x1);    
    helix(le,x1,x2,y1,a,m,4001,5,17,len,fileID);
    % First Helix bottom
    x1 = 0;
    y1 = h;
    x2 = l/8;
    y2 = t;
    % Transform
    y1 = y1 - yoff - yoff2;
    x2 = x2 - xoff;
    y2 = y2 + yoff;
    y1 = y1 - le/2;
    y2 = y2 - le/2;
    m = (y2-y1)/(x2-x1);    
    helix(le,x1,x2,y1,a,m,4002,4,16,len,fileID);   
    % Second Helix top
    %Points
    x1 = l/8;
    y1 = t;
    x2 = l/4;
    y2 = h;
    % Transform
    y1 = y1 + yoff + yoff2;
    x2 = x2 - xoff;
    y2 = y2 - yoff;
    y1 = y1 - le/2;
    y2 = y2 - le/2;
    m = (y2-y1)/(x2-x1);    
    helix(le,x1,x2,y1,a,m,4003,4,16,len,fileID);
    % Second Helix bottom
    x1 = l/8;
    y1 = t;
    x2 = l/4;
    y2 = h;
    % Transform
    x1 = x1 + xoff;
    y1 = y1 + yoff;
    y2 = y2 - yoff - yoff2;
    y1 = y1 - le/2;
    y2 = y2 - le/2;
    m = (y2-y1)/(x2-x1);    
    helix(le,x1,x2,y1,a,m,4004,5,17,len,fileID);
    % Third Helix top
    %Points
    x1 = 0+l/4;
    y1 = h;
    x2 = l/8+l/4;
    y2 = t;
    % Transform
    x1 = x1+ xoff;
    y1 = y1 - yoff;
    y2 = y2 + yoff + yoff2;
    y1 = y1 - le/2;
    y2 = y2 - le/2;
    m = (y2-y1)/(x2-x1);    
    helix(le,x1,x2,y1,a,m,4005,5,17,len,fileID);
    % Third Helix bottom
    x1 = 0+l/4;
    y1 = h;
    x2 = l/8+l/4;
    y2 = t;
    % Transform
    y1 = y1 - yoff - yoff2;
    x2 = x2 - xoff;
    y2 = y2 + yoff;
    y1 = y1 - le/2;
    y2 = y2 - le/2;
    m = (y2-y1)/(x2-x1);    
    helix(le,x1,x2,y1,a,m,4006,4,16,len,fileID);  
    % Fourth Helix top
    %Points
    x1 = l/8+l/4;
    y1 = t;
    x2 = l/4+l/4;
    y2 = h;
    % Transform
    y1 = y1 + yoff + yoff2;
    x2 = x2 - xoff;
    y2 = y2 - yoff;
    y1 = y1 - le/2;
    y2 = y2 - le/2;
    m = (y2-y1)/(x2-x1);    
    helix(le,x1,x2,y1,a,m,4007,4,16,len,fileID);    
    % Fourth Helix bottom
    x1 = l/8+l/4;
    y1 = t;
    x2 = l/4+l/4;
    y2 = h;
    % Transform
    x1 = x1 + xoff;
    y1 = y1 + yoff;
    y2 = y2 - yoff - yoff2;
    y1 = y1 - le/2;
    y2 = y2 - le/2;
    m = (y2-y1)/(x2-x1);    
    helix(le,x1,x2,y1,a,m,4008,5,17,len,fileID);
    % Fifth Helix top
    %Points
    x1 = 0+l/2;
    y1 = h;
    x2 = l/8+l/2;
    y2 = t;
    % Transform
    x1 = x1+ xoff;
    y1 = y1 - yoff;
    y2 = y2 + yoff + yoff2;
    y1 = y1 - le/2;
    y2 = y2 - le/2;
    m = (y2-y1)/(x2-x1);    
    helix(le,x1,x2,y1,a,m,4009,5,17,len,fileID);    
    % Fifth Helix bottom
    x1 = 0+l/2;
    y1 = h;
    x2 = l/8+l/2;
    y2 = t;
    % Transform
    y1 = y1 - yoff - yoff2;
    x2 = x2 - xoff;
    y2 = y2 + yoff;
    y1 = y1 - le/2;
    y2 = y2 - le/2;
    m = (y2-y1)/(x2-x1);    
    helix(le,x1,x2,y1,a,m,4010,4,16,len,fileID);        
    % Sixth Helix top
    %Points
    x1 = l/8+l/2;
    y1 = t;
    x2 = l/4+l/2;
    y2 = h;
    % Transform
    y1 = y1 + yoff + yoff2;
    x2 = x2 - xoff;
    y2 = y2 - yoff;
    y1 = y1 - le/2;
    y2 = y2 - le/2;
    m = (y2-y1)/(x2-x1);    
    helix(le,x1,x2,y1,a,m,4011,4,16,len,fileID);    
    % Sixth Helix bottom
    x1 = l/8+l/2;
    y1 = t;
    x2 = l/4+l/2;
    y2 = h;
    % Transform
    x1 = x1 + xoff;
    y1 = y1 + yoff;
    y2 = y2 - yoff - yoff2;
    y1 = y1 - le/2;
    y2 = y2 - le/2;
    m = (y2-y1)/(x2-x1);    
    helix(le,x1,x2,y1,a,m,4012,5,17,len,fileID);
    % Seven Helix top
    %Points
    x1 = 0+3*l/4;
    y1 = h;
    x2 = l/8+3*l/4;
    y2 = t;
    % Transform
    x1 = x1+ xoff;
    y1 = y1 - yoff;
    y2 = y2 + yoff + yoff2;
    y1 = y1 - le/2;
    y2 = y2 - le/2;
    m = (y2-y1)/(x2-x1);    
    helix(le,x1,x2,y1,a,m,4013,5,17,len,fileID);    
    % Seven Helix bottom
    x1 = 0+3*l/4;
    y1 = h;
    x2 = l/8+3*l/4;
    y2 = t;
    % Transform
    y1 = y1 - yoff - yoff2;
    x2 = x2 - xoff;
    y2 = y2 + yoff;
    y1 = y1 - le/2;
    y2 = y2 - le/2;
    m = (y2-y1)/(x2-x1);    
    helix(le,x1,x2,y1,a,m,4014,4,16,len,fileID);        
    % Eight Helix top
    %Points
    x1 = l/8+3*l/4;
    y1 = t;
    x2 = l/4+3*l/4;
    y2 = h;
    % Transform
    y1 = y1 + yoff + yoff2;
    x2 = x2 - xoff;
    y2 = y2 - yoff;
    y1 = y1 - le/2;
    y2 = y2 - le/2;
    m = (y2-y1)/(x2-x1);    
    helix(le,x1,x2,y1,a,m,4015,4,16,len,fileID);
    % Eight Helix bottom
    x1 = l/8+3*l/4;
    y1 = t;
    x2 = l/4+3*l/4;
    y2 = h;
    % Transform
    x1 = x1 + xoff;
    y1 = y1 + yoff;
    y2 = y2 - yoff - yoff2;
    y1 = y1 - le/2;
    y2 = y2 - le/2;
    m = (y2-y1)/(x2-x1);    
    helix(le,x1,x2,y1,a,m,4016,5,17,len,fileID);
    % Connect the Lines 
    fprintf(fileID,'Circle( 10123 ) = {2008, 1011, 6005}; \n');
    fprintf(fileID,'Circle( 10163 ) = {3008, 1016, 7005}; \n');
    fprintf(fileID,'Circle( 10093 ) = {2005, 1011, 6006}; \n');
    fprintf(fileID,'Circle( 10133 ) = {3005, 1016, 7006}; \n');
    fprintf(fileID,'Circle( 10103 ) = {2006, 1011, 6007}; \n');
    fprintf(fileID,'Circle( 10143 ) = {3006, 1016, 7007}; \n');
    fprintf(fileID,'Circle( 10113 ) = {2007, 1011, 6008}; \n');
    fprintf(fileID,'Circle( 10153 ) = {3007, 1016, 7008}; \n');
    fprintf(fileID,'Circle( 10043 ) = {3001, 1006, 7001}; \n');
    fprintf(fileID,'Circle( 10083 ) = {2001, 1001, 6001}; \n');
    fprintf(fileID,'Circle( 10033 ) = {3004, 1006, 7002}; \n');
    fprintf(fileID,'Circle( 10073 ) = {2004, 1001, 6002}; \n');
    fprintf(fileID,'Circle( 10023 ) = {2003, 1001, 6003}; \n');
    fprintf(fileID,'Circle( 10063 ) = {3003, 1006, 7003}; \n');  
    fprintf(fileID,'Circle( 10013 ) = {2002, 1001, 6004}; \n');
    fprintf(fileID,'Circle( 10053 ) = {3002, 1006, 7004}; \n');  
    fprintf(fileID,'Circle( 10122 ) = {4008, 1011, 2008}; \n'); 
    fprintf(fileID,'Circle( 10162 ) = {5008, 1016, 3008}; \n'); 
    fprintf(fileID,'Circle( 10102 ) = {4006, 1011, 2006}; \n');
    fprintf(fileID,'Circle( 10142 ) = {5006, 1016, 3006}; \n');
    fprintf(fileID,'Circle( 10112 ) = {4007, 1011, 2007}; \n');
    fprintf(fileID,'Circle( 10152 ) = {5007, 1016, 3007}; \n');
    fprintf(fileID,'Circle( 10092 ) = {4005, 1011, 2005}; \n'); 
    fprintf(fileID,'Circle( 10132 ) = {5005, 1016, 3005}; \n');
    fprintf(fileID,'Circle( 10042 ) = {5004, 1006, 3001}; \n');
    fprintf(fileID,'Circle( 10082 ) = {4004, 1001, 2001}; \n');
    fprintf(fileID,'Circle( 10032 ) = {5003, 1006, 3004}; \n');
    fprintf(fileID,'Circle( 10072 ) = {4003, 1001, 2004}; \n');    
    fprintf(fileID,'Circle( 10022 ) = {4002, 1001, 2003}; \n');
    fprintf(fileID,'Circle( 10062 ) = {5002, 1006, 3003}; \n');    
    fprintf(fileID,'Circle( 10012 ) = {4001, 1001, 2002}; \n');
    fprintf(fileID,'Circle( 10052 ) = {5001, 1006, 3002}; \n');
    fprintf(fileID,'Circle( 10011 ) = {10061, 1001, 4001}; \n');
    fprintf(fileID,'Circle( 10014 ) = {6004, 1001, 10120}; \n');
    fprintf(fileID,'Circle( 10021 ) = {10141, 1001, 4002}; \n');
    fprintf(fileID,'Circle( 10024 ) = {6003, 1001, 10200}; \n');
    fprintf(fileID,'Circle( 10031 ) = {221, 1006, 5003}; \n');
    fprintf(fileID,'Circle( 10034 ) = {7002, 1006, 280}; \n');
    fprintf(fileID,'Circle( 10041 ) = {301, 1006, 5004}; \n');
    fprintf(fileID,'Circle( 10044 ) = {7001, 1006, 40}; \n');
    fprintf(fileID,'Circle( 10051 ) = {61, 1006, 5001}; \n');
    fprintf(fileID,'Circle( 10054 ) = {7004, 1006, 120}; \n');
    fprintf(fileID,'Circle( 10061 ) = {141, 1006, 5002}; \n');
    fprintf(fileID,'Circle( 10064 ) = {7003, 1006, 200}; \n');
    fprintf(fileID,'Circle( 10071 ) = {10221, 1001, 4003}; \n');
    fprintf(fileID,'Circle( 10074 ) = {6002, 1001, 10280}; \n');
    fprintf(fileID,'Circle( 10081 ) = {10301, 1001, 4004}; \n');
    fprintf(fileID,'Circle( 10084 ) = {6001, 1001, 10040}; \n');    
    fprintf(fileID,'Circle( 10091 ) = {1, 1011, 4005}; \n'); 
    fprintf(fileID,'Circle( 10094 ) = {6006, 1011, 60}; \n');     
    fprintf(fileID,'Circle( 10101 ) = {81, 1011, 4006}; \n');
    fprintf(fileID,'Circle( 10104 ) = {6007, 1011, 140}; \n'); 
    fprintf(fileID,'Circle( 10111 ) = {161, 1011, 4007}; \n');
    fprintf(fileID,'Circle( 10114 ) = {6008, 1011, 220}; \n');
    fprintf(fileID,'Circle( 10121 ) = {241, 1011, 4008}; \n'); 
    fprintf(fileID,'Circle( 10124 ) = {6005, 1011, 300}; \n'); 
    fprintf(fileID,'Circle( 10131 ) = {10001, 1016, 5005}; \n');
    fprintf(fileID,'Circle( 10134 ) = {7006, 1016, 10060}; \n');    
    fprintf(fileID,'Circle( 10141 ) = {10081, 1016, 5006}; \n');
    fprintf(fileID,'Circle( 10144 ) = {7007, 1016, 10140}; \n');
    fprintf(fileID,'Circle( 10151 ) = {10161, 1016, 5007}; \n');
    fprintf(fileID,'Circle( 10154 ) = {7008, 1016, 10220}; \n');
    fprintf(fileID,'Circle( 10161 ) = {10241, 1016, 5008}; \n'); 
    fprintf(fileID,'Circle( 10164 ) = {7005, 1016, 10300}; \n');   
    %Middle Triangle Lines
    fprintf(fileID,'Circle( 172 ) = {1029, 1011, 1012}; \n');  
    fprintf(fileID,'Circle( 182 ) = {1025,1006, 1007}; \n'); 
    fprintf(fileID,'Circle( 192 ) = {1030, 1011, 1013}; \n'); 
    fprintf(fileID,'Circle( 202 ) = {1028, 1006, 1008}; \n');   
    fprintf(fileID,'Circle( 212 ) = {1031, 1011, 1014}; \n');    
    fprintf(fileID,'Circle( 222 ) = {1027, 1006, 1009}; \n');    
    fprintf(fileID,'Circle( 232 ) = {1032, 1011, 1015}; \n');  
    fprintf(fileID,'Circle( 242 ) = {1026, 1006, 1010}; \n');  
    %Vertical Lines  
    fprintf(fileID,'Line( 1017 ) = {1020,1015}; \n');    
    fprintf(fileID,'Line( 2017 ) = {241  ,10241}; \n');   
    fprintf(fileID,'Line( 1018 ) = {1019,1014}; \n');  
    fprintf(fileID,'Line( 2018 ) = {161   ,10161}; \n');   
    fprintf(fileID,'Line( 1019 ) = {1018,1013}; \n');  
    fprintf(fileID,'Line( 2019 ) = {81  ,10081}; \n');     
    fprintf(fileID,'Line( 1020 ) = {1017,1012}; \n');
    fprintf(fileID,'Line( 2020 ) = {1     ,10001}; \n'); 
    fprintf(fileID,'Line( 1021 ) = {1009,1004}; \n');  
    fprintf(fileID,'Line( 2021 ) = {221 ,10221}; \n');    
    fprintf(fileID,'Line( 1022 ) = {1010,1005}; \n');   
    fprintf(fileID,'Line( 2022 ) = {301,10301}; \n');   
    fprintf(fileID,'Line( 1023 ) = {1007,1002}; \n');    
    fprintf(fileID,'Line( 2023 ) = {61    ,10061}; \n');  
    fprintf(fileID,'Line( 1024 ) = {1008,1003}; \n');     
    fprintf(fileID,'Line( 2024 ) = {141  ,10141}; \n');     
    fprintf(fileID,'Line( 1025 ) = {1036,1032}; \n');
    fprintf(fileID,'Line( 2025 ) = {220  ,10220}; \n');       
    fprintf(fileID,'Line( 1026 ) = {1031,1035}; \n');     
    fprintf(fileID,'Line( 2026 ) = {140   ,10140}; \n');
    fprintf(fileID,'Line( 1027 ) = {1034,1030}; \n');    
    fprintf(fileID,'Line( 2027 ) = {60  ,10060}; \n'); 
    fprintf(fileID,'Line( 1028 ) = {1033,1029}; \n');    
    fprintf(fileID,'Line( 2028 ) = {300  ,10300}; \n');  
    fprintf(fileID,'Line( 1029 ) = {1027,1023}; \n');    
    fprintf(fileID,'Line( 2029 ) = {200  ,10200}; \n');   
    fprintf(fileID,'Line( 1030 ) = {1028,1024}; \n');    
    fprintf(fileID,'Line( 2030 ) = {120  ,10120}; \n'); 
    fprintf(fileID,'Line( 1031 ) = {1025,1021}; \n');   
    fprintf(fileID,'Line( 2031 ) = {40     ,10040}; \n'); 
    fprintf(fileID,'Line( 1032 ) = {1026,1022}; \n');
    fprintf(fileID,'Line( 2032 ) = {280,10280}; \n');    
    %Connection Lines with Helices  
    fprintf(fileID,'Line( 25 ) = {21, 1}; \n');    
    fprintf(fileID,'Line( 26 ) = {20, 61}; \n');    
    fprintf(fileID,'Line( 27 ) = {80, 81}; \n');    
    fprintf(fileID,'Line( 28 ) = {100, 141}; \n');    
    fprintf(fileID,'Line( 29 ) = {160, 161}; \n');    
    fprintf(fileID,'Line( 30 ) = {180, 221}; \n');    
    fprintf(fileID,'Line( 31 ) = {240, 241}; \n');    
    fprintf(fileID,'Line( 32 ) = {260, 301}; \n');
    fprintf(fileID,'Line( 33 ) = {300, 21}; \n');    
    fprintf(fileID,'Line( 34 ) = {40, 20}; \n');    
    fprintf(fileID,'Line( 35 ) = {60, 80}; \n');    
    fprintf(fileID,'Line( 36 ) = {120, 100}; \n');    
    fprintf(fileID,'Line( 37 ) = {140, 160}; \n');    
    fprintf(fileID,'Line( 38 ) = {200, 180}; \n');    
    fprintf(fileID,'Line( 39 ) = {220, 240}; \n');    
    fprintf(fileID,'Line( 40 ) = {280, 260}; \n');
    %Lines for ports
    fprintf(fileID,'Line( 412 ) = {2001, 3001}; \n');    
    fprintf(fileID,'Line( 422 ) = {2002, 3002}; \n');      
    fprintf(fileID,'Line( 432 ) = {2003, 3003}; \n');       
    fprintf(fileID,'Line( 442 ) = {2004, 3004}; \n');        
    fprintf(fileID,'Line( 452 ) = {2005, 3005}; \n');    
    fprintf(fileID,'Line( 462 ) = {2006, 3006}; \n');       
    fprintf(fileID,'Line( 472 ) = {2007, 3007}; \n');        
    fprintf(fileID,'Line( 482 ) = {2008, 3008}; \n');
    %Lines for the capacitors on the cylinder
    fprintf(fileID,'Line( 49 ) = {4001, 5001}; \n');    
    fprintf(fileID,'Line( 50 ) = {4002, 5002}; \n');      
    fprintf(fileID,'Line( 51 ) = {4003, 5003}; \n');       
    fprintf(fileID,'Line( 52 ) = {4004, 5004}; \n');        
    fprintf(fileID,'Line( 53 ) = {4005, 5005}; \n');    
    fprintf(fileID,'Line( 54 ) = {4006, 5006}; \n');       
    fprintf(fileID,'Line( 55 ) = {4007, 5007}; \n');        
    fprintf(fileID,'Line( 56 ) = {4008, 5008}; \n');    
    fprintf(fileID,'Line( 57 ) = {6001, 7001}; \n');    
    fprintf(fileID,'Line( 58 ) = {6002, 7002}; \n');      
    fprintf(fileID,'Line( 59 ) = {6003, 7003}; \n');       
    fprintf(fileID,'Line( 60 ) = {6004, 7004}; \n');        
    fprintf(fileID,'Line( 61 ) = {6005, 7005}; \n');    
    fprintf(fileID,'Line( 62 ) = {6006, 7006}; \n');       
    fprintf(fileID,'Line( 63 ) = {6007, 7007}; \n');        
    fprintf(fileID,'Line( 64 ) = {6008, 7008}; \n');    
    fprintf(fileID,'Circle( 171 ) = {300, 1011, 1029}; \n');       
    fprintf(fileID,'Circle( 173 ) = {1012, 1011, 1}; \n');   
    fprintf(fileID,'Circle( 174 ) = {10300, 1016, 1033}; \n');       
    fprintf(fileID,'Circle( 175 ) = {1017, 1016, 10001}; \n');   
    fprintf(fileID,'Circle( 181 ) = {40,1006, 1025}; \n');     
    fprintf(fileID,'Circle( 183 ) = {1007,1006, 61}; \n');  
    fprintf(fileID,'Circle( 184 ) = {10040,1001, 1021}; \n');     
    fprintf(fileID,'Circle( 185 ) = {1002,1001, 10061}; \n');
    fprintf(fileID,'Circle( 191 ) = {60, 1011, 1030}; \n');     
    fprintf(fileID,'Circle( 193 ) = {1013, 1011, 81}; \n');  
    fprintf(fileID,'Circle( 194 ) = {10060, 1016, 1034}; \n');     
    fprintf(fileID,'Circle( 195 ) = {1018, 1016, 10081}; \n');  
    fprintf(fileID,'Circle( 201 ) = {120, 1006, 1028}; \n'); 
    fprintf(fileID,'Circle( 203 ) = {1008, 1006, 141}; \n');   
    fprintf(fileID,'Circle( 204 ) = {10120, 1001, 1024}; \n'); 
    fprintf(fileID,'Circle( 205 ) = {1003, 1001, 10141}; \n');   
    fprintf(fileID,'Circle( 211 ) = {140, 1011, 1031}; \n');    
    fprintf(fileID,'Circle( 213 ) = {1014, 1011, 161}; \n');
    fprintf(fileID,'Circle( 214 ) = {10140, 1016, 1035}; \n');    
    fprintf(fileID,'Circle( 215 ) = {1019, 1016, 10161}; \n');
    fprintf(fileID,'Circle( 221 ) = {200, 1006, 1027}; \n'); 
    fprintf(fileID,'Circle( 223 ) = {1009, 1006, 221}; \n');   
    fprintf(fileID,'Circle( 224 ) = {10200, 1001, 1023}; \n'); 
    fprintf(fileID,'Circle( 225 ) = {1004, 1001, 10221}; \n'); 
    fprintf(fileID,'Circle( 231 ) = {220, 1011, 1032}; \n');     
    fprintf(fileID,'Circle( 233 ) = {1015, 1011, 241}; \n'); 
    fprintf(fileID,'Circle( 234 ) = {10220, 1016, 1036}; \n');     
    fprintf(fileID,'Circle( 235 ) = {1020, 1016, 10241}; \n'); 
    fprintf(fileID,'Circle( 241 ) = {280, 1006, 1026}; \n');  
    fprintf(fileID,'Circle( 243 ) = {1010, 1006, 301}; \n'); 
    fprintf(fileID,'Circle( 244 ) = {10280, 1001, 1022}; \n');  
    fprintf(fileID,'Circle( 245 ) = {1005, 1001, 10301}; \n');
    fprintf(fileID,'Line( 3001 ) = {5,24}; \n');
    fprintf(fileID,'Line( 3002 ) = {17,36}; \n'); 
    fprintf(fileID,'Line( 3003 ) = {65,44}; \n');
    fprintf(fileID,'Line( 3004 ) = {56,77}; \n');
    fprintf(fileID,'Line( 3005 ) = {97,116}; \n');
    fprintf(fileID,'Line( 3006 ) = {104,85}; \n');   
    fprintf(fileID,'Line( 3007 ) = {145,124}; \n'); 
    fprintf(fileID,'Line( 3008 ) = {136,157}; \n');
    fprintf(fileID,'Line( 3009 ) = {196,177}; \n'); 
    fprintf(fileID,'Line( 3010 ) = {165,184}; \n'); 
    fprintf(fileID,'Line( 3011 ) = {216,237}; \n');
    fprintf(fileID,'Line( 3012 ) = {225,204}; \n');
    fprintf(fileID,'Line( 3013 ) = {245,264}; \n'); 
    fprintf(fileID,'Line( 3014 ) = {276,257}; \n');
    fprintf(fileID,'Line( 3015 ) = {296,317}; \n'); 
    fprintf(fileID,'Line( 3016 ) = {305,284}; \n');   
    %Create Line Loops    
    %Triangles
    fprintf(fileID,'Line Loop( 101 ) = {-38,-30,221,222,223}; \n');
    fprintf(fileID,'Line Loop( 102 ) = {39,-231,-232,-233,31}; \n');
    fprintf(fileID,'Line Loop( 103 ) = {-40,241,242,243,-32}; \n');
    fprintf(fileID,'Line Loop( 104 ) = {33,-171,-172,-173,25}; \n');
    fprintf(fileID,'Line Loop( 105 ) = {-34,181,182,183,-26}; \n');
    fprintf(fileID,'Line Loop( 106 ) = {35,-191,-192,-193,27}; \n');
    fprintf(fileID,'Line Loop( 107 ) = {-36,-28,201,202,203}; \n');
    fprintf(fileID,'Line Loop( 108 ) = {37,-211,-212,-213,29}; \n');
    %Helices, broken to three surfaces to include capacitors
    fprintf(fileID,'Line Loop( 109 ) = {-4007,3007,4008,28}; \n');
    fprintf(fileID,'Line Loop( 110 ) = {-4037,-3007,4038,-3008}; \n');
    fprintf(fileID,'Line Loop( 111 ) = {-4067,3008,4068,-37}; \n');
    fprintf(fileID,'Line Loop( 112 ) = {4010,-3010,-4009,-29}; \n');
    fprintf(fileID,'Line Loop( 113 ) = {4040,3010,-4039,3009}; \n');
    fprintf(fileID,'Line Loop( 114 ) = {4070,-3009,-4069,38}; \n');
    fprintf(fileID,'Line Loop( 115 ) = {-4011,3012,4012,30}; \n');  
    fprintf(fileID,'Line Loop( 116 ) = {-4041,-3012,4042,-3011}; \n');    
    fprintf(fileID,'Line Loop( 117 ) = {-4071,3011,4072,-39}; \n');  
    fprintf(fileID,'Line Loop( 118 ) = {-31,-4013,-3013,4014}; \n');
    fprintf(fileID,'Line Loop( 119 ) = {3013,-4043,3014,4044}; \n');
    fprintf(fileID,'Line Loop( 120 ) = {-3014,-4073,40,4074}; \n');
    fprintf(fileID,'Line Loop( 121 ) = {32,4016,3016,-4015}; \n');    
    fprintf(fileID,'Line Loop( 122 ) = {-3016,4046,-3015,-4045}; \n');    
    fprintf(fileID,'Line Loop( 123 ) = {3015,4076,-33,-4075}; \n'); 
    fprintf(fileID,'Line Loop( 124 ) = {-25,-4001,-3001,4002}; \n');
    fprintf(fileID,'Line Loop( 125 ) = {3001,-4031,-3002,4032}; \n');
    fprintf(fileID,'Line Loop( 126 ) = {3002,-4061,34,4062}; \n');
    fprintf(fileID,'Line Loop( 127 ) = {26,4004,3003,-4003}; \n');
    fprintf(fileID,'Line Loop( 128 ) = {-3003,4034,-3004,-4033}; \n');
    fprintf(fileID,'Line Loop( 129 ) = {3004,4064,-35,-4063}; \n');
    fprintf(fileID,'Line Loop( 130 ) = {-27,-4005,3006,4006}; \n');
    fprintf(fileID,'Line Loop( 131 ) = {-3006,-4035,-3005,4036}; \n');
    fprintf(fileID,'Line Loop( 132 ) = {3005,-4065,36,4066}; \n');
    %Circles, broken to five surfaces to include capacitors and ports 
    fprintf(fileID,'Line Loop( 133 ) = {10091,53,-10131,-2020}; \n');
    fprintf(fileID,'Line Loop( 134 ) = {10092,452,-10132,-53}; \n');
    fprintf(fileID,'Line Loop( 135 ) = {10093,62,-10133,-452}; \n');
    fprintf(fileID,'Line Loop( 136 ) = {10094,2027,-10134,-62}; \n');
    fprintf(fileID,'Line Loop( 137 ) = {-10141,-2019,10101,54}; \n');
    fprintf(fileID,'Line Loop( 138 ) = {-10142,-54,10102,462}; \n');
    fprintf(fileID,'Line Loop( 139 ) = {-10143,-462,10103,63}; \n');
    fprintf(fileID,'Line Loop( 140 ) = {-10144,-63,10104,2026}; \n');
    fprintf(fileID,'Line Loop( 141 ) = {10124,2028,-10164,-61}; \n');
    fprintf(fileID,'Line Loop( 142 ) = {10123,61,-10163,-482}; \n');
    fprintf(fileID,'Line Loop( 143 ) = {10122,482,-10162,-56}; \n');
    fprintf(fileID,'Line Loop( 144 ) = {10121,56,-10161,-2017}; \n');
    fprintf(fileID,'Line Loop( 145 ) = {-10151,-2018,10111,55}; \n');
    fprintf(fileID,'Line Loop( 146 ) = {-10152,-55,10112,472}; \n');
    fprintf(fileID,'Line Loop( 147 ) = {-10153,-472,10113,64}; \n');
    fprintf(fileID,'Line Loop( 148 ) = {-10154,-64,10114,2025}; \n');
    fprintf(fileID,'Line Loop( 149 ) = {-10031,2021,10071,51}; \n');
    fprintf(fileID,'Line Loop( 150 ) = {-10032,-51,10072,442}; \n');
    fprintf(fileID,'Line Loop( 151 ) = {-10033,-442,10073,58}; \n');
    fprintf(fileID,'Line Loop( 152 ) = {-10034,-58,10074,-2032}; \n');
    fprintf(fileID,'Line Loop( 153 ) = {-10044,-2031,10084,-57}; \n');   
    fprintf(fileID,'Line Loop( 154 ) = {-10043,57,10083,-412}; \n');
    fprintf(fileID,'Line Loop( 155 ) = {-10042,412,10082,-52}; \n');
    fprintf(fileID,'Line Loop( 156 ) = {-10041,52,10081,2022}; \n');
    fprintf(fileID,'Line Loop( 157 ) = {2023,-10051,49,10011}; \n');
    fprintf(fileID,'Line Loop( 158 ) = {-49,-10052,422,10012}; \n');
    fprintf(fileID,'Line Loop( 159 ) = {-422,-10053,60,10013}; \n'); 
    fprintf(fileID,'Line Loop( 160 ) = {-60,-10054,-2030,10014}; \n');
    fprintf(fileID,'Line Loop( 161 ) = {2024,-10061,50,10021}; \n');
    fprintf(fileID,'Line Loop( 162 ) = {-50,-10062,432,10022}; \n');
    fprintf(fileID,'Line Loop( 163 ) = {-432,-10063,59,10023}; \n');
    fprintf(fileID,'Line Loop( 164 ) = {-59,-10064,-2029,10024}; \n');
    %Othogonals
    fprintf(fileID,'Line Loop( 165 ) = {1021,225,-2021,-223}; \n');
    fprintf(fileID,'Line Loop( 166 ) = {1022,245,-2022,-243}; \n');
    fprintf(fileID,'Line Loop( 167 ) = {1023,185,-2023,-183}; \n');   
    fprintf(fileID,'Line Loop( 168 ) = {1024,205,-2024,-203}; \n'); 
    fprintf(fileID,'Line Loop( 169 ) = {-1029,224,2029,-221}; \n');
    fprintf(fileID,'Line Loop( 170 ) = {-1032,244,2032,-241}; \n');
    fprintf(fileID,'Line Loop( 171 ) = {-1031,184,2031,-181}; \n');
    fprintf(fileID,'Line Loop( 172 ) = {-1030,204,2030,-201}; \n');  
    fprintf(fileID,'Line Loop( 173 ) = {-1028,-174,-2028,171}; \n');
    fprintf(fileID,'Line Loop( 174 ) = {-1027,-194,-2027,191}; \n');
    fprintf(fileID,'Line Loop( 175 ) = {1026,-214,-2026,211}; \n');
    fprintf(fileID,'Line Loop( 176 ) = {-1025,-234,-2025,231}; \n'); 
    fprintf(fileID,'Line Loop( 177 ) = {1017,233, 2017,-235}; \n');
    fprintf(fileID,'Line Loop( 178 ) = {1020,173, 2020,-175}; \n');
    fprintf(fileID,'Line Loop( 179 ) = {1019,193, 2019,-195}; \n');
    fprintf(fileID,'Line Loop( 180 ) = {1018,213, 2018,-215}; \n');
    fprintf(fileID,'\n\n');  
    % Create Surfaces   
    fprintf(fileID,'Plane Surface( 101 ) = {101}; \n');
    fprintf(fileID,'Plane Surface( 102 ) = {102}; \n');
    fprintf(fileID,'Plane Surface( 103 ) = {103}; \n');
    fprintf(fileID,'Plane Surface( 104 ) = {104}; \n');
    fprintf(fileID,'Plane Surface( 105 ) = {105}; \n');
    fprintf(fileID,'Plane Surface( 106 ) = {106}; \n');
    fprintf(fileID,'Plane Surface( 107 ) = {107}; \n');
    fprintf(fileID,'Plane Surface( 108 ) = {108}; \n');
    fprintf(fileID,'Ruled Surface( 109 ) = {109}; \n');
    fprintf(fileID,'Ruled Surface( 110 ) = {110}; \n');
    fprintf(fileID,'Ruled Surface( 111 ) = {111}; \n');
    fprintf(fileID,'Ruled Surface( 112 ) = {112}; \n');
    fprintf(fileID,'Ruled Surface( 113 ) = {113}; \n');
    fprintf(fileID,'Ruled Surface( 114 ) = {114}; \n');
    fprintf(fileID,'Ruled Surface( 115 ) = {115}; \n');
    fprintf(fileID,'Ruled Surface( 116 ) = {116}; \n');
    fprintf(fileID,'Ruled Surface( 117 ) = {117}; \n');
    fprintf(fileID,'Ruled Surface( 118 ) = {118}; \n');
    fprintf(fileID,'Ruled Surface( 119 ) = {119}; \n');
    fprintf(fileID,'Ruled Surface( 120 ) = {120}; \n');
    fprintf(fileID,'Ruled Surface( 121 ) = {121}; \n');
    fprintf(fileID,'Ruled Surface( 122 ) = {122}; \n');
    fprintf(fileID,'Ruled Surface( 123 ) = {123}; \n');
    fprintf(fileID,'Ruled Surface( 124 ) = {124}; \n');
    fprintf(fileID,'Ruled Surface( 125 ) = {125}; \n');
    fprintf(fileID,'Ruled Surface( 126 ) = {126}; \n');
    fprintf(fileID,'Ruled Surface( 127 ) = {127}; \n');
    fprintf(fileID,'Ruled Surface( 128 ) = {128}; \n');
    fprintf(fileID,'Ruled Surface( 129 ) = {129}; \n');
    fprintf(fileID,'Ruled Surface( 130 ) = {130}; \n');
    fprintf(fileID,'Ruled Surface( 131 ) = {131}; \n');
    fprintf(fileID,'Ruled Surface( 132 ) = {132}; \n');
    fprintf(fileID,'Ruled Surface( 133 ) = {133}; \n');
    fprintf(fileID,'Ruled Surface( 134 ) = {134}; \n');
    fprintf(fileID,'Ruled Surface( 135 ) = {135}; \n');
    fprintf(fileID,'Ruled Surface( 136 ) = {136}; \n');
    fprintf(fileID,'Ruled Surface( 137 ) = {137}; \n');
    fprintf(fileID,'Ruled Surface( 138 ) = {138}; \n');
    fprintf(fileID,'Ruled Surface( 139 ) = {139}; \n');
    fprintf(fileID,'Ruled Surface( 140 ) = {140}; \n');
    fprintf(fileID,'Ruled Surface( 141 ) = {141}; \n');
    fprintf(fileID,'Ruled Surface( 142 ) = {142}; \n');
    fprintf(fileID,'Ruled Surface( 143 ) = {143}; \n');
    fprintf(fileID,'Ruled Surface( 144 ) = {144}; \n');
    fprintf(fileID,'Ruled Surface( 145 ) = {145}; \n');
    fprintf(fileID,'Ruled Surface( 146 ) = {146}; \n');
    fprintf(fileID,'Ruled Surface( 147 ) = {147}; \n');
    fprintf(fileID,'Ruled Surface( 148 ) = {148}; \n');
    fprintf(fileID,'Ruled Surface( 149 ) = {149}; \n');
    fprintf(fileID,'Ruled Surface( 150 ) = {150}; \n');
    fprintf(fileID,'Ruled Surface( 151 ) = {151}; \n');
    fprintf(fileID,'Ruled Surface( 152 ) = {152}; \n');
    fprintf(fileID,'Ruled Surface( 153 ) = {153}; \n');
    fprintf(fileID,'Ruled Surface( 154 ) = {154}; \n');
    fprintf(fileID,'Ruled Surface( 155 ) = {155}; \n');
    fprintf(fileID,'Ruled Surface( 156 ) = {156}; \n');
    fprintf(fileID,'Ruled Surface( 157 ) = {157}; \n');
    fprintf(fileID,'Ruled Surface( 158 ) = {158}; \n');
    fprintf(fileID,'Ruled Surface( 159 ) = {159}; \n');
    fprintf(fileID,'Ruled Surface( 160 ) = {160}; \n');
    fprintf(fileID,'Ruled Surface( 161 ) = {161}; \n');
    fprintf(fileID,'Ruled Surface( 162 ) = {162}; \n');
    fprintf(fileID,'Ruled Surface( 163 ) = {163}; \n');
    fprintf(fileID,'Ruled Surface( 164 ) = {164}; \n');
    fprintf(fileID,'Ruled Surface( 165 ) = {165}; \n');
    fprintf(fileID,'Ruled Surface( 166 ) = {166}; \n');
    fprintf(fileID,'Ruled Surface( 167 ) = {167}; \n');
    fprintf(fileID,'Ruled Surface( 168 ) = {168}; \n');
    fprintf(fileID,'Ruled Surface( 169 ) = {169}; \n');
    fprintf(fileID,'Ruled Surface( 170 ) = {170}; \n');
    fprintf(fileID,'Ruled Surface( 171 ) = {171}; \n');
    fprintf(fileID,'Ruled Surface( 172 ) = {172}; \n');
    fprintf(fileID,'Ruled Surface( 173 ) = {173}; \n');
    fprintf(fileID,'Ruled Surface( 174 ) = {174}; \n');
    fprintf(fileID,'Ruled Surface( 175 ) = {175}; \n');
    fprintf(fileID,'Ruled Surface( 176 ) = {176}; \n');
    fprintf(fileID,'Ruled Surface( 177 ) = {177}; \n');
    fprintf(fileID,'Ruled Surface( 178 ) = {178}; \n');
    fprintf(fileID,'Ruled Surface( 179 ) = {179}; \n');
    fprintf(fileID,'Ruled Surface( 180 ) = {180}; \n');
    %% Physical Surface/Lines for Mesh  
    % Sources
    fprintf(fileID,'Physical Line (1) = {412}; \n');    
    fprintf(fileID,'Physical Line (2) = {452}; \n');      
    fprintf(fileID,'Physical Line (3) = {422}; \n');       
    fprintf(fileID,'Physical Line (4) = {462}; \n');      
    fprintf(fileID,'Physical Line (5) = {432}; \n');    
    fprintf(fileID,'Physical Line (6) = {472}; \n');       
    fprintf(fileID,'Physical Line (7) = {442}; \n');        
    fprintf(fileID,'Physical Line (8) = {482}; \n');
    % Capacitors on Cylinders
    fprintf(fileID,'Physical Line (9) = {49}; \n');    
    fprintf(fileID,'Physical Line (10) = {50}; \n');      
    fprintf(fileID,'Physical Line (11) = {51}; \n');       
    fprintf(fileID,'Physical Line (12) = {52}; \n');        
    fprintf(fileID,'Physical Line (13) = {53}; \n');    
    fprintf(fileID,'Physical Line (14) = {54}; \n');       
    fprintf(fileID,'Physical Line (15) = {55}; \n');        
    fprintf(fileID,'Physical Line (16) = {56}; \n');    
    fprintf(fileID,'Physical Line (17) = {57}; \n');    
    fprintf(fileID,'Physical Line (18) = {58}; \n');      
    fprintf(fileID,'Physical Line (19) = {59}; \n');       
    fprintf(fileID,'Physical Line (20) = {60}; \n');        
    fprintf(fileID,'Physical Line (21) = {61}; \n');    
    fprintf(fileID,'Physical Line (22) = {62}; \n');       
    fprintf(fileID,'Physical Line (23) = {63}; \n');
    fprintf(fileID,'Physical Line (24) = {64}; \n');    
    % Inductors
    fprintf(fileID,'Physical Line (25) = {171}; \n');       
    fprintf(fileID,'Physical Line (26) = {173}; \n');   
    fprintf(fileID,'Physical Line (27) = {181}; \n');     
    fprintf(fileID,'Physical Line (28) = {183}; \n');  
    fprintf(fileID,'Physical Line (29) = {191}; \n');     
    fprintf(fileID,'Physical Line (30) = {193}; \n');  
    fprintf(fileID,'Physical Line (31) = {201}; \n'); 
    fprintf(fileID,'Physical Line (32) = {203}; \n');   
    fprintf(fileID,'Physical Line (33) = {211}; \n');    
    fprintf(fileID,'Physical Line (34) = {213}; \n');
    fprintf(fileID,'Physical Line (35) = {221}; \n'); 
    fprintf(fileID,'Physical Line (36) = {223}; \n');   
    fprintf(fileID,'Physical Line (37) = {231}; \n');     
    fprintf(fileID,'Physical Line (38) = {233}; \n'); 
    fprintf(fileID,'Physical Line (39) = {241}; \n');  
    fprintf(fileID,'Physical Line (40) = {243}; \n'); 
    % Capacitors on Helices
    fprintf(fileID,'Physical Line (41) = {3001}; \n');
    fprintf(fileID,'Physical Line (42) = {3002}; \n'); 
    fprintf(fileID,'Physical Line (43) = {3003}; \n');
    fprintf(fileID,'Physical Line (44) = {3004}; \n');
    fprintf(fileID,'Physical Line (45) = {3005}; \n');
    fprintf(fileID,'Physical Line (46) = {3006}; \n');   
    fprintf(fileID,'Physical Line (47) = {3007}; \n'); 
    fprintf(fileID,'Physical Line (48) = {3008}; \n');
    fprintf(fileID,'Physical Line (49) = {3009}; \n'); 
    fprintf(fileID,'Physical Line (50) = {3010}; \n'); 
    fprintf(fileID,'Physical Line (51) = {3011}; \n');
    fprintf(fileID,'Physical Line (52) = {3012}; \n');
    fprintf(fileID,'Physical Line (53) = {3013}; \n'); 
    fprintf(fileID,'Physical Line (54) = {3014}; \n');
    fprintf(fileID,'Physical Line (55) = {3015}; \n'); 
    fprintf(fileID,'Physical Line (56) = {3016}; \n');   
    
    
    if shield.flag == 1
        fprintf(fileID,'lc_s = %f; \n\n',shield.meshing);   
        
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[99990+1  ; 0    ; 0    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[99990+2  ; 0    ; 0    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[99990+3  ; shield.radius    ; 0    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[99990+4  ; 0    ; shield.radius    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[99990+5  ; -shield.radius    ; 0    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[99990+6  ; 0    ; -shield.radius    ; shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[99990+7  ; shield.radius    ; 0    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[99990+8  ; 0    ; shield.radius    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[99990+9  ; -shield.radius    ; 0    ; -shield.length/2]);
        fprintf(fileID,'Point( %d )   = {%f, %f, %f, lc_s}; \n',[99990+10  ; 0    ; -shield.radius    ; -shield.length/2]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[99990+1; 99990+3;99990+7 ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[99990+2; 99990+4;99990+8 ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[99990+3; 99990+5;99990+9 ]);
        fprintf(fileID,'Line( %d ) = {%d,%d}; \n',[99990+4; 99990+6;99990+10 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[99990+5; 99990+3;99990+1;99990+4 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[99990+6; 99990+4;99990+1;99990+5 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[99990+7; 99990+5;99990+1;99990+6 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[99990+8; 99990+6;99990+1;99990+3 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[99990+9; 99990+7;99990+2;99990+8 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[99990+10; 99990+8;99990+2;99990+9 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[99990+11; 99990+9;99990+2;99990+10 ]);
        fprintf(fileID,'Circle( %d ) = {%d,%d,%d}; \n',[99990+12; 99990+10;99990+2;99990+7 ]);
        
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[99990+1;   99990+1 ;   99990+9    ; -(99990+2); -(99990+5)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[99990+2;   99990+2 ;   99990+10    ; -(99990+3); -(99990+6)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[99990+3;   99990+3 ;   99990+11    ; -(99990+4); -(99990+7)]);
        fprintf(fileID,'Line Loop( %d ) = {%d,%d,%d,%d}; \n',[99990+4;   99990+4 ;   99990+12    ; -(99990+1); -(99990+8)]);
        
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[99990+1  ;99990+1  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[99990+2  ;99990+2  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[99990+3  ;99990+3  ]);
        fprintf(fileID,'Ruled Surface( %d ) = {%d}; \n',[99990+4  ;99990+4  ]);
        
        fprintf(fileID,'Physical Surface(2)  = {99990+1,99990+2,99990+3,99990+4};\n\n');
        
    end
    
    % Physical Surfaces
    fprintf(fileID,'Physical Surface( 1 )  = ');
    fprintf(fileID,'{101,102,103,104,105,106,107,108,109,110,');
    fprintf(fileID, '111,112,113,114,115,116,117,118,119,120,');
    fprintf(fileID, '121,122,123,124,125,126,127,128,129,130,');
    fprintf(fileID, '131,132,133,134,135,136,137,138,139,140,');
    fprintf(fileID, '141,142,143,144,145,146,147,148,149,150,');
    fprintf(fileID, '151,152,153,154,155,156,157,158,159,160,');
    fprintf(fileID, '161,162,163,164,165,166,167,168,169,170,');
    fprintf(fileID, '171,172,173,174,175,176,177,178,179,180}; \n\n');
    fprintf(fileID,'Coherence Mesh; \n\n');
    
    %% Add comments and close
    fprintf(fileID,'View "comments" { \n');
    fprintf(fileID,'T2(10, -10, 0){ "Copyright (C) Ilias Giannakopoulos" }; \n');
    fprintf(fileID,'T2(10, 15, 0){ StrCat("File created on ", Today) }; }; \n\n');
    fclose(fileID);  
    
%     command = sprintf('./data/coils/GMSH/bin/gmsh %s', geofile);
    command = sprintf('./data/coils/GMSH/bin/gmsh %s -1 -2 -o %s', geofile,mshfile);
    [status,cmdout] = system(command); 
end

function[] = helix(le,x1,x2,y1,a,m,spline,step1,step2,len,fileID)

    x = linspace(x1,x2,len);
    y = m*(x-x1)+y1;
    xp = a*cos(x/a);
    yp = a*sin(x/a);
    zp = y;
    fprintf(fileID,'\n\n');
    for i = 1:length(x)
        fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i+(spline-4000-1)*length(x);xp(i);yp(i);zp(i)]);
        
         if i+(spline-4000-1)*length(x) == 301
            fprintf(fileID,'Point( 10301 ) = {%f, %f, %f, lc}; \n',[xp(i);yp(i);0-le/2]);
        elseif i+(spline-4000-1)*length(x) == 81
            fprintf(fileID,'Point( 10081 ) = {%f, %f, %f, lc}; \n',[xp(i);yp(i);le-le/2]);
        elseif i+(spline-4000-1)*length(x) == 280
            fprintf(fileID,'Point( 10280 ) = {%f, %f, %f, lc}; \n',[xp(i);yp(i);0-le/2]);
        elseif i+(spline-4000-1)*length(x) == 60 
            fprintf(fileID,'Point( 10060 ) = {%f, %f, %f, lc}; \n',[xp(i);yp(i);le-le/2]);
        elseif i+(spline-4000-1)*length(x) ==  221 
            fprintf(fileID,'Point( 10221 ) = {%f, %f, %f, lc}; \n',[xp(i);yp(i);0-le/2]);
        elseif i+(spline-4000-1)*length(x) ==  1   
            fprintf(fileID,'Point( 10001 ) = {%f, %f, %f, lc}; \n',[xp(i);yp(i);le-le/2]);
        elseif i+(spline-4000-1)*length(x) == 200    
            fprintf(fileID,'Point( 10200 ) = {%f, %f, %f, lc}; \n',[xp(i);yp(i);0-le/2]);
        elseif i+(spline-4000-1)*length(x) ==  300   
            fprintf(fileID,'Point( 10300 ) = {%f, %f, %f, lc}; \n',[xp(i);yp(i);le-le/2]);
        elseif i+(spline-4000-1)*length(x) ==  141   
            fprintf(fileID,'Point( 10141 ) = {%f, %f, %f, lc}; \n',[xp(i);yp(i);0-le/2]);
        elseif i+(spline-4000-1)*length(x) == 241   
            fprintf(fileID,'Point( 10241 ) = {%f, %f, %f, lc}; \n',[xp(i);yp(i);le-le/2]);
        elseif i+(spline-4000-1)*length(x) == 120    
            fprintf(fileID,'Point( 10120 ) = {%f, %f, %f, lc}; \n',[xp(i);yp(i);0-le/2]);
        elseif i+(spline-4000-1)*length(x) == 220    
            fprintf(fileID,'Point( 10220 ) = {%f, %f, %f, lc}; \n',[xp(i);yp(i);le-le/2]);
        elseif i+(spline-4000-1)*length(x) == 61    
            fprintf(fileID,'Point( 10061 ) = {%f, %f, %f, lc}; \n',[xp(i);yp(i);0-le/2]);
        elseif  i+(spline-4000-1)*length(x) == 161               
            fprintf(fileID,'Point( 10161 ) = {%f, %f, %f, lc}; \n',[xp(i);yp(i);le-le/2]);
        elseif i+(spline-4000-1)*length(x) == 40   
            fprintf(fileID,'Point( 10040 ) = {%f, %f, %f, lc}; \n',[xp(i);yp(i);0-le/2]);
        elseif i+(spline-4000-1)*length(x) == 140    
            fprintf(fileID,'Point( 10140 ) = {%f, %f, %f, lc}; \n',[xp(i);yp(i);le-le/2]);   
        end  
        
        if i+(spline-4000-1)*length(x) == 301 || i+(spline-4000-1)*length(x) == 81 || i+(spline-4000-1)*length(x) == 280 ||...
                i+(spline-4000-1)*length(x) == 60 || i+(spline-4000-1)*length(x) == 221 ||  i+(spline-4000-1)*length(x) == 1 ||...
                i+(spline-4000-1)*length(x) == 200 ||  i+(spline-4000-1)*length(x) == 300 ||  i+(spline-4000-1)*length(x) == 141 ||...
                i+(spline-4000-1)*length(x) == 241 ||  i+(spline-4000-1)*length(x) == 120 ||  i+(spline-4000-1)*length(x) == 220 ||...
                i+(spline-4000-1)*length(x) == 61 ||  i+(spline-4000-1)*length(x) == 161 ||...
                i+(spline-4000-1)*length(x) == 40 ||  i+(spline-4000-1)*length(x) == 140
                
%                 if zp(i) > 0
%                     fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i+(spline-4000-1)*length(x)+10000;xp(i);yp(i);le-le/2]);
%                 else
%                     fprintf(fileID,'Point( %i ) = {%f, %f, %f, lc}; \n',[i+(spline-4000-1)*length(x)+10000;xp(i);yp(i);0-le/2]);
%                 end
        
        end
    end
    fprintf(fileID,'\n');
    fprintf(fileID,'Spline( %d ) = {',spline);            
    %% First part
    switch spline-4000
        case 3
            fprintf(fileID,'%i, ',20);
            for i = 2:step1-1
                fprintf(fileID,'%i, ',i+(spline-4000-1)*length(x));
            end
        case 6
            fprintf(fileID,'%i, ',80);
            for i = 2:step1-1
                fprintf(fileID,'%i, ',i+(spline-4000-1)*length(x));
            end
        case 7
            fprintf(fileID,'%i, ',100);
            for i = 2:step1-1
                fprintf(fileID,'%i, ',i+(spline-4000-1)*length(x));
            end
        case 10
            fprintf(fileID,'%i, ',160);
            for i = 2:step1-1
                fprintf(fileID,'%i, ',i+(spline-4000-1)*length(x));
            end
        case 11
            fprintf(fileID,'%i, ',180);
            for i = 2:step1-1
                fprintf(fileID,'%i, ',i+(spline-4000-1)*length(x));
            end
        case 14
            fprintf(fileID,'%i, ',240);
            for i = 2:step1-1
                fprintf(fileID,'%i, ',i+(spline-4000-1)*length(x));
            end
        case 15
            fprintf(fileID,'%i, ',260);
            for i = 2:step1-1
                fprintf(fileID,'%i, ',i+(spline-4000-1)*length(x));
            end        
        otherwise
            for i = 1:step1-1
                fprintf(fileID,'%i, ',i+(spline-4000-1)*length(x));
            end
    end
    fprintf(fileID,'%i}; ',step1+(spline-4000-1)*length(x));
    fprintf(fileID,'\n');
    %% Second Part
    fprintf(fileID,'Spline( %d ) = {',spline+30);        
    for i = step1:step2-1
        fprintf(fileID,'%i, ',i+(spline-4000-1)*length(x));
    end
    fprintf(fileID,'%i};',step2+(spline-4000-1)*length(x));
    fprintf(fileID,'\n');
    %% Third Part
    fprintf(fileID,'Spline( %d ) = {',spline+60);        
    for i = step2:len-1
        fprintf(fileID,'%i, ',i+(spline-4000-1)*length(x));
    end
    if spline == 4016
        fprintf(fileID,'%i};',21);
    else
        fprintf(fileID,'%i};',length(x)+(spline-4000-1)*length(x));
    end    
    
    fprintf(fileID,'\n');
end

function[x1,y1,x2,y2] = Points_on_circle(x0,y0,a)    

    if sign(x0) == 0
        signx0 = 1;
    else
        signx0 = sign(x0);
    end
    
    if sign(y0) == 0
        signy0 = 1;
    else
        signy0 = sign(y0);
    end
    
    t = 0.01315;
    pos = atan(y0/x0);
    x1 = signx0*a*cos(pos+t);
    y1 = signx0*a*sin(pos+t);
    x2 = signy0*a*cos(pos-t);
    y2 = signy0*a*sin(pos-t);   
    
end

function[x1,y1,x2,y2] = Position_of_capacitors(x0,y0,a)    

    if sign(x0) == 0
        signx0 = 1;
    else
        signx0 = sign(x0);
    end
    
    if sign(y0) == 0
        signy0 = 1;
    else
        signy0 = sign(y0);
    end
    
    t = 0.263158;

    pos = atan(y0/x0);
    x1 = signx0*a*cos(pos+t);
    y1 = signx0*a*sin(pos+t);
    x2 = signy0*a*cos(pos-t);
    y2 = signy0*a*sin(pos-t);   
    
end