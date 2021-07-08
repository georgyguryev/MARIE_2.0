function COIL_TwoPortCoil(y_shift, opt1, opt2)
%clear all;
clc;

refine = false;
second_refinement = false;
third_refinement = false;

fw = 0;

name = '2port_asym_dg';
filename = strcat('./data/coils/', name);
geofilename = sprintf('%s.geo', filename);
scoilname = sprintf('%s.smm', filename);
meshname = sprintf('%s.msh', filename);
portspec_file = sprintf('%s.txt', filename);


%% set geometry

if(nargin>2)
    
    order_of_ref = opt1;
    length = opt2.L;
    wire_width = opt2.WW;
    width = opt2.W;
    shift = opt2.shift;
    res = wire_width * opt2.mesh_size;
    
elseif (nargin > 1)
    
    order_of_ref = opt1;
    gap = 0;
    length = 0.35;                   % length of a sub-coil
    wire_width = 10e-3;             % wire width of sub-coil
    width = 12*wire_width + gap;     % width of a sub-coil
    shift = 5e-3;
    res = wire_width;
    
else
    order_of_ref = 0;
    length = 0.350;                   % length of a sub-coil
    wire_width = 10e-3;             % wire width of sub-coil
    width = 0.12;     % width of a sub-coil
    shift = 5e-3;
    gap = 0;
    res = wire_width;
end

%% 

rho = 0.0;

Points = [];
Lines = [];
Patch = [];
numP = 2; 


N_sc = 2;  %number of sub-coils in the coil

pcount = 0;
lcount = 0;

L = length;
w = width;
ww = wire_width;



pshift = pcount;

% point 1
x = - w/2; y = y_shift; z = -L/2;
Points = [Points; x, y, z, res];
pcount = pcount + 1;

% point 2
x = -w/2; y = y_shift; z = L/2;
Points = [Points; x, y, z, res];
pcount = pcount + 1;

%point 3

x = -gap/2 + shift; y = y_shift; z = L/2;
Points = [Points; x, y, z, res];
pcount = pcount + 1;

%point 4

x = -gap/2 + shift; y = y_shift; z = L/2 - ww;
Points = [Points; x, y, z, res];
pcount = pcount + 1;

% point 5

x = -w/2 + ww; y = y_shift; z = L/2 - ww;
Points = [Points; x, y, z, res];
pcount = pcount + 1;

% point 6

x = -w/2 + ww; y = y_shift; z = -L/2 + ww;
Points = [Points; x, y, z, res];
pcount = pcount + 1;

% point 7

x = w/2 - ww; y = y_shift; z =  -L/2 + ww;
Points = [Points; x, y, z, res];
pcount = pcount + 1;

% point 8

x = w/2 - ww; y = y_shift; z = - shift;
Points = [Points; x, y, z, res];
pcount = pcount + 1;

% point 9

x = w/2; y = y_shift; z = -shift;
Points = [Points; x, y, z, res];
pcount = pcount + 1;

% point 10

x = w/2; y = y_shift; z = -L/2;
Points = [Points; x, y, z, res];
pcount = pcount + 1;


% point 11
x = w/2 -ww; y = y_shift; z = L/2- ww;
Points = [Points; x, y, z, res];
pcount = pcount + 1;

% point 12

x = w/2; y = y_shift; z = L/2;
Points = [Points; x, y, z, res];
pcount = pcount + 1;



idxp1 = pshift + 1;
idxp2 = pshift + 2;
idxp3 = pshift + 3;
idxp4 = pshift + 4;
idxp5 = pshift + 5;
idxp6 = pshift + 6;
idxp7 = pshift + 7;
idxp8 = pshift + 8;
idxp9 = pshift + 9;
idxp10 = pshift + 10;
idxp11 = pshift + 11;
idxp12 = pshift + 12;



%--------------------------------------------
% Line
%--------------------------------------------
% define lines

% Line 1
Lines = [Lines; idxp1, idxp2, -100];
lcount = lcount + 1;
idxl1 = lcount;

%Line 2
Lines = [Lines; idxp2, idxp3, 1];
lcount = lcount + 1;
idxl2 = lcount;

Lines = [Lines; idxp3, idxp4, -100];
lcount = lcount + 1;
idxl3 = lcount;

Lines = [Lines; idxp4, idxp5, 2];
lcount = lcount + 1;
idxl4 = lcount;

Lines = [Lines; idxp5, idxp6, -100];
lcount = lcount + 1;
idxl5 = lcount;

Lines = [Lines; idxp6, idxp7, 3];
lcount = lcount + 1;
idxl6 = lcount;

Lines = [Lines; idxp7, idxp8, -100];
lcount = lcount + 1;
idxl7 = lcount;

Lines = [Lines; idxp8, idxp9, 4];
lcount = lcount + 1;
idxl8 = lcount;

Lines = [Lines; idxp9, idxp10, -100];
lcount = lcount + 1;
idxl9 = lcount;

Lines = [Lines; idxp10, idxp1, 5];
lcount = lcount + 1;
idxl10 = lcount;

% Line 11
Lines = [Lines; idxp8, idxp11, -100];
lcount = lcount + 1;
idxl11 = lcount;

%Line 12
Lines = [Lines; idxp11, idxp4, 1];
lcount = lcount + 1;
idxl12 = lcount;


idxl13 = -3;

Lines = [Lines; idxp3, idxp12, 2];
lcount = lcount + 1;
idxl14 = lcount;

Lines = [Lines; idxp12, idxp9, -100];
lcount = lcount + 1;
idxl15 = lcount;

idxl16 = -8;

if fw
    
    Lines = [Lines; idxp3, idxp14, 1];
    lcount = lcount + 1;
    idxl17 = lcount;

    Lines = [Lines; idxp13, idxp4, 2];
    lcount = lcount + 1;
    idxl18 = lcount;
    
    Lines = [Lines; idxp8, idxp11, -100];
    lcount = lcount + 1;
    idxl19 = lcount;

    Lines = [Lines; idxp16, idxp9, -100];
    lcount = lcount + 1;
    idxl20 = lcount;
    
end


% Lines = [Lines; idxp8, idxp11, 3];
% lcount = lcount + 1;
% idxl17 = lcount;
% 
% Lines = [Lines; idxp16, idxp9, -100];
% lcount = lcount + 1;
% idxl18 = lcount;

Patch = cell(2,1);

Patch{1, :} = [idxl1, idxl2, idxl3, idxl4, idxl5, idxl6, idxl7, idxl8, idxl9, idxl10 ]; 
Patch{2, :} = [idxl11, idxl12, idxl13, idxl14, idxl15, idxl16];

if fw
    Patch{3, :} = [-idxl3, idxl17, -idxl13, idxl18];
    Patch{4, :} = [idxl19, -idxl16, idxl20, -idxl8];
end


Ports = zeros(numP,2);

% for i = 1:numP
%     Ports(i,1) = 3 + (i-1)*4;
%     Ports(i,2) = 19 + (i-1)*4;    
% end;

Ports = [3;
         8];



%% prepare 
% number of coils 
Cnum = 1;

% preparing individual indexes for geometry entities

Nstart = Cnum * 100;         % indexes for physical lines (ports) and physical surfaces(i.e. coils)
Pstart = Nstart + 100;       % indexes for elementary points
Lstart = Pstart + 100;       % indexes for elementary lines
Ostart = Lstart + 100;       % indexes for surfaces
Sstart = Ostart + 100;

%% Open file for output

fid = fopen(geofilename, 'w');

fprintf(fid, '\n //--------------------------------------------//');
fprintf(fid, '\n //--------------------------------------------//');
fprintf(fid, '\n //------Simple coil with %d sub-coils----------//', N_sc);
fprintf(fid, '\n //--------Resolution = %f---------------//', res);
fprintf(fid, '\n //--------------------------------------------//');

fprintf(fid,'\n\n');

%% write geometry

% write points to a file
for i = 1:pcount
    fprintf(fid, 'P%d = %d;\n', Pstart + i, Pstart +i);
    fprintf(fid, 'Point(P%d) = {%2.16g, %2.16g, %2.16g, %2.16g}; \n', Pstart+i, Points(i,1), Points(i,2), Points(i,3), Points(i,4));
end

fprintf(fid, '\n \n');

%write lines

for i = 1:lcount
    fprintf(fid, 'L%d = %d; \n', Lstart + i, Lstart +i);
    fprintf(fid, 'Line(L%d) = {P%d, P%d}; \n', Lstart + i, Pstart + Lines(i,1), Pstart + Lines(i,2));
end

fprintf(fid, '\n \n');

% write Line loops

fprintf(fid, '\n \n');
% write Line loops
for i = 1:size(Patch,1)
    fprintf(fid, 'O%d = %d;\n', Ostart + i, Ostart + i);
    fprintf(fid, 'Line Loop(O%d) = {', Ostart + i); 
    
    for  j = 1:size(Patch{i},2)-1 
        if ((Patch{i}(j)) < 0)
            fprintf(fid, '-L%d,', Lstart + abs(Patch{i}(j)));
        else
            fprintf(fid, 'L%d,', Lstart + abs(Patch{i}(j)));
        end
    end
    
    if (Patch{i}(size(Patch{i},2)) < 0)
        fprintf(fid, '-L%d}; \n', Lstart + abs(Patch{i}(size(Patch{i},2))));        
    else
        fprintf(fid, 'L%d}; \n', Lstart + abs(Patch{i}(size(Patch{i},2))));
    end
end
fprintf(fid, '\n \n');



%write Ports
for i = 1:size(Ports,1)
    fprintf(fid, 'Port%d = %d;\n', Nstart + 2*i-1, Nstart + 2*i -1);
    fprintf(fid, 'Physical Line (Port%d) = {L%d}; \n', Nstart + 2*i -1, Lstart + Ports(i,1));
end
fprintf(fid, '\n \n');

%write Ruled surfaces
for i = 1:size(Patch,1)
    fprintf(fid, 'S%d = %d;\n', Sstart + i, Sstart + i);
    fprintf(fid, 'Plane Surface(S%d) = {O%d}; \n', Sstart + i, Ostart + i);
end

fprintf(fid, '\n \n');

% write fw port physical surface

if fw
% write fw port physical surface

count = 0;

for i = size(Patch,1) - size(Ports,1) + 1 : size(Patch,1)
    
    fprintf(fid, 'Coil%d = %d;\n', Cnum + count, Cnum + count);
    fprintf(fid, 'Physical Surface(Coil%d) = {', Cnum + count);
    fprintf(fid, 'S%d}; \n ',Sstart +i);
    count = count + 1;
    
end

    %write physical surfaces
    fprintf(fid, 'Coil%d = %d;\n', Cnum + count, Cnum + count);
    fprintf(fid, 'Physical Surface(Coil%d) = {', Cnum + count);

    for i = 1:size(Patch,1)- (size(Ports,1) + 1) 
        fprintf(fid, 'S%d, ',Sstart +i);
    end

    fprintf(fid, 'S%d}; \n', Sstart + size(Patch,1) - size(Ports,1));

else
    count = 0;
    %write physical surfaces
    fprintf(fid, 'Coil%d = %d;\n', Cnum + count, Cnum + count);
    fprintf(fid, 'Physical Surface(Coil%d) = {', Cnum + count);

    for i = 1:size(Patch,1)-1 
        fprintf(fid, 'S%d, ',Sstart +i);
    end

    fprintf(fid, 'S%d}; \n', Sstart + size(Patch,1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %write physical surfaces
% fprintf(fid, 'Coil%d = %d;\n', Cnum, Nstart);
% fprintf(fid, 'Physical Surface(Coil%d) = {', Cnum);
% 
% for i = 1:size(Patch,1)-1
%     fprintf(fid, 'S%d, ',Sstart +i);
% end;
% 
% fprintf(fid, 'S%d}; \n', Sstart + size(Patch,1));


fclose(fid);

if ispc
    
    command = sprintf('.\\src_generate\\src_scoil\\gmsh %s', geofilename);
    
elseif ismac
    
    command1 = sprintf('./src_generate/src_scoil/gmsh.app/Contents/MacOS/gmsh  -2  %s ', geofilename);
    command2 = sprintf('./src_generate/src_scoil/gmsh.app/Contents/MacOS/gmsh  -refine  %s ', meshname);

else
    
%     command1 = sprintf('./data/coils/GMSH/bin/gmsh %s -1 -2 -o %s', geofilename, meshname);
    command2 = sprintf('./data/coils/GMSH/bin/gmsh  -refine  %s ', meshname);
    command1 = sprintf('/usr/bin/gmsh %s -1 -2 -o %s', geofilename, meshname);

end

system('pwd');
system(command1);

for i = 1:order_of_ref
    system(command2);
end

%% write scoilfile

clear fid;


fid = fopen(scoilname, 'w');
fprintf(fid, strcat(name, '\n'));
fprintf(fid,'%f  %f \n', rho, res);
fprintf(fid, strcat(name,'.msh'));
fclose(fid);

%% create a template with port specification

N_ports = size(Ports,1);

% create file port specs file
fid = fopen(portspec_file, 'w');

% generate port specs file
src_geo.generate_PortSpecs(fid, N_ports);  

% close file
fclose(fid);


%% create a template for 

%clear all;
clc;

fprintf('Coil was successfully generated!\n');
 
