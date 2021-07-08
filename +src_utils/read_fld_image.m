

function [img xs ys zs]=read_fld_image(filepath)


file=fopen(filepath,'r');

line1=fgets(file);
line2=fgets(file);

% decode first line
pos1=strfind(line1,'[');
pos2=strfind(line1,']');

if size(pos1,2)~=3 || size(pos2,2)~=3
    error('Error reading first line.');
end

tmp=line1( pos1(1)+1:pos2(1)-1 );
[a b c ua ub uc]=decode_line(tmp);
x_min=convert_mm(a,ua);
y_min=convert_mm(b,ub);
z_min=convert_mm(c,uc);

tmp=line1( pos1(2)+1:pos2(2)-1 );
[a b c ua ub uc]=decode_line(tmp);
x_max=convert_mm(a,ua);
y_max=convert_mm(b,ub);
z_max=convert_mm(c,uc);

tmp=line1( pos1(3)+1:pos2(3)-1 );
[a b c ua ub uc]=decode_line(tmp);
dx=convert_mm(a,ua);
dy=convert_mm(b,ub);
dz=convert_mm(c,uc);

% read input file


xs = x_min : dx : x_max;
ys = y_min : dy : y_max;
zs = z_min : dz : z_max;

nx = numel(xs);
ny = numel(ys);
nz = numel(zs);

fprintf('\tx varies from  %f mm  to  %f mm (dx=%f mm  nx=%d)\n',x_min,x_max,dx,nx);
fprintf('\ty varies from  %f mm  to  %f mm (dy=%f mm  ny=%d)\n',y_min,y_max,dy,ny);
fprintf('\tz varies from  %f mm  to  %f mm (dz=%f mm  nz=%d)\n',z_min,z_max,dz,nz);

img=zeros(nx,ny,nz);

i=1; j=1; k=1;

Nleft = nx * ny * nz; 
Bsize = 10000;  % block size (=num lines in a block)

sbuf=[];
while Nleft > 0
    
    % READ THE INPUT FILE BLOCK BY BLOCK (faster than line by line...)
    
    if Nleft >= Bsize  % if enough entries are left, read an entire block
        BLOCK = reshape( fscanf( file , '%f' , Bsize*4 ) , 4 , Bsize )';
        nBLOCK = Bsize;
        Nleft = Nleft - nBLOCK;
        
    else  % otherwise, read whatever is left in the file
        
        BLOCK = reshape( fscanf( file , '%f' , Nleft*4 ) , 4 , Nleft )';
        nBLOCK = Nleft;
        Nleft = Nleft - nBLOCK;
    end
    
    
    % DECODE THE CURRENT BLOCK
    
    for ii = 1 : nBLOCK
    
        n=k + (j-1)*nz + (i-1)*nz*ny;

        if mod( n , 10000 )==0
            for kk=1:size(sbuf,2)
                fprintf('\b');
            end
            % sbuf=sprintf('[%d/%d]',n,nx*ny*nz);
            sbuf=sprintf('[%f percent]',double(n)/double(nx*ny*nz)*100);
            fprintf(sbuf);
        end

        tmp = BLOCK( ii,: );

        test_nan=findstr(tmp,'Nan');
        if ~isempty(test_nan)
            tmp=str2num( tmp(1:test_nan(1)-1) );
            tmp=[tmp 0 0 0 0 0 0];
        end

        if abs(tmp(1)*1000-xs(i))>1e-6
            error( sprintf('x dimension error [read=%e  expected=%e].',tmp(1)*1000,xs(i)) );
        end

        if abs(tmp(2)*1000-ys(j))>1e-6
            error( sprintf('y dimension error [read=%e  expected=%e].',tmp(2)*1000,ys(j)) );
        end

        if abs(tmp(3)*1000-zs(k))>1e-6
            error( sprintf('z dimension error [read=%e  expected=%e].',tmp(3)*1000,zs(k)) );
        end

        img(i,j,k)=tmp(4);

        k=k+1;

        if k>nz
            k=1;
            j=j+1;
        end

        if j>ny
            j=1;
            i=i+1;
        end

    %     if i>nx
    %         error('x index out of bound.');
    %     end

    
    end  % end for loop over the lines in the current block   
    
    
end

fprintf('\n');

fclose(file);




function [a b c ua ub uc]=decode_line(line)

pos=strfind(line,' ');

a=str2num( line( 1:pos(1)-3 ) );
ua=line( pos(1)-2:pos(1)-1 );

b=str2num( line( pos(1)+1:pos(2)-3 ) );
ub=line( pos(2)-2:pos(2)-1 );

c=str2num( line( pos(2)+1:end-2 ) );
uc=line( end-1:end );




function b=convert_mm(a,u)

if strcmpi('mm',u)==1
    fac=1;
elseif strcmpi('cm',u)==1
    fac=10;
else
    error('Unknown input unit.');
end

b=a*fac;

