function [verts, faces, normals] = read_off(filename)

% read_off - read data from OFF file.
%
%   [verts, faces, normals] = read_off(filename);
%
%   'verts' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'faces' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%   'normals' is a 'nb.normal * 3' array specifying the normal of each vertex.
%   Copyright (c) 2003 Gabriel Peyr?
%   Changed by (2009) JJCAO

fid = fopen(filename,'r');
if( fid==-1 )
    warning('Can''t open the file.');
    verts = [];faces=[];normals=[];
    return;
end

fgetl(fid); 
str = fgetl(fid);
tmp = sscanf(str,'%d %d');
nverts = tmp(1);
nfaces = tmp(2);
str = fgetl(fid);
tmp = sscanf(str,'%f %f %f');

frewind(fid);
fgetl(fid); fgetl(fid);
if length(tmp) < 4 % only x y z
    [A,cnt] = fscanf(fid,'%f %f %f', 3*nverts);
    if cnt~=3*nverts
        warning('Problem in reading vertices.');
    end
    A = reshape(A, 3, cnt/3);
    verts = A';
    normals = [];
else % x y z nx ny nz
    [A,cnt] = fscanf(fid,'%f %f %f %f %f %f', 6*nverts);
    if cnt~=6*nverts
        warning('Problem in reading vertices.');
    end
    
    A = reshape(A, 6, cnt/6);
    A = A';
    verts = A(:,1:3);
    normals = A(:,4:6);
end

if nfaces > 0 % read Face 1  1088 480 1022
    [A,cnt] = fscanf(fid,'%d %d %d %d\n', 4*nfaces);
    if cnt~=4*nfaces
        warning('Problem in reading faces.');
    end
    A = reshape(A, 4, cnt/4);
    faces = A(2:4,:)+1;
    faces = faces';
else
    faces = [];
end

fclose(fid);
return;