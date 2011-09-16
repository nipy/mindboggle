%% FreeSurfer annotation to VTK
%
% This software is licensed under MIT license. 
% 
% Copyright (C) 2011 by Forrest Sheng Bao <http://fsbao.net>
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
%
% v. 1.2 Last update 2011-09-02 by Forrest Sheng Bao

%% Dependencies
% FreeSurfer's MATLAB library: read_surf
% 

%% Variables
% AnnotFile - string, input FreeSurfer annotation file name
% vertices - 1-D double array, vertex IDs
% label - 1-D double array, vertex labels 
% colortable - a struct, content not studied yet. 
% VTKFile - string, output VTK filename
% Surf - 3-column double array, 3-D coordinates of vertexes on surface
% SurfFile - string, input FreeSurfer surface file name
% MapTo - string, location to map annotations (Default: 'Face'; otherwise, 'Vertex')
%
%% Usage
% On MATLAB prompt >> annot2vtk(AnnotFile, SurfFile)
% Example: annot2vtk('lh.aparcNMMjt.annot','lh.inflated')
% The output file is in the same directory with AnnotFile
% For complete information, please check wiki:
% http://code.google.com/p/mindboggle-utils/wiki/annot2vtk

function annot2vtk(AnnotFile, SurfFile, MapTo)
%% Processing input arguments
if nargin < 3
    MapTo = 'Face';
elseif ~strcmp(MapTo, 'Face') && ~strcmp(MapTo, 'Vertex')
    sprintf('Wrong value for MapTo. Check documentation for usage.\n')
    return 
end

%% Load FreeSurfer files
[Vertexes,Label,Colortable]=read_annotation(AnnotFile);
Label = int32(Label);
VTKFile = [AnnotFile(1:length(AnnotFile)-6) SurfFile(strfind(SurfFile, '.'):length(SurfFile)) '.vtk'];
[Surf, Faces] = read_surf(SurfFile);
NumFace = length(Faces);

%% write VTK header 
NumVertex = length(Vertexes);
Fid = fopen(VTKFile, 'w');
fprintf(Fid, '# vtk DataFile Version 2.0\nBy annot2vtk of Forrest Sheng Bao \nASCII\n');
fprintf(Fid, 'DATASET POLYDATA\n');

%% write vertex coordinates
fprintf(Fid, 'POINTS %d float\n', NumVertex);
fclose(Fid);
dlmwrite(VTKFile, Surf, '-append', 'delimiter', ' ');

%% write faces into VTK file 
if strcmp(MapTo, 'Face')
    Fid = fopen(VTKFile, 'a');
    fprintf(Fid, 'POLYGONS %d %d\n', NumFace, NumFace*4);
    fclose(Fid);
    Faces = horzcat(3*ones(length(Faces),1),Faces);
    dlmwrite(VTKFile, Faces, '-append', 'delimiter', ' ', 'precision', '%d');
elseif strcmp(MapTo, 'Vertex')    
%% write vertex IDs as VERTICES in one line
    Fid = fopen(VTKFile, 'a');
    fprintf(Fid, 'VERTICES %d %d\n%d ', NumVertex, NumVertex+1, NumVertex+1);
    for i = 1: NumVertex
         fprintf(Fid, '%d ', i);
    end
    fprintf(Fid, '\n');
    fclose(Fid);    
else
    sprintf('Unrecognized mapping destination. Say Face or Vertex. Check documentation for help.\n')
    return
end

%% write vertex labels as an LUT
Fid = fopen(VTKFile, 'a');
fprintf(Fid, 'POINT_DATA %d\nSCALARS annot float\n', NumVertex);
fprintf(Fid, 'LOOKUP_TABLE annot\n');

fclose(Fid);

dlmwrite(VTKFile, Label, '-append', 'delimiter', ' ', 'precision', '%d');

%% write an LUT for coloring purpose
% Fid = fopen(VTKFile, 'a');
% fprintf(Fid, 'LOOKUP_TABLE annot %d\n', NumVertex);
% fclose(Fid);
% 
% % convert FreeSurfer color table contents into RGBA values
% LUT = Colortable.table(:,1:3)/255;
% Color = ones(NumVertex, 4);
% for i = 1: NumVertex
% %    disp(i);
%     [Row, Col] = find(Colortable.table(:,5) == Label(i));
%     if isempty(Row) % this might be a bug of FreeSurfer
%         Color(i, 1:3) = [0 0 0];
%     else
%         Color(i, 1:3) = LUT(Row, 1:3);
%     end
% end
% 
% dlmwrite(VTKFile, Color, '-append', 'delimiter', ' ', 'precision', '%1.5f');