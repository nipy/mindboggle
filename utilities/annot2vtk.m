%% FreeSurfer annotation to VTK
% Forrest Sheng Bao 2011-07-30 
% v. 1.1 Last update 2011-08-02 
% This software is free software, as defined in GNU GPL v.3.0
% There is ABSOLUTELY NO WARRANTY 

%% Dependencies
% FreeSurfer's MATLAB library: read_surf, 
% 

%% Variables
% AnnotFile - string, input FreeSurfer annotation file name
% vertices - 1-D double array, vertex IDs
% label - 1-D double array, vertex labels 
% colortable - a struct, content not studied yet. 
% VTKFile - string, output VTK filename
% Surf - 3-column double array, 3-D coordinates of vertexes on surface
% SurfFile - string, input FreeSurfer surface file name

function annot2vtk(AnnotFile, SurfFile)
%% Load FreeSurfer files
[Vertexes,Label,Colortable]=read_annotation(AnnotFile);
Label = int32(Label);
VTKFile = [AnnotFile(1:length(AnnotFile)-6) SurfFile(strfind(SurfFile, '.'):length(SurfFile)) '.vtk'];
[Surf, Faces] = read_surf(SurfFile);
NumFace = length(Faces);

%% write VTK header 
NumVertex = length(Vertexes);
Fid = fopen(VTKFile, 'w');
fprintf(Fid, '# vtk DataFile Version 2.0\nBy annot2VTK of Forrest Sheng Bao http://fsbao.net\nASCII\n');
fprintf(Fid, 'DATASET POLYDATA\n');

%% write vertex coordinates
fprintf(Fid, 'POINTS %d float\n', NumVertex);
fclose(Fid);
dlmwrite(VTKFile, Surf, '-append', 'delimiter', ' ');

%% write vertex IDs as VERTICES in one line
% This block is currently disabled because we wanna see faces instead of
% vertexes.
% Fid = fopen(VTKFile, 'a');
% fprintf(Fid, 'VERTICES %d %d\n%d ', NumVertex, NumVertex+1, NumVertex+1);
% for i = 1: NumVertex
%     fprintf(Fid, '%d ', i);
% end
% fprintf(Fid, '\n');
% fclose(Fid);

%% write faces into VTK file 
Fid = fopen(VTKFile, 'a');
fprintf(Fid, 'POLYGONS %d %d\n', NumFace, NumFace*4);
fclose(Fid);
Faces = horzcat(3*ones(length(Faces),1),Faces);
dlmwrite(VTKFile, Faces, '-append', 'delimiter', ' ', 'precision', '%d');

%% write vertex labels as an LUT
Fid = fopen(VTKFile, 'a');
fprintf(Fid, 'POINT_DATA %d\nSCALARS annot float\n', NumVertex);
fprintf(Fid, 'LOOKUP_TABLE annot\n');

% for i = 1:NumVertex
%     fprintf(Fid, '%d\n', i);
% end


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
