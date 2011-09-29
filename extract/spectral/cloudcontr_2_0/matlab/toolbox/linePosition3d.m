function d = linePosition3d(point, line)
%LINEPOSITION3D return position of a 3D point on a 3D line
%
%   L = linePosition3d(POINT, LINE)
%   compute position of point POINT on the line LINE, relative to origin
%   point and direction vector of the line.
%   LINE has the form [x0 y0 z0 dx dy dy],
%   POINT has the form [x y z], and is assumed to belong to line.
%   If POINT does not belong to LINE, the position of its orthogonal
%   projection is computed instead.
%
%   L = linePosition3d(POINT, LINES)
%   if LINES is an array of NL lines, return NL positions, corresponding to
%   each line.
%
%   L = linePosition3d(POINTS, LINE)
%   if POINTS is an array of NP points, return NP positions, corresponding
%   to each point.
%
%   L = linePosition3d(POINTS, LINES)
%   if POINTS is an array of NP points and LINES is an array of NL lines,
%   return an array of [NP NL] position, corresponding to each couple
%   point-line.
%
%   See also:
%   lines3d, points3d, createLine3d
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 17/02/2005.
%

%   HISTORY :
%   05/01/2007 update doc

if size(point, 1)== 1 && size(line, 1)>1
    point = repmat(point, [size(line, 1) 1]);
end

if size(point, 1)>1 && size(line, 1)==1
    line = repmat(line, [size(point, 1) 1]);
end

% vector from line origin to point
dp = point - line(:,1:3);

% direction vector of the line
dl = line(:,4:6);

% compute position using dot product normalized with norm of line vector.
d = dot(dp, dl, 2)./dot(dl, dl, 2);
