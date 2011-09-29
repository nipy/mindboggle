function h = plot_embedded_graph( vertex,edge,varargin)
i = edge(:,1);
j = edge(:,2);

% Create a long, NaN-separated list of line segments,
% rather than individual segments.

X = [ vertex(i,1) vertex(j,1) repmat(NaN,size(i))]';
Y = [ vertex(i,2) vertex(j,2) repmat(NaN,size(i))]';
Z = [ vertex(i,3) vertex(j,3) repmat(NaN,size(i))]';
X = X(:);
Y = Y(:);
Z = Z(:);

h = plot3(X, Y, Z, varargin{:} );
