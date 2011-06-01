function  [lambda,ecc,radius,diameter] = charpath(D)
%CHARPATH       Characteristic path length and related statistics
%
%   lambda = charpath(D);
%   [lambda,ecc,radius,diameter] = charpath(D);
%
%   The characteristic path length is the average shortest path length in 
%   the network.
%
%   Input:      D,          distance matrix
%
%   Outputs:    lambda,     characteristic path length
%               ecc,        eccentricity (for each vertex)
%               radius,     radius of graph
%               diameter,   diameter of graph
%
%   Note: Characteristic path length is calculated as the global mean of 
%   the distance matrix D, excludings any 'Infs' but including distances on
%   the main diagonal.
%
%
%   Olaf Sporns, Indiana University, 2002/2007/2008

% Mean of finite entries of D(G)
lambda = sum(sum(D(D~=Inf)))/length(nonzeros(D~=Inf));

% Eccentricity for each vertex (note: ignore 'Inf') 
ecc = max(D.*(D~=Inf),[],2);

% Radius of graph
radius = min(ecc);  % but what about zeros?

% Diameter of graph
diameter = max(ecc);
