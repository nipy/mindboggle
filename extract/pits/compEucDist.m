function d = compEucDist(A,B)
d = A-B;
d = d.^2;
d = sqrt(sum(d));
