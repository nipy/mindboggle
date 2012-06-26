function d = compEucDist(A,B)
d = 0;

d = A-B;
d = d.^2;
d = sqrt(sum(d));
