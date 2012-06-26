function D = compDirectionalDist(P1,P2,U)
dirV = dot(P1-P2,U);
D= sqrt(sum(dirV.^2));
