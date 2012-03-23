function p = compProbability(q,v,eD,QP,cI)

wd = 17;
%%%%%%%%%%%
p = log(q*v + (1-q)*(1-v));
p = p + wd * ((q^2) + sum(eD(cI,:)' .* (-1.* (q.*QP)))); 



