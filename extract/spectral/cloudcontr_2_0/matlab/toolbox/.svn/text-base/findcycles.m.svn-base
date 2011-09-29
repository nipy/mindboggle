function [cycles]=findcycles(G)
% for directed graph
% 
% example:
% G = sparse([1 2 3 4 5 5 5 6 7 8 8 9],[2 4 1 5 3 6 7 9 8 1 10 2],true,10,10);
% view(biograph(G));
% findcycles(G);

numNodes = size(G,1); 
cycles = cell(0);
for n = 1:numNodes
   [D,P]=graphtraverse(G,n);
   for d = D
       if G(d,n)
           cycles{end+1} = graphpred2path(P,d);
       end
   end
   G(n,:)=0; 
end
