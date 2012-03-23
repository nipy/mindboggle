function gP = compGlobalProb(Q,L,D,wd)
%wd = 20;

m = length(Q);

gP = log(Q.*L + (1-Q).*(1-L));

gP = gP + wd.*(((D * ((-1).*Q)).*Q) + Q.^2);

% for i = 1:m        
%     gP(i) = gP(i) + wd * ((Q(i).^2) + sum(D(i,:)' .* (-1.* (Q(i).*Q)))); 
%     
%     %gP = gP + compProbability(Q(i),L(i),D,Q,i);
% end

gP = sum(gP);
