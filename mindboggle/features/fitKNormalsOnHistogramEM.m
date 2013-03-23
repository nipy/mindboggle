function [clMeans clSigmas W] = fitKNormalsOnHistogramEM(data,x)

Kclasses = 3;

xrange = max(x) - min(x);

Ndata= size(data,1);

clMeans = zeros(Kclasses,1);
clSigmas = zeros(Kclasses,1);

W = zeros(Kclasses,1);

for i = 1:Kclasses
    clMeans(i) = max(x) - xrange/2 - (0.2*(i - Kclasses/2)*xrange);
    clSigmas(i) = 0.2;
end

iter = 0;

probs = zeros(Ndata,Kclasses);
weights = zeros(Ndata,Kclasses);

while (iter < 25)
    iter = iter + 1;
    
    for i = 1:Kclasses                
        probs(:,i) = (1/(clSigmas(i)*sqrt(2*pi))) .* exp((-1/(2*(clSigmas(i)^2))) .* ((data-clMeans(i)).^2));
    end
    
    for i = 1:Kclasses
        weights(:,i) = probs(:,i) ./ sum(probs,2);
    end
    for i = 1:Kclasses
        
        clSigmas(i) = sqrt(sum(weights(:,i).*((data - clMeans(i)).^2)) / sum(weights(:,i)));
        clMeans(i) = sum(weights(:,i) .* data) / sum(weights(:,i));
    end        
    disp([clMeans clSigmas]);        
end

for i = 1:Kclasses
    W(i) = sum(weights(:,i))/sum(weights(:));
end



