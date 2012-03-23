function [nlh nrh minDistlh minDistrh] = analyzeResults(datapath,lambda)

load([datapath '/miccaiData/pitResultsv6LH' num2str(1000*lambda) '.mat']);
nlh = sum(Q>0.5);
minDistlh = findMinDistOfExtractedPits(D,Q,pitInds);

load([datapath '/miccaiData/pitResultsv6RH' num2str(1000*lambda) '.mat']);
nrh = sum(Q>0.5);
minDistrh = findMinDistOfExtractedPits(D,Q,pitInds);
