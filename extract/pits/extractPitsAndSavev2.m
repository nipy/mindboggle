function a = extractPitsAndSavev2(datapath)
%
% for lambda = 0:0.2:1
%
%     load([datapath '/miccaiData/inputLH.mat']);
%     Q = compOptimalHMMF(L,D,pitInds,lambda);
%     save([datapath '/miccaiData/pitResultsLH' num2str(1000*lambda)],'Q','D','pitInds');
%
%     load([datapath '/miccaiData/inputRH.mat']);
%     Q = compOptimalHMMF(L,D,pitInds,lambda);
%     save([datapath '/miccaiData/pitResultsRH' num2str(1000*lambda)],'Q','D','pitInds');
%
% end

for lambda = 0:1:5
    
    
    
    load([datapath '/miccaiData/inputv5LH.mat']);
    Q = compOptimalHMMF(L,D,pitInds,lambda);
    save([datapath '/miccaiData/pitResultsv5LH' num2str(1000*lambda)],'Q','D','pitInds');
    
    load([datapath '/miccaiData/inputv5RH.mat']);
    Q = compOptimalHMMF(L,D,pitInds,lambda);
    save([datapath '/miccaiData/pitResultsv5RH' num2str(1000*lambda)],'Q','D','pitInds');
    
    
end

a = 0;
