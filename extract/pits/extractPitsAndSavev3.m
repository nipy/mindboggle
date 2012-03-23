function a = extractPitsAndSavev3(datapath,subjectID)
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

disp('v3');
for lambda = 0:1:5
        
    load([datapath '/miccaiData/inputv6LH.mat']);  
    load(['/Users/yrjo/Dropbox/oma/miccai2012/transfer/' subjectID '/inputdistv6LH.mat']);            
    Q = compOptimalHMMF(L,D,pitInds,lambda);
    save([datapath '/miccaiData/pitResultsv6LH' num2str(1000*lambda)],'Q','D','pitInds');
    
    load([datapath '/miccaiData/inputv6RH.mat']);
    load(['/Users/yrjo/Dropbox/oma/miccai2012/transfer/' subjectID '/inputdistv6RH.mat']);
    Q = compOptimalHMMF(L,D,pitInds,lambda);
    save([datapath '/miccaiData/pitResultsv6RH' num2str(1000*lambda)],'Q','D','pitInds');
        
end

a = 0;
