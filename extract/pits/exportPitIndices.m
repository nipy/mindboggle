function a = exportPitIndices(datapath,subjectName,lambda)

disp('v6');

load([datapath '/miccaiData/pitResultsv6LH' num2str(1000*lambda) '.mat']);
indsLH = find(Q>.5);

load([datapath '/miccaiData/pitResultsv6RH' num2str(1000*lambda) '.mat']);
indsRH = find(Q>.5);

% fidLH = fopen([datapath '/miccaiData/pitIdxLH_' subjectName '.txt'], 'w+t');
% fidRH = fopen([datapath '/miccaiData/pitIdxRH_' subjectName '.txt'], 'w+t');

fidLH = fopen(['/Users/yrjo/yrjo_work/mindboggle2/miccai2012/data/export/pitIdxLH_' subjectName '.txt'], 'w+t');
fidRH = fopen(['/Users/yrjo/yrjo_work/mindboggle2/miccai2012/data/export/pitIdxRH_' subjectName '.txt'], 'w+t');

for i = 1:size(indsLH,1)
    fprintf(fidLH,num2str(indsLH(i,1)));
    fprintf(fidLH,'\n');    
end
for i = 1:size(indsRH,1)
    fprintf(fidRH,num2str(indsRH(i,1)));
    fprintf(fidRH,'\n');    
end

fclose(fidLH);
fclose(fidRH);

a = 0;
