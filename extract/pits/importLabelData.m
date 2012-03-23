function [inds labels distances] = importLabelData(subjectID,lh)

% pathMain = '/Users/yrjo/yrjo_work/mindboggle2/miccai2012/data/KKI/';
% datapath = ([pathMain subjectID]);

if(lh)
    sideString = 'LH';
else
    sideString = 'RH';
end
data = load(['/Users/yrjo/yrjo_work/mindboggle2/miccai2012/data/seg/pitIdx' sideString '_' subjectID '.txt.seg.tsv']);
%%%%%%%%%%%%%%%%%
data(:,1) = data(:,1) + 1;
data(:,9:11) = data(:,9:11) + 1;
%%%%%%%%%%%%%%%%%
inds = data(:,1);
indsLabel = data(:,9:11);
%%%%%%%%%%%%%%%%%
labels = data(:,2:5);
hops = data(:,6:8);

%%%%%%%%%%%%%%%%%

load(['/Users/yrjo/Dropbox/oma/miccai2012/transfer/' subjectID '/pitDistv6' sideString '.mat']);
%load([datapath '/miccaiData/inputv5' sideString '.mat']);
%%%%%%%%%%%%%%%%%

distances = zeros(size(data,1),3);
%%%%%%%%%%%%%%%%%%%


for i = 1:size(distances,1)
    dColumn = find(pitInds == inds(i));
    if (D(inds(i),dColumn) ~= 0)
       disp('ERROR!'); 
    end       
    distances(i,:) = D(indsLabel(i,:)',dColumn)';
    distances(i,:) = distances(i,:) - distances(i,:)./hops(i,:);        
end

