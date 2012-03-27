function [allLabels allDistances subjectList] = importLabelDataBatch()

allLabels = -1*ones(2000,4);
allDistances = -1*ones(2000,3);
subjectList = -1*ones(2000,1);

counter = 0;
subjectCounter = 1;

subjectID = 'KKI2009-11';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

subjectID = 'KKI2009-14';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

subjectID = 'KKI2009-15';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
disp(inds(29));
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

subjectID = 'KKI2009-16';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
disp(inds(14));
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

% subjectID = 'KKI2009-17';
% [inds labels distance] = importLabelData(subjectID,lh);

subjectID = 'KKI2009-18';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

subjectID = 'KKI2009-19';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

subjectID = 'KKI2009-20';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

subjectID = 'KKI2009-21';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

subjectID = 'KKI2009-22';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

subjectID = 'KKI2009-24';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

subjectID = 'KKI2009-25';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

subjectID = 'KKI2009-27';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

subjectID = 'KKI2009-29';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

subjectID = 'KKI2009-31';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

subjectID = 'KKI2009-33';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

subjectID = 'KKI2009-34';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

% subjectID = 'KKI2009-36';
% [inds labels distance] = importLabelData(subjectID,lh);

subjectID = 'KKI2009-37';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

subjectID = 'KKI2009-40';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

subjectID = 'KKI2009-41';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

subjectID = 'KKI2009-42';
[inds labels distance] = importLabelData(subjectID,1);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;

[inds labels distance] = importLabelData(subjectID,0);
len = length(inds);
allLabels(counter+1:counter+len,:) = labels;
allDistances(counter+1:counter+len,:) = distance;
subjectList(counter+1:counter+len,:) = subjectCounter;
subjectCounter = subjectCounter +1;
counter = counter + len;


allLabels = allLabels(1:counter,:);
allDistances = allDistances(1:counter,:);

disp(subjectCounter);