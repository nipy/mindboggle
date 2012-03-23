function newLabels = modifyLabels(allLabels)

newLabels = allLabels;

newLabels(newLabels == 10) = 2;
newLabels(newLabels == 23) = 2;
newLabels(newLabels == 26) = 2;

newLabels(newLabels == 27) = 3;

newLabels(newLabels == 19) = 18;
newLabels(newLabels == 20) = 18;
