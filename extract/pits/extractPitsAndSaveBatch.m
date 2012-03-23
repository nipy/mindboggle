function a = extractPitsAndSaveBatch()

pathMain = '/Users/yrjo/yrjo_work/mindboggle2/miccai2012/data/KKI/';

for lambda = 1:.5:5
    disp(['current lambda: ' num2str(lambda)]);
    %lambda = 3.5;
    
    path = ([pathMain 'KKI2009-11']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-14']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-15']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-16']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-17']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-18']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-19']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-20']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-21']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-22']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-24']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-25']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-27']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-29']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-31']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-33']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-34']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-36']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-37']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-40']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-41']);
    extractPitsAndSave(path,lambda);
    path = ([pathMain 'KKI2009-42']);
    extractPitsAndSave(path,lambda);
    
end

a = 0;