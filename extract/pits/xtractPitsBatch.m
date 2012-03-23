function a = xtractPitsBatch()

for i = 1:12

    subjectId = i;
    
    surf_path = ['/Users/yrjo/Public/2011_HBM/CUMC12_' num2str(subjectId) '/surf/lh.pial'];
    curv_path = ['/Users/yrjo/Public/2011_HBM/CUMC12_' num2str(subjectId) '/surf/lh.curv'];
    
    [vertices faces curv checkList] = xtractPitsSurf(surf_path,curv_path,subjectId,1);

end

for i = 1:12

    subjectId = i;
    
    surf_path = ['/Users/yrjo/Public/2011_HBM/CUMC12_' num2str(subjectId) '/surf/rh.pial'];
    curv_path = ['/Users/yrjo/Public/2011_HBM/CUMC12_' num2str(subjectId) '/surf/rh.curv'];
    
    [vertices faces curv checkList] = xtractPitsSurf(surf_path,curv_path,subjectId,2);

end


a = 0;