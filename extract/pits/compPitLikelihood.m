function L = compPitLikelihood(Cgauss,Cmean,Depth)

visualize = 0;
testSurface = 0;

adaptive = 1;


if (adaptive)
    
    depthSearch = 0;
    massLeft = 1;
    while(massLeft > .15)
        depthSearch = depthSearch + .01;
        massLeft = sum(Depth > depthSearch) / length(Depth);
    end
    halfValDepth = depthSearch;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    while(massLeft > .05)
        depthSearch = depthSearch + .01;
        massLeft = sum(Depth > depthSearch) / length(Depth);
    end
    slopeDepth = (-1/(depthSearch - halfValDepth)) * log((1/.80)-1);
    
    disp([halfValDepth slopeDepth]);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    massLeft = 1;
    depthSearch = 0;
    
    while(massLeft > .15)
        depthSearch = depthSearch + .0001;
        massLeft = sum(Cgauss > depthSearch) / length(Cgauss);
    end
    slopeCgauss = (-1/(depthSearch)) * log((1/.90)-1);
    
    
   % slopeCgauss = 1300;        
    halfValCgauss = .000;   
    disp([halfValCgauss slopeCgauss]);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     slopeDepth = 3;
%     slopeCgauss = 1100;
    
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    halfValCmean = .0;
    slopeCmean = 100;
    
    disp('ok');
    
else
    if (testSurface)
        
        halfValCgauss = .002;
        halfValCmean = .05;
        halfValDepth = 2.0;
        
        slopeCgauss = 1000;
        slopeCmean = 100;
        slopeDepth = 4;
        
    else
        
        halfValCgauss = .000;
        halfValCmean = .0;
        halfValDepth = .61;
        %halfValDepth = .0;
        
        slopeCgauss = 1300;
        slopeCmean = 100;
        slopeDepth = 3.8;
        %slopeDepth = 1.0;
    end
                
end

if (visualize == 1)
    
    test1 = -.005:.001:.01;
    st1 = 1./(1+exp(-slopeCgauss*(test1 - halfValCgauss)));
    %figure;plot(testi,st1);
    
    test2 = -.05:.001:.1;
    st2 = 1./(1+exp(-slopeCmean*(test2 - halfValCmean)));
    %figure;plot(testi,st2);
    
    test3 = -1:.1:4;
    st3 = 1./(1+exp(-slopeDepth*(test3 - halfValDepth)));
    %figure;plot(testi,st3);
    
    figure;
    subplot(1,3,1);
    plot(test1,st1);title('Cgauss');
    subplot(1,3,2);
    plot(test2,st2);title('Cmean');
    subplot(1,3,3);
    plot(test3,st3);title('Depth');
end

% %%%%%%%%
% Cmean = Cmean*-1;
% CgaussTemp = Cgauss;
% CgaussTemp(Cgauss < 0) = 0;
% CgaussTemp(Cmean < 0) = CgaussTemp(Cmean < 0) *-1;
% Cgauss(Cgauss > 0) = CgaussTemp(Cgauss > 0);
% %%%%%%%%




st1 = 1./(1+exp(-slopeCgauss*(Cgauss - halfValCgauss)));
%st2 = 1./(1+exp(-slopeCmean*(Cmean - halfValCmean)));
st3 = 1./(1+exp(-slopeDepth*(Depth - halfValDepth)));

% disp(max(st1))
% disp(max(st3))
% 
% disp(size(st1));
% disp(size(st2));
% disp(size(st3));

%L = st1 .* st2 .* st3;


L = st1 .* st3;
%L = st2 .* st3;

