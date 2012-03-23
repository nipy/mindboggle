function Q = compOptimalHMMF(L,D,pitInds,wd)

stopCrit =  .003;
stepSize = .01;
likelihoodLimit = .5;
gradModifier = .1;
maxIter = 5000;
dspread = 75;
%wd = 3.5;

visualize = 0;

histxaxis = 0:0.05:1;

Q = zeros(length(L),1);

m = length(pitInds);

for i= 1:m
    D(i,:) = exp(-1.*((D(i,:)').^2)./(2*dspread));
    
    if (L(pitInds(i)) > likelihoodLimit)
        Q(pitInds(i)) = L(pitInds(i));
    end
end

Q = Q(pitInds);
%%%%%%%%%
% Q(Q>.5) = Q(Q>.5) - .25 + randn(length(Q(Q>.5)),1)/13;
% Q(Q<0) = 0;
% Q(Q>1) = 1;
%Q(Q>0.5) = rand(length(Q(Q>0.5)),1);
%%%%%%%%%

Lorig = L;
L = L(pitInds);

disp('Initial candidates:');
disp(sum(Q>0));

globalProbPrev = compGlobalProb(Q,L,D,wd);
globalProb = globalProbPrev;

diff = stopCrit + 1;

iterationCount = 0;

QNEW = Q;

tic

qchange = 1;

while(qchange == 1)
    while (diff > stopCrit && iterationCount < maxIter)
        
        if (iterationCount > 250)
            
            stepSize = .02;
        elseif (iterationCount > 500)
            stepSize = .05;
        elseif (iterationCount > 750)
            stepSize = .08;
        elseif (iterationCount > 1000)
            stepSize = 0.1;
        end
        
        iterationCount = iterationCount + 1;
        
        downCounter = 0;
        upCounter = 0;
        
        for i= 1:m
            
            %if (iterationCount < 10 || (Q(i) > 0.01 && Q(i) < .99))
            %if (1)
            %if (Q(i) > 0.01 && (iterationCount < 10 || Q(i) < .99))
            if (L(i) > 0.5)
                QUp = Q;
                QDown = Q;
                
                q = Q(i);
                qUp = min(q + stepSize, 1);
                qDown = max(q - stepSize, 0);
                
                QUp(i) = qUp;
                QDown(i) = qDown;
                
                gPUp = compGlobalProb(QUp,L,D,wd);
                gPDown = compGlobalProb(QDown,L,D,wd);
                
                
                if (gPUp > globalProb && gPUp > gPDown)
                    
                    mult = (gPUp - globalProb)/gradModifier;
                    mult = mult*(0.75 + .25*rand);
                    mult = max(mult,.2);
                    
                    QNEW(i) = min(q + mult*stepSize, 1);
                    
                    upCounter = upCounter + 1;
                    
                elseif (gPDown > globalProb && gPDown > gPUp)
                    
                    mult = (gPDown - globalProb)/gradModifier;
                    mult = mult*(0.75 + .25*rand);
                    mult = max(mult,.2);
                    
                    
                    QNEW(i) = max(q - mult*stepSize, 0);
                    
                    downCounter = downCounter + 1;
                    
                end
                
            end
        end
        
        Q = QNEW;
        
        %disp('**');
        %disp(Q(Q>0));
        
        
        %     disp('Iteration:');
        %     disp(iterationCount);
        %
        %     disp('Up and down movement');
        %     disp(upCounter);
        %     disp(downCounter);
        %
        %     disp('Current pits (q > .5):');
        %     disp(sum(Q>0.5));
        
        globalProb = compGlobalProb(Q,L,D,wd);
        
        diff = globalProb - globalProbPrev;
        
        globalProbPrev = globalProb;
        
        %     disp('Current global probability:');
        %    disp(globalProb);
        
        if (visualize)
            disp(['Iteration: ' num2str(iterationCount) ' ** Up and down movement: ' num2str(upCounter) ', ' num2str(downCounter) ' ** Current pits (q > .5): ' num2str(sum(Q>0.5)) ' ** Current global probability: ' num2str(globalProb)]);
            figure(11);hist(Q(L>0.5),histxaxis);
        end
        %pause(.01);
    end
    
    qchange = 0;
    
    %disp('L to 0');
    L(Q < 0.05) = 0.001;
    
    QTEST = Q;
    QTEST(QTEST < 0.05) = 2;
    QTEST(QTEST > 0.5) = 2;
    [qtestminval qtestminInd] = min(QTEST);
    
    if (qtestminval < 2 && iterationCount < maxIter)
        L(qtestminInd) = 0.001;
        Q(qtestminInd) = 0.001;
        QNEW(qtestminInd) = 0.001;
        if (visualize)
            disp(qtestminval);
            disp(qtestminInd);
        end
        qchange = 1;
        stopCrit =  .00001;
        diff = stopCrit + 1;
        globalProb = compGlobalProb(Q,L,D,wd);
        globalProbPrev = globalProb;
    end
    
    
    
    
    
    
    % counter = 0;
    %     while(qchange == 0 && counter < m)
    %         counter = counter + 1;
    %
    %
    %
    %         if (Q(counter)<.5 && L(counter)>.5)
    %             %             QTEST = Q;
    %             %             QTEST(counter) = .8;
    %             %             QTEST2 = QTEST;
    %             %             QTEST2(QTEST2 < .5) = 0;
    %             %             LTEST = L;
    %             %             LTEST(QTEST2 < .5) = 0;
    %             %
    %             %
    %             %             testProb = compGlobalProb(QTEST2,LTEST,D,wd);
    %             %             if (testProb > globalProb)
    %             %                 Q = QTEST;
    %             %                 qchange = 1;
    %             %                 diff = stopCrit + 1;
    %             %
    %             %                 globalProb = compGlobalProb(Q,L,D,wd);
    %             %                 iterationCount = 0;
    %             %                 stepSize = 0.01;
    %             %                 disp('BOGEY FOUND!!');
    %             %             end
    %
    %             disp('L to 0');
    %             L(Q == 0) = 0;
    %             L(counter) = 0;
    %             globalProbPrev = compGlobalProb(Q,L,D,wd);
    %             qchange = 1;
    %             diff = stopCrit + 1;
    %         end
    %     end
end



disp(['** Found pits (q > .5): ' num2str(sum(Q>0.5)) ', lambda value: ' num2str(wd)]);
toc
QOut = zeros(length(Lorig),1);
QOut(pitInds) = Q;
Q = QOut;
