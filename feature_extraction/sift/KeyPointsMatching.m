function [MatchPersentage  MatchingPoints  NearestDistance] = KeyPointsMatching(Points_1, Points_2, MatchMargin)

[M1 N1 R1] =  size(Points_1);
[M2 N2 R2] =  size(Points_2);
Counts = 0;
MatchingPoints = 0;

if M1 > M2
    Mmax = M1;
    Mmin = M2;
    MaxPoints = Points_1;
    MinPoints = Points_2;
else
    Mmax = M2;
    Mmin = M1;
    MaxPoints = Points_2;
    MinPoints = Points_1;
end

if N1 ~= N2
    fprintf('This keypoints are not comparable');
    return;
end



if N1 == 3
    
for i = 1: Mmin
    MoreThanOneMatchNumber = 0;
    MatchIndex = 0;
    MatchDifference = 0;
    for j=1:Mmax
     
         f = MinPoints(i,:) - MaxPoints(j,:);
         if ((sqrt(f(1)^2+f(2)^2)) < (MatchMargin))
            if abs(f(3)) < (2^.5) 
                 MoreThanOneMatchNumber = MoreThanOneMatchNumber +1;
                 MatchIndex(MoreThanOneMatchNumber) = j;
                 MatchDifference(MoreThanOneMatchNumber) = sqrt(f(1)^2+f(2)^2);
                 
            end
         end
    end
    NearestDistance(i) = MatchMargin;
    if (MoreThanOneMatchNumber >0)
        Counts = Counts + 1;
        [m n] = min(MatchDifference);
        MatchingPoints(Counts, 1:3) = MinPoints(i,:);
        MaxPoints(MatchIndex(n),:) = [0 0 0];
        NearestDistance(i) = MatchDifference(n);
    end
end

end

if N1 == 4
for i = 1: Mmin
    MoreThanOneMatchNumber = 0;
    MatchIndex = 0;
    MatchDifference = 0;
    for j=1:Mmax
     
         f = MinPoints(i,:) - MaxPoints(j,:);
         if ((sqrt(f(1)^2+f(2)^2+f(3)^2)) < (MatchMargin))
             if abs(f(3)) < (2^.5) 
                 MoreThanOneMatchNumber = MoreThanOneMatchNumber +1;
                 MatchIndex(MoreThanOneMatchNumber) = j;
                 MatchDifference(MoreThanOneMatchNumber) = sqrt(f(1)^2+f(2)^2+f(3)^2);
                 
            end
         end
    end
    NearestDistance(i) = MatchMargin;
    if (MoreThanOneMatchNumber >0)
        Counts = Counts + 1;
        [m n] = min(MatchDifference);
        MatchingPoints(Counts, 1:4) = MinPoints(i,:);
        MaxPoints(MatchIndex(n),:) = [0 0 0 0];
        NearestDistance(i) = MatchDifference(n);
    end
end

end

    
    
    
    

MatchPersentage = (Counts/Mmin)*100;


