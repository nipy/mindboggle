function [fv A Depth] = makeTestSurf()

% asize = [26 50];
% 
% cent1 =1.1*[asize(1)/2 asize(2)/4];
% cent2 =1.1*[asize(1)/2 asize(2)*3/4];
% 
% spread= 4;
% 
% 
% A = ones(asize);
% 
% for i = 1:asize(1)
%     for j = 1:asize(2)
%         A(i,j) = A(i,j) - exp(-1*compEucDist(cent1,[i j])/spread);
%         A(i,j) = A(i,j) - exp(-1*compEucDist(cent2,[i j])/spread);
%     end
% end
% 
% figure; surf(A);
% 

asize = [50 40 20];
spread= 2;

depth1 = 7;
depth2 = 10;
depth3 = 16;
depth4 = 11;


A = zeros(asize);

coordMap = zeros(asize(1),asize(2));
cent1 =1*[asize(1)/2 asize(2)/4];
cent2 =1*[asize(1)/2 asize(2)/4+15];
cent3 =1*[asize(1)/2 asize(2)/4+7];
cent4 =1*[asize(1)/2+9 asize(2)/4+5];


disp([cent1;cent2;cent3;cent4])

%cent2 =1.1*[asize(1)/2 asize(2)*3/4];

for i = 1:asize(1)
    for j = 1:asize(2)
        coordMap(i,j) = asize(3)/2 - depth1*exp(-1*compEucDist(cent1,[i j])/spread);
        coordMap(i,j) = coordMap(i,j) - depth2*exp(-1*compEucDist(cent2,[i j])/spread);
        coordMap(i,j) = coordMap(i,j) - depth3*exp(-1*compEucDist(cent3,[i j])/spread);
        coordMap(i,j) = coordMap(i,j) - depth4*exp(-1*compEucDist(cent4,[i j])/spread);
        
        A(i,j,1:round(coordMap(i,j))) = 1;
        
    end
end


%A(asize(1)/2:asize(1),:,:) = 1;

%cent1 = [asize(1)/2 asize(2)/2 asize(3)/4];
%cent1 =1.1*[asize(1)/2 asize(2)/4];
%cent2 =1.1*[asize(1)/2 asize(2)*3/4];

% for i = 1:asize(1)
%     for j = 1:asize(2)
%         for k = 1:asize(3)
%             if (compEucDist(cent1,[i j k]) < spread)
%                 A(i,j,k) = 0;
%                 
%             end
%             %A(i,j) = A(i,j) - exp(-1*compEucDist(cent1,[i j])/spread);
%         %A(i,j) = A(i,j) - exp(-1*compEucDist(cent2,[i j])/spread);
%         end
%     end
% end

A = smooth3(A,'gaussian',5,3);

% surf2D = zeros(asize(1),asize(2));
% for i = 1:asize(1)
%     for j = 1:asize(2)
%         temp = A(i,j,:);
%         surf2D(i,j) = max(find(temp(:)>0));
%     end
% end
% figure;surf(surf2D);

fv = isosurface(A,.5);

figure;
p = patch(isosurface(A,.5));
isonormals(A,p)
set(p,'FaceColor','red','EdgeColor','none');
daspect([1 1 1])
view(3); axis tight
camlight 
lighting gouraud

Depth = abs(fv.vertices(:,3)-10.5);


