clear all
clc
format long

%%% Drawing a typical third order curve and computing it derivatives in
%%% continous domain

x = -10:1:10;
y = x.^3 - 5*x.^2 + x + 4;

dy  = 3*x.^2 - 10*x +1;
ddy = 6*x -10;

figure(1)
plot(x,y, 'c')
hold on
grid on
plot(x, dy, 'b')
hold on 
plot(x, ddy, 'r')

%%% Computing the Taylor expansion estimation of the curve in the points
%%% with Delta deviation from the actual points

Delta = -.1;
OriginalEstimation = y + dy*Delta + ddy*(Delta^2);

%%% Computing the original curve at the point with Delta deviation

Dx = (-10+Delta):1:(10+Delta);
ActualValues = Dx.^3 - 5*Dx.^2 + Dx + 4;
OriginalEstimationError = norm(ActualValues(2:end-2) - OriginalEstimation(2:end-2))

%%% Compute the derivatives in the descrete domain and draw them in compare
%%% with the ones computed in continues domains

Dy = y(2:end)-y(1:end-1);
DDy = y(3:end)-2*y(2:end-1)+y(1:end-2);
DDy(21) = DDy(19);
DDy(20) = DDy(19);
Dy(21) = Dy(20);

figure(2)
plot(x, dy, 'b')
hold on
plot(x, ddy, 'b')
hold on 
plot(x, Dy, 'r')
hold on
plot(x, DDy, 'r')
grid on

%%% Compute the Taylor expansion of the curve in the desceret domain for
%%% the theoritical derivatives

DescreteDomainEstimation = y + Dy*Delta + DDy*(Delta^2);
ErrorFromActualValue1 = norm(DescreteDomainEstimation(2:end-2) - ActualValues(2:end-2))
ErrorFromEstimatedValue1 = norm(DescreteDomainEstimation(2:end-2) - OriginalEstimation(2:end-2))

%%% Compute the Taylor expansion of the curve in the desceret domain for
%%% the conventional second derivatives

DDy(2:20) = y(1:end-2)-2*y(2:end-1)+y(3:end);
DDy(21) = DDy(20);
DDy(1) = DDy(2);

DescreteDomainEstimation = y + Dy*Delta + DDy*(Delta^2);
ErrorFromActualValue2 = norm(DescreteDomainEstimation(2:end-2) - ActualValues(2:end-2))
ErrorFromEstimatedValue2 = norm(DescreteDomainEstimation(2:end-2) - OriginalEstimation(2:end-2))

%%% Compute the Taylor expansion of the curve in the desceret domain for
%%% the conventional first and  derivatives

Dy(2:20) = (y(1:end-2)-y(3:end))/2;
Dy(21) = Dy(20);
Dy(1) = Dy(2);

DescreteDomainEstimation = y + Dy*Delta + DDy*(Delta^2);
ErrorFromActualValue3 = norm(DescreteDomainEstimation(2:end-2) - ActualValues(2:end-2))
ErrorFromEstimatedValue3 = norm(DescreteDomainEstimation(2:end-2) - OriginalEstimation(2:end-2))



