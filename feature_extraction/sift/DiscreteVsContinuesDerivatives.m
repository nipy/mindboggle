clear all
clc
format long

%%% Drawing a typical third order curve and computing it derivatives in
%%% continous domain

x = -10:1:10;

% y = x.^3 - 5*x.^2 + x + 4;
% dy  = 3*x.^2 - 10*x +1;
% ddy = 6*x -10;


y = exp(-x.^2);
dy = -2*x.*exp(-x.^2);
ddy = 4*(x.^2).*exp(-x.^2);

figure(1)
plot(x,y, 'k')
hold on
grid on

%%% Computing the Taylor expansion estimation of the curve in the points
%%% with Delta deviation from the actual points

Delta = .1;
OriginalEstimation = y + dy*Delta + ddy*(Delta^2);

%%% Computing the original curve at the point with Delta deviation

Dx = (-10+Delta):1:(10+Delta);
%ActualValues = Dx.^3 - 5*Dx.^2 + Dx + 4;
ActualValues = exp(-Dx.^2);

OriginalEstimationError = norm(ActualValues(2:end-2) - OriginalEstimation(2:end-2))

%%% Compute the derivatives in the discrete domain and draw them in compare
%%% with the ones computed in continues domains
%%% f' = x(n+1) - x(n)
%%% f" = x(n+2) -2x(n+1) + x(n)

Dy = y(2:end)-y(1:end-1);
Dy(21) = Dy(20);
DDy = y(3:end)-2*y(2:end-1)+y(1:end-2);
DDy(21) = DDy(19);
DDy(20) = DDy(19);



plot(x, dy, '-b')
hold on
plot(x, ddy, '--b')
hold on 
plot(x, Dy, '-rd')
hold on
plot(x, DDy, '--rd')
grid on

%%% Compute the Taylor expansion of the curve in the discrete domain for
%%% the theoritical derivatives


DescreteDomainEstimation = y + Dy*Delta + DDy*(Delta^2);
ErrorFromActualValue1 = norm(DescreteDomainEstimation(2:end-2) - ActualValues(2:end-2))
ErrorFromEstimatedValue1 = norm(DescreteDomainEstimation(2:end-2) - OriginalEstimation(2:end-2))

%%% Compute the Taylor expansion of the curve in the discrete domain for
%%% the conventional second derivatives
%%% f' = x(n+1) - x(n)
%%% f" = x(n+1) -2x(n) + x(n-1)


DDy(2:20) = y(1:end-2)-2*y(2:end-1)+y(3:end);
DDy(21) = DDy(20);
DDy(1) = DDy(2);

hold on 
plot(x, Dy, '-cp')
hold on
plot(x, DDy, '--cp')
grid on


DescreteDomainEstimation = y + Dy*Delta + DDy*(Delta^2);
ErrorFromActualValue2 = norm(DescreteDomainEstimation(2:end-2) - ActualValues(2:end-2))
ErrorFromEstimatedValue2 = norm(DescreteDomainEstimation(2:end-2) - OriginalEstimation(2:end-2))

%%% Compute the Taylor expansion of the curve in the discrete domain for
%%% the conventional first and  derivatives
%%% f' = (x(n+1) - x(n-1))/2
%%% f" = x(n+1) -2x(n) + x(n-1)

Dy(2:20) = (y(3:end) - y(1:end-2))/2;
Dy(21) = Dy(20);
Dy(1) = Dy(2);

hold on 
plot(x, Dy, '-ms')
hold on
plot(x, DDy, '--ms')
grid on


DescreteDomainEstimation = y + Dy*Delta + DDy*(Delta^2);
ErrorFromActualValue3 = norm(DescreteDomainEstimation(2:end-2) - ActualValues(2:end-2))
ErrorFromEstimatedValue3 = norm(DescreteDomainEstimation(2:end-2) - OriginalEstimation(2:end-2))


%%% Compute the Taylor expansion of the curve in the discrete domain for
%%% the conventional first and  derivatives
%%% f' = (x(n+1) - x(n-1))/2
%%% f" = x(n+2)/2 -x(n) + x(n-2)/2


DDy(3:19) = (y(5:end)/2 - y(3:end-2) + y(1:end-4)/2)/2;
DDy(1) = DDy(3);
DDy(2) = DDy(3);
DDy(21) = DDy(19);
DDy(20) = DDy(19);

hold on 
plot(x, Dy, '-gp')
hold on
plot(x, DDy, '--gp')
grid on


DescreteDomainEstimation = y + Dy*Delta + DDy*(Delta^2);
ErrorFromActualValue4 = norm(DescreteDomainEstimation(2:end-2) - ActualValues(2:end-2))
ErrorFromEstimatedValue4 = norm(DescreteDomainEstimation(2:end-2) - OriginalEstimation(2:end-2))

