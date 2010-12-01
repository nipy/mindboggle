n = -pi:pi

N=10;
q=0;

T = 2*pi*q/N

for k=1:2:10
    x=sin(2*pi*n*k/N +T);
    plot(n,x)
    hold on
end