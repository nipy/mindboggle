function P=participation_coef(A,Ci)
%PARTICIPATION_COEF     Participation coefficient
%
%   P = participation_coef(A,Ci);
%
%   Participation coefficient is a measure of diversity of intermodular
%   connections of individual nodes.
%
%   Inputs:     A,      binary/weighted, directed/undirected connection matrix
%               Ci,     community affiliation vector
%
%   Output:     P,      participation coefficient
%
%   Note: The output for directed graphs is the "out-neighbor"
%         participation coefficient.
%
%   Reference: Guimera R, Amaral L. Nature (2005) 433:895-900.
%
%
%   Mika Rubinov, UNSW, 2008-2010

n=length(A);                        %number of vertices
Ko=sum(A,2);                        %(out)degree
Ko(~Ko)=inf;                        %p_ind=0 if no (out)neighbors
Gc=A*diag(Ci);                      %neighbor community affiliation
Kc2=zeros(n,1);                     %community-specific neighbors

for i=1:max(Ci);
   Kc2=Kc2+(sum(Gc==i,2).^2);
end

P=ones(n,1)-Kc2./(Ko.^2);