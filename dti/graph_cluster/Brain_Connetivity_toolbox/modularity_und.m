function [Ci Q NodeStrength NormNodeStrength]=modularity_und(A)
%MODULARITY_UND     Optimal community structure and modularity
%
%   Ci = modularity_und(W);
%   [Ci Q Nc NNc] = modularity_und(W);
%
%   The optimal community structure is a subdivision of the network into
%   nonoverlapping groups of nodes in a way that maximizes the number of
%   within-group edges, and minimizes the number of between-group edges. 
%   The modularity is a statistic that quantifies the degree to which the
%   network may be subdivided into such clearly delineated groups. 
%
%   Input:      W,      undirected (weighted or binary) connection matrix.
%
%   Outputs:    Ci,     optimal community structure
%               Q,      maximized modularity
%               Nc,     nodewise community centrality
%               NNc,    nodewise normalized community centrality
%                       (normalized to [0,1] for each community)
%
%   Note: See Newman (2006) Phys Rev E 74:036104 for the description of
%           community centrality.
%
%   Reference: Newman (2006) Phys Rev E 74:036104; PNAS 23:8577-8582.
%
%
%   2008-2010
%   Mika Rubinov, UNSW
%   Jonathan Power, WUSTL
%   Alexandros Goulas, Maastricht University


%   Modification History:
%   Jul 2008: Original
%   Oct 2008: Positive eigenvalues made insufficient for division (Jonathan Power, WUSTL)
%   Dec 2008: Fine-tuning made consistent with Newman's description (Jonathan Power)
%   Dec 2008: Fine-tuning vectorized (Mika Rubinov)
%   Jul 2010: Incorporated community centrality (Alexandros Goulas)

K=sum(A);                               %degree
N=length(A);                            %number of vertices
m=sum(K);                               %number of edges
B=A-(K.'*K)/m;                          %modularity matrix
Ci=ones(N,1);                           %community indices
cn=1;                                   %number of communities
U=[1 0];                                %array of unexamined communites

%%community centrality
NodeStrength=zeros(N,1);                %community centrality
NormNodeStrength=zeros(N,1);            %normalized community centrality

ind=1:N;
Bg=B;
Ng=N;

while U(1)                              %examine community U(1)
    [V D]=eig(Bg);
    [d1 i1]=max(diag(D));               %most positive eigenvalue of Bg
    v1=V(:,i1);                         %corresponding eigenvector

    S=ones(Ng,1);
    S(v1<0)=-1;
    q=S.'*Bg*S;                         %contribution to modularity

    if q>1e-10                       	%contribution positive: U(1) is divisible
        qmax=q;                         %maximal contribution to modularity
        Bg(logical(eye(Ng)))=0;      	%Bg is modified, to enable fine-tuning
        indg=ones(Ng,1);                %array of unmoved indices
        Sit=S;
        
        %%community centrality
        NodeStrength(ind)=abs(v1);      %community centrality
        
        while any(indg);                %iterative fine-tuning
            Qit=qmax-4*Sit.*(Bg*Sit); 	%this line is equivalent to:
            qmax=max(Qit.*indg);        %for i=1:Ng
            imax=(Qit==qmax);           %	Sit(i)=-Sit(i);
            Sit(imax)=-Sit(imax);       %	Qit(i)=Sit.'*Bg*Sit;
            indg(imax)=nan;             %	Sit(i)=-Sit(i);
            if qmax>q;                  %end
                q=qmax;
                S=Sit;
            end
        end

        if abs(sum(S))==Ng              %unsuccessful splitting of U(1)
            U(1)=[];
        else
            cn=cn+1;
            Ci(ind(S==1))=U(1);         %split old U(1) into new U(1) and into cn
            Ci(ind(S==-1))=cn;
            U=[cn U];
        end
    else                                %contribution nonpositive: U(1) is indivisible
        U(1)=[];
    end

    ind=find(Ci==U(1));                 %indices of unexamined community U(1)
    bg=B(ind,ind);
    Bg=bg-diag(sum(bg));                %modularity matrix for U(1)
    Ng=length(ind);                     %number of vertices in U(1)
end

s=Ci(:,ones(1,N));                      %compute modularity
Q=~(s-s.').*B/m;
Q=sum(Q(:));

%%normalize community centrality to 0-1 for each community.
unique_ind=unique(Ci);
for i=1:size(unique_ind,1)
    ind=find(Ci==unique_ind(i));
    NormNodeStrength(ind) = ...
        (NodeStrength(ind) - repmat(min(NodeStrength(ind)),size(NodeStrength(ind),1),1)) ./ ...
        repmat(max(NodeStrength(ind))-min(NodeStrength(ind)),size(NodeStrength(ind),1),1);
end