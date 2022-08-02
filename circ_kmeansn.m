function [c X]=circ_kmeansN(x,K)

%%%%% USAGE:   [c X]=circ_kmeansN(x,K)
%%%%%  where  x   are input vectors arranged in columns [x(t1) x(t2) .... x(tN)].
%%%%%         K   is the number of desired clusters.
%%%%%         c   contains the cluster centres at each column.
%%%%%         X   is a 3D array containing the input vectors segmented into K clusters.
%%%%%               Contains 0 where the corresponding point has been
%%%%%               attributed to another cluster.
%%%%%
%%%%%  Implements the paper:
%%%%%  Mitianoudis N., "A Generalised Directional Laplacian Distribution: 
%%%%%  Estimation, Mixture Models and Audio Source Separation" , IEEE Trans.
%%%%%  on Audio, Speech and Language Process., Vol. 20, No. 9, pp. 2397- 2408, 
%%%%%  Nov. 2012.
%%%%%
%%%%%  Code by Nikolaos Mitianoudis, Democritus University of Thrace,
%%%%%  Greece, 2013

maxIter=100;
[DD N]=size(x);
c=randn(DD,K);
c=c./(ones(DD,1)*sqrt(sum(c.^2)));

dff=9;
DF=9;
j=0;
while j<maxIter
    j=j+1;
    for i=1:K
        D(i,:)=sqrt(1-(c(:,i)'*x).^2);
    end
    
    for i=1:N
        tmp=min(D(:,i));
        tmp2=find(tmp==D(:,i));
        INDX(i)=tmp2(1);
    end
        
    X=zeros(DD,size(x,2),K);
    for i=1:K
        t=find(INDX==i);
        X(:,t,i)=x(:,t);
%
        c(:,i)=sum(x(:,t),2)/length(t);
        c(:,i)=c(:,i)./sqrt(sum(c(:,i).^2));
    end
    
%     dff=abs(DF-mean(mean(D,2)));
%     DF=mean(mean(D,2));
end
    
    