function [a m k]=ml_mixdirlaplace_fun(x,K)

%%%%% USAGE:   [a m k]=ml_mixcirclaplace_fun(x,K)
%%%%%  where  x   are input vectors arranged in columns [x(t1) x(t2) .... x(tN)].
%%%%%             They must be directional vectors, i.e. ||x||=1    
%%%%%         a   is a 1xK matrix containing the mixture weights of the K individual DLDs.
%%%%%         m   is a DxK matrix contains the K individual DLDs' mean vectors (D-dim).
%%%%%         k   is a 1xK matrix containing the mixture weights of the K
%%%%%         individual DLDs.
%%%%%
%%%%%  Implements the paper:
%%%%%  Mitianoudis N., "A Generalised Directional Laplacian Distribution: 
%%%%%  Estimation, Mixture Models and Audio Source Separation" , IEEE Trans.
%%%%%  on Audio, Speech and Language Process., Vol. 20, No. 9, pp. 2397- 2408, 
%%%%%  Nov. 2012.
%%%%%
%%%%%  Code by Nikolaos Mitianoudis, Democritus University of Thrace,
%%%%%  Greece, 2013



%%%% Loads some pre-calculated bessel integrals for the estimation of k.
load I1I0_integ;

[Q,M]=size(x);
if Q>11 | Q==1
    disp('The algorithm can handle from 2d up to 11-d directional vectors. You need to numerically evaluate the integral for greater dimensions')
    exit;
end

I1Io=[G(1,:);G(Q+1,:)];

maxIter=100;

m=[];

[m X]=circ_kmeansn(x,K); %%%% Initialise the algorithm using Circular k_means.
X=[];

k=ones(1,K)*5;
a=ones(1,K)/K;

pa=[];pdmi=[];pk=[];



for i=1:maxIter
i
 %%% Estimate a   
    for j=1:K
       Io(j)=interp1(I1Io(1,:),I0(1,:),k(j),'spline');
       prb(j,:)=a(j)*exp(-k(j)*sqrt(1-(m(:,j)'*x).^2))./Io(j);      
    end
    prb=prb./ (ones(K,1)*sum(prb));
    a=mean(prb,2);
     

%%% Estimate m
    for j=1:K
      thm=m(:,j)'*x;
      tmp=thm./sqrt(1-thm.^2);
      tmp=tmp.*prb(j,:);
      dJdm(:,j)=k(j)*mean((ones(Q,1)*tmp).*x,2);
    end
    m2=m+0.01*dJdm;

    for j=1:K
        m2(:,j)=m2(:,j)./norm(m2(:,j));
    end
    dm=sum(abs(m-m2));
    m=m2;
 %%% Estimate k
     for j=1:K
        I=sum(sqrt(1-(m(:,j)'*x).^2).*prb(j,:))/sum(prb(j,:)); 
        k(j)= interp1(I1Io(2,:),I1Io(1,:),I,'spline');
        if k(j)<1 
            k(j)=2;
        end
        if k(j)>25
            k(j)=25;
        end
     end


  %%%%% Convergence check variables  
    pa=[pa ;a'];pk=[pk ;k];pdmi=[pdmi ;dm]; 
    
  
end
 figure
  plot(pdmi);
figure
plot(pa);
figure
plot(pk);

