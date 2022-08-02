function xp=reduce_points(x,ratio,frame)

X1=reshape(x(1,:),frame,size(x,2)/frame);
X2=reshape(x(2,:),frame,size(x,2)/frame);
X1p=[];
X2p=[];
for i=1:size(X1,2)
    pp=[X1(:,i) X2(:,i)]';
    Cp=pp*pp'/1024;
    [V, D]=eig(Cp);
%     tmp(i)=D(2,2)/D(1,1);
    if D(2,2)/D(1,1)>ratio
        X1p=[X1p X1(:,i)];
        X2p=[X2p X2(:,i)];
%     else
%         X1p=[X1p zeros(frame,1)];
%         X2p=[X2p zeros(frame,1)];
    end
end
xp=[reshape(X1p,1,prod(size(X1p)));reshape(X2p,1,prod(size(X2p)))];