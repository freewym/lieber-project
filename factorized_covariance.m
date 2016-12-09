function [K,Cxy,invCxx05_diag,invCyy05]=factorized_covariance(X,Y)
    n=size(X,1);
    p=size(X,2);
    q=size(Y,2);
    X=X-repmat(mean(X),n,1);
    Y=Y-repmat(mean(Y),n,1);
    K=zeros(p,q);
    %for i=1:n
    %    K=K+X(i,:)'*Y(i,:);
    %end
    invCxx05_diag=zeros(p,1);
    for i=1:p
        cxx=(X(:,i)'*X(:,i))/(n-1);
        if cxx <= 0
            cxx = 1e-3;
        end
        invCxx05_diag(i)=cxx^(-0.5);
    end
    %Cxx=(X'*X)./(n-1);invCxx05=diag(diag(Cxx))^(-0.5);
    Cyy=(Y'*Y)./(n-1);diagCyy=diag(Cyy);diagCyy(diagCyy<=0)=1e-3;invCyy05=diag(diagCyy)^(-0.5);
    Cxy=(X'*Y)./(n-1);
    for i=1:p
        K(i,:)=Cxy(i,:).*invCxx05_diag(i);
    end
    K=K*invCyy05;
end
