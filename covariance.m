function [K,invCxx05,invCyy05]=covariance(X,Y)
    p=size(X,2);
    q=size(Y,2); 
    z=[X';Y'];
    try
        C=cov(z.');
    catch
        disp('Start to compute covariance in least memory requirement...');
        C=zeros(p+q,p+q);
        for i = 1:p+q
            currentVector = z(i, :);
            inds = (i:n);
            % Loop over nInd
            for nInd = inds
                temp = currentVector*z(nInd, :)';
                %% Utilise the symmetric nature of the covariance matrix
                C(i, nInd) = temp;
                C(nInd, i) = temp;
            end
        end
        C=C./(size(X,1)-1);
        disp('Finished.');
    end
    C=(C+C')/2;
    Cxx=C(1:p,1:p);invCxx05=diag(diag(Cxx))^(-0.5);
    Cyy=C(p+1:p+q,p+1:p+q);diagCyy=diag(Cyy);diagCyy(find(diagCyy<=0))=1e-3;invCyy05=diag(diagCyy)^(-0.5);
    Cxy=C(1:p,p+1:p+q);
    %invCxx05=eye(p);invCyy05=eye(q);
    K=invCxx05*Cxy*invCyy05;
end
