% scca_lasso n=500 p=50000 q=1000 init=svd lambda_u=0.002 lambda_v=0.06
% scca_lasso n=500 p=10000 q=1000 init=avg lambda_u=0.01 lambda_v=0.06
% scca_lasso n=200 p=200 q=200 init=svd lambda_u=0.186 lambda_v=0.186
% scca       n=500 p=50000 q=1000 init=avg lambda_u=0.007 lambda_v=0.06
% scca       n=500 p=10000 q=1000 init=avg lambda_u=0.015 lambda_v=0.06
% scca       n=200 p=200 q=200 init=avg lambda_u=0.12  lambda_v=0.12
% pmi        n=500 p=50000 q=1000 lambda_u=0.3 lambda_v=0.3
% pmi        n=500 p=10000 q=1000 lambda_u=0.3 lambda_v=0.3
% pmi        n=200 p=200   q=200 lambda_u=0.3 lambda_v=0.3
% pmi        n=500 p=50000 q=1000 discete lambda_u=0.3 lambda_v=0.3
% pmi        n=500 p=10000 q=1000 discete lambda_u=0.8 lambda_v=0.8
n=500;
p=50000;
q=1000;
alg='scca_lasso';%'scca_lasso'
init='svd';%'svd'
criterion='perm';%'cv'
mode='normal';%'write' 'normal'

[X,Y,U,V]=simulation(n,p,q);
if strcmpi(mode,'write')
    dlmwrite('X.txt',X,'delimiter',' ');
    dlmwrite('Y.txt',Y,'delimiter',' ');
    return;
elseif strcmpi(mode,'read')
    U_hat=dlmread('U_hat');
    V_hat=dlmread('V_hat');
    U_hat=U_hat/norm(U_hat);V_hat=V_hat/norm(V_hat);
    fprintf('num of non-zero elements in U_hat=%d\n',sum(U_hat ~= 0));
    fprintf('num of non-zero elements in V_hat=%d\n',sum(V_hat ~= 0));
    if (U'*U_hat) < 0, U_hat=-U_hat;end;if (V'*V_hat) < 0, V_hat=-V_hat;end
    fprintf('cos(U,U_hat)=%f\n',(U'*U_hat)/(norm(U)*norm(U_hat)));
    fprintf('cos(V,V_hat)=%f\n',(V'*V_hat)/(norm(V)*norm(V_hat)));
    fprintf('recall U=%d,V=%d\n',sum((U~=0) & (U_hat~=0))/sum(U~=0),sum((V~=0) & (V_hat~=0))/sum(V~=0));
    fprintf('precision U=%d,V=%d\n',sum((U~=0) & (U_hat~=0))/sum(U_hat~=0),sum((V~=0) & (V_hat~=0))/sum(V_hat~=0));
    fprintf('max_corr on whole data=%f\n',abs(corr(X*U_hat,Y*V_hat)));

    figure
    subplot(2,2,1);
    plot (U,'r');
    subplot(2,2,2);
    plot(U_hat,'b');
    subplot(2,2,3);
    plot(V,'r');
    subplot(2,2,4);
    plot(V_hat,'b');
    return
end

lambda_u_seq=0.002:0.002:0.002;%0.186:0.02:0.186;%0.12:0.02:0.12;
lambda_v_seq=0.06:0.02:0.06;%0.186:0.02:0.186;%0.12:0.02:0.12;
num_lambdas_u=length(lambda_u_seq);
num_lambdas_v=length(lambda_v_seq);
grid_mat=zeros(num_lambdas_u,num_lambdas_v);
if strcmpi(criterion,'cv')
nfolds=5;
for f=1:nfolds
    train_idx=true(n,1);
    train_idx((f-1)*n/nfolds+1:f*n/nfolds)=false;
    test_idx = ~train_idx;
    X_train=X(train_idx,:);X_test=X(test_idx,:);
    Y_train=Y(train_idx,:);Y_test=Y(test_idx,:);
    [K,Cxy,invCxx05_diag,invCyy05]=factorized_covariance(X_train,Y_train);
    if strcmpi(alg,'scca_lasso')
        U_init=Cxy*(ones(q,1)/q);V_init=Cxy'*(ones(p,1)/p);
        if strcmpi(init,'svd')
            [U_init,~,V_init]=svds(Cxy,1);
        end
        U_init=U_init / norm(invCxx05_diag.^(-1) .* U_init);V_init=V_init/norm(diag(invCyy05).^(-1).*V_init);
    elseif strcmpi(alg,'scca')
        U_init=K*(ones(q,1)/q);V_init=K'*(ones(p,1)/p);
        U_init = U_init /norm(U_init);V_init = V_init /norm(V_init);     
    end
    for i =1:num_lambdas_u
        flag_nan = 0;
        for j=1:num_lambdas_v
            lambda_u = lambda_u_seq(i);
            lambda_v = lambda_v_seq(j);
            if flag_nan == 0
                if strcmpi(alg,'scca_lasso')
                    [ U_hat,V_hat ] = scca_lasso( Cxy,invCxx05_diag,invCyy05,U_init,V_init,lambda_u,lambda_v );
                elseif strcmpi(alg,'scca')
                    [U_hat,V_hat]=scca(K,U_init,V_init,lambda_u,lambda_v);
                    U_hat=invCxx05_diag.*U_hat;V_hat=invCyy05*V_hat;
                end
                U_hat=U_hat/norm(U_hat);V_hat=V_hat/norm(V_hat);
                %corr(X_test*U_hat,Y_test*V_hat)
                grid_mat(i,j)=grid_mat(i,j)+ abs(corr(X_test*U_hat,Y_test*V_hat));
                if isnan(grid_mat(i,j)), flag_nan = 1;end
            else
                grid_mat(i,j:end)=grid_mat(i,j:end)+NaN;
                break;
            end
        end
    end
end

elseif strcmpi(criterion,'perm')
rng(1);
T=1;
for i =1:num_lambdas_u
    flag_nan = 0;
    for j=1:num_lambdas_v
        lambda_u = lambda_u_seq(i);
        lambda_v = lambda_v_seq(j);
        if flag_nan == 0
            [K,Cxy,invCxx05_diag,invCyy05]=factorized_covariance(X,Y);
            if strcmpi(alg,'scca_lasso')
                U_init=Cxy*(ones(q,1)/q);V_init=Cxy'*(ones(p,1)/p);
                if strcmpi(init,'svd')
                    [U_init,~,V_init]=svds(Cxy,1);
                end
                U_init=U_init / norm(invCxx05_diag.^(-1) .* U_init);V_init=V_init/norm(diag(invCyy05).^(-1).*V_init);
            elseif strcmpi(alg,'scca')
                U_init=K*(ones(q,1)/q);V_init=K'*(ones(p,1)/p);
                U_init = U_init /norm(U_init);V_init = V_init /norm(V_init);
            end
            if strcmpi(alg,'scca_lasso')
                [ U_hat,V_hat ] = scca_lasso( Cxy,invCxx05_diag,invCyy05,U_init,V_init,lambda_u,lambda_v );
            elseif strcmpi(alg,'scca')
                [U_hat,V_hat]=scca(K,U_init,V_init,lambda_u,lambda_v);
                %U_hat=invCxx05_diag.*U_hat;V_hat=invCyy05*V_hat;
            end
            U_hat=U_hat/norm(U_hat);V_hat=V_hat/norm(V_hat);
            grid_mat(i,j)=abs(corr(X*U_hat,Y*V_hat));
            if isnan(grid_mat(i,j)), flag_nan = 1; end
        else
            grid_mat(i,j:end)=grid_mat(i,j:end)+NaN;
        end
        tmp=zeros(T,1);
        for t=1:T
            perm=randperm(n);
            Y_test=Y(perm,:);
            [K,Cxy,invCxx05_diag,invCyy05]=factorized_covariance(X,Y_test);
            if strcmpi(alg,'scca_lasso')
                U_init=Cxy*(ones(q,1)/q);V_init=Cxy'*(ones(p,1)/p);
                if strcmpi(init,'svd')
                    [U_init,~,V_init]=svds(Cxy,1);
                end
                U_init=U_init / norm(invCxx05_diag.^(-1) .* U_init);V_init=V_init/norm(diag(invCyy05).^(-1).*V_init);
            elseif strcmpi(alg,'scca')
                U_init=K*(ones(q,1)/q);V_init=K'*(ones(p,1)/p);
                U_init = U_init /norm(U_init);V_init = V_init /norm(V_init);
            end
            if strcmpi(alg,'scca_lasso')
                [ U_hat,V_hat ] = scca_lasso( Cxy,invCxx05_diag,invCyy05,U_init,V_init,lambda_u,lambda_v );
            elseif strcmpi(alg,'scca')
                [U_hat,V_hat]=scca(K,U_init,V_init,lambda_u,lambda_v);
                %U_hat=invCxx05_diag.*U_hat;V_hat=invCyy05*V_hat;
            end
            U_hat=U_hat/norm(U_hat);V_hat=V_hat/norm(V_hat);
            tmp(t)=abs(corr(X*U_hat,Y_test*V_hat));
        end
        grid_mat(i,j) = (grid_mat(i,j) - mean(tmp));% / std(tmp);
    end
end
end

grid_mat(isnan(grid_mat))=0;
if strcmpi(criterion,'cv')
    grid_mat=grid_mat./(nfolds);
end
max_corr=max(max(abs(grid_mat)));
[best_i,best_j] = find(abs(grid_mat)==max_corr,1);
best_lambda_u=lambda_u_seq(best_i);
best_lambda_v=lambda_v_seq(best_j);
fprintf('max_corr on test=%f,lambda_u=%f,lambda_v=%f\n',max_corr,best_lambda_u,best_lambda_v);

[K,Cxy,invCxx05_diag,invCyy05]=factorized_covariance(X,Y);
if strcmpi(alg,'scca_lasso')
    U_init=Cxy*(ones(q,1)/q);V_init=Cxy'*(ones(p,1)/p);
    if strcmpi(init,'svd')
        [U_init,~,V_init]=svds(Cxy,1);
    end
    U_init=U_init / norm(invCxx05_diag.^(-1) .* U_init);V_init=V_init/norm(diag(invCyy05).^(-1).*V_init);
elseif strcmpi(alg,'scca')
    U_init=K*(ones(q,1)/q);V_init=K'*(ones(p,1)/p);
    U_init = U_init /norm(U_init);V_init = V_init /norm(V_init);
end
if strcmpi(alg,'scca_lasso')
    [ U_hat,V_hat ] = scca_lasso( Cxy,invCxx05_diag,invCyy05,U_init,V_init,best_lambda_u,best_lambda_v );
elseif strcmpi(alg,'scca')
    [U_hat,V_hat]=scca(K,U_init,V_init,best_lambda_u,best_lambda_v);
    %U_hat=invCxx05_diag.*U_hat;V_hat=invCyy05*V_hat;
end
U_hat=U_hat/norm(U_hat);V_hat=V_hat/norm(V_hat);
fprintf('num of non-zero elements in U_hat=%d\n',sum(U_hat ~= 0));
fprintf('num of non-zero elements in V_hat=%d\n',sum(V_hat ~= 0));
if (U'*U_hat) < 0, U_hat=-U_hat;end;if (V'*V_hat) < 0, V_hat=-V_hat;end
if strcmpi(alg,'scca_lasso')
    [~,idx]=sort(abs(U_hat));
    U_hat(idx(1:p*0.9)) = 0;
    [~,idx]=sort(abs(V_hat));
    V_hat(idx(1:q*0.9)) = 0;
    fprintf('num of non-zero elements after thesholding in U_hat=%d\n',sum(U_hat ~= 0));
    fprintf('num of non-zero elements after thesholding in V_hat=%d\n',sum(V_hat ~= 0));
end
fprintf('cos(U,U_hat)=%f\n',(U'*U_hat)/(norm(U)*norm(U_hat)));
fprintf('cos(V,V_hat)=%f\n',(V'*V_hat)/(norm(V)*norm(V_hat)));
fprintf('recall U=%d,V=%d\n',sum((U~=0) & (U_hat~=0))/sum(U~=0),sum((V~=0) & (V_hat~=0))/sum(V~=0));
fprintf('precision U=%d,V=%d\n',sum((U~=0) & (U_hat~=0))/sum(U_hat~=0),sum((V~=0) & (V_hat~=0))/sum(V_hat~=0));
fprintf('max_corr on whole data=%f\n',abs(corr(X*U_hat,Y*V_hat)));

figure
subplot(2,2,1);
plot (U,'r');
title('the true U');
ylim([-0.03,0.03]);
subplot(2,2,2);
plot(U_hat,'b');
title('the recovered U')
ylim([-0.03,0.03]);
subplot(2,2,3);
plot(V,'r');
title('the true V')
ylim([-0.2,0.2]);
subplot(2,2,4);
plot(V_hat,'b');
title('the recovered V')
ylim([-0.2,0.2]);