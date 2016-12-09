function [X,Y,U,V] = simulation(n,p,q)
% input:  n - sample size
%         p - dim for X
%         q - dim for Y
% output: X - n x p
%         Y - n x q
rng(1);
sigma_e=0.6; % noise level of latent variable relative to 1.0
sparsity1=0.1; % percentage of non-zero elements in weight vector of X
sparsity2=0.1; % percentage of non-zero elements in weight vector of Y
sigma1=1; % variance of irrelavant variables in X
sigma2=1; % variance pf irrelant variables in Y

% latent variable E[(s+e1)(s+e2)]/sqrt(E[(s+e1)(s+e1)]E[(s+e2)(s+e2)])=1/(1+sigma_e^2)
gamma=normrnd(0,1,[n,1]);
gamma1=gamma + normrnd(0,sigma_e,[n,1]);
gamma2=gamma + normrnd(0,sigma_e,[n,1]);
corr=corrcoef(gamma1,gamma2);
fprintf('corr of latent variables=%f\n',corr(1,2));


X=zeros(n,p);Y=zeros(n,q);
irre_cols_X=false(p,1);irre_cols_Y=false(q,1);
% weight vectors
U=zeros(p,1);
for i=1:p
    if rand <= sparsity1
        U(i) = rand * 2 - 1;
    else
        X(:,i) = normrnd(0,sigma1,[n,1]);
        irre_cols_X(i) = true;
    end
end
%X(:,irre_cols_X) = mvnrnd(zeors(1,p),SIGMA1,n);
fprintf('num of non-zero elements in U=%d\n',sum(U ~= 0));
U=U/norm(U);
V=zeros(q,1);
for i=1:q
    if rand <= sparsity2
        V(i) = rand * 2 - 1;
    else
        Y(:,i) = normrnd(0, sigma2, [n,1]);
        irre_cols_Y(i) = true;
    end
end
%Y(:,irre_cols_Y) = mvnrnd(zeors(1,q),SIGMA2,n);
fprintf('num of non-zero elements in V=%d\n',sum(V~=0));
V=V/norm(V);

for i=1:n
    X(i,:) = X(i,:) + gamma1(i) * U';
    Y(i,:) = Y(i,:) + gamma2(i) * V';
end

MAF=0.2+0.2*rand(q,1);
for i=1:q
    n1 = round(MAF(i)^2 * n);
    n2 = round(2 * (1-MAF(i)) * MAF(i) * n);
    n3 = round((1-MAF(i))^2 * n);
    [~,I]=sort(Y(:,i),'descend');
    Y(I(1:n1),i) = 1;
    Y(I(n1+1:n1+n2),i) = 0;
    Y(I(n1+n2+1:end),i) = -1;
end
%}


        
