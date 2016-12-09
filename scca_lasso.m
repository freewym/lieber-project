function [ U,V ] = scca_lasso( K,invCxx05_diag,invCyy05,U_init,V_init,lambda_u,lambda_v )
% sparse cca lasso
% input: K - pxq covariance matrix
%        invCxx05_diag - diagonal of invCxx05, for memory issue
%        invCyy05 - inverse of Cyy^(-0.5)
%        U_init - px1 initial U
%        V_init - qx1 initial V
%        lambda_u - sparsity para for U
%        lambda_v - sparsity para for V
% output: U - px1
%         V - qx1 

eps = 0.0001;
max_iter = 100;
diff_u=eps*10;
diff_v=eps*10;

p=size(K,1);
q=size(K,2);
K1=zeros(p,q);
for i=1:p
    K1(i,:)=K(i,:).*(invCxx05_diag(i)^2);
end
K2=invCyy05.^2*K';
for i=1:max_iter
    if diff_u<eps && diff_v<eps, break, end
    
    % Update left singular vector
    U=K1*V_init;
    % soft shresholding
    U_sign=sign(U);
    U_sign(U_sign==0)=-1;
    U=abs(U)-lambda_u.*(invCxx05_diag.^2);
    U=(U+abs(U))/2;
    U=U.*U_sign;
    % normalize U
    U_norm=norm(invCxx05_diag.^(-1) .* U);
    if U_norm == 0, U_norm=1; end
    U=U/U_norm;
   
    % Update right singular vector
    V=K2*U;
    % soft shresholding
    V_sign=sign(V);
    V_sign(V_sign==0)=-1;
    V=abs(V)-lambda_v.*(diag(invCyy05).^2);
    V=(V+abs(V))/2;
    V=V.*V_sign;
    % normalize V
    V_norm=norm(diag(invCyy05).^(-1).*V);
    if V_norm == 0, V_norm=1; end
    V=V/V_norm;
    
    % Convergence measures 
    diff_u=max(abs(U_init - U));
    diff_v=max(abs(V_init - V));
    
    U_init = U;
    V_init = V;
end


end

