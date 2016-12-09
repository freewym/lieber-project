function [U, V] = run_scca(best_lambda_u, best_lambda_v, best_tau_u, best_tau_v, alg, data_fmri, data_snp)
p=size(data_fmri,2);
q=size(data_snp,2);

% Compute singular vectors using the optimal sparseness parameter combination for the whole data
%[K,invCxx05,invCyy05]=factorized_covariance(data_fmri,data_snp);
[K,Cxy,invCxx05_diag,invCyy05]=factorized_covariance(data_fmri,data_snp);
% Get starting values for singular vectors as column and row means from matrix K
U_init=K*(ones(q,1)/q);
U_init = U_init /norm(U_init);
V_init=K'*(ones(p,1)/p);
V_init = V_init /norm(V_init);
%[U_init,~,V_init]=svds(Cxy,1); %%%
if strcmpi(alg, 'scca')
    [U,V]=scca(K,U_init,V_init,best_lambda_u,best_lambda_v);
    U=invCxx05_diag.*U;V=invCyy05*V;
elseif strcmpi(alg, 'scca_lasso')
    [ U,V ] = scca_lasso( Cxy,invCxx05_diag,invCyy05,U_init,V_init,best_lambda_u,best_lambda_v );
elseif strcmpi(alg, 'group_scca_lasso')
    [U, V] = group_scca_lasso(Cxy, invCxx05_diag, invCyy05,U_init, V_init, best_lambda_u, best_lambda_v, best_tau_u, best_tau_v);
end
U=U/norm(U);V=V/norm(V);
coeff_fmri=data_fmri*U;coeff_snp=data_snp*V;
result=abs(corr(coeff_fmri,coeff_snp));
fprintf('the resulting correlation is %f, at lambda_u=%f,lambda_v=%f.\n',result,best_lambda_u,best_lambda_v);
end

