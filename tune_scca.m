function [best_lambda_u, best_lambda_v] = tune_scca(tuning_params, uu, vv, alg, data_fmri, data_snp, labels, num_SZ)
if strcmpi(tuning_params,'manual')
    best_lambda_u=uu;
    best_lambda_v=vv;
else
    p=size(data_fmri,2);
    q=size(data_snp,2);
    % tune lambda
    rng(1);
    % sparsity para of u and v, since u and v will be normalized, the upper
    % bound can be 2
    lambda_u_seq=0.00:0.002:0.005;%0.050:0.002:0.060;%0.1:0.02:0.15;%0.003:0.01:0.003;%0.2:0.01:0.4;
    lambda_v_seq=1.5:0.05:1.7;%0.8:0.05:1.5;%0.2:0.02:0.4;%0.00:0.001:0.09;%0.2:0.01:0.5;%0.05:0.002:0.09;
    % num of paras to be tuned
    num_lambdas_u=length(lambda_u_seq);
    num_lambdas_v=length(lambda_v_seq);

    % tunning on grid matrix
    grid_mat=zeros(num_lambdas_u,num_lambdas_v);
    if strcmpi(tuning_params, 'cv')
        % Cross-validation to select optimal combination of sparseness parameters
        iters=5; % number of iters   
        k=5; % K-folder cross validation
        cp = cvpartition(labels,'kfold',k); % stratified K-folds
        for iter=1:iters
            for kk=1:k
                testing_idx=test(cp,kk);
                training_idx=training(cp,kk);
                fmri_training=data_fmri(training_idx,:);snp_training=data_snp(training_idx,:);
                [K,Cxy,invCxx05_diag,invCyy05]=factorized_covariance(fmri_training,snp_training);
                %[K,invCxx05,invCyy05]=factorized_covariance(fmri_training,snp_training);         
                % Get starting values for singular vectors as column and row means from matrix K
                U_init=K*(ones(q,1)/q);
                U_init = U_init /norm(U_init);
                V_init=K'*(ones(p,1)/p);
                V_init = V_init /norm(V_init);
    
                fmri_test=data_fmri(testing_idx,:);snp_test=data_snp(testing_idx,:);
                
                % Loops for sparseness parameter combinations
                for i =1:num_lambdas_u
                    flag_nan = 0;
                    for j=1:num_lambdas_v
                        lambda_u = lambda_u_seq(i);	% sparseness parameter for X
                        lambda_v = lambda_v_seq(j);	% sparseness parameter for Y
                        if flag_nan == 0
                            if strcmpi(alg, 'scca')
                                [U,V]=scca(K,U_init,V_init,lambda_u,lambda_v);
                            elseif strcmpi(alg, 'scca_lasso')
                                [ U,V ] = scca_lasso( Cxy,invCxx05_diag,invCyy05,U_init,V_init,lambda_u,lambda_v );
                            end
                            if norm(U)==0 || norm(V)==0, grid_mat(i,j)=NaN; flag_nan = 1;continue; end
                            if strcmpi(alg, 'scca')
                                 U=invCxx05_diag.*U;V=invCyy05*V;
                            end
                            U=U/norm(U);V=V/norm(V);
                            coeff_fmri=fmri_test*U;
                            coeff_fmri_SZ=coeff_fmri(1:num_SZ,:);coeff_fmri_NC=coeff_fmri(num_SZ+1:end,:);
                            [h,pval]=ttest2(coeff_fmri_NC,coeff_fmri_SZ);
                            %corr_mat(i,j)=corr_mat(i,j)+ abs(abs(corr(fmri_training*U,snp_training*V))-abs(corr(fmri_test*U,snp_test*V)));if abs(corr(fmri_test*U,snp_test*V))>max_corr, max_corr=abs(corr(fmri_test*U,snp_test*V));end
                            grid_mat(i,j)=grid_mat(i,j)+pval;
                            if isnan(grid_mat(i,j)), flag_nan = 1;end
                        else grid_mat(i,j:end)=NaN;break;end
                    end
                end 
            end
            cp=repartition(cp);
        end

        % Identify optimal sparseness parameter combination
        grid_mat(isnan(grid_mat))=Inf;
        grid_mat=grid_mat./(k*iters);
        min_pval= min(min(abs(grid_mat)));
        [best_i,best_j] = find(abs(grid_mat)==min_pval,1);
        best_lambda_u=lambda_u_seq(best_i);
        best_lambda_v=lambda_v_seq(best_j);
        fprintf('best lambda_u=%f, best lambda_v=%f, min pval=%f\n',best_lambda_u,best_lambda_v,min_pval);
    elseif strcmpi(tuning_params, 'perm')
        iters=10; % number of iters   
        max_corr=-inf; 
        for i =1:num_lambdas_u
            flag_nan = 0;
            for j=1:num_lambdas_v
                lambda_u = lambda_u_seq(i);	% sparseness parameter for X
                lambda_v = lambda_v_seq(j);	% sparseness parameter for Y
                if flag_nan == 0
                    [K,invCxx05,invCyy05]=covariance(data_fmri,data_snp);
                    U_init=K*(ones(q,1)/q);
                    U_init = U_init /norm(U_init);
                    V_init=K'*(ones(p,1)/p);
                    V_init = V_init /norm(V_init);
                    [U,V]=scca(K,U_init,V_init,lambda_u,lambda_v);
                    U=invCxx05*U;U=U/norm(U);V=invCyy05*V;V=V/norm(V);
                    corr_data=abs(corr(data_fmri*U,data_snp*V));
                    perm_corr_data=zeros(iters,1);
                    for iter=1:iters
                        perm=randperm(size(data_fmri,1));
                        perm_data_fmri=data_fmri(perm,:);
                        [K,invCxx05,invCyy05]=covariance(perm_data_fmri,data_snp);
                        U_init=K*(ones(q,1)/q);
                        U_init = U_init /norm(U_init);
                        V_init=K'*(ones(p,1)/p);
                        V_init = V_init /norm(V_init);
                        [U,V]=scca(K,U_init,V_init,lambda_u,lambda_v);
                        U=invCxx05*U;U=U/norm(U);V=invCyy05*V;V=V/norm(V);
                        perm_corr_data(iter)=abs(corr(perm_data_fmri*U,data_snp*V));if perm_corr_data(iter)>max_corr, max_corr=perm_corr_data(iter);end
                    end
                    if sum(isnan(perm_corr_data)) || isnan(corr_data), grid_mat(i,j)=NaN;flag_nan = 1;
                    else [h,pval]=ttest(perm_corr_data,corr_data,'Tail','left');grid_mat(i,j)=pval;end
                    %if lambda_u==0.24 && lambda_v==0.20 %%%%%%%
                    %    perm_corr_data,corr_data %%%%%%%%%%
                    %end
                else grid_mat(i,j:end)=NaN;break;end
            end
        end
        fprintf('max perm_corr=%f\n',max_corr);%%%%%%%%%%%%%%%%
        grid_mat(isnan(grid_mat))=Inf;
        min_pval=min(min(abs(grid_mat)));
        [best_i,best_j] = find(abs(grid_mat)==min_pval,1);
        best_lambda_u=lambda_u_seq(best_i);
        best_lambda_v=lambda_v_seq(best_j);
    elseif strcmpi(tuning_params, 'case_control')
        min_pval=inf;
        %[K,invCxx05,invCyy05]=covariance(data_fmri,data_snp);
        [K,Cxy,invCxx05_diag,invCyy05]=factorized_covariance(data_fmri,data_snp);
        for i =1:num_lambdas_u
            flag_nan = 0;
            for j=1:num_lambdas_v
                lambda_u = lambda_u_seq(i);	% sparseness parameter for X
                lambda_v = lambda_v_seq(j);	% sparseness parameter for Y
                if flag_nan == 0
                    U_init=K*(ones(q,1)/q);
                    U_init = U_init /norm(U_init);
                    V_init=K'*(ones(p,1)/p);
                    V_init = V_init /norm(V_init);
                    %[U_init,~,V_init]=svds(Cxy,1); %%%
                    if strcmpi(alg, 'scca')
                        [U,V]=scca(K,U_init,V_init,lambda_u,lambda_v);
                    elseif strcmpi(alg, 'scca_lasso')
                        [ U,V ] = scca_lasso( Cxy,invCxx05_diag,invCyy05,U_init,V_init,lambda_u,lambda_v );
                    end
                    if norm(U)==0 || norm(V)==0, grid_mat(i,j)=NaN; flag_nan = 1;continue; end
                    if strcmpi(alg, 'scca')
                        U=invCxx05_diag.*U;V=invCyy05*V;
                    end
                    U=U/norm(U);V=V/norm(V);
                    coeff_fmri=data_fmri*U;
                    coeff_fmri_SZ=coeff_fmri(1:num_SZ,:);coeff_fmri_NC=coeff_fmri(num_SZ+1:end,:);
                    %[b,~,stats]=mnrfit([coeff_fmri_SZ;coeff_fmri_NC],[ones(size(coeff_fmri_SZ,1),1);2*ones(size(coeff_fmri_NC,1),1)]);
                    %b=glmfit([coeff_fmri_SZ;coeff_fmri_NC],[ones(size(coeff_fmri_SZ,1),1);zeros(size(coeff_fmri_NC,1),1)],'binomial');
                    [h,pval]=ttest2(coeff_fmri_NC,coeff_fmri_SZ);
                    grid_mat(i,j)=pval;%stats.p(2);
                else grid_mat(i,j:end)=NaN;break;end
            end
        end
        grid_mat(isnan(grid_mat))=Inf;
        min_pval=min(min(abs(grid_mat)));
        [best_i,best_j] = find(abs(grid_mat)==min_pval,1);
        best_lambda_u=lambda_u_seq(best_i);
        best_lambda_v=lambda_v_seq(best_j);
        fprintf('best lambda_u=%f, best lambda_v=%f, min pval=%f\n',best_lambda_u,best_lambda_v,min_pval);
    elseif strcmpi(tuning_params,'num_non_zeros')
        for i =1:num_lambdas_u
            flag_nan = 0;
            for j=1:num_lambdas_v
                lambda_u = lambda_u_seq(i);	% sparseness parameter for X
                lambda_v = lambda_v_seq(j);	% sparseness parameter for Y
                if flag_nan == 0
                    [K,invCxx05,invCyy05]=covariance(data_fmri,data_snp);
                    U_init=K*(ones(q,1)/q);
                    U_init = U_init /norm(U_init);
                    V_init=K'*(ones(p,1)/p);
                    V_init = V_init /norm(V_init);
                    [U,V]=scca(K,U_init,V_init,lambda_u,lambda_v);
                    if norm(U)==0 || norm(V)==0, flag_nan = 1, continue; end
                    U=invCxx05*U;U=U/norm(U);V=invCyy05*V;V=V/norm(V);
                    if sum(U~=0) <= uu && sum(V~=0) <= vv
                        best_lambda_u = uu;
                        best_lambda_v = vv;
                        found = 1;
                        break;
                    end
                else break; end
            end
            if found, break;end
        end
    end
end

end


