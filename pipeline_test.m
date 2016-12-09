%function pipeline_test(s)
addpath('~/spm8');

s.type_fmri = 'nback';
s.type_mask = 'ROIs_cerebrum_act';
s.type_summary = 'avg';
s.use_demong = 'yes';
s.tuning_params = 'case_control';
s.alg = 'scca_lasso';
s.uu = 0.005;%0.005
s.vv = 0.075;%0.076
s.tau_uu = 0.001;
s.tau_vv = 0.001;
s.num_comps = 10;

p = inputParser;
addParameter(p, 'type_fmri', 'nback', @(x) any(validatestring(x,{'nback','flanker','mtl'})));
addParameter(p, 'type_mask', 'ROIs_cerebrum', @(x) any(validatestring(x,{'ROIs_cerebrum','ROIs_cerebrum_grey','ROIs_cerebrum_act', 'ROIs_all','none','grey_matter'})));
addParameter(p, 'type_summary', 'avg', @(x) any(validatestring(x,{'avg','eig'})));
addParameter(p, 'use_demong', 'yes', @(x) any(validatestring(x,{'yes','no'})));
addParameter(p, 'tuning_params', 'cv', @(x) any(validatestring(x,{'manual','case_control','perm','cv','num_non_zeros'})));
addParameter(p, 'alg', 'scca', @(x) any(validatestring(x,{'scca', 'scca_lasso','group_scca_lasso'})));
addParameter(p, 'uu', 0.2, @isnumeric);
addParameter(p, 'vv', 0.2, @isnumeric);
addParameter(p, 'tau_uu', 0.2, @isnumeric);
addParameter(p, 'tau_vv', 0.2, @isnumeric);
addParameter(p, 'num_comps', 1, @isnumeric)
addParameter(p, 'mul_var_hyp_test', 0, @isnumeric)
parse(p, s);
params=p.Results;
%{
[data_fmri, data_snp, dx_csv, csv_entries, num_NC, num_SB, num_SZ]=extract_data(params.type_fmri);
[data_fmri, ROI, AAL_ROI_VI] = load_and_apply_mask(params.type_fmri, params.type_mask, params.type_summary, data_fmri);
data_fmri = extract_and_use_demong(params.type_fmri, params.use_demong, csv_entries, data_fmri);

% only use SZ and NC
data_fmri=[data_fmri(num_NC+num_SB+1:end,:);data_fmri(1:num_NC,:)];
data_snp=[data_snp(num_NC+num_SB+1:end,:);data_snp(1:num_NC,:)];
labels=[dx_csv(num_NC+num_SB+1:end,1);dx_csv(1:num_NC,1)];

p=size(data_fmri,2);
q=size(data_snp,2);
n=size(data_fmri,1);
assert(n==size(data_snp,1),'sample size should match.');

%group_diff_per_entry(data_fmri, num_SZ, 'snp');
%compute_MAF(data_snp,0.10);
%}

% sparse CCA
% allocate space for U and V and best_lambda_U and best_lambda_V
U=zeros(size(data_fmri,2),params.num_comps);V=zeros(size(data_snp,2),params.num_comps);
best_lambda_U=zeros(1,params.num_comps);best_lambda_V=zeros(1,params.num_comps);
best_tau_U=zeros(1, params.num_comps); best_tau_V=zeros(1,params.num_comps);
% initalize the residual of data as the original one
data_fmri_res=data_fmri;data_snp_res=data_snp;
%data_snp_res=data_snp_res(randperm(size(data_snp_res,1)),:);%%%%%%%%%%%%%%%%%%%%%%
for i=1:params.num_comps
    [best_lambda_u, best_lambda_v, best_tau_u, best_tau_v] = tune_scca(params.tuning_params, params.uu, params.vv, params.tau_uu, params.tau_vv, params.alg, data_fmri_res, data_snp_res, labels, num_SZ);
    [u, v] = run_scca(best_lambda_u, best_lambda_v, best_tau_u, best_tau_v, params.alg, data_fmri_res, data_snp_res);
    U(:,i) = u;V(:,i) = v;
    best_lambda_U(i) = best_lambda_u;best_lambda_V(i) = best_lambda_v;
    best_tau_U(i) = best_tau_u;best_tau_V(i)=best_tau_v;
    coeff_fmri=data_fmri_res*u;coeff_snp=data_snp_res*v;
    coeff_fmri_SZ=coeff_fmri(1:num_SZ,:);coeff_snp_SZ=coeff_snp(1:num_SZ,:);
    coeff_fmri_NC=coeff_fmri(num_SZ+1:end,:);coeff_snp_NC=coeff_snp(num_SZ+1:end,:);
    test_group_diff(i, coeff_fmri_SZ, coeff_fmri_NC, coeff_snp_SZ, coeff_snp_NC, params.mul_var_hyp_test);
    data_fmri_res=data_fmri_res-data_fmri_res*u*u';data_snp_res=data_snp_res-data_snp_res*v*v';
end
coeff_fmri=data_fmri*U;coeff_snp=data_snp*V;
coeff_fmri_SZ=coeff_fmri(1:num_SZ,:);coeff_snp_SZ=coeff_snp(1:num_SZ,:);
coeff_fmri_NC=coeff_fmri(num_SZ+1:end,:);coeff_snp_NC=coeff_snp(num_SZ+1:end,:);
result=zeros(params.num_comps,1);
for i=1:params.num_comps
    %{
    for j=1:10%%%%%
        coeff_snp_perm=coeff_snp(randperm(size(coeff_snp,1)),i);%%%%%%%%%%
        co=corr(coeff_fmri(:,i),coeff_snp_perm);%%%%%
        fprintf('%f ',co);%%%%%
    end%%%%%
    fprintf('\n');%%%%
    %}
    [rho,pval]=corr(coeff_fmri(:,i),coeff_snp(:,i));
    result(i)=abs(rho);
    fprintf('the %d-th resulting correlation is %f(p=%f), at lambda_u=%f(%d non-0s),tau_u=%f,lambda_v=%f(%d non-0s),tau_v=%f\n.',i,result(i),pval,best_lambda_U(i), sum(U(:,i)~=0),best_tau_U(i),best_lambda_V(i),sum(V(:,i)~=0),best_tau_V(i));
end

% evaluation and output
test_group_diff(params.num_comps, coeff_fmri_SZ, coeff_fmri_NC, coeff_snp_SZ, coeff_snp_NC, params.mul_var_hyp_test);
%plot_loadings(params.num_comps, coeff_fmri_SZ, coeff_fmri_NC, coeff_snp_SZ, coeff_snp_NC, U, V);
write_components(params.type_fmri, params.type_mask, params.num_comps, U, V, best_lambda_U, best_lambda_V, best_tau_U, best_tau_V, ROI);
%}
labels_int=zeros(size(labels));
for i = 1:length(labels_int)
    if strcmpi(labels{i},'SZ')
        labels_int(i)=1;
    end
end
fprintf('discriminative training on raw snp...\n');
loss = discriminative_training(data_snp, labels);
fprintf('loss=%f\n', loss);
fprintf('discriminative training on pca snp...\n');
[~,score,~]=pca(data_snp, 'NumComponents', s.num_comps-1);
loss = discriminative_training(score,labels);
fprintf('loss=%f\n', loss);
fprintf('discriminative training on cca snp...\n');
loss = discriminative_training(coeff_snp(:,1:6-1), labels);
fprintf('loss=%f\n', loss);

fprintf('discriminative training on raw fmri...\n');
loss = discriminative_training(data_fmri, labels);
fprintf('loss=%f\n', loss);
fprintf('discriminative training on pca fmri...\n');
[~,score,~]=pca(data_fmri, 'NumComponents', s.num_comps-1);
loss = discriminative_training(score,labels);
fprintf('loss=%f\n', loss);
fprintf('discriminative training on cca fmri...\n');
loss = discriminative_training(coeff_fmri(:,1:6-1), labels);
fprintf('loss=%f\n', loss);
%end



