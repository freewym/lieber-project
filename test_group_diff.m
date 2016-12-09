function test_group_diff(i, coeff_fmri_SZ, coeff_fmri_NC, coeff_snp_SZ, coeff_snp_NC, is_mul_var)
if ~is_mul_var
    p1=zeros(1,1);p2=zeros(1,1);
    %for i=1:num_comps
        [h,pval]=ttest2(coeff_fmri_NC(:,1),coeff_fmri_SZ(:,1));p1(1)=pval;
        fprintf('ttest2 for %d-th component of fmri:h=%d,p=%f,mean_SZ=%f,mean_NC=%f\n',i,h,p1(1),mean(coeff_fmri_SZ(:,1)),mean(coeff_fmri_NC(:,1)));
        [h,pval]=ttest2(coeff_snp_NC(:,1),coeff_snp_SZ(:,1));p2(1)=pval;
        fprintf('ttest2 for %d-th component of snp:h=%d,p=%f,mean_SZ=%f,mean_NC=%f\n',i,h,p2(1),mean(coeff_snp_SZ(:,1)),mean(coeff_snp_NC(:,1)));
    %end
else
    [h,pval]=ttest2_mul_var(num_comps,coeff_fmri_SZ,coeff_fmri_NC);
    fprintf('ttest2_mul_var of fmri:h=%d,p=%f\n',h,pval);
    [h,pval]=ttest2_mul_var(num_comps,coeff_snp_SZ,coeff_snp_NC);
    fprintf('ttest2_mul_var of snp:h=%d,p=%f\n',h,pval);
end
end

function [h,p]=ttest2_mul_var(num_comps,coeff_SZ,coeff_NC)
    n=size(coeff_SZ,1);m=size(coeff_NC,1);
    mean_coeff_SZ=repmat(mean(coeff_SZ),[n,1]);mean_coeff_NC=repmat(mean(coeff_NC),[m,1]);
    coeff_SZ=coeff_SZ-mean_coeff_SZ;coeff_NC=coeff_NC-mean_coeff_NC;
    abs_coeff_SZ=sqrt(sum(coeff_SZ.^2,2));abs_coeff_NC=sqrt(sum(coeff_NC.^2,2));
    var_coeff_SZ=nanvar(abs_coeff_SZ);var_coeff_NC=nanvar(abs_coeff_NC);
    dfe = n+m-2;
    s=sqrt(((n-1).*var_coeff_SZ+(m-1).*var_coeff_NC)./dfe);
    t=norm(mean_coeff_SZ(1,:)-mean_coeff_NC(1,:))./(s.*sqrt(1/n+1/m));
    p = 2 * tcdf(-abs(t),dfe);
    h=0;
    if p<=0.05
        h=1;
    end
end
