function test_group_diff(i, coeff_fmri_SZ, coeff_fmri_NC, coeff_snp_SZ, coeff_snp_NC)
    p1=zeros(1,1);p2=zeros(1,1);
    %for i=1:num_comps
        [h,pval]=ttest2(coeff_fmri_NC(:,1),coeff_fmri_SZ(:,1));p1(1)=pval;
        fprintf('ttest2 for %d-th component of fmri:h=%d,p=%f,mean_SZ=%f,mean_NC=%f\n',i,h,p1(1),mean(coeff_fmri_SZ(:,1)),mean(coeff_fmri_NC(:,1)));
        [h,pval]=ttest2(coeff_snp_NC(:,1),coeff_snp_SZ(:,1));p2(1)=pval;
        fprintf('ttest2 for %d-th component of snp:h=%d,p=%f,mean_SZ=%f,mean_NC=%f\n',i,h,p2(1),mean(coeff_snp_SZ(:,1)),mean(coeff_snp_NC(:,1)));
    %end
end

