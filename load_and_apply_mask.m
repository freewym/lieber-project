function [data_fmri, ROI, AAL_ROI_VI] = load_and_apply_mask(type_fmri, type_mask, type_summary, data_fmri)
if strcmpi(type_fmri,'nback') || strcmpi(type_fmri,'mtl')
    ROI_file='raal_MNI_V4.img';
elseif strcmpi(type_fmri,'flanker')
    ROI_file='raal_MNI_V4_flanker.img';
end
AAL_ROI_VI=spm_vol(ROI_file);
ROI=spm_read_vols(AAL_ROI_VI);
nROIs=max(max(max(ROI)));
ROI_ori=ROI;

Grey_file='rgrey_for_nback.nii';
Grey_VI=spm_vol(Grey_file);
if strcmpi(type_mask, 'ROIs_cerebrum_grey')
  ROI=spm_read_vols(Grey_VI);
end

act_file='nb02_wb_p05.nii';
act_VI=spm_vol(act_file);
if strcmpi(type_mask, 'ROIs_cerebrum_act')
  ROI=spm_read_vols(act_VI);
end

if strcmpi(type_mask,'none')
    % only exclude margins
    for i=1:nROIs
        if i>1
            ROI(ROI==i)=1;
        end
    end
    nROIs=sum(sum(sum(ROI)));n=1;
    for k=1:size(ROI,3)
        for j=1:size(ROI,2)
            for i=1:size(ROI,1)
                if ROI(i,j,k) ~= 0
                    ROI(i,j,k) = n;
                    n=n+1;
                end
            end
        end
    end
end

if strcmpi(type_mask,'grey_matter')
    ROI(ROI>=0.20)=1;
    ROI(ROI<0.20)=0;
    nROIs=sum(sum(sum(ROI)));n=1;
    for k=1:size(ROI,3)
        for j=1:size(ROI,2)
            for i=1:size(ROI,1)
                if ROI(i,j,k) ~= 0
                    ROI(i,j,k) = n;
                    n=n+1;
                end
            end
        end
    end
end

if strcmpi(type_mask,'ROIs_cerebrum')
    % select the first 90 masks corresponding to frontal brain region
    for i=1:nROIs
        if i>90
            ROI(ROI==i)=0;
        end
    end
    nROIs=90;
elseif strcmpi(type_mask, 'ROIs_cerebrum_grey') || strcmpi(type_mask, 'ROIs_cerebrum_act')
    for i=1:nROIs
        if i>90
            ROI_ori(ROI_ori==i)=0;
        end
    end
    ROI(ROI>=0.20)=1;
    ROI(ROI<0.20)=0;
    %ROI= ROI.*(ROI_ori~=0);
    nROIs=sum(sum(sum(ROI)));n=1;
    for k=1:size(ROI,3)
        for j=1:size(ROI,2)
            for i=1:size(ROI,1)
                if ROI(i,j,k) ~= 0
                    ROI(i,j,k) = n;
                    n=n+1;
                end
            end
        end
    end
elseif strcmpi(type_mask,'ROIs_all')
    %do nothing
end
fprintf('There are %d ROIs in total\n',nROIs);

% apply mask
data_ROI=zeros(size(data_fmri,1),nROIs);
if strcmpi(type_summary,'avg')
    for i=1:nROIs
        region_i=(ROI==i);
        region_i=region_i(:)';
        for j=1:size(data_fmri,1)
            data_ROI(j,i)=mean(data_fmri(j,region_i));
        end
    end
    % for group scca
    if strcmpi(type_mask, 'ROIs_cerebrum_grey')
        voxel_membership=zeros(nROIs,1);
        for i=1:90
            pos_in_vec=ROI(ROI_orig==i);
            voxel_membership(pos_in_vec)=i;
        end
        assert(sum(voxel_membership==0)==0,'every element should belong to a group.');
        ROICount=zeros(90,1);
        for i=1:90
            ROICount(i) = sum(voxel_membership==i);
        end
        save('/home/yiming.wang/fmri_snp/voxel_membership.mat', 'voxel_membership');
        save('/home/yiming.wang/fmri_snp/voxel_group_count.mat', 'ROICount');
    end
elseif strcmpi(type_summary,'eig')
    % test the pca variances for each region, and use the 1st PC as the summary of the region
    s_sum=[];
    for i=1:nROIs
        region_i=(ROI==i);
        region_i=region_i(:)';
        fmri_within_region_i=data_fmri(:,region_i);
        [~,score,~,~,explained,~]=pca(fmri_within_region_i,'NumComponents',10);
        data_ROI(:,i)=score(:,1);
        s_sum=[s_sum explained(1:10)];
    end
end
data_fmri=data_ROI; 
end
