function write_components(type_fmri, type_mask, num_comps, U, V, best_lambda_U, best_lambda_V, best_tau_U, best_tau_V, ROI)
% load header of each fmri file for writing out
load('csv_data.mat');
if strcmpi(type_fmri,'nback')
    path=strtrim(csv_entries_nback{2,end});
elseif strcmpi(type_fmri,'flanker')
    path=strtrim(csv_entries_flanker{2,end});
elseif strcmpi(type_fmri,'mtl')
    path=strtrim(csv_entries_mtl{2,end});
end
load(path);fname=SPM.xCon(4).Vspm.fname;fmri_file=strrep(path,'SPM.mat',fname);VI = spm_vol(fmri_file); % header of the fmri file

if strcmpi(type_mask,'none')
    nROIs = max(max(max(ROI)));
    U_ROI=zeros(nROIs,size(U,2));
    ROI_vec=ROI(:);
    for i=1:nROIs
        region_i=(ROI_vec==i);
        tmp=U(region_i,:);
        U_ROI(i,:)=mean(tmp);
    end
    U=U_ROI;
end
%}
if strcmpi(type_mask,'ROIs_all') || strcmpi(type_mask,'ROIs_cerebrum') || strcmpi(type_mask,'none') || strcmpi(type_mask,'grey_matter') || strcmpi(type_mask,'ROIs_cerebrum_grey') || strcmpi(type_mask,'ROIs_cerebrum_act')
    nROIs = max(max(max(ROI)));
    nROIs = max(max(max(ROI)));
    if strcmpi(type_mask,'ROIs_all') || strcmpi(type_mask,'ROIs_cerebrum')
        %read ROI mapping from AAL files 
        fname='/home/yiming.wang/fmri_snp/aal.txt';
        fid=fopen(fname);
        aal=textscan(fid,'%d %s %d');
        fclose(fid);
        aal_mapping=cell(nROIs,1);
        for i=1:nROIs
            aal_mapping{i}=aal{2}{i};
        end
    end

    for i=1:num_comps
        % output region weights for each component
        vis_CCA_brain=zeros(size(ROI));
    
        if strcmpi(type_mask,'ROIs_all') || strcmpi(type_mask,'ROIs_cerebrum')
            fname=['/home/yiming.wang/fmri_snp/results_folder_',type_fmri,'_raw9/','vis_cca_brain',num2str(i,'%02d'),num2str(best_lambda_U(i),'%1.2f'),'_',num2str(best_lambda_V(i),'%1.2f'),'_aal.txt'];
            fid=fopen(fname,'w');
            for j=1:nROIs
                if U(j,i) ~= 0
                    vis_CCA_brain(ROI==j)=U(j,i);
                    fprintf(fid,'%.4g %d %s\n',U(j,i),j,aal_mapping{j});
                end
            end
            fclose(fid);strcmpi(type_mask,'ROIs_cerebrum')
        elseif strcmpi(type_mask,'none') || strcmpi(type_mask,'grey_matter') || strcmpi(type_mask,'ROIs_cerebrum_grey') || strcmpi(type_mask,'ROIs_cerebrum_act')
            % in this case ROI is actually a mapping from the pos in U to the pos in brain
            for j=1:nROIs
                if U(j,i) ~= 0
                    vis_CCA_brain(ROI==j)=U(j,i);
                end
            end
            max_U=max(abs(U(:,i)));

            %read ROI mask
            ROI_file='raal_MNI_V4.img';
            AAL_ROI_VI=spm_vol(ROI_file);
            ROI_mask=spm_read_vols(AAL_ROI_VI);
            nROIs_mask=max(max(max(ROI_mask)));

            %read ROI mapping from AAL files 
            fname='/home/yiming.wang/fmri_snp/aal.txt';
            fid=fopen(fname);
            aal=textscan(fid,'%d %s %d');
            fclose(fid);
            aal_mapping=cell(nROIs_mask,1);
            for j=1:nROIs_mask
                aal_mapping{j}=aal{2}{j};
            end

            fname=['/home/yiming.wang/fmri_snp/results_folder_',type_fmri,'_raw9/','vis_cca_brain',num2str(i,'%02d'),num2str(best_lambda_U(i),'%1.4f'),'_',num2str(best_lambda_V(i),'%1.4f'),'_aal.txt'];
            fid=fopen(fname,'w');
            weights_ROI=zeros(nROIs_mask,1);
            for j=1:nROIs_mask
                %abc=(ROI_mask==j & ROI);
                %tot_vox=sum(sum(sum(abc)));
                %weights_ROI(j)=sum(abs(vis_CCA_brain(ROI_mask==j)))/(tot_vox*max_U);
                weights_ROI(j)=sum(abs(vis_CCA_brain(ROI_mask==j)))/(sum(sum(sum(ROI_mask==j)))*max_U);
                fprintf(fid,'%.10f %d %s\n',weights_ROI(j),j,aal_mapping{j});
            end
            fclose(fid);
            fname=['/home/yiming.wang/fmri_snp/results_folder_',type_fmri,'_raw9/','vis_cca_brain',num2str(i,'%02d'),num2str(best_lambda_U(i),'%1.4f'),'_',num2str(best_lambda_V(i),'%1.4f'),'_aal_pos.txt'];
            fid=fopen(fname,'w');
            weights_ROI=zeros(nROIs_mask,1);
            for j=1:nROIs_mask
                weights_ROI(j)=0.5*sum(abs(vis_CCA_brain(ROI_mask==j))+vis_CCA_brain(ROI_mask==j))/(sum(sum(sum(ROI_mask==j)))*max_U);
                fprintf(fid,'%.10f %d %s\n',weights_ROI(j),j,aal_mapping{j});
            end
            fclose(fid);
            fname=['/home/yiming.wang/fmri_snp/results_folder_',type_fmri,'_raw9/','vis_cca_brain',num2str(i,'%02d'),num2str(best_lambda_U(i),'%1.4f'),'_',num2str(best_lambda_V(i),'%1.4f'),'_aal_neg.txt'];
            fid=fopen(fname,'w');
            weights_ROI=zeros(nROIs_mask,1);
            for j=1:nROIs_mask
                weights_ROI(j)=0.5*sum(abs(vis_CCA_brain(ROI_mask==j))-vis_CCA_brain(ROI_mask==j))/(sum(sum(sum(ROI_mask==j)))*max_U);
                fprintf(fid,'%.10f %d %s\n',weights_ROI(j),j,aal_mapping{j});
            end
            fclose(fid);
            fname=['/home/yiming.wang/fmri_snp/results_folder_',type_fmri,'_raw9/','vis_cca_brain','_ROI_aal.txt'];
            fid=fopen(fname,'w');
            weights_ROI=zeros(nROIs_mask,1);
            for j=1:nROIs_mask
                weights_ROI(j)=sum(ROI(ROI_mask==j)~=0)/sum(sum(sum(ROI_mask==j)));
                fprintf(fid,'%.10f %d %s\n',weights_ROI(j),j,aal_mapping{j});
            end
            fclose(fid);
        end

        % output fmri weights on brain volumn for each component
        fname=['/home/yiming.wang/fmri_snp/results_folder_',type_fmri,'_raw9/','vis_cca_brain',num2str(i,'%02d'),num2str(best_lambda_U(i),'%1.4f'),'_',num2str(best_lambda_V(i),'%1.4f'),'.nii'];
        VI.fname=fname;VI.private.dat.fname=VI.fname;
        spm_write_vol(VI,vis_CCA_brain.*2000);

        % output snp weights for each component
        fname=['/home/yiming.wang/fmri_snp/results_folder_',type_fmri,'_raw9/','snp_vector',num2str(i,'%02d'),sprintf('_%1.4f_%1.4f.txt',best_lambda_U(i),best_lambda_V(i))];
        fid=fopen(fname,'w');
        for j=1:size(V,1)
            fprintf(fid,'%f\n',V(j,i));
        end
        fclose(fid);
    end
end
end
