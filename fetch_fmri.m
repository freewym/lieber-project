%% this script is used to fetch all fmri data of subjects specified
%% in csv_entries, save it as a subjects by voxels matrix "fmri" into 
%% fmri_data.mat. 
%% SPM8 toolkit should be installed in the home directory first.
addpath('~/spm8')
load('csv_data.mat');
%{
paths=csv_entries_nback(2:end,end); %% the last col
for i=1:length(paths)
    load(paths{i}); % SPM
    % data corresponding to nback2-nback0
    fname=SPM.xCon(4).Vspm.fname;
    fmri_file=strrep(paths{i},'SPM.mat',fname);
    VI = spm_vol(fmri_file); % header of the fmri file
    if i==1
        % allocate space
        fmri_nback=zeros(size(paths,1),prod(VI.dim));
    end
    img = spm_read_vols(VI);
    fmri_nback(i,:) = img(:)';
end
save('fmri_data.mat','fmri_nback','-v7.3');
%}
fmri_nback=[];fmri_flanker=[];fmri_mtl=[];
for i=1:length(union_subjects)
    idx=find(strcmpi(union_subjects{i},strtrim(csv_entries_nback(:,2))));
    if length(idx)==1
        path=strtrim(csv_entries_nback{idx,end});
        if exist(path,'file') ~= 2
            fprintf('A nback fMRI file(idx=%d) does not exist.\n',idx);
	else
        load(path);
        fname=SPM.xCon(4).Vspm.fname;
        fmri_file=strrep(path,'SPM.mat',fname);
        VI = spm_vol(fmri_file); % header of the fmri file
        if isempty(fmri_nback)
            fmri_nback=zeros(length(union_subjects),prod(VI.dim));
	    fprintf('fMRI filename of nback is %s. dim=%d x %d x %d\n',fname,VI.dim(1),VI.dim(2),VI.dim(3));
        end
        img = spm_read_vols(VI);
        fmri_nback(i,:) = img(:)';
	end
    end
    idx=find(strcmpi(union_subjects{i},strtrim(csv_entries_flanker(:,2))));
    if length(idx)==1
        path=strtrim(csv_entries_flanker{idx,end});
        if exist(path,'file') ~= 2
            fprintf('A flanker fMRI file(idx=%d) does not exist.\n',idx);
	    return;
	else
        load(path);
        fname=SPM.xCon(10).Vspm.fname; %nogo
        fmri_file=strrep(path,'SPM.mat',fname);
        VI = spm_vol(fmri_file); % header of the fmri file
        if isempty(fmri_flanker)
            fmri_flanker=zeros(length(union_subjects),prod(VI.dim));
	    fprintf('fMRI filename of flanker is %s. dim=%d x %d x %d\n',fname,VI.dim(1),VI.dim(2),VI.dim(3));
        end
        img = spm_read_vols(VI);
        fmri_flanker(i,:) = img(:)';
	end
    end
    idx=find(strcmpi(union_subjects{i},strtrim(csv_entries_mtl(:,2))));
    if length(idx)==1
        path=strtrim(csv_entries_mtl{idx,end});
        if exist(path,'file') ~= 2
            fprintf('An mtl fMRI file(idx=%d) does not exist.\n',idx);
	else
        load(path);
        fname=SPM.xCon(51).Vspm.fname;%%%%% neural encoding
        fmri_file=strrep(path,'SPM.mat',fname);
        VI = spm_vol(fmri_file); % header of the fmri file
        if isempty(fmri_mtl)
            fmri_mtl=zeros(length(union_subjects),prod(VI.dim));
	    fprintf('fMRI filename of mtl is %s. dim=%d x %d x %d\n',fname,VI.dim(1),VI.dim(2),VI.dim(3));
        end
        img = spm_read_vols(VI);
        fmri_mtl(i,:) = img(:)';
        end
    end
end
save('fmri_data.mat','fmri_nback','fmri_flanker','fmri_mtl','-v7.3');
