function [data_fmri,data_snp,dx_csv,csv_entries,num_NC,num_SB,num_SZ]=extract_data(type_fmri)
load('fmri_data.mat');
load('csv_data.mat');
load('raw_data.mat');

% extract labels of the same order as in union_subjects
dx_csv=cell(length(union_subjects),1);
indice_in_nback=zeros(length(union_subjects),1);
indice_in_flanker=zeros(length(union_subjects),1);
indice_in_mtl=zeros(length(union_subjects),1);
for i=1:length(union_subjects)
    idx_in_nback=find(strcmpi(union_subjects{i},strtrim(csv_entries_nback(:,2))));
    idx_in_flanker=find(strcmpi(union_subjects{i},strtrim(csv_entries_flanker(:,2))));
    idx_in_mtl=find(strcmpi(union_subjects{i},strtrim(csv_entries_mtl(:,2))));
    % sanity check
    if ~isempty(idx_in_nback) && ~isempty(idx_in_flanker)
        assert(strcmpi(strtrim(csv_entries_nback{idx_in_nback,3}),strtrim(csv_entries_flanker{idx_in_flanker,3})),'label between nback and flanker not consistent');
    end
    if ~isempty(idx_in_nback) && ~isempty(idx_in_mtl)
        assert(strcmpi(strtrim(csv_entries_nback{idx_in_nback,3}),strtrim(csv_entries_mtl{idx_in_mtl,3})),'label between nback and mtl not consistent');
    end
    if ~isempty(idx_in_flanker) && ~isempty(idx_in_mtl)
        assert(strcmpi(strtrim(csv_entries_flanker{idx_in_flanker,3}),strtrim(csv_entries_mtl{idx_in_mtl,3})),'label between flanker and mtl not consistent');
    end
    % store label into dx_csv
    if ~isempty(idx_in_nback)
        dx_csv{i}=strtrim(csv_entries_nback{idx_in_nback,3});
    elseif ~isempty(idx_in_flanker)
        dx_csv{i}=strtrim(csv_entries_flanker{idx_in_flanker,3});
    else
        dx_csv{i}=strtrim(csv_entries_mtl{idx_in_mtl,3});
    end
    % store index of i-th entry in each fmri , 0 if not appear 
    if ~isempty(idx_in_nback)
        indice_in_nback(i)=idx_in_nback;
    else indice_in_nback(i)=0;
    end
    if ~isempty(idx_in_flanker)
        indice_in_flanker(i)=idx_in_flanker;
    else indice_in_flanker(i)=0;
    end
    if ~isempty(idx_in_mtl)
        indice_in_mtl(i)=idx_in_mtl;
    else indice_in_mtl(i)=0;
    end
end
[dx_csv,I]=sort(dx_csv);% order: NC<SB<SZ
fmri_nback=fmri_nback(I,:);fmri_flanker=fmri_flanker(I,:);fmri_mtl=fmri_mtl(I,:);
K_nback=diag(logical(K_nback));K_nback=K_nback(I);
K_flanker=diag(logical(K_flanker));K_flanker=K_flanker(I);
K_mtl=diag(logical(K_mtl));K_mtl=K_mtl(I);
fmri_nback=fmri_nback(K_nback,:);
fmri_flanker=fmri_flanker(K_flanker,:);
fmri_mtl=fmri_mtl(K_mtl,:);
indice_in_nback=indice_in_nback(I);indice_in_nback(indice_in_nback==0)=[];
indice_in_flanker=indice_in_flanker(I);indice_in_flanker(indice_in_flanker==0)=[];
indice_in_mtl=indice_in_mtl(I);indice_in_mtl(indice_in_mtl==0)=[];
csv_entries_nback=csv_entries_nback([1;indice_in_nback],:);
csv_entries_flanker=csv_entries_flanker([1;indice_in_flanker],:);
csv_entries_mtl=csv_entries_mtl([1;indice_in_mtl],:);
if strcmpi(type_fmri,'nback')
    dx_csv=dx_csv(K_nback==1,:);
    iid=strtrim(csv_entries_nback(:,2)); % 2nd col
elseif strcmpi(type_fmri,'flanker')
    dx_csv=dx_csv(K_flanker==1,:);
    iid=strtrim(csv_entries_flanker(:,2)); % 2nd col
elseif strcmpi(type_fmri,'mtl')
    dx_csv=dx_csv(K_mtl==1,:);
    iid=strtrim(csv_entries_mtl(:,2)); % 2nd col
end
% retain those iids that also appear in raw_entries
iid_raw=raw_entries(:,2);
for i=2:length(iid_raw)
    iid_raw{i}=int2str(iid_raw{i});
end
to_be_retained=[];
for i=1:length(iid)
    if length(find(strcmpi(iid{i},iid_raw)))~=0
        to_be_retained=[to_be_retained i];
    end
end

% for each entry in csv_entries, find its index in iid_raw
indices=[];
if strcmpi(type_fmri,'nback')
    (length(to_be_retained)-1)/(size(csv_entries_nback,1)-1) %%%%%%%%%%%%
    csv_entries_nback=csv_entries_nback(to_be_retained,:);
    dx_csv=dx_csv(to_be_retained(1,2:end)-1,:);
    fmri_nback=fmri_nback(to_be_retained(1,2:end)-1,:);
    for i=1:size(csv_entries_nback,1)
        idx=find(strcmpi(strtrim(csv_entries_nback(i,2)),iid_raw));
        assert(length(idx)==1,'some entries in csv_entries_nback couldnt been found in raw');
        indices=[indices;idx];
    end
    data_fmri=fmri_nback;
    csv_entries=csv_entries_nback;
elseif strcmpi(type_fmri,'flanker')
    (length(to_be_retained)-1)/(size(csv_entries_flanker,1)-1) %%%%%%%%%%%%
    csv_entries_flanker=csv_entries_flanker(to_be_retained,:);
    dx_csv=dx_csv(to_be_retained(1,2:end)-1,:);
    fmri_flanker=fmri_flanker(to_be_retained(1,2:end)-1,:);
    for i=1:size(csv_entries_flanker,1)
        idx=find(strcmpi(strtrim(csv_entries_flanker(i,2)),iid_raw));
        assert(length(idx)==1,'some entries in csv_entries_flanker couldnt been found in raw');
        indices=[indices;idx];
    end
    data_fmri=fmri_flanker;
    csv_entries=csv_entries_flanker;
elseif strcmpi(type_fmri,'mtl')
    (length(to_be_retained)-1)/(size(csv_entries_mtl,1)-1) %%%%%%%%%%%%
    csv_entries_mtl=csv_entries_mtl(to_be_retained,:);
    dx_csv=dx_csv(to_be_retained(1,2:end)-1,:);
    fmri_mtl=fmri_mtl(to_be_retained(1,2:end)-1,:);
    for i=1:size(csv_entries_mtl,1)
        idx=find(strcmpi(strtrim(csv_entries_mtl(i,2)),iid_raw));
        assert(length(idx)==1,'some entries in csv_entries_mtl couldnt been found in raw');
        indices=[indices;idx];
    end
    data_fmri=fmri_mtl;
    csv_entries=csv_entries_mtl;
end
assert(to_be_retained(1,1)==1,'wrong');%%%%%%%%%%%%%

%raw_entries(:,933+6)%%%%%%%%%%%%%%%%%

raw_entries=raw_entries(indices,:);
data_snp=cell2mat(raw_entries(2:end,7:end));
% impute NaN with 1 (the middle of the three possible values)
data_snp(isnan(data_snp))=1;

num_NC=length(find(strcmpi('NC',dx_csv)));
num_SB=length(find(strcmpi('SB',dx_csv)));
num_SZ=length(find(strcmpi('SZ',dx_csv)));
fprintf('Total num of NC:%d, total num of SB:%d,total num of SZ:%d\n',num_NC,num_SB,num_SZ);
end

function data=impute_by_average(data)
    avg=zeros(1,size(data,2));
    for i=1:size(data,2)
        s=0;count=0;
        for j=1:size(data,1)
            if ~isnan(data(j,i)),s=s+data(j,i);count=count+1;end
            avg(i)=s/count;
        end
    end
    for i=1:size(data,2)
        for j=1:size(data,1)
            if isnan(data(j,i)), data(j,i)=avg(i);end
        end
    end
end
