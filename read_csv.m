%% this script is used to parse the csv file and store
%% it as a cell matrix "csv_entries" into csv_data.mat

filename='nb02.sz_sb_nc.rhand.20140601.cr0.5.csv';
fid = fopen(filename);
data = fread(fid, '*char')';
fclose(fid);
csv_entries_nback = regexp(data, ',|\n', 'split');
csv_entries_nback=csv_entries_nback(1:end-1);
csv_entries_nback=reshape(csv_entries_nback',14,[])';
%save('csv_data.mat','csv_entries_nback');

filename='flanker.sz_sb_nc.rhand.20140601.csv';
fid = fopen(filename);
data = fread(fid, '*char')';
fclose(fid);
csv_entries_flanker = regexp(data, ',|\n', 'split');
csv_entries_flanker=csv_entries_flanker(1:end-1);
csv_entries_flanker=reshape(csv_entries_flanker',17,[])';%save('csv_flanker.mat','csv_entries_flanker');

filename='mtl_enc.sz_sb_nc.rhand0.20141208.csv';
fid = fopen(filename);
data = fread(fid, '*char')';
fclose(fid);
csv_entries_mtl = regexp(data, ',|\n', 'split');
csv_entries_mtl=csv_entries_mtl(1:end-1);
csv_entries_mtl=reshape(csv_entries_mtl',14,[])';%save('csv_mtl.mat','csv_entries_mtl');

iid_nback=strtrim(csv_entries_nback(2:end,2)); % 2nd col
iid_flanker=strtrim(csv_entries_flanker(2:end,2));
iid_mtl=strtrim(csv_entries_mtl(2:end,2));

%% get some stats
fprintf('length of iid_nback=%d\n',length(iid_nback)); %548
fprintf('length of iid_flanker=%d\n',length(iid_flanker)); %361
fprintf('length of iid_mtl=%d\n',length(iid_mtl)); %433
fprintf('length of intersect of iid_nback & iid_flanker=%d\n',length(intersect(iid_flanker,iid_nback))); %277
fprintf('length of intersect of iid_nback & iid_mtl=%d\n',length(intersect(iid_mtl,iid_nback))); %250
fprintf('length of intersect of iid_flanker & iid_mtl=%d\n',length(intersect(iid_mtl,iid_flanker))); %205
fprintf('length of intersect of iid_nback & iid_flanker & iid_mtl=%d\n',length(intersect(iid_nback,intersect(iid_mtl,iid_flanker)))); %167
%inters=intersect(iid_nback,intersect(iid_mtl,iid_flanker)); %SZ:19 SB:27 NC:121
%inters=intersect(iid_flanker,iid_nback); %SZ:33 SB:43 NC:201
%inters=intersect(iid_mtl,iid_nback); %SZ:33 SB:42 NC:175
inters=intersect(iid_mtl,iid_flanker); %SZ:33 SB:31 NC:141
SZ_count=0;NC_count=0;SB_count=0;
for i=1:length(inters)
    idx=find(strcmpi(inters{i},iid_mtl))+1;
    assert(length(idx)<2,'duplicate entries in csv.');
    if length(idx)==1
        if strcmpi(csv_entries_mtl(idx,3),'SZ')
            SZ_count=SZ_count+1;
        elseif strcmpi(csv_entries_mtl(idx,3),'NC')
            NC_count=NC_count+1;
        elseif strcmpi(csv_entries_mtl(idx,3),'SB')
            SB_count=SB_count+1;
        else fprintf('unknown sub %s\n',csv_entries_mtl{idx,3});
        end
    end
end

%% optionally exclude SB
exclude_SB=0;
if exclude_SB
    for i=1:length(iid_nback)
        if strcmpi(csv_entries_nback(i+1,3),'SB')
            iid_nback(i)=[];
        end
    end
    for i=1:length(iid_flanker)
        if strcmpi(csv_entries_flanker(i+1,3),'SB')
            iid_flanker(i)=[];
        end
    end
    for i=1:length(iid_mtl)
        if strcmpi(csv_entries_mtl(i+1,3),'SB')
            iid_mtl(i)=[];
        end
    end
end

%% calculate the group size for the union of all the subjects
union_subjects=union(iid_nback,union(iid_mtl,iid_flanker)); %SZ:148 SB:135 NC:494
SZ_count=0;NC_count=0;SB_count=0;
for i=1:length(union_subjects)
    idx=find(strcmpi(union_subjects{i},strtrim(csv_entries_nback(:,2))));
    assert(length(idx)<2,'duplicate entries in csv\n');
    if length(idx)==1
        if strcmpi(csv_entries_nback(idx,3),'SZ')
            SZ_count=SZ_count+1;
        elseif strcmpi(csv_entries_nback(idx,3),'NC')
            NC_count=NC_count+1;
        elseif strcmpi(csv_entries_nback(idx,3),'SB')
            SB_count=SB_count+1;
        else fprintf('unknown sub %s\n',csv_entries_nback{idx,3});
        end
        continue;
    end
    idx=find(strcmpi(union_subjects{i},strtrim(csv_entries_flanker(:,2))));
    assert(length(idx)<2,'duplicate entries in csv\n');
    if length(idx)==1
        if strcmpi(csv_entries_flanker(idx,3),'SZ')
            SZ_count=SZ_count+1;
        elseif strcmpi(csv_entries_flanker(idx,3),'NC')
            NC_count=NC_count+1;
        elseif strcmpi(csv_entries_flanker(idx,3),'SB')
            SB_count=SB_count+1;
        else fprintf('unknown sub %s\n',csv_entries_flanker{idx,3});
        end
        continue;
    end
    idx=find(strcmpi(union_subjects{i},strtrim(csv_entries_mtl(:,2))));
    assert(length(idx)<2,'duplicate entries in csv\n');
    if length(idx)==1
        if strcmpi(csv_entries_mtl(idx,3),'SZ')
            SZ_count=SZ_count+1;
        elseif strcmpi(csv_entries_mtl(idx,3),'NC')
            NC_count=NC_count+1;
        elseif strcmpi(csv_entries_mtl(idx,3),'SB')
            SB_count=SB_count+1;
        else fprintf('unknown sub %s\n',csv_entries_mtl{idx,3});
        end
        continue;
    end
end

%% compute the subject occurence matrix for each experiment
K_nback=zeros([length(union_subjects),length(union_subjects)]);K_flanker=zeros([length(union_subjects),length(union_subjects)]);K_mtl=zeros([length(union_subjects),length(union_subjects)]);
for i=1:length(union_subjects)
    if length(find(strcmpi(union_subjects{i},iid_nback)))==1
        K_nback(i,i)=1;
    end
    if length(find(strcmpi(union_subjects{i},iid_flanker)))==1
        K_flanker(i,i)=1;
    end
    if length(find(strcmpi(union_subjects{i},iid_mtl)))==1
        K_mtl(i,i)=1;
    end
end
assert(sum(diag(K_nback))==length(iid_nback) && sum(diag(K_flanker))==length(iid_flanker) && sum(diag(K_mtl))==length(iid_mtl),'occurence matrix wrong.');
save('csv_data.mat','csv_entries_nback','K_nback','csv_entries_flanker','K_flanker','csv_entries_mtl','K_mtl','union_subjects');

