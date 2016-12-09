%% this script is used to parse the raw file, find its
%% matched entries with csv_entries_nback, and save them as 
%% a cell matrix "raw_entries" into raw_data.mat. csv_data.mat
%% is also been overwritten to keep only those matched entries

load('csv_data.mat');
filename='~/fmri_snp/libd_round123456.qced.phased.imputed.20140915.info0.1.call0.8.no_dup.geno0.02.pgs2.p1e-04.raw';
fid=fopen(filename);
header=fgetl(fid);
header=strsplit(header);
% columns
length(header)
% rows
lc=0;
while 1
	ln=fgetl(fid);
	if ~ischar(ln), break;end;
	lc=lc+1;
end
% allocate space
raw_entries=cell(size(csv_entries_nback,1),length(header));
% reset the cursor
fseek(fid,0,'bof');
% IID (individual ids) in csv
iid=csv_entries_nback(:,2); % 2nd col
% count the matched subjects
count=0;
% the indices of matched subjects in csv
indices=[];
for i=1:lc
	row=strsplit(fgetl(fid));
    % filter out IID of all SBs at the same time
	%idx=find((strcmpi(row{2},strtrim(iid)) & ~strcmpi('SB',strtrim(csv_entries_nback(:,3)))));
	idx=find(strcmpi(row{2},strtrim(iid)));
	assert(length(idx)<2,'duplicate entries in csv\n');
	if length(idx)==1
		raw_entries(idx,:)=row;
		indices=[indices;idx];
		count=count+1;
	end
end
fclose(fid);
indices=sort(indices);
fprintf('%f of csv_entries_nback have their SNP information.\n',(length(indices)-1)/(length(iid)-1));
% filter out mismatched subjects
csv_entries_nback=csv_entries_nback(indices,:);
raw_entries=raw_entries(indices,:);
fprintf('total number of matched entries:%d\n',count);
fprintf('starting converting...\n');
raw_entries(2:end,:)=num2cell(round(str2double(raw_entries(2:end,:))));
fprintf('starting saving...\n');
%save('csv_data.mat','csv_entries_nback');
save('raw_data.mat','raw_entries','-v7.3');
