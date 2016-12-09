function data_fmri_new = extract_and_use_demong(type_fmri, use_demong, csv_entries, data_fmri)
if strcmpi(use_demong,'no')
    return;
end
%% extract demographics
age_csv=csv_entries(2:end,5);
age=str2double(age_csv);
sex_csv=csv_entries(2:end,6);
sex=str2double(sex_csv);
if strcmpi(type_fmri,'nback') || strcmpi(type_fmri,'flanker')
    edu_csv=csv_entries(2:end,9);
elseif strcmpi(type_fmri,'mtl')
    edu_csv=csv_entries(2:end,8);
end
edu=str2double(edu_csv);
if strcmpi(type_fmri,'nback') || strcmpi(type_fmri,'flanker')
    IQ_csv=csv_entries(2:end,10);
elseif strcmpi(type_fmri,'mtl')
    IQ_csv=csv_entries(2:end,9);
end
IQ=str2double(IQ_csv);
if strcmpi(type_fmri,'nback') || strcmpi(type_fmri,'flanker')
    cr_csv=csv_entries(2:end,12); % crNB2 for nback, cr_nogo for flanker
elseif strcmpi(type_fmri,'mtl')
    cr_csv=[];   % empty for mtl
end
if ~isempty(cr_csv)
    cr=str2double(cr_csv);
else cr=[];
end
if strcmpi(type_fmri,'nback')
    rt_csv=csv_entries(2:end,13); % rtNB2 for nback
    snav_csv=csv_entries(2:end,11); % SNAV for nback
elseif strcmpi(type_fmri,'flanker')
    rt_csv=[]; % empty for flanker nogo
elseif strcmpi(type_fmri,'mtl')
    rt_csv=[]; % empty for mtl
end
if ~isempty(rt_csv)
    rt=str2double(rt_csv);
else rt=[];
end
if ~isempty(snav_csv)
    snav=str2double(snav_csv);
else snav=[];
end
%dmg=[age sex edu IQ cr rt];
dmg=[age sex snav];
dmg=impute_by_average(dmg); % only some IQs are missing

%dmg=zscore(dmg);
data_fmri_new = zeros(size(data_fmri));

% remove all-zero cols, record their locations
zero_cols_idx = sum((data_fmri == 0)) >= 0.99 * size(data_fmri,1);
fprintf('remove %d cols\n',sum(zero_cols_idx));
data_fmri(:,zero_cols_idx)=[];

% regress out to obtain residuals
data_fmri_temp=zeros(size(data_fmri));
part_size=ceil(size(data_fmri,2)/500);
fprintf('part_size=%d\n',part_size);
for i=1:(size(data_fmri,2)/part_size+1)
    [beta,~,data_part,~,~]=mvregress([ones(size(dmg,1),1),dmg],data_fmri(:,((i-1)*part_size+1):min(i*part_size,size(data_fmri,2))));
    data_fmri_temp(:,((i-1)*part_size+1):min(i*part_size,size(data_fmri,2)))=data_part;
end

% add back all zero cols
data_fmri_new(:,~zero_cols_idx) = data_fmri_temp;

%or equivalently create object fitlm
%for i=1:size(data_ROI,2)
%    mdl=fitlm(dmg,data_ROI(:,i));data_fmri=mdl.Residuals{:,'Raw'});beta=mdl.Coefficients{:,'Estimate'});
%end
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
