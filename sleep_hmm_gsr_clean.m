cd /raid/common/sleep1/derivatives/nils/
addpath(genpath('/misc/imeel/dezwart/matlab/'))
addpath(genpath('/misc/imeel/yangf7/matlab/HMM-MAR-master'))
addpath(genpath('/misc/imeel/yangf7/matlab/BCT'))
addpath(genpath('/misc/imeel/yangf7/matlab/matlab_code_Nils'))


load('/misc/imeel/yangf7/matlab/SW_detection/label_table_0705.mat')
load('/raid/common/sleep1/derivatives/nils/ROI_list.mat');
load('/raid/common/sleep1/derivatives/nils/HMM_results/GSR_mean_variance.mat');
load('/raid/common/sleep1/derivatives/nils/HMM_results/GSR_mean_variance_all.mat');

demo_table=readtable('/raid/common/sleep1/derivatives/nils/HMM_results/output.csv');

demo_table=demo_table(demo_table.include_analysis==1,:);


ROI_list.Network_Label_v2=nan(300,1);
for ii=1:height(ROI_list)
    if ROI_list.Network_Label(ii)<18
        ROI_list.Network_Label_v2(ii)=ROI_list.Network_Label(ii);
    elseif ROI_list.Network_Label(ii)>=18&&ROI_list.Network_Label(ii)<39
        ROI_list.Network_Label_v2(ii)=18;
    elseif ROI_list.Network_Label(ii)>=39 &&ROI_list.Network_Label(ii)<60
        ROI_list.Network_Label_v2(ii)=19;
    end

end


ROI_list.Network_Label_v3=nan(300,1);
for ii=1:height(ROI_list)
    if ROI_list.Network_Label(ii)<18
        ROI_list.Network_Label_v3(ii)=ROI_list.Network_Label(ii);
    elseif ROI_list.Network_Label_v2(ii)==18 && ROI_list.Network_Label(ii)==37
        ROI_list.Network_Label_v3(ii)=17;
    elseif ROI_list.Network_Label_v2(ii)==18 && ROI_list.Anatomical_lable(ii)==3
        ROI_list.Network_Label_v3(ii)=18;
    elseif ROI_list.Network_Label_v2(ii)==18 && (ROI_list.Anatomical_lable(ii)==5||ROI_list.Anatomical_lable(ii)==4)
        ROI_list.Network_Label_v3(ii)=19;
    elseif ROI_list.Network_Label_v2(ii)==18 && ROI_list.Anatomical_lable(ii)==6
        ROI_list.Network_Label_v3(ii)=20;
    elseif ROI_list.Network_Label(ii)>=39 &&ROI_list.Network_Label(ii)<60
        ROI_list.Network_Label_v3(ii)=21;
    end

end







net_n_sep={'DMN','VIS','FPN','REW','DAN','VAN','SAL','CON','dSMN','lSMN','AUD','PMN','MTL',...
    'Sub_DMN','Sub_VIS','Sub_FPN','Sub_REW','Sub_VAN','Sub_SAL','Sub_CON','Sub_dSMN','Sub_lSMN','Sub_MTL',...
    'CerebellarDMN','CerebellarVIS','CerebellarFPN','CerebellarDAN','CerebellarSAL','CerebellarCON','CerebellardSMN','CerebellarlSMN'}';
prepro={
          'procbasic_ric_gsrwb' % with RETROICOR and whole-brain GS
        };
% 
% prepro={
%         'procbasic' %basic preprocessing w/o physio correction
%         'procbasic_ric' % with RETROICOR only
%         'procbasic_ric_gsrwb' % with RETROICOR and whole-brain GS
%         'procbasic_ric_rvt_ppga'};%with RETROICOR +RVT +PPGA






label_table{:,'subs'}=ones(height(label_table),1);
label_table{:,'sess'}=ones(height(label_table),1);
subs=unique(label_table(:,'subject'));
subs=cellstr(subs{:,1});






id=1:12;
for i = 1:height(label_table)
    label_table{i,'subs'}=id(strcmp(label_table{i,'subject'},subs));
    label_table{i,'sess'}=double(strcmp(label_table{i,'session'},'ses-2'))+1;
end





label_table.uni_id=label_table.sess+2*(label_table.subs-1);


labels=unique(cellstr(label_table.label_fmri),'stable');


for ii = 1: height(label_table)% per label_fmri
    try
        label_table{ii,'Dream'}=dream_table{ismember(cellstr(dream_table.RunName),{label_table{ii,'label_fmri'}}),'DreamType'};
    catch
        label_table{ii,'Dream'}=nan;
    end

end




sleep_score=[];
for i=1:length(labels)
temp_table=volume_table_all(ismember(cellstr(volume_table_all.label_fmri),labels{i}),:);

for j=1:height(temp_table)
    sleep_score{i}(temp_table.volume_start(j):temp_table.volume_end(j))=temp_table.Score(j);

end
end


test_mat=cell(height(label_table),length(prepro)  );
FC_mat=test_mat;
sleep_score_censor=sleep_score;
sleep_score{22}(755)=6;%missing one TR of sleep score
test_mat=[];
for i = 1: height(label_table)% per label_fmri
    for j=1    % per preprocessing type
%         try
%
%         catch
%             continue
%         end

        if label_table{i,'diff_count'}==1|label_table{i,'EEG_missing'}==1
            continue
        end


          censor_data=readtable(fullfile('/raid/common/projects/OLD/sleep1/dante/1601_Ap5a6_AlwReg/',label_table{i,'subject'},...
            label_table{i,'session'},['run-',label_table{i,'label_fmri'}],[label_table{i,'label_fmri'},'.results'],...
            ['censor_',label_table{i,'label_fmri'},'_combined_2.1D']),'FileType','text');
            censor_array=censor_data{4:end,2};%first 3 volumes were removed
            censor_ind=1:length(censor_array);
            censor_ind=censor_ind(censor_array>0);

            
         test_mat{i,1}=readtable(['/raid/common/sleep1/derivatives/nils/dfc_results/output_', label_table{i,'label_fmri'},'_',char(prepro{j}),'_300_000.netts'],'FileType','text');
         test_mat{i,1}=table2array(test_mat{i,1}(1:300,censor_ind));%with censor
         sleep_score_censor{i}=sleep_score{i}(censor_ind+3);% first 3 volumes were removed


    end
end








%% slow-wave amplitute time course
subjects = {'00003','00003','00006','00006','00014','00014','00029','00029','00047','00047','00056','00056','00060','00060'...
    ,'00063','00063','00065','00065','00078','00078','00086','00086','00105','00105'};
sessione = {'1','2','1','2','1','2','1','2','1','2','1','2','1','2','1','2','1','2','1','2','1','2','1','2'};


sw_amp_mat=[];sw_time_mat=[];sw_freq_mat=[];
for i = 1: height(label_table)% per label_fmri
    censor_data=[];
    try
    censor_data=readtable(fullfile('/raid/common/projects/OLD/sleep1/dante/1601_Ap5a6_AlwReg/',label_table{i,'subject'},...
        label_table{i,'session'},['run-',label_table{i,'label_fmri'}],[label_table{i,'label_fmri'},'.results'],...
        ['censor_',label_table{i,'label_fmri'},'_combined_2.1D']),'FileType','text');

%     censor_data=readtable(fullfile('/raid/common/projects/OLD/sleep1/dante/1601_Ap5a6_AlwReg/',label_table_clean{i,'subject'}{:},...
%         label_table_clean{i,'session'}{:},['run-',label_table_clean{i,'label_fmri'}{:}],[label_table_clean{i,'label_fmri'}{:},'.results'],...
%         ['censor_',label_table_clean{i,'label_fmri'}{:},'_combined_2.1D']),'FileType','text');
% 
    catch
        continue
    end
    censor_array=censor_data{4:end,2};%first 3 volumes were removed
    censor_ind=1:length(censor_array);
    censor_ind=censor_ind(censor_array>0);
  
    sleepscore_table=readtable(fullfile('/raid/common/sleep1/derivatives/sleep_scoring_mr/',label_table{i,'subject'},...
        label_table{i,'session'},[label_table{i,'subject'},'_',label_table{i,'session'},'_task-sleep_run-',...
        label_table{i,'label_fmri'},'.txt']));
    fmri_start=(find(sleepscore_table{:,4}>0,1)-sleepscore_table{find(sleepscore_table{:,4}>0,1),4})*30*250;
    fmri_start=fmri_start+3*3*250;%first 3 volumes discarded.
    fmri_end=(find(sleepscore_table{:,4}>0,1,'last')-1+sleepscore_table{find(sleepscore_table{:,4}>0,1,'last'),4})*30*250;
    load(fullfile('/misc/imeel/yangf7/matlab/SW_detection/derivatives_icarej025mod_altbadch',label_table{i,'subject'},...
        ['report_',label_table{i,'session'}],[label_table{i,'subject'},'_',label_table{i,'session'},'_clean_swa_results_thr0.mat']),'swa_results')
    load(fullfile('/misc/imeel/yangf7/matlab/SW_detection/derivatives_icarej025mod_altbadch',label_table{i,'subject'},...
        ['report_',label_table{i,'session'}],[label_table{i,'subject'},'_',label_table{i,'session'},'.mat']),'N_EPOCHE')
    sw_amp=nan(1,length(censor_array));sw_time=nan(1,length(censor_array));sw_freq=nan(1,length(censor_array));
    if label_table{i,'EEG_sess'}>1
    fmri_start=fmri_start+sum(N_EPOCHE(1:label_table{i,'EEG_sess'}-1))*250;
    end
    for j=1:length(censor_array)
        t1=fmri_start+(j-1)*250*3;
        t2=t1+250*3;
        ind=double(swa_results.maxnegpk)>=t1&double(swa_results.maxnegpk)<=t2;
        if sum(ind)==0
            sw_amp(j)=nan;
        else
            sw_amp(j)=mean(double(swa_results.maxnegpkamp(ind)));
            sw_freq(j)=length(swa_results.maxnegpkamp(ind));
        end
        sw_time(j)=t1;
    end

            sw_amp_mat{i}=sw_amp(censor_ind);
            sw_time_mat{i}=sw_time(censor_ind);
            sw_freq_mat{i}=sw_freq(censor_ind);
        
        



    clearvars N_EPOCHE swa_results
    
    
end
%% PPG analysis Catie matlab code

ppg_mat=[];
for ii=1:height(label_table)
    if label_table{ii,'diff_count'}==1
        continue
    end
    censor_data=readtable(fullfile('/raid/common/projects/OLD/sleep1/dante/1601_Ap5a6_AlwReg/',label_table{ii,'subject'},...
    label_table{ii,'session'},['run-',label_table{ii,'label_fmri'}],[label_table{ii,'label_fmri'},'.results'],...
    ['censor_',label_table{ii,'label_fmri'},'_combined_2.1D']),'FileType','text');
    d=read_data(fullfile('/raid/common/sleep1/derivatives/fmriprepro_025/',[label_table{ii,'subject'}],...
        [label_table{ii,'session'}],[label_table{ii,'subject'},'_',label_table{ii,'session'},'_task-sleep_run-'...
        label_table{ii,'label_fmri'},'_nils','.svd']));
%     ppg_paras=d.(['scan_',label_table{ii,'label_fmri'}]);
%     
%     ppg_paras=ppg_paras(1:height(censor_data),:);
    ppg_mat{ii}=d;
end






fs=1000;
tr=3;
CardStd=[];CardMean=[];power_idx=[];

for ii=1:height(label_table)
    if isempty(ppg_mat{ii})||~isstruct(ppg_mat{ii})
        continue
    end
    
    d=ppg_mat{ii};
    d.ibi(end+1)=nan;
    d.pks=filloutliers(d.pks,"clip","movmedian",[30 30]);
    d.ibi=filloutliers(d.ibi,"clip","movmedian",[30 30]);
    censor_data=readtable(fullfile('/raid/common/projects/OLD/sleep1/dante/1601_Ap5a6_AlwReg/',label_table{ii,'subject'},...
    label_table{ii,'session'},['run-',label_table{ii,'label_fmri'}],[label_table{ii,'label_fmri'},'.results'],...
    ['censor_',label_table{ii,'label_fmri'},'_combined_2.1D']),'FileType','text');
    censor_array=censor_data{4:end,2};%first 3 volumes were removed
    censor_ind=1:length(censor_array);
    censor_ind=censor_ind(censor_array>0);
    censor_ind=censor_ind+3;


    for n_tr=1:ceil(height(censor_data)/10)
        indx_b=max(((n_tr-1)*fs*tr*10+1)-0,1);%included 0 tr before current tr;
        indx_f=min((n_tr*fs*tr*10)+0,height(censor_data)*fs*tr*10);%included 0 tr after current tr;
        indx=(d.locs>indx_b&d.locs<indx_f);
        start_indx=(n_tr-1)*10+1;
        end_indx=min((n_tr)*10,height(censor_data));
        CardStd{ii}(start_indx:end_indx,1)=std(d.pks(indx));
        CardStd{ii}(start_indx:end_indx,2)=std(d.ibi(indx),'omitnan');
        CardStd{ii}(start_indx:end_indx,3)=std(60./d.ibi(indx),'omitnan');


    end
    CardStd{ii}=CardStd{ii}(censor_ind,:);


end

%% PPG/RVT analysis hendrik python

ppg_rvt_mat=[];
for ii=1:height(label_table)
    if label_table{ii,'diff_count'}==1
        continue
    end
    CardAmp=[];RespRvt=[];
    censor_data=readtable(fullfile('/raid/common/projects/OLD/sleep1/dante/1601_Ap5a6_AlwReg/',label_table{ii,'subject'},...
    label_table{ii,'session'},['run-',label_table{ii,'label_fmri'}],[label_table{ii,'label_fmri'},'.results'],...
    ['censor_',label_table{ii,'label_fmri'},'_combined_2.1D']),'FileType','text');
    load(fullfile('/raid/common/sleep1/derivatives/1602_Ap5a6_Ric2AlwReg/',[label_table{ii,'subject'}],...
        [label_table{ii,'session'}],['run-',label_table{ii,'label_fmri'}],[...
        label_table{ii,'label_fmri'},'_phy','.mat']),'CardAmp','RespRvt');
%     ppg_paras=d.(['scan_',label_table{ii,'label_fmri'}]);
%     
%     ppg_paras=ppg_paras(1:height(censor_data),:);

    ppg_rvt_mat{ii}=[CardAmp;RespRvt];
end






fs=20;
tr=3;
ppg_rvt_Std=[];
for ii=1:height(label_table)
    if isempty(ppg_rvt_mat{ii})
        continue
    end
    % 
    % d=ppg_mat{ii};
    % d.ibi(end+1)=nan;
    % d.pks=filloutliers(d.pks,"clip","movmedian",[30 30]);
    % d.ibi=filloutliers(d.ibi,"clip","movmedian",[30 30]);
    censor_data=readtable(fullfile('/raid/common/projects/OLD/sleep1/dante/1601_Ap5a6_AlwReg/',label_table{ii,'subject'},...
    label_table{ii,'session'},['run-',label_table{ii,'label_fmri'}],[label_table{ii,'label_fmri'},'.results'],...
    ['censor_',label_table{ii,'label_fmri'},'_combined_2.1D']),'FileType','text');
    censor_array=censor_data{4:end,2};%first 3 volumes were removed
    censor_ind=1:length(censor_array);
    censor_ind=censor_ind(censor_array>0);
    censor_ind=censor_ind+3;


    for n_tr=1:ceil(height(censor_data)/10)
        indx_b=max(((n_tr-1)*fs*tr*10+1)-0,1);%included 0 tr before current tr;
        indx_f=min((n_tr*fs*tr*10)+0,size(ppg_rvt_mat{ii},2));%included 0 tr after current tr;
        % indx=(d.locs>indx_b&d.locs<indx_f);
        start_indx=(n_tr-1)*10+1;
        end_indx=min((n_tr)*10,height(censor_data));
        ppg_rvt_Std{ii}(start_indx:end_indx,1)=std(ppg_rvt_mat{ii}(1,indx_b:indx_f),'omitnan');
        ppg_rvt_Std{ii}(start_indx:end_indx,2)=std(ppg_rvt_mat{ii}(2,indx_b:indx_f),'omitnan');
        % CardStd{ii}(start_indx:end_indx,3)=mean(60./d.ibi(indx),'omitnan');


    end
    ppg_rvt_Std{ii}=ppg_rvt_Std{ii}(censor_ind,:);


end

%% load eog data


load('/raid/common/sleep1/derivatives/nils/eog/eye_eeg.mat')





%% hmm raw fc pca prep


dfc_mat_sep_clean=test_mat(label_table.diff_count~=1&label_table.EEG_missing~=1,:);
label_table_clean=label_table(label_table.diff_count~=1&label_table.EEG_missing~=1,:);

sleep_score_mat_clean=sleep_score_censor(:,label_table.diff_count~=1&label_table.EEG_missing~=1&label_table.sess==2)';%night2 only

sleep_score_night2=cat(2,sleep_score_mat_clean{:})';

sleep_score_mat_clean1=sleep_score_censor(:,label_table.diff_count~=1&label_table.EEG_missing~=1&label_table.sess==1)';%night1 only

sleep_score_night1=cat(2,sleep_score_mat_clean1{:})';


CardStd_night=cat(2,CardStd(:,label_table.diff_count~=1&label_table.EEG_missing~=1)');
CardStd_night{39}=nan(548,3);
% Card_night2=cell2mat(CardStd_night2);
% CardStd_night1=cat(2,CardStd(:,label_table.diff_count~=1&label_table.EEG_missing~=1&label_table.sess==1)');
% Card_night1=cell2mat(CardStd_night1);
ppg_rvt_night=cat(2,ppg_rvt_Std(:,label_table.diff_count~=1&label_table.EEG_missing~=1)');
eog_night=cat(2,eog_std_mat(:,label_table.diff_count~=1&label_table.EEG_missing~=1)');



sw_amp_mat_clean=sw_amp_mat(:,label_table.diff_count~=1&label_table.EEG_missing~=1)';
sw_freq_mat_clean=sw_freq_mat(:,label_table.diff_count~=1&label_table.EEG_missing~=1)';




n=1;dfc_sub_mat=[];sw_amp_sub=[];sw_freq_sub=[];Card_night1=[];ppg_rvt_night1=[];ppg_rvt_night2=[];
dfc_sub_mat1=[];sw_amp_sub1=[];sw_freq_sub1=[];Card_night2=[];ppg_night1=[];ppg_night2=[];eog_night1=[];eog_night2=[];
T_runs=[];T_runs1=[];
for i=1:length(unique(label_table_clean{:,'uni_id'}))/2
for n_prepro=1
%     for i = 1:height(lable_table)
        ind= (label_table_clean{:,'uni_id'}==i*2);
        ind1=(label_table_clean{:,'uni_id'}==i*2-1);
        inx=1:172;
        inx1=inx(ind1);
        inx=inx(ind);
       
        dfc_sub_mat{i,n_prepro}=cat(2,dfc_mat_sep_clean{ind,n_prepro});
        dfc_sub_mat1{i,n_prepro}=cat(2,dfc_mat_sep_clean{label_table_clean{:,'uni_id'}==i*2-1,n_prepro});
        sw_amp_sub{i,n_prepro}=cat(2,sw_amp_mat_clean{ind,n_prepro});
        sw_amp_sub1{i,n_prepro}=cat(2,sw_amp_mat_clean{label_table_clean{:,'uni_id'}==i*2-1,n_prepro});

        Card_night1{i,n_prepro}=zscore(cat(1,CardStd_night{label_table_clean{:,'uni_id'}==i*2-1,n_prepro}));
        Card_night2{i,n_prepro}=zscore(cat(1,CardStd_night{ind,n_prepro}));
        
        ppg_rvt_night1{i,n_prepro}=zscore(cat(1,ppg_rvt_night{ind1,n_prepro}));
        ppg_rvt_night2{i,n_prepro}=zscore(cat(1,ppg_rvt_night{ind,n_prepro}));
        
        eog_night1{i,n_prepro}=zscore(cat(1,eog_night{ind1,n_prepro}));
        eog_night2{i,n_prepro}=zscore(cat(1,eog_night{ind,n_prepro}));
        


        sw_freq_sub{i,n_prepro}=cat(2,sw_freq_mat_clean{ind,n_prepro});
        sw_freq_sub1{i,n_prepro}=cat(2,sw_freq_mat_clean{label_table_clean{:,'uni_id'}==i*2-1,n_prepro});


        for j=1:sum(ind)
            T_runs{i,n_prepro}(j)=size(dfc_mat_sep_clean{inx(j),n_prepro},2);
        end
        for k=1:sum(ind1)
            T_runs1{i,n_prepro}(k)=size(dfc_mat_sep_clean{inx1(k),n_prepro},2);
        end
end
end

Card_night1=cell2mat(Card_night1);
Card_night2=cell2mat(Card_night2);
ppg_rvt_night1=cell2mat(ppg_rvt_night1);
ppg_rvt_night2=cell2mat(ppg_rvt_night2);

eog_night1=cell2mat(eog_night1);
eog_night2=cell2mat(eog_night2);


sw_amp_sub=cat(2,sw_amp_sub{:,n_prepro});
sw_freq_sub=cat(2,sw_freq_sub{:,n_prepro});
sw_amp_sub1=cat(2,sw_amp_sub1{:,n_prepro});
sw_freq_sub1=cat(2,sw_freq_sub1{:,n_prepro});


dfc_sub_mat_all=[];dfc_sub_mat_all1=[];
for n_prepro=1

    dfc_sub_mat_all{1,n_prepro}=cat(2,dfc_sub_mat{:,n_prepro})';
    dfc_sub_mat_all1{1,n_prepro}=cat(2,dfc_sub_mat1{:,n_prepro})';
end
dfc_sub_mat_all{n_prepro}(isnan(dfc_sub_mat_all{n_prepro}))=0;
dfc_sub_mat_all1{n_prepro}(isnan(dfc_sub_mat_all1{n_prepro}))=0;


T=[];T1=[];

for n_prepro=1
    for i=1:length(unique(label_table_clean{:,'uni_id'}))/2
    T{n_prepro}(i,1)=size(dfc_sub_mat{i,n_prepro},2);
    T1{n_prepro}(i,1)=size(dfc_sub_mat1{i,n_prepro},2);
    end
end

%% HMM code

options = struct();
e = explainedvar_PCA(dfc_sub_mat_all{n_prepro},T{n_prepro},options); % how much variance PCA explains on the data
pc1 = find(e>0.5,1); % get no. of PCA components to explain 50% of variance
pc2 = find(e>0.6,1); % get no. of PCA components to explain 60% of variance
pc3 = find(e>0.7,1); % get no. of PCA components to explain 70% of variance
pc4 = find(e>0.8,1); % get no. of PCA components to explain 80% of variance
pc5 = find(e>0.9,1); % get no. of PCA components to explain 90% of variance
pc6 = 0;
pc7 = find(e>0.4,1); % get no. of PCA components to explain 40% of variance

individual_e=[e(1);diff(e)];
pc_comp=[pc1,pc2,pc3];

% figure;
% yyaxis left
% plot(e)
% ylabel('Cumulative Variance Explained','FontSize',24)
% yline([ e(13)],'--')
% yyaxis right
% plot(individual_e)
% ylabel('Variance Explained','FontSize',24)
% xlabel('Number of Components','FontSize',24)
% xline([13],'--')
% set(gcf,'Position',[800 800 1600 1600])
% 


hmm_all_zscore=[];template_configuration=[];
lambda_all_zscore=[];gamma_all_zscore=[];vpath_all_zscore=[];
n_prepro=1;vpath=[];
d=[];p=[];stats=[];mean_perm=[];std_perm=[];maxFO=[];FO=[];MeanLifeTimes=[];
for n_pca=1
    for k=4:25
        template_configuration = struct();
        template_configuration.order = 0; 
        template_configuration.dropstates = 0; 
        template_configuration.verbose = 0;
        template_configuration.cyc = 1000;
        template_configuration.initcyc = 25;
        template_configuration.pca = 13;
        template_configuration.K = k;
        template_configuration.standardise = 1;
        template_configuration.zeromean = 0; 
        template_configuration.covtype = 'full';
        
        tic
        try
        [hmm,Gamma,~,vpath] = hmmmar(dfc_sub_mat_all,T,template_configuration);
        catch me
            disp(me);
            continue
        end
        toc
%         test_gamma=zeros(size(Gamma,1),1);
%         for i=1:size(Gamma,2)
%             test_gamma=test_gamma+round(Gamma(:,i))*i;
%         end
        
%         [~,test_gamma]=max(Gamma,[],2);
        [d{k-3,n_pca},p{k-3,n_pca},stats{k-3,n_pca}]=manova1(vpath,sleep_score_night2);
        
        hmm_all_zscore{k-3,n_pca}=hmm;
        gamma_all_zscore{k-3,n_pca}=Gamma;
        vpath_all_zscore{k-3,n_pca}=vpath;
        lambda_all_zscore{k-3,n_pca}=stats{k-3,n_pca}.lambda;
        lambda_perm=nan(1000,1);
%         for i=1:1000
%         [~,~,stats1]=manova1(test_gamma,sleep_score_night2(randperm(size(sleep_score_night2,1))));
%         lambda_perm(i)=stats1.lambda;
% %         end
%         mean_perm(k-3,n_pca)=mean(lambda_perm);
%         std_perm(k-3,n_pca)=std(lambda_perm);
        maxFO{k-3,n_pca} = getMaxFractionalOccupancy(Gamma,T{:},template_configuration); % useful to diagnose if the HMM 
                    % is capturing dynamics or grand between-subject 
                    % differences (see Wiki)
        FO{k-3,n_pca} = getFractionalOccupancy (Gamma,T{:},template_configuration); % state fractional occupancies per session
%         LifeTimes = getStateLifeTimes (Gamma,T,template_configuration); % state life times
        LifeTimes = getStateLifeTimes (vpath,T{:},template_configuration); % state life times
        disp(k);
        MeanLifeTimes{k-3,n_pca} = mean(cat(2,LifeTimes{:}));
                % 
        %     t = hmmtest(Gamma{i},T,Tsubject,Y,options_test,hmm);
        %     test_group{i} = t.grouplevel; % only doing group-level testing
        %     disp([num2str(i) ' of ' num2str(L)])
        % 
        % 
    end
end



for n_pca=1
    for k=4:25
         template_configuration.K = k;
        maxFO{k-3,n_pca} = getMaxFractionalOccupancy(gamma_all_zscore{k-3,n_pca},T{:},template_configuration); % useful to diagnose if the HMM 
                    % is capturing dynamics or grand between-subject 
                    % differences (see Wiki)
        FO{k-3,n_pca} = getFractionalOccupancy (gamma_all_zscore{k-3,n_pca},T{:},template_configuration); % state fractional occupancies per session
%         LifeTimes = getStateLifeTimes (Gamma,T,template_configuration); % state life times
        LifeTimes = getStateLifeTimes (vpath_all_zscore{k-3,n_pca},T{:},template_configuration); % state life times
        disp(k);
        MeanLifeTimes{k-3,n_pca} = mean(cat(2,LifeTimes{:}));
                % 
    end
end


medianFO=[];
for ii=1:size(FO,1)
    for n_pca=1
        try
            medianFO(1:12,ii,n_pca)=median(FO{ii,n_pca},2);
        catch
            medianFO(1:12,ii,n_pca)=zeros(12,1);
            continue
        end
    end


end

free_energy=ones(21,1);
for ii=1:22
    free_energy(ii)=hmmfe(dfc_sub_mat_all,T,hmm_all_zscore{ii,1},gamma_all_zscore{ii,1});
end

LifeTimes = getStateLifeTimes (vpath_all_zscore{k-3,n_pca},T{:},template_configuration);
mean_LT=cellfun(@mean, LifeTimes(1:21));
std_LT=cellfun(@std, LifeTimes(1:21));
length_LT=cellfun(@length,LifeTimes(1:21));
group_LT=[];
for i=1:21
group_LT=[group_LT;ones(length_LT(i),1)*i];
end
LifeTimes=cat(2,LifeTimes{:});
figure;
% boxplot(LifeTimes',group_LT)
bar(1:21,mean_LT(I2))
hold on

er = errorbar(1:21,mean_LT(I2)',std_LT(I2)./sqrt(length_LT(I2)),std_LT(I2)./sqrt(length_LT(I2)));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off
title('Mean LifeTime of states')
xticks(1:21)
xticklabels(I2)
set(gcf,'Position',[1600 1600 2000 2000])
xlabel('HMM states')
ylabel('TR (3s)')
% xlim([3 size(hmm_all_zscore,1)+4])

% mean(reshape(cell2mat(maxFO(:,2)),12,[]))'% mean maxFo for each k order
n_pca=1;
figure
subplot(5,1,1)
plot(4:size(hmm_all_zscore,1)+3,free_energy)
ylabel('Free Energy')
xlabel('Model Order (K)')
xlim([3 size(hmm_all_zscore,1)+4])

subplot(5,1,2)
errorbar(4:size(hmm_all_zscore,1)+3,mean(reshape(cell2mat(maxFO(:,n_pca)),12,[]))',std(reshape(cell2mat(maxFO(:,n_pca)),12,[]))'./sqrt(12))
ylabel('Max Fractional Occupancy')
xlabel('Model Order (K)')

xlim([3 size(hmm_all_zscore,1)+4])


subplot(5,1,3)
errorbar(4:size(hmm_all_zscore,1)+3,mean(medianFO(:,:,n_pca) )',std(medianFO(:,:,n_pca) )'./sqrt(12))
ylabel('Median Fractional Occupancy')
xlim([3 size(hmm_all_zscore,1)+4])
xlabel('Model Order (K)')


subplot(5,1,4)
plot(4:size(hmm_all_zscore,1)+3,cell2mat(lambda_all_zscore(:,n_pca)))
ylabel('Wilk''s lambda ')
xlim([3 size(hmm_all_zscore,1)+4])
xlabel('Model Order (K)')


subplot(5,1,5)
plot(4:size(hmm_all_zscore,1)+3,cell2mat(MeanLifeTimes(:,n_pca)))
ylabel('Mean HMM State Life Time')
xlim([3 size(hmm_all_zscore,1)+4])
xlabel('Model Order (K)')

set(gcf,'Position',[1600 1600 2000 2000])




ind_sub_sep=[0;T{1}];

for i=1:25
    ind_sub(i)=sum(ind_sub_sep(1:i));
end
for j=1:length(gamma_all)
for i=1:24
   num_of_state(i,j)=length(unique(gamma_all{j}((ind_sub(i)+1):ind_sub(i+1))));
end
end


lambda_perm=[];
for j=1:length(gamma_all_zscore)
    for i=1:1000
        [d,p,stats]=manova1(gamma_all_zscore{j},sleep_score_night2(randperm(size(sleep_score_night2,1))));
        lambda_perm{j}(i)=stats.lambda;
    end
end

for j=1:length(gamma_all_zscore)

mean_perm(j)=mean(lambda_perm{j});
std_perm(j)=std(lambda_perm{j});
end
% 
% N = 25; % subjects
% Q = 4; % sessions per subject
% ttrial = 500; % time points
% nregions = 50; % regions or voxels
% Y = zeros(N*Q,2); % design matrix with conditions
% X = randn(Q*N*ttrial,nregions); % all data concatenated
% T = ttrial * ones(N*Q,1);  % length of data for each session
% Tsubject = Q*ttrial * ones(N,1);  % length of data for each subjec





%% gsr 21 states
k=18;sleep_states21=[];n_pca=1;
% for k=1:length(gamma_all_zscore)

for i=1:k+3

sleep_states21(i)=mode(sleep_score_night2(vpath_all_zscore{k,n_pca}==i));
end

Time_night=[];Time_run=[];Time_night1=[];Time_run1=[];Time_dream=[];Time_dream1=[];
for ii=1:12
Temp1=1:T{1}(ii);
Temp1_n1=1:T1{1}(ii);
Temp2=[];Temp2_n1=[];Temp4=[];Temp4_n1=[];
    for jj=1:length(T_runs{ii})
    
        Temp3=1:T_runs{ii}(jj);
        Temp5=zeros(size(Temp3));
        Temp5(end:-1:end-59)=1;
        Temp2=[Temp2;Temp3'];
        Temp4=[Temp4;Temp5'];
    end

    for kk=1:length(T_runs1{ii})
    
        Temp3=1:T_runs1{ii}(kk);
        Temp5=zeros(size(Temp3));
        Temp5(end:-1:end-59)=1;
        Temp2_n1=[Temp2_n1;Temp3'];
        Temp4_n1=[Temp4_n1;Temp5'];

    end

Time_run=[Time_run;Temp2];
Time_night=[Time_night;Temp1'];
Time_run1=[Time_run1;Temp2_n1];
Time_night1=[Time_night1;Temp1_n1'];
Time_dream=[Time_dream;Temp4];
Time_dream1=[Time_dream1;Temp4_n1];
end

TP=getTransProbs(hmm_all_zscore{18,1});%21 states
Ci=modularity_und(TP,1.5);

indx=1:21;
indx=[indx;sleep_states21;Ci'];
[B,I]=sort(indx(2,:));
[~,I2]=sort(indx(3,:));
%  I2=[ 13 16  18  1  2  3  6  8  10   4 5  9    11    12    14  7    15    17    19 ];


TP_module=nan(21,21);
TP_rerange=nan(21,21);
for ii=1:21
temp=TP(:,I(ii));
temp_mod=TP(:,I2(ii));
TP_rerange(:,ii)=temp(I);
TP_module(:,ii)=temp_mod(I2);
end


    
figure
imagesc(TP_rerange)
set(gcf,'Position',[100 100 800 800])
set(gca, 'XTick', 1:21); % center x-axis ticks on bins
set(gca, 'YTick', 1:21); % center y-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel', indx(1,I)); % set x-axis labels
set(gca, 'YTickLabel', indx(1,I)); % set y-axis labels
c = colorbar;
c.FontSize=12;
colormap(bluewhitered)
% clim([-1 1])
c.TickLength=0;
c.Box='off';
xtickangle(45)

set(gca,'box','off') 




    
figure
imagesc(TP_module)
set(gcf,'Position',[100 100 800 800])
set(gca, 'XTick', 1:21); % center x-axis ticks on bins
set(gca, 'YTick', 1:21); % center y-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
set(gca, 'YTickLabel', indx(1,I2)); % set y-axis labels
c = colorbar;
c.FontSize=12;
c.Label.String = 'Transtion Probability';
c.FontSize=12;
colormap(bluewhitered)
c.Label.FontSize = 18;
% clim([-1 1])
c.TickLength=0;
c.Box='off';
xtickangle(45)

set(gca,'box','off') 




example_mat=[];
for ii=1:21
example_mat(ii,1)=mean(sw_amp_sub(vpath_all_zscore{18,1}==ii),'omitnan');

example_mat(ii,2)=mean(sw_freq_sub(vpath_all_zscore{18,1}==ii),'omitnan');

example_mat(ii,3)=mean(Card_night2(vpath_all_zscore{18,1}==ii,1),'omitnan');

example_mat(ii,4)=mean(Card_night2(vpath_all_zscore{18,1}==ii,3),'omitnan');


example_mat(ii,5)=std(Card_night2(vpath_all_zscore{18,1}==ii,1),'omitnan')./...
            sqrt(sum(~isnan(Card_night2(vpath_all_zscore{18,1}==ii,1))));

example_mat(ii,6)=std(Card_night2(vpath_all_zscore{18,1}==ii,3),'omitnan')./...
            sqrt(sum(~isnan(Card_night2(vpath_all_zscore{18,1}==ii,3))));
example_mat(ii,7)=sum(~isnan(sw_amp_sub(vpath_all_zscore{18,1}==ii)))./sum((vpath_all_zscore{18,1}==ii));


example_mat(ii,8)=mean(Time_run(vpath_all_zscore{18,1}==ii),'omitnan');
example_mat(ii,9)=mean(Time_night(vpath_all_zscore{18,1}==ii),'omitnan');
example_mat(ii,10)=std(Time_run(vpath_all_zscore{18,1}==ii),'omitnan')./...
            sqrt(sum(~isnan(Time_run(vpath_all_zscore{18,1}==ii))));
example_mat(ii,11)=std(Time_night(vpath_all_zscore{18,1}==ii),'omitnan')./...
            sqrt(sum(~isnan(Time_night(vpath_all_zscore{18,1}==ii))));


example_mat(ii,12)=std(sw_amp_sub(vpath_all_zscore{18,1}==ii),'omitnan')./...
            sqrt(sum(~isnan(sw_amp_sub(vpath_all_zscore{18,1}==ii))));

example_mat(ii,13)=std(sw_freq_sub(vpath_all_zscore{18,1}==ii),'omitnan')./...
            sqrt(sum(~isnan(sw_freq_sub(vpath_all_zscore{18,1}==ii))));


example_mat(ii,14)=mean(ppg_rvt_night2(vpath_all_zscore{18,1}==ii,1),'omitnan');

example_mat(ii,15)=mean(ppg_rvt_night2(vpath_all_zscore{18,1}==ii,2),'omitnan');


example_mat(ii,16)=std(ppg_rvt_night2(vpath_all_zscore{18,1}==ii,1),'omitnan')./...
            sqrt(sum(~isnan(ppg_rvt_night2(vpath_all_zscore{18,1}==ii,1))));

example_mat(ii,17)=std(ppg_rvt_night2(vpath_all_zscore{18,1}==ii,2),'omitnan')./...
            sqrt(sum(~isnan(ppg_rvt_night2(vpath_all_zscore{18,1}==ii,2))));



example_mat(ii,18)=mean(eog_night2(vpath_all_zscore{18,1}==ii,1),'omitnan');

example_mat(ii,19)=mean(eog_night2(vpath_all_zscore{18,1}==ii,2),'omitnan');


example_mat(ii,20)=std(eog_night2(vpath_all_zscore{18,1}==ii,1),'omitnan')./...
            sqrt(sum(~isnan(eog_night2(vpath_all_zscore{18,1}==ii,1))));

example_mat(ii,21)=std(eog_night2(vpath_all_zscore{18,1}==ii,2),'omitnan')./...
            sqrt(sum(~isnan(eog_night2(vpath_all_zscore{18,1}==ii,2))));



end

font_size=18;
figure;
subplot(6,1,1)
plot(1:21,example_mat(I2,7))
set(gca, 'XTick', 1:21); % center x-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
ylabel('Slow Wave Density','FontSize',font_size)
xlim([0 22])
ylim([0 0.9])



subplot(6,1,2)
errorbar(1:21,example_mat(I2,14),example_mat(I2,16))

set(gca, 'XTick', 1:21); % center x-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
ylabel('Variation in PPG AMP','FontSize',font_size)
% xlabel('HMM states')
xlim([0 22])



subplot(6,1,3)
errorbar(1:21,example_mat(I2,15),example_mat(I2,17))

set(gca, 'XTick', 1:21); % center x-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
ylabel('Variation in RespRVT','FontSize',font_size)
% xlabel('HMM states')
xlim([0 22])



% subplot(7,1,4)
% errorbar(1:21,example_mat(I2,3),example_mat(I2,5))
% set(gca, 'XTick', 1:21); % center x-axis ticks on bins
% set(gca, 'TickLength',[0 0])
% set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
% ylabel('Card Amp')
% xlabel('HMM states')
% xlim([0 22])



subplot(6,1,4)
errorbar(1:21,example_mat(I2,4),example_mat(I2,6))
set(gca, 'XTick', 1:21); % center x-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
ylabel('Variation in Heart Rates','FontSize',font_size)
% xlabel('HMM states')
xlim([0 22])

% subplot(6,1,5)
% errorbar(1:21,example_mat(I2,18),example_mat(I2,20))
% set(gca, 'XTick', 1:21); % center x-axis ticks on bins
% set(gca, 'TickLength',[0 0])
% set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
% ylabel('EOG1','FontSize',font_size)
% % xlabel('HMM states')
% xlim([0 22])
% 
% 
% 
% subplot(6,1,6)
% errorbar(1:21,example_mat(I2,19),example_mat(I2,21))
% set(gca, 'XTick', 1:21); % center x-axis ticks on bins
% set(gca, 'TickLength',[0 0])
% set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
% ylabel('EOG2','FontSize',font_size)
% % xlabel('HMM states')
% xlim([0 22])



subplot(6,1,5)
errorbar(1:21,example_mat(I2,8),example_mat(I2,10))
set(gca, 'XTick', 1:21); % center x-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
ylabel('Time of Run (TR)','FontSize',font_size)
% xlabel('HMM states')
xlim([0 22])



subplot(6,1,6)
errorbar(1:21,example_mat(I2,9),example_mat(I2,11))
set(gca, 'XTick', 1:21); % center x-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
ylabel('Time of Night (TR)','FontSize',font_size)
% xlabel('HMM states')
xlim([0 22])



set(gcf,'Position',[1600 1600 2000 2000])



%semisupervised learning
data.X=[dfc_sub_mat_all{1};dfc_sub_mat_all1{1}];
data.C=[gamma_all_zscore{18,1};nan(65231,21)];
    template_configuration = struct();
    template_configuration.order = 0; 
    template_configuration.dropstates = 0; 
    template_configuration.verbose = 0;
    template_configuration.cyc = 1000;
    template_configuration.initcyc = 25;
    template_configuration.pca = 13;
    template_configuration.K = 21;
    template_configuration.standardise = 1;
    template_configuration.zeromean = 0; 
    template_configuration.covtype = 'full';
    
    [hmm,Gamma,~,vpath1] = hmmmar(data,[T{1};T1{1}],template_configuration);




set(gcf,'Position',[1600 1600 2000 2000])

sleep_label_num={'0','2','3','4','5', '6'};
sleep_per=[];
% figure
for i=1:21
    temp_s=sleep_score_night2(vpath_all_zscore{18,1}==i);
    
    % subplot(21,1,i)
    sleep_per(i,:)=[ sum(temp_s==0) sum(temp_s==2) sum(temp_s==3) sum(temp_s==4) sum(temp_s==5)...
    sum(temp_s==6)]./sum([ sum(temp_s==0) sum(temp_s==2) sum(temp_s==3) sum(temp_s==4) sum(temp_s==5)...
    sum(temp_s==6)]);
    % pie(sleep_per(I2(i),:),{'','','','','', ''})
    % lgd = legend(labels_sleep);
    % title(['state',num2str(I2(i))])
end
labels_sleep={'Undefined','N3','N2','N1','REM','Awake'};

figure
for i=1:21
    temp_s=sleep_score_night2(vpath_all_zscore{18,1}==i);
    
    subplot(3,7,i)
   
    pie(sleep_per(i,:),{'','','','','', ''})
    % lgd = legend(labels_sleep);
    title(['State ',num2str(i)],'FontSize',20)
end



set(gcf,'Position',[800 1000 2000 1000])

sleep_per1=[];
for i=1:21
    temp_s1=sleep_score_night1(vpath1(66830:end,:)==i);
    sleep_per1(i,:)=[ sum(temp_s1==0) sum(temp_s1==2) sum(temp_s1==3) sum(temp_s1==4) sum(temp_s1==5)...
    sum(temp_s1==6)]./sum([ sum(temp_s1==0) sum(temp_s1==2) sum(temp_s1==3) sum(temp_s1==4) sum(temp_s1==5)...
    sum(temp_s1==6)]);
end


figure
for i=1:21
    subplot(3,7,i)

        pie(sleep_per1(i,:),{'','','','','', ''})

    title(['State ',num2str(i)],'FontSize',20)
    
end

set(gcf,'Position',[800 1000 2000 1000])





figure
pie([sum(sleep_score_night2==0) sum(sleep_score_night2==2) sum(sleep_score_night2==3) sum(sleep_score_night2==4) sum(sleep_score_night2==5)...
    sum(sleep_score_night2==6) ])
lgd = legend(labels_sleep);
lgd.Layout.Tile = 'east';

figure
pie([sum(sleep_score_night1==0) sum(sleep_score_night1==2) sum(sleep_score_night1==3) sum(sleep_score_night1==4) sum(sleep_score_night1==5)...
    sum(sleep_score_night1==6) ])
lgd = legend(labels_sleep);
title('Night1')







[t,p]=corr(sleep_per(:),sleep_per1(:)) %t=0.94 p<1e-58




example_mat1=[];
for ii=1:21
example_mat1(ii,1)=mean(sw_amp_sub1(vpath1(66830:end,:)==ii),'omitnan');

example_mat1(ii,2)=mean(sw_freq_sub1(vpath1(66830:end,:)==ii),'omitnan');

example_mat1(ii,3)=mean(Card_night1(vpath1(66830:end,:)==ii,1),'omitnan');

example_mat1(ii,4)=mean(Card_night1(vpath1(66830:end,:)==ii,3),'omitnan');


example_mat1(ii,5)=std(Card_night1(vpath1(66830:end,:)==ii,1),'omitnan')./...
            sqrt(sum(~isnan(Card_night1(vpath1(66830:end,:)==ii,1))));

example_mat1(ii,6)=std(Card_night1(vpath1(66830:end,:)==ii,3),'omitnan')./...
            sqrt(sum(~isnan(Card_night1(vpath1(66830:end,:)==ii,3))));

example_mat1(ii,7)=sum(~isnan(sw_amp_sub1(vpath1(66830:end,:)==ii)))./sum((vpath1(66830:end,:)==ii));


example_mat1(ii,8)=mean(Time_run1(vpath1(66830:end,:)==ii),'omitnan');
example_mat1(ii,9)=mean(Time_night1(vpath1(66830:end,:)==ii),'omitnan');
example_mat1(ii,10)=std(Time_run1(vpath1(66830:end,:)==ii),'omitnan')./...
            sqrt(sum(~isnan(Time_run1(vpath1(66830:end,:)==ii))));
example_mat1(ii,11)=std(Time_night1(vpath1(66830:end,:)==ii),'omitnan')./...
            sqrt(sum(~isnan(Time_night1(vpath1(66830:end,:)==ii))));

example_mat1(ii,12)=std(sw_amp_sub1(vpath1(66830:end,:)==ii),'omitnan')./...
            sqrt(sum(~isnan(sw_amp_sub1(vpath1(66830:end,:)==ii))));

example_mat1(ii,13)=std(sw_freq_sub1(vpath1(66830:end,:)==ii),'omitnan')./...
            sqrt(sum(~isnan(sw_freq_sub1(vpath1(66830:end,:)==ii))));




example_mat1(ii,14)=mean(ppg_rvt_night1(vpath1(66830:end,:)==ii,1),'omitnan');

example_mat1(ii,15)=mean(ppg_rvt_night1(vpath1(66830:end,:)==ii,2),'omitnan');


example_mat1(ii,16)=std(ppg_rvt_night1(vpath1(66830:end,:)==ii,1),'omitnan')./...
            sqrt(sum(~isnan(ppg_rvt_night1(vpath1(66830:end,:)==ii,1))));

example_mat1(ii,17)=std(ppg_rvt_night1(vpath1(66830:end,:)==ii,2),'omitnan')./...
            sqrt(sum(~isnan(ppg_rvt_night1(vpath1(66830:end,:)==ii,2))));





end


font_size=18;
figure;
subplot(6,1,1)
plot(1:21,example_mat1(I2,7))
set(gca, 'XTick', 1:21); % center x-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
ylabel('Slow Wave Density','FontSize',font_size)
xlim([0 22])
ylim([0 0.9])



subplot(6,1,2)
errorbar(1:21,example_mat1(I2,14),example_mat1(I2,16))

set(gca, 'XTick', 1:21); % center x-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
ylabel('Variation in PPG AMP','FontSize',font_size)
% xlabel('HMM states')
xlim([0 22])



subplot(6,1,3)
errorbar(1:21,example_mat1(I2,15),example_mat1(I2,17))

set(gca, 'XTick', 1:21); % center x-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
ylabel('Variation in RespRVT','FontSize',font_size)
% xlabel('HMM states')
xlim([0 22])



% subplot(7,1,4)
% errorbar(1:21,example_mat(I2,3),example_mat(I2,5))
% set(gca, 'XTick', 1:21); % center x-axis ticks on bins
% set(gca, 'TickLength',[0 0])
% set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
% ylabel('Card Amp')
% xlabel('HMM states')
% xlim([0 22])



subplot(6,1,4)
errorbar(1:21,example_mat1(I2,4),example_mat1(I2,6))
set(gca, 'XTick', 1:21); % center x-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
ylabel('Variation in Heart Rates','FontSize',font_size)
% xlabel('HMM states')
xlim([0 22])


subplot(6,1,5)
errorbar(1:21,example_mat1(I2,8),example_mat1(I2,10))
set(gca, 'XTick', 1:21); % center x-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
ylabel('Time of Run (TR)','FontSize',font_size)
% xlabel('HMM states')
xlim([0 22])



subplot(6,1,6)
errorbar(1:21,example_mat1(I2,9),example_mat1(I2,11))
set(gca, 'XTick', 1:21); % center x-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
ylabel('Time of Night (TR)','FontSize',font_size)
% xlabel('HMM states')
xlim([0 22])



set(gcf,'Position',[1600 1600 2000 2000])



k=21;n_pca=1;mean_act=nan(k,300);corrmat=[];corrmat_all=[];
for n_state=1:k
    mean_act(n_state,:)=getMean(hmm_all_zscore{k-3,n_pca},n_state); % activation
    [~,corrmat{n_state}]=getFuncConn(hmm_all_zscore{k-3,n_pca},n_state); % fc
    % [~,corrmat_all{n_state}]=getFuncConn(hmm,n_state); % fc
end
mean_act=mean_act-mean(mean_act,1);% activation against all states

r_act=ones(21,21);p_act=ones(21,21);
for ii = 1:21
    for jj =1:21
        [r_act(ii,jj),p_act(ii,jj)]=corr(mean_act(ii,:)',mean_act(jj,:)');

    end
end

figure
imagesc(r_act(I2,I2));
c = colorbar;

c.FontSize=12;
colormap(bluewhitered)
% clim([-1 1])
c.TickLength=0;
c.Box='off';
xtickangle(45)
set(gca, 'YTick', 1:21)
set(gca, 'YTickLabel', indx(1,I2)); % set x-axis labels
set(gca, 'XTick', 1:21)
set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
set(gca,'box','off') 
set(gcf,'Position',[8000 8000 2000 2000])
    


% mean_act_diff=nan(k,300);
% for n_state=1:k
% per_pos=prctile(mean_act(n_state,mean_act(n_state,:)>0),50);
% per_neg=prctile(mean_act(n_state,mean_act(n_state,:)<0),50);
% temp=mean_act(n_state,:);
% temp(mean_act(n_state,:)<per_pos&mean_act(n_state,:)>per_neg)=0;
% mean_act_diff(n_state,:)=temp;
% %     save(['GSR_21states_mean_',num2str(n_state),'.txt'], 'temp','-ascii');
% end 
% 

figure
imagesc(mean_act(I2,:));
% h=gca;
% h.XRuler.Axle.LineStyle='none';
% Create title
% title(['State\_',num2str(I(n_state)),'\_SleepStages\_',num2str(sleep_states21(I(n_state)))],'FontSize',14);
c = colorbar;

c.FontSize=12;
colormap(bluewhitered)
% clim([-1 1])
c.TickLength=0;
c.Box='off';
xtickangle(45)
set(gca, 'YTick', 1:21)
set(gca, 'YTickLabel', indx(1,I2)); % set x-axis labels
set(gca,'box','off') 
set(gcf,'Position',[8000 8000 2000 2000])









corr_mat=nan(300,300,21);
corr_mat_all=nan(300,300,21);
corr_diff=nan(300,300,21);

for n_state=1:k
corr_mat(:,:,n_state)=corrmat{n_state};
% corr_mat_all(:,:,n_state)=corrmat_all{n_state};
% corr_diff(:,:,n_state)=corr_mat(:,:,n_state)-corr_mat_mean;
end

corr_mat_mean=mean(corr_mat,3);

corr_diff=corr_mat-corr_mat_mean;
% corr_diff_all=corr_mat_all-mean(corr_mat_all,3);




% 
% corr_diff_threshold=nan(300,300,21);
% for n_state=1:k
%     temp=corr_diff(:,:,n_state);
% per_pos=prctile(temp(temp(:)>0),99.5,'all');
% per_neg=prctile(temp(temp(:)<0),0.5,'all');
% temp(temp(:)<per_pos&temp(:)>per_neg)=0;
% 
% corr_diff_threshold(:,:,n_state)=temp;
% %     save(['GSR_21states_mean_',num2str(n_state),'.txt'], 'temp','-ascii');
% end





figure

for n_state=1:k
subplot(4,6,n_state)
imagesc(corr_diff_all(:,:,I2(n_state)));
% h=gca;
% h.XRuler.Axle.LineStyle='none';
% Create title
title(['State\_',num2str(I2(n_state)),'\_SleepStages\_',num2str(sleep_states21(I2(n_state)))],'FontSize',14);
c = colorbar;

c.FontSize=12;
colormap(bluewhitered)
clim([-0.5 0.5])
c.TickLength=0;
c.Box='off';
xtickangle(45)

set(gca,'box','off') 
end
set(gcf,'Position',[8000 8000 2000 2000])





net_code=unique(ROI_list{:,'Network_Label_v2'});

w_wch=zeros(length(net_code),1);%
b_wch=zeros(length(net_code),1);
bw_wch=zeros(length(net_code)*(length(net_code)-1)/2,1);


net_mat=ones(length(net_code),length(net_code));
net_mat_tril=tril(net_mat,-1);
% net_mat_tril=[net_mat_tril,zeros(length(net_code),1)];
net_mat_diag=eye(length(net_code),length(net_code));
% net_mat_diag=[net_mat_diag,zeros(length(net_code),1)];
corr_mat_sep=nan(16,17,21);
temp_mat=corr_diff_all;
for k=1:size(corr_mat,3) 


    for n_net = 1:length(net_code)
        net_nodes=ROI_list{ROI_list{:,'Network_Label_v2'}==net_code(n_net),'ROI_num'};
            m_net=temp_mat(net_nodes,net_nodes,k);
            m_net=tril(m_net,-1);
    %         m_net_s(:,n_sub)=m_net(m_net>0);
            ind_mat_temp=ones(numel(net_nodes),numel(net_nodes));
            ind_mat_temp=tril(ind_mat_temp,-1);
            ind_mat_temp=ind_mat_temp.*~isnan(m_net);
            w_wch(n_net)=nansum(nansum(m_net))/sum(sum(ind_mat_temp));
        
    end


    for n_net = 1:length(net_code)
        net_nodes=ROI_list{ROI_list{:,'Network_Label_v2'}==net_code(n_net),'ROI_num'};
        net_nodes_2=ROI_list{ROI_list{:,'Network_Label_v2'}~=net_code(n_net),'ROI_num'};
            m_net=temp_mat(net_nodes,net_nodes_2,k);
            ind_mat_temp=ones(numel(net_nodes),numel(net_nodes_2));
            ind_mat_temp=ind_mat_temp.*~isnan(m_net);
%             ind_mat_temp=tril(ind_mat,-1);
            b_wch(n_net)=nansum(nansum(m_net))/sum(ind_mat_temp(:));
        
    end




  n=1;
    for n_net = 1:length(net_code)-1
        for n_net2= n_net+1:length(net_code)
            net_nodes=ROI_list{ROI_list{:,'Network_Label_v2'}==net_code(n_net),'ROI_num'};
            net_nodes_2=ROI_list{ROI_list{:,'Network_Label_v2'}==net_code(n_net2),'ROI_num'};
                m_net=temp_mat(net_nodes,net_nodes_2,k);
                ind_mat_temp=ones(numel(net_nodes),numel(net_nodes_2));
                ind_mat_temp=ind_mat_temp.*~isnan(m_net);

%                 ind_mat_temp=tril(ind_mat,-1);
    
                bw_wch(n)=nansum(nansum(m_net))/sum(ind_mat_temp(:));
              
            
         n=n+1;
        end
        
    end




       mat_test=zeros(length(net_code),length(net_code));
       mat_test(net_mat_diag>0)=w_wch;
       mat_test(net_mat_tril>0)=bw_wch;
       mat_test(:,length(net_code)+1)=b_wch;
       corr_mat_sep(:,:,k)=mat_test;

end

net_n_sep={'Un','DMN','VIS','FPN','REW','DAN','VAN','SAL','CON','dSMN','lSMN','AUD','PMN','MTL','Sub','CB'};







%% seperate subcortical to BG/HIP/THA

net_code=unique(ROI_list{:,'Network_Label_v3'});

w_wch=zeros(length(net_code),1);%
b_wch=zeros(length(net_code),1);
bw_wch=zeros(length(net_code)*(length(net_code)-1)/2,1);


net_mat=ones(length(net_code),length(net_code));
net_mat_tril=tril(net_mat,-1);
% net_mat_tril=[net_mat_tril,zeros(length(net_code),1)];
net_mat_diag=eye(length(net_code),length(net_code));
% net_mat_diag=[net_mat_diag,zeros(length(net_code),1)];
corr_mat_sep=nan(18,19,21);
temp_mat=corr_diff;
for k=1:size(temp_mat,3) 


    for n_net = 1:length(net_code)
        net_nodes=ROI_list{ROI_list{:,'Network_Label_v3'}==net_code(n_net),'ROI_num'};
            m_net=temp_mat(net_nodes,net_nodes,k);
            m_net=tril(m_net,-1);
    %         m_net_s(:,n_sub)=m_net(m_net>0);
            ind_mat_temp=ones(numel(net_nodes),numel(net_nodes));
            ind_mat_temp=tril(ind_mat_temp,-1);
            ind_mat_temp=ind_mat_temp.*~isnan(m_net);
            w_wch(n_net)=nansum(nansum(m_net))/sum(sum(ind_mat_temp));
        
    end


    for n_net = 1:length(net_code)
        net_nodes=ROI_list{ROI_list{:,'Network_Label_v3'}==net_code(n_net),'ROI_num'};
        net_nodes_2=ROI_list{ROI_list{:,'Network_Label_v3'}~=net_code(n_net),'ROI_num'};
            m_net=temp_mat(net_nodes,net_nodes_2,k);
            ind_mat_temp=ones(numel(net_nodes),numel(net_nodes_2));
            ind_mat_temp=ind_mat_temp.*~isnan(m_net);
%             ind_mat_temp=tril(ind_mat,-1);
            b_wch(n_net)=nansum(nansum(m_net))/sum(ind_mat_temp(:));
        
    end




  n=1;
    for n_net = 1:length(net_code)-1
        for n_net2= n_net+1:length(net_code)
            net_nodes=ROI_list{ROI_list{:,'Network_Label_v3'}==net_code(n_net),'ROI_num'};
            net_nodes_2=ROI_list{ROI_list{:,'Network_Label_v3'}==net_code(n_net2),'ROI_num'};
                m_net=temp_mat(net_nodes,net_nodes_2,k);
                ind_mat_temp=ones(numel(net_nodes),numel(net_nodes_2));
                ind_mat_temp=ind_mat_temp.*~isnan(m_net);

%                 ind_mat_temp=tril(ind_mat,-1);
    
                bw_wch(n)=nansum(nansum(m_net))/sum(ind_mat_temp(:));
              
            
         n=n+1;
        end
        
    end




       mat_test=zeros(length(net_code),length(net_code));
       mat_test(net_mat_diag>0)=w_wch;
       mat_test(net_mat_tril>0)=bw_wch;
       mat_test(:,length(net_code)+1)=b_wch;
       corr_mat_sep(:,:,k)=mat_test;

end












net_n_sep={'Un','DMN','VIS','FPN','REW','DAN','VAN','SAL','CON','dSMN','lSMN','AUD','PMN','MTL','pHIP','BG','THAL','CB'};

figure

for n_state=1:21
subplot(3,7,n_state)
imagesc(corr_mat_sep(1:18,1:18,I2(n_state)));
% h=gca;
% h.XRuler.Axle.LineStyle='none';
% Create title
title(['State ',num2str(I2(n_state))],'FontSize',14);
set(gca,'TickLength',[0, 0])
c = colorbar;

c.FontSize=12;
clim([-0.5 0.5])
colormap(bluewhitered)

c.TickLength=0;
c.Box='off';
xtickangle(90)
set(gca, 'XTick', 1:18); % center x-axis ticks on bins
set(gca, 'YTick', 1:18); % center y-axis ticks on bins
set(gca, 'XTickLabel', strrep([net_n_sep],'_','\_')); % set x-axis labels
set(gca, 'YTickLabel', strrep([net_n_sep],'_','\_')); % set y-axis label
set(gca,'box','off') 
end
set(gcf,'Position',[8000 8000 2000 2000])


mat_inx=ones(18,18);
mat_inx=tril(mat_inx);
mat_inx=[mat_inx,zeros(18,1)];

for  n_state=1:21
temp_val=corr_mat_sep(:,:,n_state).*mat_inx;
cor_val(n_state,:)=temp_val(temp_val~=0);

end

r=corr(cor_val');

r_vis=r(I2,I2);
r_vis=tril(r_vis,-1);






figure

imagesc(r_vis)
set(gca,'TickLength',[0, 0])
c = colorbar;

c.FontSize=12;
clim([-1 1])
colormap(bluewhitered)

c.TickLength=0;
c.Box='off';
% xtickangle(90)
set(gca, 'XTick', 1:21); % center x-axis ticks on bins
set(gca, 'YTick', 1:21); % center y-axis ticks on bins
set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
set(gca, 'YTickLabel', indx(1,I2)); % set y-axis label
set(gca,'box','off') 
set(gcf,'Position',[8000 8000 2000 2000])


cor_val=[];
mat_inx300=ones(300,300);
mat_inx300=tril(mat_inx300);
% mat_inx=[mat_inx,zeros(18,1)];

for  n_state=1:21
temp_val=corr_diff(:,:,n_state).*mat_inx300;
cor_val(n_state,:)=temp_val(temp_val~=0);

end

r=corr(cor_val');

r_vis300=r(I2,I2);
% r_vis300=tril(r_vis300,-1);


figure

imagesc(r_vis300)
set(gca,'TickLength',[0, 0])
c = colorbar;

c.FontSize=12;
clim([-1 1])
colormap(bluewhitered)

c.TickLength=0;
c.Box='off';
% xtickangle(90)
set(gca, 'XTick', 1:21); % center x-axis ticks on bins
set(gca, 'YTick', 1:21); % center y-axis ticks on bins
set(gca, 'XTickLabel', indx(1,I2)); % set x-axis labels
set(gca, 'YTickLabel', indx(1,I2)); % set y-axis label
set(gca,'box','off') 
set(gcf,'Position',[8000 8000 2000 2000])





net_act=[];

for ii=1:length(net_code)

    net_nodes=ROI_list{ROI_list{:,'Network_Label_v3'}==net_code(ii),'ROI_num'};
    net_act(:,ii)=mean(mean_act(:,net_nodes),2);


end

figure;
imagesc(net_act(I2,:));

title(['Mean Activation'],'FontSize',14);
c = colorbar;

c.FontSize=12;
% clim([-0.5 0.5])
colormap(bluewhitered)

c.TickLength=0;
c.Box='off';
xtickangle(45)
set(gca, 'YTick', 1:21); % center x-axis ticks on bins
set(gca, 'XTick', 1:18); % center y-axis ticks on bins
set(gca, 'YTickLabel', indx(1,I2)); % set x-axis labels
set(gca, 'XTickLabel', strrep([net_n_sep],'_','\_')); % set y-axis label
set(gca,'box','off') 

set(gcf,'Position',[8000 8000 2000 2000])



