%% initial parameters
% This code is for ABCD dasta release 3.0
% Yang, F. N.,1 Xie, W.,1 & Wang, Z. (2022). Effects of sleep duration on neurocognitive development in U.S. early adolescents: a propensity score matched, longitudinal, observational study. Lancet Child & adolescent health. (1 contributed equally)




father_dir='/home/zwang/Documents/Nils/ABCD';%Working directory
cd(father_dir)
beh_dir='/home/zwang/Documents/Nils/ABCD/Beh_tabulated_data/';%folder with all the tabulated data

network_table=readtable('/Beh_tabulated_data/abcd_betnet02.txt');% network conenctivity
network_table_FL2=network_table(strcmp(network_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
network_table_BL=network_table(strcmp(network_table{:,'eventname'},'baseline_year_1_arm_1'),:);
network_table_BL.Properties.RowNames=network_table_BL.src_subject_id;
network_table_FL2.Properties.RowNames=network_table_FL2.src_subject_id;
qc_table=readtable(fullfile(beh_dir, 'abcd_imgincl01.txt'));
qc_table_BL=qc_table(strcmp(qc_table{:,'eventname'},'baseline_year_1_arm_1'),:);
qc_table_FL2=qc_table(strcmp(qc_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
qc_table_BL.Properties.RowNames=qc_table_BL.src_subject_id;
qc_table_FL2.Properties.RowNames=qc_table_FL2.src_subject_id;
sub_id_BL=network_table_BL.Properties.RowNames;% define subject list based on data on network connectivity




sub_id_BL=sub_id_BL(ismember(sub_id_BL,qc_table_BL.Properties.RowNames));% all subs that have qc values

sub_id_BL=sub_id_BL(qc_table_BL{sub_id_BL,'imgincl_rsfmri_include'}==1);% all subs that passed rsfmri qc


sub_network_table=readtable('/Beh_tabulated_data/mrirscor02.txt');% network connectivity between subcortical ROIs and cortical network
sub_network_table_FL2=sub_network_table(strcmp(sub_network_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
sub_network_table_BL=sub_network_table(strcmp(sub_network_table{:,'eventname'},'baseline_year_1_arm_1'),:);
sub_network_table_BL.Properties.RowNames=sub_network_table_BL.src_subject_id;
sub_network_table_FL2.Properties.RowNames=sub_network_table_FL2.src_subject_id;

cbcl_sum=readtable('/Beh_tabulated_data/abcd_cbcls01.txt');% behavioral problems sum scores
cbcl_sum_FL2=cbcl_sum(strcmp(cbcl_sum{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
cbcl_sum_BL=cbcl_sum(strcmp(cbcl_sum{:,'eventname'},'baseline_year_1_arm_1'),:);
cbcl_sum_BL.Properties.RowNames=cbcl_sum_BL.src_subject_id;
cbcl_sum_FL2.Properties.RowNames=cbcl_sum_FL2.src_subject_id;
cbcl_sum_FL1=cbcl_sum(strcmp(cbcl_sum{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
cbcl_sum_FL1.Properties.RowNames=cbcl_sum_FL1.src_subject_id;

sleep_parent_SDSC=readtable(fullfile(beh_dir,'abcd_ssphp01.txt'));% total score (sds_p_ss_total) for Sleep Disturbance Scale for Children 
sleep_parent_SDSC_BL=sleep_parent_SDSC(strcmp(sleep_parent_SDSC{:,'eventname'},'baseline_year_1_arm_1'),:);
sleep_parent_SDSC_BL.Properties.RowNames=sleep_parent_SDSC_BL.src_subject_id;
sleep_parent_SDSC_FL2=sleep_parent_SDSC(strcmp(sleep_parent_SDSC{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
sleep_parent_SDSC_FL2.Properties.RowNames=sleep_parent_SDSC_FL2.src_subject_id;
sleep_parent_SDSC_FL1=sleep_parent_SDSC(strcmp(sleep_parent_SDSC{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
sleep_parent_SDSC_FL1.Properties.RowNames=sleep_parent_SDSC_FL1.src_subject_id;




NIH_TB=readtable(fullfile(beh_dir,'abcd_tbss01.txt'));% total score (sds_p_ss_total) for NIH cognition toolbox 
NIH_TB_BL=NIH_TB(strcmp(NIH_TB{:,'eventname'},'baseline_year_1_arm_1'),:);
NIH_TB_BL.Properties.RowNames=NIH_TB_BL.src_subject_id;
NIH_TB_FL2=NIH_TB(strcmp(NIH_TB{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
NIH_TB_FL2.Properties.RowNames=NIH_TB_FL2.src_subject_id;
NIH_TB_FL1=NIH_TB(strcmp(NIH_TB{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
NIH_TB_FL1.Properties.RowNames=NIH_TB_FL1.src_subject_id;

sub_id_all=NIH_TB_BL.Properties.RowNames;% all subjects with cognition data


sleep_parent_SDS=readtable(fullfile(beh_dir,'abcd_sds01.txt'));% total hours of sleep (sleep_1_p) for SDS 
%(Parent Sleep Disturbance Scale for Children)
sleep_parent_SDS_BL=sleep_parent_SDS(strcmp(sleep_parent_SDS{:,'eventname'},'baseline_year_1_arm_1'),:);
sleep_parent_SDS_BL.Properties.RowNames=sleep_parent_SDS_BL.src_subject_id;
sleep_parent_SDS_FL2=sleep_parent_SDS(strcmp(sleep_parent_SDS{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
sleep_parent_SDS_FL2.Properties.RowNames=sleep_parent_SDS_FL2.src_subject_id;
sleep_parent_SDS_FL1=sleep_parent_SDS(strcmp(sleep_parent_SDS{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
sleep_parent_SDS_FL1.Properties.RowNames=sleep_parent_SDS_FL1.src_subject_id;


sub_id_FL1=sub_id_BL(ismember(sub_id_BL,sleep_parent_SDS_FL1.Properties.RowNames));% 1-year follow-up subjects with sleep data

sub_id_FL2=sub_id_BL(ismember(sub_id_BL,sleep_parent_SDS_FL2.Properties.RowNames));% 2-year follow-up subjects with sleep data



dem_table=readtable(fullfile(beh_dir,'ABCD_sites_effect.xlsx'));% demographic information for subjects, downlowded from DEAP
for i=1:size(dem_table,1)
    dem_table{i,'id'}{1}=[dem_table{i,'id'}{1}(1:4),'_',dem_table{i,'id'}{1}(5:end)];

end
for i=1:size(dem_table,1)
    dem_table{i,'abcd_site'}{1}=[dem_table{i,'abcd_site'}{1}(5:end)];

end
dem_table.abcd_site=str2num(cell2mat(dem_table.abcd_site));

dem_table.Properties.RowNames=dem_table.id;
dem_table.sex=double(strcmp(dem_table.sex_at_birth,'M'));
% coding race 
dem_table.white=double(strcmp(dem_table.race_4level,'White'));
dem_table.black=double(strcmp(dem_table.race_4level,'Black'));
dem_table.mixed=double(strcmp(dem_table.race_4level,'Other/Mixed'));
dem_table.asian=double(strcmp(dem_table.race_4level,'Asian'));
dem_table.other=double(~any([dem_table.white,dem_table.black]'))';
dem_table.race=dem_table.white+2*dem_table.black+3*dem_table.other;
dem_table.race4=dem_table.white+2*dem_table.black+3*dem_table.asian+4*dem_table.mixed;
dem_table.race4(dem_table.race4==0)=4;

pdem_table=readtable(fullfile(beh_dir,'pdem02.txt'));% parent info
pdem_table.Properties.RowNames=pdem_table.src_subject_id;

dem_table.edu1=double(strcmp(dem_table.high_educ,'< HS Diploma'));
dem_table.edu2=double(strcmp(dem_table.high_educ,'HS Diploma/GED'));
dem_table.edu3=double(strcmp(dem_table.high_educ,'Some College'));
dem_table.edu4=double(strcmp(dem_table.high_educ,'Bachelor'));
dem_table.edu5=double(strcmp(dem_table.high_educ,'Post Graduate Degree'));
dem_table.edu=dem_table.edu1+2*dem_table.edu2+3*dem_table.edu3+4*dem_table.edu4+5*dem_table.edu5;%coding educational level





BMI=readtable('/Beh_tabulated_data/BMI.csv');
BMI_BL=BMI(strcmp(BMI{:,'event_name'},'baseline_year_1_arm_1'),:);
BMI_BL.Properties.RowNames=BMI_BL.src_subject_id;





mental_health=readtable('/Beh_tabulated_data/abcd_mhy02.txt');
mental_health_BL=mental_health(strcmp(mental_health{:,'eventname'},'baseline_year_1_arm_1'),:);
mental_health_BL.Properties.RowNames=mental_health_BL.src_subject_id;
mental_health_FL1=mental_health(strcmp(mental_health{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
mental_health_FL1.Properties.RowNames=mental_health_FL1.src_subject_id;
mental_health_FL2=mental_health(strcmp(mental_health{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
mental_health_FL2.Properties.RowNames=mental_health_FL2.src_subject_id;




network_conn_table=([network_table_BL(sub_id_BL,23:191) sub_network_table_BL(sub_id_BL,23:269)]);% combine network connecitivity and subcortical network connectivity
network_conn_table_FL2=([network_table_FL2(sub_id_FL2,23:191) sub_network_table_FL2(sub_id_FL2,23:269)]);

T1_table=readtable('/Beh_tabulated_data/abcd_mrisdp101.txt');% brain structural measurements
T1_table_BL=T1_table(strcmp(T1_table{:,'eventname'},'baseline_year_1_arm_1'),:);
T1_table_BL.Properties.RowNames=T1_table_BL.src_subject_id;
T1_table_FL2=T1_table(strcmp(T1_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
T1_table_FL2.Properties.RowNames=T1_table_FL2.src_subject_id;


T1_sub=readtable('/Beh_tabulated_data/abcd_smrip201.txt');%subcortical structural measurements 
T1_sub_BL=T1_sub(strcmp(T1_sub{:,'eventname'},'baseline_year_1_arm_1'),:);
T1_sub_BL.Properties.RowNames=T1_sub_BL.src_subject_id;
T1_sub_FL2=T1_sub(strcmp(T1_sub{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
T1_sub_FL2.Properties.RowNames=T1_sub_FL2.src_subject_id;



DTI_table=readtable('/Beh_tabulated_data/abcd_dmdtifp101.txt');%DTI structural connectivity
DTI_table_BL=DTI_table(strcmp(DTI_table{:,'eventname'},'baseline_year_1_arm_1'),:);
DTI_table_BL.Properties.RowNames=DTI_table_BL.src_subject_id;
DTI_table_FL2=DTI_table(strcmp(DTI_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
DTI_table_FL2.Properties.RowNames=DTI_table_FL2.src_subject_id;

total_sleep_duration=sleep_parent_SDS_BL{sub_id_BL,'sleepdisturb1_p'};%
%1 = 9-11 hours/ 9 a 11 horas; 2 = 8-9 hours /8 a 9 horas; 3 = 7-8 hours /7 a 8 horas; 4 = 5-7 hours /5 a 7 horas; 5 = Less than 5 hours/ Menos de 5 horas// Consider each question pertaining to the PAST 6 MONTHS of the child's life 
total_sleep_duration(total_sleep_duration==5)=4; % combine 4 & 5 as only 26 subjects have a value of 5
total_sleep_duration=5-total_sleep_duration; %reverse sleep duration to 4 are recommed sleep time.


pubertal_f=sleep_parent_SDSC_BL{sub_id_BL,'pds_p_ss_female_category'};% pubertal status
pubertal_m=sleep_parent_SDSC_BL{sub_id_BL,'pds_p_ss_male_category'};
pubertal=sum([pubertal_f,pubertal_m],2,'omitnan');
pubertal(pubertal==0)=nan;

House_in=pdem_table{sub_id_BL,'demo_comb_income_v2'};% household income
House_in(House_in>100)=nan;


psm_table=[network_table_BL(sub_id_BL,'interview_age'),...
        dem_table(sub_id_BL,'race4'),dem_table(sub_id_BL,'sex'),dem_table(sub_id_BL,'edu')...
        ,dem_table(sub_id_BL,'abcd_site'),dem_table(sub_id_BL,'rel_family_id'),...
        network_table_BL(sub_id_BL,'rsfmri_c_ngd_ntpoints'),network_table_BL(sub_id_BL,'rsfmri_c_ngd_meanmotion')...
        BMI_BL(sub_id_BL,'anthro_bmi_calc')];% table to be used in propensity score matching

psm_table.Properties.VariableNames={'age','race','sex','prt_edu','sites','family','nframe','FD','BMI'};
psm_table.age_sex=psm_table.age.*psm_table.sex;%  interaction between sex and age
psm_table.puberty=pubertal;
psm_table.HI=House_in;
TS=double(total_sleep_duration==4);% binary sleep duration 1=sufficient sleep 0=insiffucient sleep
psm_table.TS=TS;


% psm_table_all=[network_table_BL(sub_id_all,'interview_age'),...
%         dem_table(sub_id_all,'race4'),dem_table(sub_id_all,'sex'),dem_table(sub_id_all,'edu')...
%         ,dem_table(sub_id_all,'abcd_site'),dem_table(sub_id_all,'rel_family_id'),...
%         network_table_BL(sub_id_all,'rsfmri_c_ngd_ntpoints'),network_table_BL(sub_id_all,'rsfmri_c_ngd_meanmotion')...
%         BMI_BL(sub_id_all,'anthro_bmi_calc')];
% 
% psm_table.Properties.VariableNames={'age','race','sex','prt_edu','sites','family','nframe','FD','BMI'};
% psm_table.age_sex=psm_table.age.*psm_table.sex;
% psm_table.puberty=pubertal;
% psm_table.HI=House_in;
% TS=double(total_sleep_duration==4);
% psm_table.TS=TS;






%% run psm

psm_table_nonan=psm_table(~isnan(sum(psm_table{:,1:12},2)),:);% exclued subs with missing values; n=8323 
psm_table_nonan.num=(1:height(psm_table_nonan))';
writetable(psm_table_nonan,'psm_table_nn_30.csv')
% psm_table_nonan_cryst=psm_table_nonan;
% psm_table_nonan_cryst.crystal=NIH_TB_BL{psm_table_nonan.Properties.RowNames,'nihtbx_cryst_uncorrected'};
% writetable(psm_table_nonan_cryst,'psm_table_nn_30_cryst.csv')% added one outcome variable to run sensitivity analysis


% R code to run psm
% % % library("MatchIt")
% % % data("lalonde", package = "MatchIt")
% % % mydata <- read.csv("/home/zwang/Documents/Nils/ABCD/psm_table_nn_30.csv")
% % % attach(mydata)
% % % mydata[1:10,]
% % % m.out=matchit(TS ~ age + race + sex  +puberty + prt_edu + sites +
% % %                 FD + nframe + BMI +HI,
% % %               data = mydata, method = "nearest",estimand = "ATC",
% % %               caliper = 0.1, ratio = 1)
% % % summary(m.out)
% % % plot(m.out,type = "jitter")
% % % plot(m.out,type = "hist")
% % % 
% % % plot(summary(m.out, subclass = TRUE),
% % %      var.order = "unmatched", abs = FALSE)
% % % 
% % % 
% % % mydata1 <- match.data(m.out)
% % % write.csv(mydata1, file = "/home/zwang/Documents/Nils/ABCD/psm_table_nn_matched30.csv")
modelspec='TS ~ age + race + sex +age_sex +puberty + prt_edu + sites + FD + nframe + BMI +HI';
mdl=fitglm(psm_table_nonan,modelspec,'Distribution','binomial');

% load psm
psm_table_nn_matched=readtable('psm_table_nn_matched30.csv');
psm_table_nn_matched.Properties.RowNames=psm_table_nonan.Properties.RowNames(psm_table_nn_matched.num);
psm_sorted_table=sortrows(psm_table_nn_matched,[14 18]);
sleep_NE=psm_sorted_table.Properties.RowNames(1:height(psm_sorted_table)/2);%sub_id without enough sleep
sleep_FS=psm_sorted_table.Properties.RowNames(height(psm_sorted_table)/2+1:end);%sub_id with full sleep
sleep_psc=[sleep_FS;sleep_NE];% subjects used in PSM study

ind_FL1=ismember(sleep_NE,sub_id_FL1)&ismember(sleep_FS,sub_id_FL1);%matched subject that both have FL1 data
sleep_psc_FL1=[sleep_FS(ind_FL1);sleep_NE(ind_FL1)];
ind_FL2=ismember(sleep_NE,sub_id_FL2)&ismember(sleep_FS,sub_id_FL2);%matched subject that both have FL2 data
sleep_psc_FL2=[sleep_FS(ind_FL2);sleep_NE(ind_FL2)];



total_sleep_duration=sleep_parent_SDS_BL{sleep_psc,'sleepdisturb1_p'};%clc
%1 = 9-11 hours/ 9 a 11 horas; 2 = 8-9 hours /8 a 9 horas; 3 = 7-8 hours /7 a 8 horas; 4 = 5-7 hours /5 a 7 horas; 5 = Less than 5 hours/ Menos de 5 horas// Consider each question pertaining to the PAST 6 MONTHS of the child's life 
total_sleep_duration(total_sleep_duration==5)=4; % combine 4 & 5 as only 26 subjects have a value of 5
total_sleep_duration=5-total_sleep_duration; %reverse sleep duration to 4 are recommed sleep time.






total_sleep_duration_FL1=sleep_parent_SDS_FL1{sleep_psc_FL1,'sleepdisturb1_p'};%1 = 9-11 hours/ 9 a 11 horas; 2 = 8-9 hours /8 a 9 horas; 3 = 7-8 hours /7 a 8 horas; 4 = 5-7 hours /5 a 7 horas; 5 = Less than 5 hours/ Menos de 5 horas// Consider each question pertaining to the PAST 6 MONTHS of the child's life 
total_sleep_duration_FL1(total_sleep_duration_FL1==5)=4; % combine 4 & 5 as only 26 subjects have a value of 5
total_sleep_duration_FL1=5-total_sleep_duration_FL1; %reverse sleep duration to 4 are recommed sleep time.

total_sleep_duration_FL2=sleep_parent_SDS_FL2{sleep_psc_FL2,'sleepdisturb1_p'};%1 = 9-11 hours/ 9 a 11 horas; 2 = 8-9 hours /8 a 9 horas; 3 = 7-8 hours /7 a 8 horas; 4 = 5-7 hours /5 a 7 horas; 5 = Less than 5 hours/ Menos de 5 horas// Consider each question pertaining to the PAST 6 MONTHS of the child's life 
total_sleep_duration_FL2(total_sleep_duration_FL2==5)=4; % combine 4 & 5 as only 26 subjects have a value of 5
total_sleep_duration_FL2=5-total_sleep_duration_FL2; %reverse sleep duration to 4 are recommed sleep time.

TS=double(total_sleep_duration==4);
TS_FL1=double(total_sleep_duration_FL1==4);
TS_FL2=double(total_sleep_duration_FL2==4);



%% independt t test on network connectivities

p_stats=nan(1,size(network_conn_table,2));
cohen_d_stats=nan(1,size(network_conn_table,2));
p_value=nan(1,size(network_conn_table,2));
for i=1:size(network_conn_table,2)
    [h,p]=ttest2(network_conn_table{sleep_FS,i},network_conn_table{sleep_NE,i});
    d=computeCohen_d(network_conn_table{sleep_FS,i},network_conn_table{sleep_NE,i},'independent');
    
    p_stats(i)=-log10(p)*sign(d);
    p_value(i)=p;
    cohen_d_stats(i)=d;
end


[~,~,~,baseline_net.p_values]=fdr_bh(p_value([1:78 92:283 303:end]));% fdr correction
p_stats_shaped=[reshape(p_stats(1:169),13,13);reshape(p_stats(170:end),19,13);];
p_value=[reshape(p_value(1:169),13,13);reshape(p_value(170:end),19,13);];

cohen_d_shaped=[reshape(cohen_d_stats(1:169),13,13);reshape(cohen_d_stats(170:end),19,13);];
cohen_d_shaped=cohen_d_shaped([1:6 8:end],[1:6 8:end]);
p_test=p_value([1:6 8:end],[1:6 8:end]);% excluding none network
p_test=[tril(p_test(1:12,1:12));p_test(13:end,:)];
a=triu(ones(12,12),1)*99;
p_test(1:12,1:12)=p_test(1:12,1:12)+a;
p_test(p_test<99)=fdr_bh(p_test(p_test<99),0.05);
p_test(p_test>1)=0;
cohen_d_shaped= cohen_d_shaped.*p_test;



net_names={'auditory network', 'cingulo-opercular network', 'cingulo-parietal network', 'default network'...
    ,'dorsal attention network', 'fronto-parietal network','retrosplenial temporal network',...
    'sensorimotor hand network','sensorimotor mouth network', 'salience network',...
    'ventral attention network','visual network','left-cerebellum-cortex','left-thalamus-proper',...
    'left-caudate',' left-putamen','left-pallidum','brain-stem','left-hippocampus',...
    'left-amygdala','left-accumbens-area','left-ventraldc','right-cerebellum-cortex',...
    'right-thalamus-proper','right-caudate','right-putamen','right-pallidum', 'right-hippocampus'...
    'right-amygdala','right-accumbens-area','right-ventraldc'};

% create_matrix_ABCD(p_stats_shaped,'Total Sleep Duration',net_names,1)

create_matrix_ABCD_cohensD(cohen_d_shaped,'Network connectiviy baseline',net_names,0) % drew network matrix

baseline_net.p_value=p_value;
baseline_net.cohen_d_stats=cohen_d_stats;
baseline_net.cohen_d_threhold=cohen_d_shaped;

%% ind t test on greymatter volume 463:613 with TIV regressed out

p_stats=nan(1,151);
cohen_d_stats=nan(1,151);
p_value=nan(1,151);
ind=qc_table_BL{sleep_FS,'imgincl_t1w_include'}==1&qc_table_BL{sleep_NE,'imgincl_t1w_include'}==1;%matched subject that both passed qc
for i=1:151
    tempv=T1_table_BL{sleep_psc,i+462};
    ind_nn=~isnan(tempv)&[ind;ind];
    b=regress(tempv(ind_nn),[ones(length(tempv(ind_nn)),1) T1_sub_BL{sleep_psc(ind_nn),'smri_vol_scs_intracranialv'}]);
    Yhat=[ones(length(tempv(ind_nn)),1) T1_sub_BL{sleep_psc(ind_nn),'smri_vol_scs_intracranialv'}]*b;
    resid(ind_nn)=tempv(ind_nn)-Yhat+mean(tempv(ind_nn));
    resid(~ind_nn)=nan;
    [h,p]=ttest2(resid(ismember(sleep_psc,sleep_FS)),resid(ismember(sleep_psc,sleep_NE)));
    d=computeCohen_d(resid(ismember(sleep_psc,sleep_FS)),resid(ismember(sleep_psc,sleep_NE)),'independent');

    p_stats(i)=-log10(p)*sign(d);
    p_value(i)=p;
%     cohen_d_stats(i)=abs(d);
        cohen_d_stats(i)=d;
end


table_GMV_TIV=T1_table_BL(1,463:613);
table_GMV_TIV{1,:}=p_stats;
table_GMV_TIV{2,:}=cohen_d_stats;
table_GMV_TIV{3,:}=p_value;
% table_GMV_TIV.Properties.VariableNames(table_GMV_TIV{2,1:148}>0.09)'
% table_GMV_TIV(2,table_GMV_TIV{2,1:148}>0.09)

%% ind t test on thickness 10:160 with TIV regressed out

p_stats=nan(1,151);
cohen_d_stats=nan(1,151);
p_value=nan(1,151);
ind=qc_table_BL{sleep_FS,'imgincl_t1w_include'}==1&qc_table_BL{sleep_NE,'imgincl_t1w_include'}==1;%matched subject that both passed qc

for i=1:151
    tempv=T1_table_BL{sleep_psc,i+9};
    ind_nn=~isnan(tempv)&[ind;ind];
    b=regress(tempv(ind_nn),[ones(length(tempv(ind_nn)),1) T1_sub_BL{sleep_psc(ind_nn),'smri_vol_scs_intracranialv'}]);
    Yhat=[ones(length(tempv(ind_nn)),1) T1_sub_BL{sleep_psc(ind_nn),'smri_vol_scs_intracranialv'}]*b;
    resid(ind_nn)=tempv(ind_nn)-Yhat+mean(tempv(ind_nn));
    resid(~ind_nn)=nan;
    [~,p]=ttest2(resid(ismember(sleep_psc,sleep_FS)),resid(ismember(sleep_psc,sleep_NE)));
    d=computeCohen_d(resid(ismember(sleep_psc,sleep_FS)),resid(ismember(sleep_psc,sleep_NE)),'independent');

    
    p_stats(i)=-log10(p)*sign(d);
    p_value(i)=p;
    cohen_d_stats(i)=d;
    
end


table_thk_TIV=T1_table_BL(1,10:160);
table_thk_TIV{1,:}=p_stats;
table_thk_TIV{2,:}=cohen_d_stats;% left/right lingual gyrus 
table_thk_TIV{3,:}=p_value;%

% table_thk_TIV.Properties.VariableNames(table_thk_TIV{2,1:148}>0.07)'
% table_thk_TIV(2,table_thk_TIV{2,1:148}>0.07)

%% ind t test on subcortical volume 330:375 with TIV regressed out 

p_stats=nan(1,46);
cohen_d_stats=nan(1,46);
p_value_sub=nan(1,46);
% ind=qc_table_BL{sleep_FS,'imgincl_t1w_include'}==1&qc_table_BL{sleep_NE,'imgincl_t1w_include'}==1;%matched subject that both passed qc
for i=1:46
       try
    tempv=T1_sub_BL{sleep_psc,i+329};
    ind_nn=~isnan(tempv)&[ind;ind];
    b=regress(tempv(ind_nn),[ones(length(tempv(ind_nn)),1) T1_sub_BL{sleep_psc(ind_nn),'smri_vol_scs_intracranialv'}]);
    Yhat=[ones(length(tempv(ind_nn)),1) T1_sub_BL{sleep_psc(ind_nn),'smri_vol_scs_intracranialv'}]*b;
    resid(ind_nn)=tempv(ind_nn)-Yhat+mean(tempv(ind_nn));
    resid(~ind_nn)=nan;
    [h,p]=ttest2(resid(ismember(sleep_psc,sleep_FS)),resid(ismember(sleep_psc,sleep_NE)));
    d=computeCohen_d(resid(ismember(sleep_psc,sleep_FS)),resid(ismember(sleep_psc,sleep_NE)),'independent');

    
    p_stats(i)=-log10(p)*sign(d);
      p_value_sub(i)=p;
   cohen_d_stats(i)=d;
        catch me
        disp(me)
        p_stats(i)=0;
        cohen_d_stats(i)=0;
            p_value_sub(i)=1;
    end
end


table_subvol_TIV=T1_sub_BL(1,330:375);
table_subvol_TIV{1,:}=p_stats;
table_subvol_TIV{2,:}=cohen_d_stats;
table_subvol_TIV{3,:}=p_value_sub;


% table_subvol_TIV.Properties.VariableNames(table_subvol_TIV{2,1:46}>0.06)'
% table_subvol_TIV(2,table_subvol_TIV{2,1:46}>0.06)

%% ind t test on cortical area 312:462 with TIV regressed out

p_stats=nan(1,151);
cohen_d_stats=nan(1,151);
p_value=nan(1,151);
ind=qc_table_BL{sleep_FS,'imgincl_t1w_include'}==1&qc_table_BL{sleep_NE,'imgincl_t1w_include'}==1;%matched subject that both passed qc

for i=1:151
    tempv=T1_table_BL{sleep_psc,i+311};
    ind_nn=~isnan(tempv)&[ind;ind];
    b=regress(tempv(ind_nn),[ones(length(tempv(ind_nn)),1) T1_sub_BL{sleep_psc(ind_nn),'smri_vol_scs_intracranialv'}]);
    Yhat=[ones(length(tempv(ind_nn)),1) T1_sub_BL{sleep_psc(ind_nn),'smri_vol_scs_intracranialv'}]*b;
    resid(ind_nn)=tempv(ind_nn)-Yhat+mean(tempv(ind_nn));
    resid(~ind_nn)=nan;
    [h,p]=ttest2(resid(ismember(sleep_psc,sleep_FS)),resid(ismember(sleep_psc,sleep_NE)));
    d=computeCohen_d(resid(ismember(sleep_psc,sleep_FS)),resid(ismember(sleep_psc,sleep_NE)),'independent');

    p_value(i)=p;
    p_stats(i)=-log10(p)*sign(d);
     cohen_d_stats(i)=d;
end


table_area_TIV=T1_table_BL(1,312:462);
table_area_TIV{1,:}=p_stats;
table_area_TIV{2,:}=cohen_d_stats;% left/right temporal pole
table_area_TIV{3,:}=p_value;%
% table_area_TIV.Properties.VariableNames(table_area_TIV{2,1:148}>0.09)'
% table_area_TIV(2,table_area_TIV{2,1:148}>0.09)

%% ind t test on DTI FA 10:51

p_stats=nan(1,42);
cohen_d_stats=nan(1,42);
p_value=nan(1,42);
ind=qc_table_BL{sleep_FS,'imgincl_dmri_include'}==1&qc_table_BL{sleep_NE,'imgincl_dmri_include'}==1;%matched subject that both passed qc

for i=1:42
    try
    [h,p]=ttest2(DTI_table_BL{sleep_FS(ind),i+9},DTI_table_BL{sleep_NE(ind),i+9});
    d=computeCohen_d(DTI_table_BL{sleep_FS(ind),i+9},DTI_table_BL{sleep_NE(ind),i+9},'independent');
    
    p_stats(i)=-log10(p)*sign(d);
    p_value(i)=p;
     cohen_d_stats(i)=d;
    catch me
        disp(me)
        p_stats(i)=0;
        cohen_d_stats(i)=0;
        p_value(i)=nan;
    end
end


table_DTI_FA=DTI_table_BL(1,10:51);
table_DTI_FA{1,:}=p_stats;
table_DTI_FA{2,:}=cohen_d_stats;%  left superior corticostriate and 
table_DTI_FA{3,:}=p_value;%

% table_DTI_FA.Properties.VariableNames(table_DTI_FA{2,1:37}>0.06)'
% table_DTI_FA(2,table_DTI_FA{2,1:37}>0.07)
%% behavior baseline
cbcl_raw_BL=cbcl_sum_BL(sleep_psc,10:4:86);
cbcl_raw_BL=cbcl_raw_BL(:,sort(cbcl_raw_BL.Properties.VariableNames));

cbcl_raw_BL=cbcl_raw_BL(:,[1:18 20 19]);% put total pro to the last one

cognition_row_BL=NIH_TB_BL(sleep_psc,[12:5:42 46:4:54]);
mhy_row_BL=mental_health_BL(sleep_psc,[28 31 37 41 44 47 50 53 56:3:65]);
screentime_row_BL=screen_time_BL(sleep_psc,[10:21]);
% cbcl_raw_BL=[cbcl_raw_BL,cognition_row_BL,mhy_row_BL,screentime_row_BL];
cbcl_raw_BL=[cbcl_raw_BL,cognition_row_BL,mhy_row_BL];


cohend=nan(1,width(cbcl_raw_BL));
p_value=nan(1,width(cbcl_raw_BL));
for i=1:width(cbcl_raw_BL)
    cohend(i)=computeCohen_d(cbcl_raw_BL{sleep_FS,i},cbcl_raw_BL{sleep_NE,i},'independent');
     [h,p_value(i)]=ttest2(cbcl_raw_BL{sleep_FS,i},cbcl_raw_BL{sleep_NE,i});
    
end
% [h,p_value_dep,ci]=ttest2(cbcl_raw_BL{sleep_FS,7},cbcl_raw_BL{sleep_NE,7});
 %[h,p_value_cc,ci]=ttest2(cbcl_raw_BL{sleep_FS,29},cbcl_raw_BL{sleep_NE,29});
p_value=-log10(p_value);
[~,~,~,beh_fdr_p_BL]=fdr_bh(p_value);
beh_fdr_p_BL(beh_fdr_p_BL<0.0001)=0.0001;

figure;
b=bar(cohend(1:42),'FaceColor','flat');

for i=21:21+width(cognition_row_BL)-1
b.CData(i,:) = [0 1 0];
end

for i=21+width(cognition_row_BL):21+width(cognition_row_BL)-1+width(mhy_row_BL)
b.CData(i,:) = [1 1 0];
end
% for i=21+width(cognition_row_BL)+width(mhy_row_BL):21+width(cognition_row_BL)-1+width(mhy_row_BL)+width(screentime_row_BL)
% b.CData(i,:) = [1 0 0];
% end
hold on
yline([-0.15 0.15],'--')
set(gca, 'FontName', 'Helvetica')
xticks(1:42);
xticklabels(strrep(cbcl_raw_BL.Properties.VariableNames(1:42),'_','\_'))
ylabel({'cohen''s d'},'FontName', 'Helvetica')
title({'behavioral variables'})
set(gca,'TickLength',[0 0])
set(gcf,'units','inches','position',[0, 0, 15, 6])
box off







%% follow up 1 




cbcl_raw_FL1=cbcl_sum_FL1(sleep_psc_FL1,10:4:86);
cbcl_raw_FL1=cbcl_raw_FL1(:,sort(cbcl_raw_FL1.Properties.VariableNames));

cbcl_raw_FL1=cbcl_raw_FL1(:,[1:18 20 19]);% put total pro to the last one


% cognition_row_FL1=NIH_TB_FL1(sleep_psc_FL1,[12:5:42 46:4:54]);
% mhy_row_FL1=mental_health_FL1(sleep_psc_FL1,[28 31 37 41 44 47 50 53 56:3:68]);
screentime_row_FL1=screen_time_FL1(sleep_psc_FL1,[10:21]);



cbcl_raw_FL1=[cbcl_raw_FL1,screentime_row_FL1];



cohend_FL1=nan(1,width(cbcl_raw_FL1));
p_value=nan(1,width(cbcl_raw_FL1));
for i=1:width(cbcl_raw_FL1)
    cohend_FL1(i)=computeCohen_d(cbcl_raw_FL1{sleep_FS(ind_FL1),i},cbcl_raw_FL1{sleep_NE(ind_FL1),i},'independent');
     [h,p_value(i)]=ttest2(cbcl_raw_FL1{sleep_FS(ind_FL1),i},cbcl_raw_FL1{sleep_NE(ind_FL1),i});
    
end
p_value=-log10(p_value);



set(gca, 'FontName', 'Helvetica')
xticks(1:width(cbcl_raw_FL1));
xtickangle(90)
xticklabels(strrep(cbcl_raw_FL1.Properties.VariableNames(1:32),'_','\_'))
ylabel({'cohen''s d'},'FontName', 'Helvetica')
title({'behavioral variables FL1'})
set(gca,'TickLength',[0 0])
set(gcf,'units','inches','position',[0, 0, 15, 6])
box off




%% follow up 2 




cbcl_raw_FL2=cbcl_sum_FL2(sleep_psc_FL2,10:4:86);
cbcl_raw_FL2=cbcl_raw_FL2(:,sort(cbcl_raw_FL2.Properties.VariableNames));

cbcl_raw_FL2=cbcl_raw_FL2(:,[1:18 20 19]);% put total pro to the last one


cognition_row_FL2=NIH_TB_FL2(sleep_psc_FL2,[12:5:42 46:4:54]);
mhy_row_FL2=mental_health_FL2(sleep_psc_FL2,[28 31 37 41 44 47 50 53 56:3:65]);
% screentime_row_FL2=screen_time_FL2(sleep_psc_FL2,[10:21]);
% 
% 
% 
cbcl_raw_FL2=[cbcl_raw_FL2,cognition_row_FL2,mhy_row_FL2];



cohend_FL2=nan(1,width(cbcl_raw_FL2));  
p_value=nan(1,width(cbcl_raw_FL2));
for i=1:width(cbcl_raw_FL2)
    if mean(isnan(cbcl_raw_FL2{sleep_FS(ind_FL2),i}))>0.9
          cohend_FL2(i)=0;
        p_value(i)=1;
    else
    cohend_FL2(i)=computeCohen_d(cbcl_raw_FL2{sleep_FS(ind_FL2),i},cbcl_raw_FL2{sleep_NE(ind_FL2),i},'independent');
     [h,p_value(i)]=ttest2(cbcl_raw_FL2{sleep_FS(ind_FL2),i},cbcl_raw_FL2{sleep_NE(ind_FL2),i});
    
      
    end
end
p_value=-log10(p_value);

[~,~,~,beh_fdr_p_FL2]=fdr_bh(p_value);%FDR correction
beh_fdr_p_FL2(beh_fdr_p_FL2<0.0001)=0.0001;


cohend_combined=[cohend_FL2(1:20),cohend_FL2(21:30),cohend_FL2(31:42)];
figure;
b=bar(cohend_combined,'FaceColor','flat');
for i=21:21+width(cognition_row_BL)-1
b.CData(i,:) = [0 1 0];
end

for i=21+width(cognition_row_BL):21+width(cognition_row_BL)-1+width(mhy_row_BL)
b.CData(i,:) = [1 1 0];
end
% for i=21+width(cognition_row_BL)+width(mhy_row_BL):21+width(cognition_row_BL)-1+width(mhy_row_BL)+width(screentime_row_BL)
% b.CData(i,:) = [1 0 0];
% end
hold on
yline([-0.15 0.15],'--')

% xticks(1:width(cbcl_raw_BL));
set(gca, 'FontName', 'Helvetica')
xticks(1:width(cohend_combined));
xticklabels(strrep(cbcl_raw_BL.Properties.VariableNames(1:43),'_','\_'))
ylabel({'cohen''s d'},'FontName', 'Helvetica')


title({'behavioral variables follow-up 2'})
set(gcf,'units','inches','position',[0,0,15,6])
set(gca,'TickLength',[0 0])
xtickangle(45)
box off






%% network connectivities FL2
sleep_FS_FL2=sleep_FS(ind_FL2);
sleep_NE_FL2=sleep_NE(ind_FL2);
p_stats=nan(1,size(network_conn_table_FL2,2));
cohen_d_stats=nan(1,size(network_conn_table_FL2,2));
p_value=nan(1,size(network_conn_table_FL2,2));

ind=qc_table_FL2{sleep_FS(ind_FL2),'imgincl_rsfmri_include'}==1&qc_table_FL2{sleep_NE(ind_FL2),'imgincl_rsfmri_include'}==1;%matched subject that both passed qc

for i=1:size(network_conn_table_FL2,2)
    [h,p]=ttest2(network_conn_table_FL2{sleep_FS_FL2(ind),i},network_conn_table_FL2{sleep_NE_FL2(ind),i});
    d=computeCohen_d(network_conn_table_FL2{sleep_FS_FL2(ind),i},network_conn_table_FL2{sleep_NE_FL2(ind),i},'independent');
    
    p_stats(i)=-log10(p)*sign(d);
    p_value(i)=p;
    cohen_d_stats(i)=d;
end
p_stats_shaped=[reshape(p_stats(1:169),13,13);reshape(p_stats(170:end),19,13);];
cohen_d_shaped=[reshape(cohen_d_stats(1:169),13,13);reshape(cohen_d_stats(170:end),19,13);];
cohen_d_shaped=cohen_d_shaped([1:6 8:end],[1:6 8:end]);
p_value=[reshape(p_value(1:169),13,13);reshape(p_value(170:end),19,13);];
p_value=p_value([1:6 8:end],[1:6 8:end]);
p_value=[tril(p_value(1:12,1:12));p_value(13:end,:)];
% p_test=p_value;
% p_test(p_value>0)=fdr_bh(p_value(p_value>0),0.05);
cohen_d_shaped= cohen_d_shaped.*p_test;

net_names={'auditory network', 'cingulo-opercular network', 'cingulo-parietal network', 'default network'...
    ,'dorsal attention network', 'fronto-parietal network','retrosplenial temporal network',...
    'sensorimotor hand network','sensorimotor mouth network', 'salience network',...
    'ventral attention network','visual network','left-cerebellum-cortex','left-thalamus-proper',...
    'left-caudate',' left-putamen','left-pallidum','brain-stem','left-hippocampus',...
    'left-amygdala','left-accumbens-area','left-ventraldc','right-cerebellum-cortex',...
    'right-thalamus-proper','right-caudate','right-putamen','right-pallidum', 'right-hippocampus'...
    'right-amygdala','right-accumbens-area','right-ventraldc'};

% create_matrix_ABCD(p_stats_shaped,'Total Sleep Duration FL2',net_names,0)

create_matrix_ABCD_cohensD(cohen_d_shaped,'Total Sleep Duration Cohen''s D FL2',net_names,0) 



FL_net.p_value=p_value;
FL_net.cohen_d_stats=cohen_d_stats;
FL_net.cohen_d_threhold=cohen_d_shaped;

%% independent t test on greymatter volume 463:613 with TIV regressed out FL2

p_stats=nan(1,151);
cohen_d_stats=nan(1,151);
p_values=nan(1,151);

sleep_temp_FL2=sleep_psc_FL2(ismember(sleep_psc_FL2,T1_table_FL2.Properties.RowNames));
ind_FL2_gmv=ismember(sleep_NE,sleep_temp_FL2)&ismember(sleep_FS,sleep_temp_FL2);%matched subject that both have FL2 data
sleep_temp_FL2=[sleep_FS(ind_FL2_gmv);sleep_NE(ind_FL2_gmv)];
ind=qc_table_FL2{sleep_FS(ind_FL2_gmv),'imgincl_t1w_include'}==1&qc_table_FL2{sleep_NE(ind_FL2_gmv),'imgincl_t1w_include'}==1;%matched subject that both passed qc

for i=1:151
    tempv=T1_table_FL2{sleep_temp_FL2,i+462};
    ind_nn=~isnan(tempv)&[ind;ind];
    b=regress(tempv(ind_nn),[ones(length(tempv(ind_nn)),1) T1_sub_FL2{sleep_temp_FL2(ind_nn),'smri_vol_scs_intracranialv'}]);
    Yhat=[ones(length(tempv(ind_nn)),1) T1_sub_FL2{sleep_temp_FL2(ind_nn),'smri_vol_scs_intracranialv'}]*b;
    resid(ind_nn)=tempv(ind_nn)-Yhat+mean(tempv(ind_nn));
    resid(~ind_nn)=nan;
    [h,p]=ttest2(resid(ismember(sleep_temp_FL2,sleep_FS)),resid(ismember(sleep_temp_FL2,sleep_NE)));
    d=computeCohen_d(resid(ismember(sleep_temp_FL2,sleep_FS)),resid(ismember(sleep_temp_FL2,sleep_NE)),'independent');

    p_stats(i)=-log10(p)*sign(d);
    p_values(i)=p;
     cohen_d_stats(i)=d;
end


table_GMV_TIV_FL2=T1_table_BL(1,463:613);
table_GMV_TIV_FL2{1,:}=p_stats;
table_GMV_TIV_FL2{2,:}=cohen_d_stats;
table_GMV_TIV_FL2{3,:}=p_values;







%% independent t test on thickness 10:160 with TIV regressed out FL2

p_stats=nan(1,151);
p_values=nan(1,151);
cohen_d_stats=nan(1,151);
sleep_temp_FL2=sleep_psc_FL2(ismember(sleep_psc_FL2,T1_table_FL2.Properties.RowNames));
ind_FL2_gmv=ismember(sleep_NE,sleep_temp_FL2)&ismember(sleep_FS,sleep_temp_FL2);%matched subject that both have FL2 data
sleep_temp_FL2=[sleep_FS(ind_FL2_gmv);sleep_NE(ind_FL2_gmv)];
ind=qc_table_FL2{sleep_FS(ind_FL2_gmv),'imgincl_t1w_include'}==1&qc_table_FL2{sleep_NE(ind_FL2_gmv),'imgincl_t1w_include'}==1;%matched subject that both passed qc

for i=1:151
    tempv=T1_table_FL2{sleep_temp_FL2,i+9};
    ind_nn=~isnan(tempv)&[ind;ind];
    b=regress(tempv(ind_nn),[ones(length(tempv(ind_nn)),1) T1_sub_FL2{sleep_temp_FL2(ind_nn),'smri_vol_scs_intracranialv'}]);
    Yhat=[ones(length(tempv(ind_nn)),1) T1_sub_FL2{sleep_temp_FL2(ind_nn),'smri_vol_scs_intracranialv'}]*b;
    resid(ind_nn)=tempv(ind_nn)-Yhat+mean(tempv(ind_nn));
    resid(~ind_nn)=nan;
    [h,p]=ttest2(resid(ismember(sleep_temp_FL2,sleep_FS)),resid(ismember(sleep_temp_FL2,sleep_NE)));
    d=computeCohen_d(resid(ismember(sleep_temp_FL2,sleep_FS)),resid(ismember(sleep_temp_FL2,sleep_NE)),'independent');

    
    p_stats(i)=-log10(p)*sign(d);
    p_values(i)=p;
 cohen_d_stats(i)=d;
end


table_thk_TIV_FL2=T1_table_BL(1,10:160);
table_thk_TIV_FL2{1,:}=p_stats;
table_thk_TIV_FL2{2,:}=cohen_d_stats;% left/right lingual gyrus 
table_thk_TIV_FL2{3,:}=p_values;

%% independent t test on subcortical volume 330:375 with TIV regressed out FL2

p_stats=nan(1,46);
cohen_d_stats=nan(1,46);
p_values=nan(1,46);
sleep_temp_FL2=sleep_psc_FL2(ismember(sleep_psc_FL2,T1_table_FL2.Properties.RowNames));
ind_FL2_gmv=ismember(sleep_NE,sleep_temp_FL2)&ismember(sleep_FS,sleep_temp_FL2);%matched subject that both have FL2 data
sleep_temp_FL2=[sleep_FS(ind_FL2_gmv);sleep_NE(ind_FL2_gmv)];
ind=qc_table_FL2{sleep_FS(ind_FL2_gmv),'imgincl_t1w_include'}==1&qc_table_FL2{sleep_NE(ind_FL2_gmv),'imgincl_t1w_include'}==1;%matched subject that both passed qc

for i=1:46
       try
    tempv=T1_sub_FL2{sleep_temp_FL2,i+329};
    ind_nn=~isnan(tempv)&[ind;ind];
    b=regress(tempv(ind_nn),[ones(length(tempv(ind_nn)),1) T1_sub_FL2{sleep_temp_FL2(ind_nn),'smri_vol_scs_intracranialv'}]);
    Yhat=[ones(length(tempv(ind_nn)),1) T1_sub_FL2{sleep_temp_FL2(ind_nn),'smri_vol_scs_intracranialv'}]*b;
    resid(ind_nn)=tempv(ind_nn)-Yhat+mean(tempv(ind_nn));
    resid(~ind_nn)=nan;
    [h,p]=ttest2(resid(ismember(sleep_temp_FL2,sleep_FS)),resid(ismember(sleep_temp_FL2,sleep_NE)));
    d=computeCohen_d(resid(ismember(sleep_temp_FL2,sleep_FS)),resid(ismember(sleep_temp_FL2,sleep_NE)),'independent');

    
    p_stats(i)=-log10(p)*sign(d);
    p_values(i)=p;
    
     cohen_d_stats(i)=d;
        catch me
        disp(me)
        p_stats(i)=0;
        cohen_d_stats(i)=0;
        p_values(i)=1;
    end
end


table_subvol_TIV_FL2=T1_sub_BL(1,330:375);
table_subvol_TIV_FL2{1,:}=p_stats;
table_subvol_TIV_FL2{2,:}=cohen_d_stats;
table_subvol_TIV_FL2{3,:}=p_values;
% table_subvol_TIV_FL2.Properties.VariableNames(table_subvol_TIV_FL2{2,1:40}>0.08)'
% table_thk_TIV(2,table_thk_TIV{2,1:148}>0.09








%% independent t test on cortical area 312:462 with TIV regressed out FL2

p_stats=nan(1,151);
cohen_d_stats=nan(1,151);
p_values=nan(1,151);
sleep_temp_FL2=sleep_psc_FL2(ismember(sleep_psc_FL2,T1_table_FL2.Properties.RowNames));
ind_FL2_gmv=ismember(sleep_NE,sleep_temp_FL2)&ismember(sleep_FS,sleep_temp_FL2);%matched subject that both have FL2 data
sleep_temp_FL2=[sleep_FS(ind_FL2_gmv);sleep_NE(ind_FL2_gmv)];
ind=qc_table_FL2{sleep_FS(ind_FL2_gmv),'imgincl_t1w_include'}==1&qc_table_FL2{sleep_NE(ind_FL2_gmv),'imgincl_t1w_include'}==1;%matched subject that both passed qc

for i=1:151
    tempv=T1_table_FL2{sleep_temp_FL2,i+311};
    ind_nn=~isnan(tempv)&[ind;ind];
    b=regress(tempv(ind_nn),[ones(length(tempv(ind_nn)),1) T1_sub_FL2{sleep_temp_FL2(ind_nn),'smri_vol_scs_intracranialv'}]);
    Yhat=[ones(length(tempv(ind_nn)),1) T1_sub_FL2{sleep_temp_FL2(ind_nn),'smri_vol_scs_intracranialv'}]*b;
    resid(ind_nn)=tempv(ind_nn)-Yhat+mean(tempv(ind_nn));
    resid(~ind_nn)=nan;
    [h,p]=ttest2(resid(ismember(sleep_temp_FL2,sleep_FS)),resid(ismember(sleep_temp_FL2,sleep_NE)));
    d=computeCohen_d(resid(ismember(sleep_temp_FL2,sleep_FS)),resid(ismember(sleep_temp_FL2,sleep_NE)),'independent');

    
    p_stats(i)=-log10(p)*sign(d);
    p_values(i)=p;
     cohen_d_stats(i)=d;;
end


table_area_TIV_FL2=T1_table_BL(1,312:462);
table_area_TIV_FL2{1,:}=p_stats;
table_area_TIV_FL2{2,:}=cohen_d_stats;% left/right temporal pole
table_area_TIV_FL2{3,:}=p_values;


%% independent t test on DTI FA 10:51 FL2

p_stats=nan(1,42);
cohen_d_stats=nan(1,42);
p_values=nan(1,42);
sleep_temp_FL2=sleep_psc_FL2(ismember(sleep_psc_FL2,DTI_table_FL2.Properties.RowNames));
ind_FL2_gmv=ismember(sleep_NE,sleep_temp_FL2)&ismember(sleep_FS,sleep_temp_FL2);%matched subject that both have FL2 data
sleep_temp_FL2=[sleep_FS(ind_FL2_gmv);sleep_NE(ind_FL2_gmv)];
ind=qc_table_FL2{sleep_FS(ind_FL2_gmv),'imgincl_dmri_include'}==1&qc_table_FL2{sleep_NE(ind_FL2_gmv),'imgincl_dmri_include'}==1;%matched subject that both passed qc
sleep_temp_FL2=sleep_temp_FL2([ind;ind]);

for i=1:42
    try
    [h,p]=ttest2(DTI_table_FL2{sleep_FS(ismember(sleep_FS,sleep_temp_FL2)),i+9},DTI_table_FL2{sleep_NE(ismember(sleep_FS,sleep_temp_FL2)),i+9});
    d=computeCohen_d(DTI_table_FL2{sleep_FS(ismember(sleep_FS,sleep_temp_FL2)),i+9},DTI_table_FL2{sleep_NE(ismember(sleep_FS,sleep_temp_FL2)),i+9},'independent');
    
    p_stats(i)=-log10(p)*sign(d);
    p_values(i)=p;
     cohen_d_stats(i)=d;
    catch me
        disp(me)
        p_stats(i)=0;
        cohen_d_stats(i)=0;
        p_values(i)=1;
    end
end


table_DTI_FA_FL2=DTI_table_FL2(1,10:51);
table_DTI_FA_FL2{1,:}=p_stats;
table_DTI_FA_FL2{2,:}=cohen_d_stats;%  left superior corticostriate and 
table_DTI_FA_FL2{3,:}=p_values;




%%  mediation analysis
clearvars temp
temp(:,1:3)=T1_table_BL{sleep_psc,[496 533 570]+9};%temporal poles +ACC volume
temp(:,4:6)=T1_table_BL{sleep_psc,[345 382 419]+9};% temporal pole +ACC cortical area
temp(:,7:9)=T1_table_BL{sleep_psc,[19 22 81]+9};% mofc lingual gyrus mcc
temp(:,10:12)=DTI_table_BL{sleep_psc,[14 27 31]};%cortio-striatal track left superior corticostriate-parietal cortex only left inferior longitudinal fasiculus
for i=1:9 % regress out TIV
    tempv=temp(:,i);
ind_nn=~isnan(tempv);
b=regress(tempv(ind_nn),[ones(length(tempv(ind_nn)),1) T1_sub_BL{sleep_psc(ind_nn),'smri_vol_scs_intracranialv'}]);
Yhat=[ones(length(tempv(ind_nn)),1) T1_sub_BL{sleep_psc(ind_nn),'smri_vol_scs_intracranialv'}]*b;
resid(ind_nn)=tempv(ind_nn)-Yhat+mean(tempv(ind_nn));
resid(~ind_nn)=nan;
temp(:,i)=resid;
end


ind=1:416;
ind=ind(baseline_net.cohen_d_stats>0.15);
ind=ind([1:3 5:6]);% network measurements that have effect size higher than 0.15, excluded one none network





mediators=[network_conn_table{sleep_psc,ind},temp];%mediators





% full behavior mediation
ind_t1=qc_table_BL{sleep_FS,'imgincl_t1w_include'}==1&qc_table_BL{sleep_NE,'imgincl_t1w_include'}==1;%matched subject that both passed qc

p_med=nan(length(ind)+6,42);
for i=1:length(ind)
    for j=1:42
    [paths, stats_med]=mediation(psm_table_nn_matched{sleep_psc,'TS'},cbcl_raw_BL{sleep_psc,j},...
    mediators(:,i)...
    ,'cov', psm_table_nn_matched{sleep_psc,2:13},'verbose','boot','bootsamples',100000);
    p_med(i,j)=stats_med.p(5);
    end
end
for i=length(ind)+1:length(ind)+6 % sturctural variables only run in subs passed t1w qc
    for j=1:42
    [paths, stats_med]=mediation(psm_table_nn_matched{sleep_psc([ind_t1;ind_t1]),'TS'},cbcl_raw_BL{sleep_psc([ind_t1;ind_t1]),j},...
    mediators([ind_t1;ind_t1],i)...
    ,'cov', psm_table_nn_matched{sleep_psc([ind_t1;ind_t1]),2:13},'verbose','boot','bootsamples',100000);
    p_med(i,j)=stats_med.p(5);
    end
end

p_med5=p_med;
p_med5(p_med5>0.05)=1;
p_med_log=-log10(p_med5);
imagesc(p_med_log(1:11,:))
title('Baseline','FontSize',18);
set(gcf,'Position',[100 100 800 800])
set(gca, 'XTick', 1:42); % center x-axis ticks on bins
set(gca, 'YTick', 1:26); % center y-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel',strrep(cbcl_raw_BL.Properties.VariableNames(1:42),'_','\_')); % set x-axis labels
set(gca, 'YTickLabel', [strrep(network_conn_table.Properties.VariableNames(ind),'_','\_')...
    ,{'Volume-left temporal pole'},{'Volume-right ACC'},{'Volume-right temporal pole'},...
    {'Area-left temporal pole'},{'Area-right ACC'},{'Area-right temporal pole'}...
    ]); % set y-axis labels
set(gcf,'units','inches','position',[0,0,20,15])
c = colorbar;
c.Limits= [0 5];
c.TickLength = [0];
c.Label.String = '-log10(p) of mediation effect';
c.FontSize=16;
colormap(bluewhitered)

% baseline 5 mediation 
p_med5=p_med;
p_med5(p_med5>0.05)=1;
p_med_log5=-log10(p_med5);
imagesc(p_med_log5(1:11,[ 7 18 21 29]))
title('Baseline','FontSize',18);
set(gcf,'Position',[100 100 800 800])
set(gca, 'XTick', 1:5); % center x-axis ticks on bins
set(gca, 'YTick', 1:26); % center y-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel',{'thought','depress','stress','picture vocabulary','crystal cognition'}); % set x-axis labels
set(gca, 'YTickLabel', [strrep(network_conn_table.Properties.VariableNames(ind),'_','\_')...
    ,{'Volume - left temporal pole'},{'Volume - right ACC'},{'Volume - right temporal pole'},...
    {'Area - left temporal pole'},{'Area - right ACC'},{'Area - right temporal pole'}...
    ]); % set y-axis labels
set(gcf,'units','inches','position',[0,0,9,9])
c = colorbar;
c.Limits= [0 5];
c.TickLength = [0];
c.Label.String = '-log10(p) of mediation effect';
c.FontSize=16;
colormap(bluewhitered)


% full behavior mediation fl2
p_med_FL=nan(length(ind)+12,42);%
ind_fmri_fl=qc_table_FL2{sleep_FS_FL2,'imgincl_rsfmri_include'}==1&qc_table_FL2{sleep_NE_FL2,'imgincl_rsfmri_include'}==1;%matched subject that both passed qc

for i=1:length(ind)
    for j=[1:22 25:27 29 31:42]
     
    [paths, stats_med]=mediation(psm_table_nn_matched{sleep_psc_FL2([ind_fmri_fl;ind_fmri_fl]),'TS'},...
        cbcl_raw_FL2{sleep_psc_FL2([ind_fmri_fl;ind_fmri_fl]),j},...
    mediators(ismember(sleep_psc,sleep_psc_FL2([ind_fmri_fl;ind_fmri_fl])),i)...
    ,'cov', [psm_table_nn_matched{sleep_psc_FL2([ind_fmri_fl;ind_fmri_fl]),2:13},cbcl_raw_BL{sleep_psc_FL2([ind_fmri_fl;ind_fmri_fl]),j}]...
    ,'verbose','boot','bootsamples',100000);
    p_med_FL(i,j)=stats_med.p(5);
    end
end
sleep_temp_FL2=sleep_psc_FL2(ismember(sleep_psc_FL2,T1_table_FL2.Properties.RowNames));
ind_FL2_gmv=ismember(sleep_NE,sleep_temp_FL2)&ismember(sleep_FS,sleep_temp_FL2);%matched subject that both have FL2 data
sleep_temp_FL2=[sleep_FS(ind_FL2_gmv);sleep_NE(ind_FL2_gmv)];
ind_t1_fl=qc_table_FL2{sleep_FS(ind_FL2_gmv),'imgincl_t1w_include'}==1&qc_table_FL2{sleep_NE(ind_FL2_gmv),'imgincl_t1w_include'}==1;%matched subject that both passed qc

for i=length(ind)+1:length(ind)+6
    for j=[1:22 25:27 29 31:42]
     
    [paths, stats_med]=mediation(psm_table_nn_matched{sleep_temp_FL2([ind_t1_fl;ind_t1_fl]),'TS'},...
        cbcl_raw_FL2{sleep_temp_FL2([ind_t1_fl;ind_t1_fl]),j},...
    mediators(ismember(sleep_psc,sleep_temp_FL2([ind_t1_fl;ind_t1_fl])),i)...
    ,'cov', [psm_table_nn_matched{sleep_temp_FL2([ind_t1_fl;ind_t1_fl]),2:13},cbcl_raw_BL{sleep_temp_FL2([ind_t1_fl;ind_t1_fl]),j}]...
    ,'verbose','boot','bootsamples',100000);
    p_med_FL(i,j)=stats_med.p(5);
    end
end

p_med5=p_med_FL;
p_med5(p_med5>0.05)=1;
p_med_log=-log10(p_med5);
imagesc(p_med_log(1:11,:))
title('FL2','FontSize',18);
set(gcf,'Position',[100 100 800 800])
set(gca, 'XTick', 1:42); % center x-axis ticks on bins
set(gca, 'YTick', 1:26); % center y-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel',strrep(cbcl_raw_BL.Properties.VariableNames(1:42),'_','\_')); % set x-axis labels
set(gca, 'YTickLabel', [strrep(network_conn_table.Properties.VariableNames(ind),'_','\_')...
    ,{'Volume - left temporal pole'},{'Volume - right ACC'},{'Volume - right temporal pole'},...
    {'Area - left temporal pole'},{'Area - right ACC'},{'Area - right temporal pole'}...
    ]); % set y-axis labels
set(gcf,'units','inches','position',[0,0,20,15])
c = colorbar;
c.Limits= [0 5];
c.TickLength = [0];
c.Label.String = '-log10(p) of mediation effect';
c.FontSize=16;
colormap(bluewhitered)


% baseline 5 mediation fl2
p_med5=p_med_FL;
p_med5(p_med5>0.05)=1;
p_med_log5=-log10(p_med5);
imagesc(p_med_log5(1:11,[7 18 21 29]))
title('FL2','FontSize',18);
set(gcf,'Position',[100 100 800 800])
set(gca, 'XTick', 1:5); % center x-axis ticks on bins
set(gca, 'YTick', 1:26); % center y-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel',{'thought','depress','stress','picture vocabulary','crystal cognition'}); % set x-axis labels
set(gca, 'YTickLabel', [strrep(network_conn_table.Properties.VariableNames(ind),'_','\_')...
    ,{'Volume - left temporal pole'},{'Volume - right ACC'},{'Volume - right temporal pole'},...
    {'Area - left temporal pole'},{'Area - right ACC'},{'Area - right temporal pole'}...
    ]); % set y-axis labels
set(gcf,'units','inches','position',[0,0,9,9])
c = colorbar;
c.Limits= [0 5];
c.TickLength = [0];
c.Label.String = '-log10(p) of mediation effect';
c.FontSize=16;
colormap(bluewhitered)



%% mediation sensitivity analysis
ind_t1=qc_table_BL{sleep_FS,'imgincl_t1w_include'}==1&qc_table_BL{sleep_NE,'imgincl_t1w_include'}==1;%matched subject that both passed qc

p_med=nan(length(ind)+6,42);
n=1;ef=[];
for i=[2 5]
    for j=[29 ]
            
     X_med=psm_table_nn_matched{sleep_psc,'TS'};
     M_med=mediators(:,i);
     Y_med=cbcl_raw_BL{sleep_psc,j};
     cov_med= psm_table_nn_matched{sleep_psc,2:13};
     
     
     b_x=regress(X_med,[ones(length(X_med),1) cov_med]);
     X_med_r=X_med-[ones(length(X_med),1) cov_med]*b_x;
     
     b_m=regress(M_med,[ones(length(M_med),1) cov_med]);
     M_med_r=M_med-[ones(length(M_med),1) cov_med]*b_m;
     b_y=regress(Y_med,[ones(length(Y_med),1) cov_med]);
     Y_med_r=Y_med-[ones(length(Y_med),1) cov_med]*b_y;
     b_m_x=regress(M_med,[ones(length(M_med),1) cov_med X_med]);
     M_med_r_x=M_med-[ones(length(M_med),1) cov_med  X_med]*b_m_x;
     b_y_x=regress(Y_med,[ones(length(Y_med),1) cov_med X_med]);
     Y_med_r_x=Y_med-[ones(length(Y_med),1) cov_med X_med]*b_y_x;
     r_my_x=corr(M_med_r_x,Y_med_r_x,'rows','pairwise');
     
     
     
     ef(n,1)=sqrt(r_my_x*nanstd(M_med_r_x)*nanstd(Y_med_r_x)/(nanstd(M_med_r)*nanstd(Y_med_r)));
     
        
    
        
        
    [paths, stats_med]=mediation(psm_table_nn_matched{sleep_psc,'TS'},cbcl_raw_BL{sleep_psc,j},...
    mediators(:,i)...
    ,'cov', psm_table_nn_matched{sleep_psc,2:13},'verbose','boot','bootsamples',10000);
%     p_med(i,j)=stats_med.p(5);

    ef(n,2)=stats_med.paths(1)*nanstd(X_med_r)/nanstd(M_med_r);
    ef(n,3)=stats_med.paths(4)*nanstd(X_med_r)/nanstd(Y_med_r);	

 n=n+1;  
    end
end





for i=[6 9]
    for j=[29 ]
            
     X_med=psm_table_nn_matched{sleep_psc([ind_t1;ind_t1]),'TS'};
     M_med=mediators([ind_t1;ind_t1],i);
     Y_med=cbcl_raw_BL{sleep_psc([ind_t1;ind_t1]),j};
     cov_med= psm_table_nn_matched{sleep_psc([ind_t1;ind_t1]),2:13};
     
     
     b_x=regress(X_med,[ones(length(X_med),1) cov_med]);
     X_med_r=X_med-[ones(length(X_med),1) cov_med]*b_x;
     
     b_m=regress(M_med,[ones(length(M_med),1) cov_med]);
     M_med_r=M_med-[ones(length(M_med),1) cov_med]*b_m;
     b_y=regress(Y_med,[ones(length(Y_med),1) cov_med]);
     Y_med_r=Y_med-[ones(length(Y_med),1) cov_med]*b_y;
     b_m_x=regress(M_med,[ones(length(M_med),1) cov_med X_med]);
     M_med_r_x=M_med-[ones(length(M_med),1) cov_med  X_med]*b_m_x;
     b_y_x=regress(Y_med,[ones(length(Y_med),1) cov_med X_med]);
     Y_med_r_x=Y_med-[ones(length(Y_med),1) cov_med X_med]*b_y_x;
     r_my_x=corr(M_med_r_x,Y_med_r_x,'rows','pairwise');
     
     
     
     ef(n,1)=sqrt(r_my_x*nanstd(M_med_r_x)*nanstd(Y_med_r_x)/(nanstd(M_med_r)*nanstd(Y_med_r)));

%         
    [paths, stats_med]=mediation(psm_table_nn_matched{sleep_psc([ind_t1;ind_t1]),'TS'},cbcl_raw_BL{sleep_psc([ind_t1;ind_t1]),j},...
    mediators([ind_t1;ind_t1],i)...
    ,'cov', psm_table_nn_matched{sleep_psc([ind_t1;ind_t1]),2:13},'verbose','boot','bootsamples',10000);
%     p_med(i,j)=stats_med.p(5);
    ef(n,2)=stats_med.paths(1)*nanstd(X_med_r)/nanstd(M_med_r);
    ef(n,3)=stats_med.paths(4)*nanstd(X_med_r)/nanstd(Y_med_r)	;
;
      n=n+1;  
    end
end














ind_fmri_fl=qc_table_FL2{sleep_FS_FL2,'imgincl_rsfmri_include'}==1&qc_table_FL2{sleep_NE_FL2,'imgincl_rsfmri_include'}==1;%matched subject that both passed qc
% n=1;
for i=[2 5]
    for j=[29 ]
     

    [paths, stats_med]=mediation(psm_table_nn_matched{sleep_psc_FL2([ind_fmri_fl;ind_fmri_fl]),'TS'},...
        cbcl_raw_FL2{sleep_psc_FL2([ind_fmri_fl;ind_fmri_fl]),j},...
    mediators(ismember(sleep_psc,sleep_psc_FL2([ind_fmri_fl;ind_fmri_fl])),i)...
    ,'cov', [psm_table_nn_matched{sleep_psc_FL2([ind_fmri_fl;ind_fmri_fl]),2:13},cbcl_raw_BL{sleep_psc_FL2([ind_fmri_fl;ind_fmri_fl]),j}]...
    ,'verbose','boot','bootsamples',10000);
    
     X_med=psm_table_nn_matched{sleep_psc_FL2([ind_fmri_fl;ind_fmri_fl]),'TS'};
     M_med=mediators(ismember(sleep_psc,sleep_psc_FL2([ind_fmri_fl;ind_fmri_fl])),i);
     Y_med=cbcl_raw_FL2{sleep_psc_FL2([ind_fmri_fl;ind_fmri_fl]),j};
     cov_med=[psm_table_nn_matched{sleep_psc_FL2([ind_fmri_fl;ind_fmri_fl]),2:13},cbcl_raw_BL{sleep_psc_FL2([ind_fmri_fl;ind_fmri_fl]),j}];
     
     b_m=regress(M_med,[ones(length(M_med),1) cov_med]);
     M_med_r=M_med-[ones(length(M_med),1) cov_med]*b_m;
     b_y=regress(Y_med,[ones(length(Y_med),1) cov_med]);
     Y_med_r=Y_med-[ones(length(Y_med),1) cov_med]*b_y;
     b_m_x=regress(M_med,[ones(length(M_med),1) cov_med X_med]);
     M_med_r_x=M_med-[ones(length(M_med),1) cov_med  X_med]*b_m_x;
     b_y_x=regress(Y_med,[ones(length(Y_med),1) cov_med X_med]);
     Y_med_r_x=Y_med-[ones(length(Y_med),1) cov_med X_med]*b_y_x;
     r_my_x=corr(M_med_r_x,Y_med_r_x,'rows','pairwise');
     
     
          
    ef(n,1)=sqrt(r_my_x*nanstd(M_med_r_x)*nanstd(Y_med_r_x)/(nanstd(M_med_r)*nanstd(Y_med_r)));
    ef(n,2)=stats_med.paths(1)*nanstd(X_med_r)/nanstd(M_med_r);
    ef(n,3)=stats_med.paths(4)*nanstd(X_med_r)/nanstd(Y_med_r)	;
      n=n+1;  
     

        
%     p_med_FL(i,j)=stats_med.p(5);
    end
end





for i=[6 9]
    for j=[29 ]
     
    [paths, stats_med]=mediation(psm_table_nn_matched{sleep_temp_FL2([ind_t1_fl;ind_t1_fl]),'TS'},...
        cbcl_raw_FL2{sleep_temp_FL2([ind_t1_fl;ind_t1_fl]),j},...
    mediators(ismember(sleep_psc,sleep_temp_FL2([ind_t1_fl;ind_t1_fl])),i)...
    ,'cov', [psm_table_nn_matched{sleep_temp_FL2([ind_t1_fl;ind_t1_fl]),2:13},cbcl_raw_BL{sleep_temp_FL2([ind_t1_fl;ind_t1_fl]),j}]...
    ,'verbose','boot','bootsamples',10000);
%     p_med_FL(i,j)=stats_med.p(5);
     X_med=psm_table_nn_matched{sleep_temp_FL2([ind_t1_fl;ind_t1_fl]),'TS'};
     M_med=mediators(ismember(sleep_psc,sleep_temp_FL2([ind_t1_fl;ind_t1_fl])),i);
     Y_med=cbcl_raw_FL2{sleep_temp_FL2([ind_t1_fl;ind_t1_fl]),j};
     cov_med=[psm_table_nn_matched{sleep_temp_FL2([ind_t1_fl;ind_t1_fl]),2:13},cbcl_raw_BL{sleep_temp_FL2([ind_t1_fl;ind_t1_fl]),j}];
     
     b_m=regress(M_med,[ones(length(M_med),1) cov_med]);
     M_med_r=M_med-[ones(length(M_med),1) cov_med]*b_m;
     b_y=regress(Y_med,[ones(length(Y_med),1) cov_med]);
     Y_med_r=Y_med-[ones(length(Y_med),1) cov_med]*b_y;
     b_m_x=regress(M_med,[ones(length(M_med),1) cov_med X_med]);
     M_med_r_x=M_med-[ones(length(M_med),1) cov_med  X_med]*b_m_x;
     b_y_x=regress(Y_med,[ones(length(Y_med),1) cov_med X_med]);
     Y_med_r_x=Y_med-[ones(length(Y_med),1) cov_med X_med]*b_y_x;
     r_my_x=corr(M_med_r_x,Y_med_r_x,'rows','pairwise');




     
    ef(n,1)=sqrt(r_my_x*nanstd(M_med_r_x)*nanstd(Y_med_r_x)/(nanstd(M_med_r)*nanstd(Y_med_r)));
    ef(n,2)=stats_med.paths(1)*nanstd(X_med_r)/nanstd(M_med_r);
    ef(n,3)=stats_med.paths(4)*nanstd(X_med_r)/nanstd(Y_med_r)	;
      n=n+1;  
    end
end








%% cohend correlation behavioral plot
cohend_combined(cohend_combined==0)=nan;
figure1=figure;
axes1 = axes('Parent',figure1);
hold(axes1,'all');
 plot1=plot(cohend(1:20)',cohend_combined(1:20)','Parent',axes1,'LineStyle','none','Marker','o','Markersize',9,'MarkerFaceColor','b','DisplayName','cbcl')
 hold on
 plot(cohend(21:30)',cohend_combined(21:30)','Parent',axes1,'LineStyle','none','Marker','o','Markersize',9,'MarkerFaceColor','g','DisplayName','nihtbx')

  hold on
 plot(cohend(31:42)',cohend_combined(31:42)','Parent',axes1,'LineStyle','none','Marker','o','Markersize',9,'MarkerFaceColor','y','DisplayName','mental health')
%   hold on
%  plot(cohend(44:55)',cohend_combined(44:55)','LineStyle','none','Marker','o','MarkerFaceColor','r','DisplayName','screen time')
%  
 set(gca,'TickLength',[0 0])
box off
xdata1 = get(plot1, 'xdata');
% 获取图中的 ydata
ydata1 = get(plot1, 'ydata');
 ylabel('Cohens'' d of behavioral assessemnts at FL2')
 xlabel('Cohens'' d of behavioral assessemnts at baseline')
 [r,p]=corr(cohend(1:42)',cohend_combined(1:42)','rows','pairwise');
txt = ['r = ',num2str(r,'%5.3f'),', p = ',num2str(p,'%5.3f')];
t=text(min(xdata1)*0.7,max(ydata1)*0.7,txt);
t.Color=[0.5 0.5 0.5];
t.FontSize=24;
set(axes1,'FontSize',24,'LineWidth',2);
set(gcf,'units','inches','position',[0,0,12,9])
% 创建 legend
legend(axes1,'show','FontSize',28,'LineWidth',2);

legend('location','best')
 [r,p,rl,ru]=corrcoef(cohend(1:42)',cohend_combined(1:42)','rows','pairwise');
%% cohend correlation behavioral plot fl1

figure1=figure;
axes1 = axes('Parent',figure1);
hold(axes1,'all');
 plot1=plot(cohend(1:20)',cohend_FL1(1:20)','Parent',axes1,'LineStyle','none','Marker','o','MarkerFaceColor','b','DisplayName','cbcl')
 hold on
%  plot(cohend(21:30)',cohend_combined(21:30)','Parent',axes1,'LineStyle','none','Marker','o','MarkerFaceColor','g','DisplayName','nihtbx')
% 
%   hold on
%  plot(cohend(31:42)',cohend_combined(31:42)','Parent',axes1,'LineStyle','none','Marker','o','MarkerFaceColor','y','DisplayName','mental health')
%   hold on
 plot(cohend(43:54)',cohend_FL1(21:32)','LineStyle','none','Marker','o','MarkerFaceColor','r','DisplayName','screen time')
 
 set(gca,'TickLength',[0 0])
box off
xdata1 = get(plot1, 'xdata');
% 获取图中的 ydata
ydata1 = get(plot1, 'ydata');
 ylabel('Cohens'' d of behavioral assessemnts at FL2')
 xlabel('Cohens'' d of behavioral assessemnts at baseline')
 [r,p]=corr(cohend([1:20 43:54])',cohend_FL1(1:32)','rows','pairwise');
txt = ['r = ',num2str(r,'%5.3f'),', p = ',num2str(p,'%5.3f')];
t=text(min(xdata1)*0.7,max(ydata1)*0.7,txt);
t.Color=[0.5 0.5 0.5];
t.FontSize=24;
set(axes1,'FontSize',24,'LineWidth',2);
set(gcf,'units','inches','position',[0,0,12,9])
% 创建 legend
legend(axes1,'show','FontSize',28,'LineWidth',2);

legend('location','best')

[r,p,rl,ru]=corrcoef(cohend([1:20 43:54])',cohend_FL1(1:32)','rows','pairwise');

%% cohend correlation network


cohen_d_shaped=[reshape(baseline_net.cohen_d_stats(1:169),13,13);reshape(baseline_net.cohen_d_stats(170:end),19,13);];
cohen_d_shaped=cohen_d_shaped([1:6 8:end],[1:6 8:end]);
cohen_d_shaped=[tril(cohen_d_shaped(1:12,1:12));cohen_d_shaped(13:end,1:12)];
cohensd_baseline=cohen_d_shaped(cohen_d_shaped(:)~=0);

cohen_d_shaped=[reshape(FL_net.cohen_d_stats(1:169),13,13);reshape(FL_net.cohen_d_stats(170:end),19,13);];
cohen_d_shaped=cohen_d_shaped([1:6 8:end],[1:6 8:end]);
cohen_d_shaped=[tril(cohen_d_shaped(1:12,1:12));cohen_d_shaped(13:end,1:12)];
cohensd_fl=cohen_d_shaped(cohen_d_shaped(:)~=0);

figure1=figure;


 
axes1 = axes('Parent',figure1);

hold(axes1,'all');
plot1=plot(cohensd_baseline,cohensd_fl,'Parent',axes1,'LineStyle','none','Marker','o','MarkerSize',10,'MarkerEdgeColor','b','DisplayName','Data');
ylabel('Cohen''s d of network connectivity at FL2')
xlabel('Cohen''s d of network connectivity at baseline')
set(gca,'TickLength',[0 0])
box off
xdata1 = get(plot1, 'xdata');
% 获取图中的 ydata
ydata1 = get(plot1, 'ydata');

xplot1 = linspace(min(xdata1)*1.05, max(xdata1))*1.05;
% ylimits=ylim(axes1);


fitResults1 = polyfit(xdata1, ydata1, 1);
% 计算多项式
yplot1 = polyval(fitResults1, xplot1);
% 绘制拟合图
fitLine1 = plot(xplot1,yplot1,'DisplayName','linear fit','Parent',axes1,...
    'Tag','linear',...
    'Color',[0.5 0.5 0.5],'LineWidth',2);

% 将新线设置到正确位置


[r,p,ru,ul]=corrcoef(cohensd_baseline,cohensd_fl);
txt = ['r = ',num2str(r,'%5.3f'),', p = ',num2str(p,'%5.3f')];
t=text(min(xdata1)*0.7,max(ydata1)*0.7,txt);
t.Color=[0.5 0.5 0.5];
t.FontSize=24;
set(axes1,'FontSize',26,'LineWidth',2);
set(gcf,'units','inches','position',[0,0,12,9])
% 创建 legend
legend(axes1,'show','FontSize',28,'LineWidth',2);

% legend('location','best')


%% cohend correlation cortical area

figure1=figure;


 
axes1 = axes('Parent',figure1);

hold(axes1,'all');
plot1=plot(table_area_TIV{2,1:148}',table_area_TIV_FL2{2,1:148}','Parent',axes1,...
    'LineStyle','none','Marker','o','MarkerSize',10,'MarkerEdgeColor','b','DisplayName','Data');
ylabel('Cohen''s d of cortical area at FL2')
xlabel('Cohen''s d of cortical area at baseline')
set(gca,'TickLength',[0 0])
box off
xdata1 = get(plot1, 'xdata');
% 获取图中的 ydata
ydata1 = get(plot1, 'ydata');

xplot1 = linspace(min(xdata1)*1.05, max(xdata1))*1.05;
% ylimits=ylim(axes1);


fitResults1 = polyfit(xdata1, ydata1, 1);
% 计算多项式
yplot1 = polyval(fitResults1, xplot1);
% 绘制拟合图
fitLine1 = plot(xplot1,yplot1,'DisplayName','linear fit','Parent',axes1,...
    'Tag','linear',...
    'Color',[0.5 0.5 0.5],'LineWidth',2);

% 将新线设置到正确位置

[r,p]=corr(table_area_TIV{2,1:148}',table_area_TIV_FL2{2,1:148}');
txt = ['r = ',num2str(r,'%5.3f'),', p = ',num2str(p,'%5.3f')];
t=text(min(xdata1)*0.7,max(ydata1)*0.7,txt);
t.Color=[0.5 0.5 0.5];
t.FontSize=24;
set(axes1,'FontSize',26,'LineWidth',2);
set(gcf,'units','inches','position',[0,0,12,9])
% 创建 legend
legend(axes1,'show','FontSize',28,'LineWidth',2);

% legend('location','best')


[r,p,ru,rl]=corrcoef(table_area_TIV{2,1:148}',table_area_TIV_FL2{2,1:148}');

%% cohend correlation gmv

figure1=figure;


 
axes1 = axes('Parent',figure1);

hold(axes1,'all');
plot1=plot([table_GMV_TIV{2,1:148}'; table_subvol_TIV{2,[1:15 17:29 31:33 36:40]}'],...
    [table_GMV_TIV_FL2{2,1:148}';table_subvol_TIV_FL2{2,[1:15 17:29 31:33 36:40]}'],'Parent',axes1,'LineStyle','none',...
    'Marker','o','MarkerSize',10,'MarkerEdgeColor','b','DisplayName','Data');
ylabel('Cohen''s d of GMV at FL2')
xlabel('Cohen''s d of GMV at baseline')
set(gca,'TickLength',[0 0])
box off
xdata1 = get(plot1, 'xdata');
% 获取图中的 ydata
ydata1 = get(plot1, 'ydata');

xplot1 = linspace(min(xdata1)*1.05, max(xdata1))*1.05;
% ylimits=ylim(axes1);


fitResults1 = polyfit(xdata1, ydata1, 1);
% 计算多项式
yplot1 = polyval(fitResults1, xplot1);
% 绘制拟合图
fitLine1 = plot(xplot1,yplot1,'DisplayName','linear fit','Parent',axes1,...
    'Tag','linear',...
    'Color',[0.5 0.5 0.5],'LineWidth',2);

% 将新线设置到正确位置

[r,p]=corr([table_GMV_TIV{2,1:148}'; table_subvol_TIV{2,1:40}'],...
    [table_GMV_TIV_FL2{2,1:148}';table_subvol_TIV_FL2{2,1:40}'],'rows','pairwise');
txt = ['r = ',num2str(r,'%5.3f'),', p = ',num2str(p,'%5.3f')];
t=text(min(xdata1)*0.7,max(ydata1)*0.7,txt);
t.Color=[0.5 0.5 0.5];
t.FontSize=24;
set(axes1,'FontSize',26,'LineWidth',2);
set(gcf,'units','inches','position',[0,0,12,9])
% 创建 legend
legend(axes1,'show','FontSize',28,'LineWidth',2);

% legend('location','best')

[r,p,rl,ru]=corrcoef([table_GMV_TIV{2,1:148}'; table_subvol_TIV{2,1:40}'],...
    [table_GMV_TIV_FL2{2,1:148}';table_subvol_TIV_FL2{2,1:40}'],'rows','pairwise');

%% cohend correlation thk

figure1=figure;


 
axes1 = axes('Parent',figure1);

hold(axes1,'all');
plot1=plot(table_thk_TIV{2,1:148}',table_thk_TIV_FL2{2,1:148}','Parent',axes1,'LineStyle','none',...
    'Marker','o','MarkerSize',10,'MarkerEdgeColor','b','DisplayName','Data');
ylabel('Cohen''s d of cortical thickness at FL2')
xlabel('Cohen''s d of cortical thickness at baseline')
set(gca,'TickLength',[0 0])
box off
xdata1 = get(plot1, 'xdata');
% 获取图中的 ydata
ydata1 = get(plot1, 'ydata');

xplot1 = linspace(min(xdata1)*1.05, max(xdata1))*1.05;
% ylimits=ylim(axes1);


fitResults1 = polyfit(xdata1, ydata1, 1);
% 计算多项式
yplot1 = polyval(fitResults1, xplot1);
% 绘制拟合图
fitLine1 = plot(xplot1,yplot1,'DisplayName','linear fit','Parent',axes1,...
    'Tag','linear',...
    'Color',[0.5 0.5 0.5],'LineWidth',2);

% 将新线设置到正确位置

[r,p]=corr(table_thk_TIV{2,1:148}',table_thk_TIV_FL2{2,1:148}');
txt = ['r = ',num2str(r,'%5.3f'),', p = ',num2str(p,'%5.5f')];
t=text(min(xdata1)*0.7,max(ydata1)*0.7,txt);
t.Color=[0.5 0.5 0.5];
t.FontSize=24;
set(axes1,'FontSize',26,'LineWidth',2);
set(gcf,'units','inches','position',[0,0,12,9])
% 创建 legend
legend(axes1,'show','FontSize',28,'LineWidth',2);

% legend('location','best')
[r,p,ru,rl]=corrcoef(table_thk_TIV{2,1:148}',table_thk_TIV_FL2{2,1:148}');

%% cohend correlation FA

figure1=figure;


 
axes1 = axes('Parent',figure1);

hold(axes1,'all');
plot1=plot(table_DTI_FA{2,1:37}',table_DTI_FA_FL2{2,1:37}','Parent',axes1,'LineStyle','none','Marker','o',...
    'MarkerSize',10,'MarkerEdgeColor','b','DisplayName','Data');
ylabel('Cohen''s d of FA at FL2')
xlabel('Cohen''s d of FA at baseline')
set(gca,'TickLength',[0 0])
box off
xdata1 = get(plot1, 'xdata');
% 获取图中的 ydata
ydata1 = get(plot1, 'ydata');

xplot1 = linspace(min(xdata1)*1.05, max(xdata1))*1.05;
% ylimits=ylim(axes1);


fitResults1 = polyfit(xdata1, ydata1, 1);
% 计算多项式
yplot1 = polyval(fitResults1, xplot1);
% 绘制拟合图
fitLine1 = plot(xplot1,yplot1,'DisplayName','linear fit','Parent',axes1,...
    'Tag','linear',...
    'Color',[0.5 0.5 0.5],'LineWidth',2);

% 将新线设置到正确位置

[r,p]=corr(table_DTI_FA{2,1:37}',table_DTI_FA_FL2{2,1:37}');
txt = ['r = ',num2str(r,'%5.2f'),', p = ',num2str(p,'%5.3f')];
t=text(min(xdata1)*0.7,max(ydata1)*0.7,txt);
t.Color=[0.5 0.5 0.5];
t.FontSize=24;
set(axes1,'FontSize',26,'LineWidth',2);
set(gcf,'units','inches','position',[0,0,12,9])
% 创建 legend
legend(axes1,'show','FontSize',28,'LineWidth',2);

% legend('location','best')
[r,p,rl,ru]=corrcoef(table_DTI_FA{2,1:37}',table_DTI_FA_FL2{2,1:37}');
%% sankey alluvialflow plot


input_data=zeros(length(sleep_psc),3);
input_data(:,1)=TS+1;


total_sleep_duration_FL1=sleep_parent_SDS_FL1{sleep_psc(ismember(sleep_psc,sleep_psc_FL1)),'sleepdisturb1_p'};%1 = 9-11 hours/ 9 a 11 horas; 2 = 8-9 hours /8 a 9 horas; 3 = 7-8 hours /7 a 8 horas; 4 = 5-7 hours /5 a 7 horas; 5 = Less than 5 hours/ Menos de 5 horas// Consider each question pertaining to the PAST 6 MONTHS of the child's life 
total_sleep_duration_FL1(total_sleep_duration_FL1==5)=4; % combine 4 & 5 as only 26 subjects have a value of 5
total_sleep_duration_FL1=5-total_sleep_duration_FL1; %reverse sleep duration to 4 are recommed sleep time.
TS_FL1=double(total_sleep_duration_FL1==4);

input_data(ismember(sleep_psc,sleep_psc_FL1),2)=TS_FL1+1;

total_sleep_duration_FL2=sleep_parent_SDS_FL2{sleep_psc(ismember(sleep_psc,sleep_psc_FL2)),'sleepdisturb1_p'};%1 = 9-11 hours/ 9 a 11 horas; 2 = 8-9 hours /8 a 9 horas; 3 = 7-8 hours /7 a 8 horas; 4 = 5-7 hours /5 a 7 horas; 5 = Less than 5 hours/ Menos de 5 horas// Consider each question pertaining to the PAST 6 MONTHS of the child's life 
total_sleep_duration_FL2(total_sleep_duration_FL2==5)=4; % combine 4 & 5 as only 26 subjects have a value of 5
total_sleep_duration_FL2=5-total_sleep_duration_FL2; %reverse sleep duration to 4 are recommed sleep time.
TS_FL2=double(total_sleep_duration_FL2==4);

input_data(ismember(sleep_psc,sleep_psc_FL2),3)=TS_FL2+1;
input_data(input_data(:,3)==2,3)=4;
% for i=1:size(input_data,1)
%    switch input_data(i,1)
%        case 1
%            if input_data(i,2)==1
%                
       
    test_str=input_data(:,1)+input_data(:,2)*2;
    
    
    
%     end

test_str=test_str-2;
test_str(test_str<=0)=0;
    input_data(input_data(:,1)==2,1)=4;
CreateSankeyPlot([input_data(:,1) test_str input_data(:,3)])
set(gcf,'units','inches','position',[0,0,12,9])

% combine groups
TS_1_always=sleep_psc(input_data(:,1)==2&input_data(:,3)==2);
TS_1_changed=sleep_psc(input_data(:,1)==2&input_data(:,3)==1);
TS_2_always=sleep_psc(input_data(:,1)==1&input_data(:,3)==1);
TS_2_changed=sleep_psc(input_data(:,1)==1&input_data(:,3)==2);

boxplot([cognition_row_FL2{TS_1_always,9};cognition_row_FL2{TS_1_changed,9};...
    cognition_row_FL2{TS_2_always,9};cognition_row_FL2{TS_2_changed,9}]...
    ,[ones(length(TS_1_always),1);ones(length(TS_1_changed),1)+1;...
    ones(length(TS_2_always),1)+2;ones(length(TS_2_changed),1)+3;])

data_ts=[nanmedian(cognition_row_FL2{TS_1_always,9}), nanmedian(cognition_row_FL2{TS_1_changed,9}),...
    nanmedian(cognition_row_FL2{TS_2_always,9}),nanmedian(cognition_row_FL2{TS_2_changed,9})];
err=[iqr(cognition_row_FL2{TS_1_always,9}), iqr(cognition_row_FL2{TS_1_changed,9}),...
    iqr(cognition_row_FL2{TS_2_always,9}),iqr(cognition_row_FL2{TS_2_changed,9})];


figure;
bar(1:4,data_ts)
% hold on
% er = errorbar(1:4,data_ts,err,err);
% er.Color=[0 0 0];
% er.LineStyle='None';
% hold off
ylim([85 95])
set(gca, 'FontName', 'Helvetica')
xticks(1:4);
xticklabels({'SS-SS','SS-IS','IS-IS','IS-SS'})
ylabel({'Median crystalized intelligence'},'FontName', 'Helvetica')
set(gca,'TickLength',[0 0])
set(gcf,'units','inches','position',[0, 0, 10, 6])
box off


%% cohen d table  supplemental
cohend_brain=table;

cohend_baseline_brain=[baseline_net.cohen_d_stats(ind)';table_GMV_TIV{2,[43 80 117]}'...
    ;table_area_TIV{2,[43 80 117]}'];
cohend_FL2_brain=[FL_net.cohen_d_stats(ind)';table_GMV_TIV_FL2{2,[43 80 117]}'...
    ;table_area_TIV_FL2{2,[43 80 117]}'];
cohend_brain{:,1:2}=[cohend_baseline_brain,cohend_FL2_brain];
cohend_brain.Properties.RowNames=[network_conn_table.Properties.VariableNames(ind)...
    ,{'Volume - left temporal pole'},{'Volume - right ACC'},{'Volume - right temporal pole'},...
    {'Area - left temporal pole'},{'Area - right ACC'},{'Area - right temporal pole'}...
    ];
cohend_brain.Properties.VariableNames=[{'cohens''d BL'},{'cohens''d FL2'}];
writetable(cohend_brain, 'cohend_brain.xlsx')


comp_table=table;
comp_table{:,1:2}=[mean(psm_table_nonan{sleep_psc,[1:5 7:9 11:12]})',...
    nanmean(psm_table_nonan{psm_table_nonan.Properties.RowNames(~ismember(psm_table_nonan.Properties.RowNames,sleep_psc)),[1:5 7:9 11:12]})'];
comp_table.Properties.RowNames=psm_table_nonan.Properties.VariableNames([1:5 7:9 11:12]);
n=1;
p=[];
for i=[1:5 7:9 11:12]
    temp=psm_table_nonan{psm_table_nonan.Properties.RowNames(~ismember(psm_table_nonan.Properties.RowNames,sleep_psc)),i};
    temp=temp(~isnan(temp));
    [t,p(n)]=ttest2(psm_table_nonan{sleep_psc,i},...
       temp);
   n=n+1;
end
comp_table{:,3}=p';
comp_table.Properties.VariableNames={'matched sample','unmatched sample','p_value'};
writetable(comp_table,'dropped_test.xlsx')



%% pvalue table supplemental material

data_table=table('Size',[42 3],'VariableTypes',["string" "string"  "string"],...
    'VariableNames',{'Behaviors','Baseline','2 years follow-up'});

% summary_table_SDS{summary_table_SDS{sub_id_BL,6}==0,6}=nan;

beh_fdr_p_BL=round(beh_fdr_p_BL*10000)/10000;

for i = 1:42
    
    
    
    data_table{i,2}={sprintf('%.2g',beh_fdr_p_BL(i))};
    data_table{i,3}={sprintf('%.2g',beh_fdr_p_FL2(i))};
    
   
end

 
data_table{:,1}=cbcl_raw_FL2.Properties.VariableNames';
writetable(data_table,'beh_p.xlsx')



%%[~,~,~,baseline_net.p_values]=fdr_bh(p_value([1:78 92:283 303:end]));

%[~,~,~,FL_net.p_values]=fdr_bh(p_value([1:78 92:283 303:end]));

data_table=table('Size',[384 3],'VariableTypes',["string" "string"  "string"],...
    'VariableNames',{'Resting-state functional connectivity','Baseline','2 years follow-up'});

% summary_table_SDS{summary_table_SDS{sub_id_BL,6}==0,6}=nan;

baseline_net.p_values=round(baseline_net.p_values*10000)/10000;
baseline_net.p_values(baseline_net.p_values<0.0001)=0.0001;
FL_net.p_values=round(FL_net.p_values*10000)/10000;
for i = 1:384
    
    
    
    data_table{i,2}={sprintf('%.2g',baseline_net.p_values(i))};
    data_table{i,3}={sprintf('%.2g',FL_net.p_values(i))};
    
   
end

 
data_table{:,1}=network_conn_table.Properties.VariableNames([1:78 92:283 303:end])';
writetable(data_table,'net_p.xlsx')



[~,~,~,CA_BL_p]=fdr_bh(table_area_TIV{3,1:148});

[~,~,~,CA_FL2_p]=fdr_bh(table_area_TIV_FL2{3,1:148});

data_table=table('Size',[148 3],'VariableTypes',["string" "string"  "string"],...
    'VariableNames',{'Cortical Area','Baseline','2 years follow-up'});

% summary_table_SDS{summary_table_SDS{sub_id_BL,6}==0,6}=nan;

CA_BL_p=round(CA_BL_p*10000)/10000;
CA_BL_p(CA_BL_p<0.0001)=0.0001;
CA_FL2_p=round(CA_FL2_p*10000)/10000;
CA_FL2_p(CA_FL2_p<0.0001)=0.0001;
for i = 1:148
    
    
    
    data_table{i,2}={sprintf('%.2g',CA_BL_p(i))};
    data_table{i,3}={sprintf('%.2g',CA_FL2_p(i))};
    
   
end

 
data_table{:,1}=table_area_TIV_FL2.Properties.VariableNames(1:148)';
writetable(data_table,'CA_p.xlsx')



[~,~,~,GMV_BL_p]=fdr_bh([table_GMV_TIV{3,1:148} table_subvol_TIV{3,[1:15 17:29 31:33 36:40]}]);

[~,~,~,GMV_FL2_p]=fdr_bh([table_GMV_TIV_FL2{3,1:148} table_subvol_TIV_FL2{2,[1:15 17:29 31:33 36:40]}]);

data_table=table('Size',[184 3],'VariableTypes',["string" "string"  "string"],...
    'VariableNames',{'GMV','Baseline','2 years follow-up'});

% summary_table_SDS{summary_table_SDS{sub_id_BL,6}==0,6}=nan;

GMV_BL_p=round(GMV_BL_p*10000)/10000;
GMV_BL_p(GMV_BL_p<0.0001)=0.0001;
GMV_FL2_p=round(GMV_FL2_p*10000)/10000;
GMV_FL2_p(GMV_FL2_p<0.0001)=0.0001;
for i = 1:184
    
    
    
    data_table{i,2}={sprintf('%.2g',GMV_BL_p(i))};
    data_table{i,3}={sprintf('%.2g',GMV_FL2_p(i))};
    
   
end

 
data_table{:,1}=[table_GMV_TIV_FL2.Properties.VariableNames(1:148),table_subvol_TIV_FL2.Properties.VariableNames([1:15 17:29 31:33 36:40])]';
writetable(data_table,'GMV_p.xlsx')



[~,~,~,thk_BL_p]=fdr_bh(table_thk_TIV{3,1:148});

[~,~,~,thk_FL2_p]=fdr_bh(table_thk_TIV_FL2{3,1:148});

data_table=table('Size',[148 3],'VariableTypes',["string" "string"  "string"],...
    'VariableNames',{'Cortical Thickness','Baseline','2 years follow-up'});

% summary_table_SDS{summary_table_SDS{sub_id_BL,6}==0,6}=nan;

thk_BL_p=round(thk_BL_p*10000)/10000;
thk_BL_p(thk_BL_p<0.0001)=0.0001;
thk_FL2_p=round(thk_FL2_p*10000)/10000;
thk_FL2_p(thk_FL2_p<0.0001)=0.0001;
for i = 1:148
    
    
    
    data_table{i,2}={sprintf('%.2g',thk_BL_p(i))};
    data_table{i,3}={sprintf('%.2g',thk_FL2_p(i))};
    
   
end

 
data_table{:,1}=table_thk_TIV_FL2.Properties.VariableNames(1:148)';
writetable(data_table,'thk_p_test.xlsx')





[~,~,~,FA_BL_p]=fdr_bh(table_DTI_FA{3,1:37});

[~,~,~,FA_FL2_p]=fdr_bh(table_DTI_FA_FL2{3,1:37});

data_table=table('Size',[37 3],'VariableTypes',["string" "string"  "string"],...
    'VariableNames',{'FA','Baseline','2 years follow-up'});

% summary_table_SDS{summary_table_SDS{sub_id_BL,6}==0,6}=nan;

FA_BL_p=round(FA_BL_p*10000)/10000;
FA_BL_p(FA_BL_p<0.0001)=0.0001;
FA_FL2_p=round(FA_FL2_p*10000)/10000;
FA_FL2_p(FA_FL2_p<0.0001)=0.0001;
for i = 1:37
    
    
    
    data_table{i,2}={sprintf('%.2g',FA_BL_p(i))};
    data_table{i,3}={sprintf('%.2g',FA_FL2_p(i))};
    
   
end

 
data_table{:,1}=table_DTI_FA_FL2.Properties.VariableNames(1:37)';
writetable(data_table,'FA_p.xlsx')



p_values_med=p_med(1:11,[ 7 18 21 29]);
p_values_med_FL=p_med_FL(1:11,[ 7 18 21 29]);
p_values_med=round(p_values_med*10000)/10000;
p_values_med(p_values_med<0.0001)=0.0001;
p_values_med_FL=round(p_values_med_FL*10000)/10000;
p_values_med_FL(p_values_med_FL<0.0001)=0.0001;



