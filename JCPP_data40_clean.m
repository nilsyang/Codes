%% initial parameters

father_dir='/home/yangf7/Documents/Nils/ABCD';
cd(father_dir)
beh_dir='/home/yangf7/Documents/Nils/ABCD/Beh_tabulated_data_40/';

network_table=readtable(fullfile(beh_dir,'abcd_betnet02.txt'));
network_table_FL2=network_table(strcmp(network_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
network_table_BL=network_table(strcmp(network_table{:,'eventname'},'baseline_year_1_arm_1'),:);
network_table_BL.Properties.RowNames=network_table_BL.src_subject_id;
network_table_FL2.Properties.RowNames=network_table_FL2.src_subject_id;
qc_table=readtable(fullfile(beh_dir, 'abcd_fmriqc01.txt'));%new qc table for ABCD40
qc_table_BL=qc_table(strcmp(qc_table{:,'eventname'},'baseline_year_1_arm_1'),:);
qc_table_FL2=qc_table(strcmp(qc_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
qc_table_BL.Properties.RowNames=qc_table_BL.src_subject_id;
qc_table_FL2.Properties.RowNames=qc_table_FL2.src_subject_id;
DTI_qc_table=readtable(fullfile(beh_dir, 'dmriqc01.txt'));%new qc table for ABCD40
DTI_qc_table_BL=DTI_qc_table(strcmp(DTI_qc_table{:,'eventname'},'baseline_year_1_arm_1'),:);
DTI_qc_table_FL2=DTI_qc_table(strcmp(DTI_qc_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
DTI_qc_table_BL.Properties.RowNames=DTI_qc_table_BL.src_subject_id;
DTI_qc_table_FL2.Properties.RowNames=DTI_qc_table_FL2.src_subject_id;
DTI_table=readtable(fullfile(beh_dir,'abcd_dmdtifp101.txt'));%10:51
DTI_table_BL=DTI_table(strcmp(DTI_table{:,'eventname'},'baseline_year_1_arm_1'),:);
DTI_table_BL.Properties.RowNames=DTI_table_BL.src_subject_id;
DTI_table_FL2=DTI_table(strcmp(DTI_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
DTI_table_FL2.Properties.RowNames=DTI_table_FL2.src_subject_id;
qc_inclusion=readtable(fullfile(beh_dir, 'abcd_imgincl01.txt'));
qc_inclusion_BL=qc_inclusion(strcmp(qc_inclusion{:,'eventname'},'baseline_year_1_arm_1'),:);
qc_inclusion_FL2=qc_inclusion(strcmp(qc_inclusion{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
qc_inclusion_BL.Properties.RowNames=qc_inclusion_BL.src_subject_id;
qc_inclusion_FL2.Properties.RowNames=qc_inclusion_FL2.src_subject_id;



sub_id_BL=DTI_table_BL.Properties.RowNames;

sub_id_BL=sub_id_BL(ismember(sub_id_BL,network_table_BL.Properties.RowNames));% all subs that have network value
sub_id_BL=sub_id_BL(sum(isnan(network_table_BL{sub_id_BL,23:191}'))==0); %excluding ppl have missing values
sub_id_BL=sub_id_BL(sum(isnan(DTI_table_BL{sub_id_BL,10:93}'))==0); %excluding ppl have missing values

sub_id_BL=sub_id_BL(ismember(sub_id_BL,qc_inclusion_BL.Properties.RowNames));% all subs that have qc values


sub_id_BL=sub_id_BL(qc_inclusion_BL{sub_id_BL,'imgincl_rsfmri_include'}==1&qc_inclusion_BL{sub_id_BL,'imgincl_dmri_include'}==1);% all subs that passed dti & fmri mri qc



sub_network_table=readtable(fullfile(beh_dir,'mrirscor02.txt'));
sub_network_table_FL2=sub_network_table(strcmp(sub_network_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
sub_network_table_BL=sub_network_table(strcmp(sub_network_table{:,'eventname'},'baseline_year_1_arm_1'),:);
sub_network_table_BL.Properties.RowNames=sub_network_table_BL.src_subject_id;
sub_network_table_FL2.Properties.RowNames=sub_network_table_FL2.src_subject_id;

cbcl_sum=readtable(fullfile(beh_dir,'abcd_cbcls01.txt'));
cbcl_sum_FL2=cbcl_sum(strcmp(cbcl_sum{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
cbcl_sum_BL=cbcl_sum(strcmp(cbcl_sum{:,'eventname'},'baseline_year_1_arm_1'),:);
cbcl_sum_BL.Properties.RowNames=cbcl_sum_BL.src_subject_id;
cbcl_sum_FL2.Properties.RowNames=cbcl_sum_FL2.src_subject_id;
cbcl_sum_FL1=cbcl_sum(strcmp(cbcl_sum{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
cbcl_sum_FL1.Properties.RowNames=cbcl_sum_FL1.src_subject_id;

sleep_parent_SDSC=readtable(fullfile(beh_dir,'abcd_ssphp01.txt'));% total score (sds_p_ss_total) for SDSC 
sleep_parent_SDSC_BL=sleep_parent_SDSC(strcmp(sleep_parent_SDSC{:,'eventname'},'baseline_year_1_arm_1'),:);
sleep_parent_SDSC_BL.Properties.RowNames=sleep_parent_SDSC_BL.src_subject_id;
sleep_parent_SDSC_FL2=sleep_parent_SDSC(strcmp(sleep_parent_SDSC{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
sleep_parent_SDSC_FL2.Properties.RowNames=sleep_parent_SDSC_FL2.src_subject_id;
sleep_parent_SDSC_FL1=sleep_parent_SDSC(strcmp(sleep_parent_SDSC{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
sleep_parent_SDSC_FL1.Properties.RowNames=sleep_parent_SDSC_FL1.src_subject_id;




NIH_TB=readtable(fullfile(beh_dir,'abcd_tbss01.txt'));% total score (sds_p_ss_total) for NIH toolbox 
NIH_TB_BL=NIH_TB(strcmp(NIH_TB{:,'eventname'},'baseline_year_1_arm_1'),:);
NIH_TB_BL.Properties.RowNames=NIH_TB_BL.src_subject_id;
NIH_TB_FL2=NIH_TB(strcmp(NIH_TB{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
NIH_TB_FL2.Properties.RowNames=NIH_TB_FL2.src_subject_id;
NIH_TB_FL1=NIH_TB(strcmp(NIH_TB{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
NIH_TB_FL1.Properties.RowNames=NIH_TB_FL1.src_subject_id;



sleep_parent_SDS=readtable(fullfile(beh_dir,'abcd_sds01.txt'));% total hours of sleep (sleep_1_p) for SDS 
%(Parent Sleep Disturbance Scale for Children)
sleep_parent_SDS_BL=sleep_parent_SDS(strcmp(sleep_parent_SDS{:,'eventname'},'baseline_year_1_arm_1'),:);
sleep_parent_SDS_BL.Properties.RowNames=sleep_parent_SDS_BL.src_subject_id;
sleep_parent_SDS_FL2=sleep_parent_SDS(strcmp(sleep_parent_SDS{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
sleep_parent_SDS_FL2.Properties.RowNames=sleep_parent_SDS_FL2.src_subject_id;
sleep_parent_SDS_FL1=sleep_parent_SDS(strcmp(sleep_parent_SDS{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
sleep_parent_SDS_FL1.Properties.RowNames=sleep_parent_SDS_FL1.src_subject_id;





dem_table=readtable(fullfile(beh_dir,'ABCD_sites_effect.xlsx'));
for i=1:size(dem_table,1)
    dem_table{i,'id'}{1}=[dem_table{i,'id'}{1}(1:4),'_',dem_table{i,'id'}{1}(5:end)];

end
for i=1:size(dem_table,1)
    dem_table{i,'abcd_site'}{1}=[dem_table{i,'abcd_site'}{1}(5:end)];

end
dem_table.abcd_site=str2num(cell2mat(dem_table.abcd_site));

dem_table.Properties.RowNames=dem_table.id;
dem_table.sex=double(strcmp(dem_table.sex_at_birth,'M'));
dem_table.white=double(strcmp(dem_table.race_4level,'White'));
dem_table.black=double(strcmp(dem_table.race_4level,'Black'));
dem_table.mixed=double(strcmp(dem_table.race_4level,'Other/Mixed'));
dem_table.asian=double(strcmp(dem_table.race_4level,'Asian'));
dem_table.other=double(~any([dem_table.white,dem_table.black]'))';
dem_table.race=dem_table.white+2*dem_table.black+3*dem_table.other;
dem_table.race4=dem_table.white+2*dem_table.black+3*dem_table.asian+4*dem_table.mixed;
dem_table.race4(dem_table.race4==0)=4;
pdem_table=readtable(fullfile(beh_dir,'pdem02.txt'));
pdem_table.Properties.RowNames=pdem_table.src_subject_id;


dem_table.edu1=double(strcmp(dem_table.high_educ,'< HS Diploma'));
dem_table.edu2=double(strcmp(dem_table.high_educ,'HS Diploma/GED'));
dem_table.edu3=double(strcmp(dem_table.high_educ,'Some College'));
dem_table.edu4=double(strcmp(dem_table.high_educ,'Bachelor'));
dem_table.edu5=double(strcmp(dem_table.high_educ,'Post Graduate Degree'));
dem_table.edu=dem_table.edu1+2*dem_table.edu2+3*dem_table.edu3+4*dem_table.edu4+5*dem_table.edu5;



% 
BMI=readtable(fullfile(beh_dir,'BMI.csv'));
BMI_BL=BMI(strcmp(BMI{:,'event_name'},'baseline_year_1_arm_1'),:);
BMI_BL.Properties.RowNames=BMI_BL.src_subject_id;
BMI_FL2=BMI(strcmp(BMI{:,'event_name'},'2_year_follow_up_y_arm_1'),:);
BMI_FL2.Properties.RowNames=BMI_FL2.src_subject_id;
% 



mental_health=readtable(fullfile(beh_dir,'abcd_mhy02.txt'));
mental_health_BL=mental_health(strcmp(mental_health{:,'eventname'},'baseline_year_1_arm_1'),:);
mental_health_BL.Properties.RowNames=mental_health_BL.src_subject_id;
mental_health_FL1=mental_health(strcmp(mental_health{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
mental_health_FL1.Properties.RowNames=mental_health_FL1.src_subject_id;
mental_health_FL2=mental_health(strcmp(mental_health{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
mental_health_FL2.Properties.RowNames=mental_health_FL2.src_subject_id;




parent_dem_table=readtable(fullfile(beh_dir,'abcd_lpds01.txt'));
parent_dem_table_BL=parent_dem_table(strcmp(parent_dem_table{:,'eventname'},'baseline_year_1_arm_1'),:);
parent_dem_table_BL.Properties.RowNames=parent_dem_table_BL.src_subject_id;
parent_dem_table_FL2=parent_dem_table(strcmp(parent_dem_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
parent_dem_table_FL2.Properties.RowNames=parent_dem_table_FL2.src_subject_id;

% network_conn_table=([network_table_BL(sub_id_BL,23:191) sub_network_table_BL(sub_id_BL,23:269)]);


T1_table=readtable(fullfile(beh_dir,'abcd_mrisdp10201.txt'));
T1_table_BL=T1_table(strcmp(T1_table{:,'eventname'},'baseline_year_1_arm_1'),:);
T1_table_BL.Properties.RowNames=T1_table_BL.src_subject_id;
T1_table_FL2=T1_table(strcmp(T1_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
T1_table_FL2.Properties.RowNames=T1_table_FL2.src_subject_id;


T1_sub=readtable(fullfile(beh_dir,'abcd_smrip201.txt'));%330:375
T1_sub_BL=T1_sub(strcmp(T1_sub{:,'eventname'},'baseline_year_1_arm_1'),:);
T1_sub_BL.Properties.RowNames=T1_sub_BL.src_subject_id;
T1_sub_FL2=T1_sub(strcmp(T1_sub{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
T1_sub_FL2.Properties.RowNames=T1_sub_FL2.src_subject_id;





% 
DTI_roi_table=readtable(fullfile(beh_dir,'abcd_dmdtifp202.txt'));%10:51
DTI_roi_table_BL=DTI_roi_table(strcmp(DTI_roi_table{:,'eventname'},'baseline_year_1_arm_1'),:);
DTI_roi_table_BL.Properties.RowNames=DTI_roi_table_BL.src_subject_id;
DTI_roi_table_FL2=DTI_roi_table(strcmp(DTI_roi_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
DTI_roi_table_FL2.Properties.RowNames=DTI_roi_table_FL2.src_subject_id;



alcohol_table=readtable(fullfile(beh_dir,'abcd_ysu02.txt'));
alcohol_table_BL=alcohol_table(strcmp(alcohol_table{:,'eventname'},'baseline_year_1_arm_1'),:);
alcohol_table_BL.Properties.RowNames=alcohol_table_BL.src_subject_id;
alcohol_table_FL2=alcohol_table(strcmp(alcohol_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
alcohol_table_FL2.Properties.RowNames=alcohol_table_FL2.src_subject_id;



substance_table=readtable(fullfile(beh_dir,'abcd_tlfb01.txt'));
substance_table_BL=substance_table(strcmp(substance_table{:,'eventname'},'baseline_year_1_arm_1'),:);
substance_table_BL.Properties.RowNames=substance_table_BL.src_subject_id;
substance_table_FL1=substance_table(strcmp(substance_table{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
substance_table_FL1.Properties.RowNames=substance_table_FL1.src_subject_id;
substance_table_FL2=substance_table(strcmp(substance_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
substance_table_FL2.Properties.RowNames=substance_table_FL2.src_subject_id;
substance_table_FL3=substance_table(strcmp(substance_table{:,'eventname'},'3_year_follow_up_y_arm_1'),:);
substance_table_FL3.Properties.RowNames=substance_table_FL3.src_subject_id;





% abcd_ysu02 'tlfb_alc' ever used alcohol
% tlfb_tob
alcohol_curiosity_table=readtable(fullfile(beh_dir,'abcd_ysua01.txt'));
alcohol_curiosity_table_FL1=alcohol_curiosity_table(strcmp(alcohol_curiosity_table{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
alcohol_curiosity_table_FL1.Properties.RowNames=alcohol_curiosity_table_FL1.src_subject_id;
alcohol_curiosity_table_FL2=alcohol_curiosity_table(strcmp(alcohol_curiosity_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
alcohol_curiosity_table_FL2.Properties.RowNames=alcohol_curiosity_table_FL2.src_subject_id;
% abcd_ysua01 'path_alc_youth2_l' curiosity of drink alcohol
%path_alc_youth1 tobacco

nicotine_use=nansum([alcohol_table_BL{sub_id_BL,'tlfb_tob_puff'},alcohol_table_BL{sub_id_BL,'tlfb_chew_use'}...
    ,alcohol_table_BL{sub_id_BL,'tlfb_cigar_use'},alcohol_table_BL{sub_id_BL,'tlfb_hookah_use'}...
    ,alcohol_table_BL{sub_id_BL,'tlfb_pipes_use'},alcohol_table_BL{sub_id_BL,'tlfb_nicotine_use'}],2);
nicotine_use=nicotine_use>0;


drug_use=any(alcohol_table_BL{sub_id_BL,36:71},2);%105 total



sub_id_BL=sub_id_BL(drug_use==0);%drug naive subject
sub_id_FL1=sub_id_BL(ismember(sub_id_BL,sleep_parent_SDS_FL1.Properties.RowNames));
sub_id_FL2=sub_id_BL(ismember(sub_id_BL,sleep_parent_SDS_FL2.Properties.RowNames));
sub_id_FL2img=sub_id_BL(ismember(sub_id_BL,network_table_FL2.Properties.RowNames));
sub_id_FL3=sub_id_BL(ismember(sub_id_BL,substance_table_FL3.Properties.RowNames));
    




% p_edu=parent_dem_table_FL2{sub_id_BL,'demo_prnt_ed_v2_2yr_l'};
% p2_edu=parent_dem_table_FL2{sub_id_BL,'demo_prtnr_ed_v2_2yr_l'};
% p1_edu=[p_edu,p2_edu];
pubertal_f=sleep_parent_SDSC_BL{sub_id_BL,'pds_p_ss_female_category'};
pubertal_m=sleep_parent_SDSC_BL{sub_id_BL,'pds_p_ss_male_category'};
pubertal=sum([pubertal_f,pubertal_m],2,'omitnan');
pubertal(pubertal==0)=nan;






House_in=pdem_table{sub_id_BL,'demo_comb_income_v2'};
House_in(House_in>100)=nan;

total_sleep_duration=sleep_parent_SDS_BL{sub_id_BL,'sleepdisturb1_p'};%
%1 = 9-11 hours/ 9 a 11 horas; 2 = 8-9 hours /8 a 9 horas; 3 = 7-8 hours /7 a 8 horas; 4 = 5-7 hours /5 a 7 horas; 5 = Less than 5 hours/ Menos de 5 horas// Consider each question pertaining to the PAST 6 MONTHS of the child's life 
total_sleep_duration(total_sleep_duration==5)=4; % combine 4 & 5 as only 26 subjects have a value of 5
total_sleep_duration=5-total_sleep_duration; %reverse sleep duration to 4 are recommed sleep time.


total_sleep_duration_FL2=sleep_parent_SDS_FL2{sub_id_FL2img,'sleepdisturb1_p'};%
%1 = 9-11 hours/ 9 a 11 horas; 2 = 8-9 hours /8 a 9 horas; 3 = 7-8 hours /7 a 8 horas; 4 = 5-7 hours /5 a 7 horas; 5 = Less than 5 hours/ Menos de 5 horas// Consider each question pertaining to the PAST 6 MONTHS of the child's life 
total_sleep_duration_FL2(total_sleep_duration_FL2==5)=4; % combine 4 & 5 as only 26 subjects have a value of 5
total_sleep_duration_FL2=5-total_sleep_duration_FL2; %reverse sleep duration to 4 are recommed sleep time.




psm_table=[network_table_BL(sub_id_BL,'interview_age'),...
        dem_table(sub_id_BL,'race4'),dem_table(sub_id_BL,'sex'),dem_table(sub_id_BL,'edu')...
        ,dem_table(sub_id_BL,'abcd_site'),dem_table(sub_id_BL,'rel_family_id'),...
        network_table_BL(sub_id_BL,'rsfmri_c_ngd_ntpoints'),network_table_BL(sub_id_BL,'rsfmri_c_ngd_meanmotion')...
        BMI_BL(sub_id_BL,'anthro_bmi_calc')];

psm_table.Properties.VariableNames={'age','race','sex','prt_edu','sites','family','nframe','FD','BMI'};
psm_table.age_sex=psm_table.age.*psm_table.sex;
psm_table.puberty=pubertal;
psm_table.HI=House_in;
TS=double(total_sleep_duration);
psm_table.TS=TS;




dh_table=readtable(fullfile(beh_dir,'dhx01'));% developmental history aka prenatal substance use exposure
dh_table.Properties.RowNames=dh_table.src_subject_id;


bef_drugs_pre=dh_table{sub_id_BL,[49 52 59 62 65 68 71]};%used any drug during pregenancy
aft_drugs_pre=dh_table{sub_id_BL,[147 150 157 160 163 166 169]};
bef_drugs_pre(bef_drugs_pre==999)=nan;
aft_drugs_pre(aft_drugs_pre==999)=nan;
all_drug_pre=any([bef_drugs_pre,aft_drugs_pre],2);


SU_env_table=readtable(fullfile(beh_dir,'abcd_pssudse01'));% susbtance use environment 
SU_env_table_BL=SU_env_table(strcmp(SU_env_table{:,'eventname'},'baseline_year_1_arm_1'),:);% no baseline data

FM_his_table=readtable(fullfile(beh_dir,'abcd_fhxssp01'));% family history
FM_his_table.Properties.RowNames=FM_his_table.src_subject_id;

overall_alc=FM_his_table{sub_id_BL,'famhx_ss_parent_alc_p'};
overall_drug=FM_his_table{sub_id_BL,'famhx_ss_parent_dg_p'};

overall_alc(overall_alc==-1)=2;
overall_alc(overall_alc==-2)=1;
overall_alc(overall_alc==1|overall_alc==2)=0.5;
overall_alc(overall_alc==3)=1;% recoding overal_alc;

overall_drug(overall_drug==-1)=2;
overall_drug(overall_drug==-2)=1;
overall_drug(overall_drug==1|overall_drug==2)=0.5;
overall_drug(overall_drug==3)=1;% recoding overal_drug;

overall_SU=sum([overall_drug,overall_alc],2);


mid_beh_table=readtable(fullfile(beh_dir,'abcd_mid02.txt'));
mid_beh_table_BL=mid_beh_table(strcmp(mid_beh_table{:,'eventname'},'baseline_year_1_arm_1'),:);
mid_beh_table_BL.Properties.RowNames=mid_beh_table_BL.src_subject_id;
mid_beh_table_FL2=mid_beh_table(strcmp(mid_beh_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
mid_beh_table_FL2.Properties.RowNames=mid_beh_table_FL2.src_subject_id;

% su_tlfb_cal_scr_nic_days_yr  total days of any type of nicotine use

nicotine_use_FL1=substance_table_FL1{sub_id_FL1, 'su_tlfb_cal_scr_nic_days_yr'};
nicotine_use_FL2=substance_table_FL2{sub_id_FL2, 'su_tlfb_cal_scr_nic_days_yr'}; %31 tried
nicotine_use_FL3=substance_table_FL3{sub_id_FL3, 'su_tlfb_cal_scr_nic_days_yr'}; %55 tried
nicotine_use_FL23=substance_table_FL3{sub_id_FL3(ismember(sub_id_FL3,sub_id_FL2)), 'su_tlfb_cal_scr_nic_days_yr'}>0 ...
|substance_table_FL2{sub_id_FL3(ismember(sub_id_FL3,sub_id_FL2)), 'su_tlfb_cal_scr_nic_days_yr'}>0;%64 tried

drugs_use_FL1=any([substance_table_FL1{sub_id_FL1, 'tlfb_cal_scr_alc_td_cum'},substance_table_FL1{sub_id_FL1, 'su_tlfb_cal_scr_nic_days_cum'},...
    substance_table_FL1{sub_id_FL1, 'su_tlfb_cal_scr_mj_days_cum'},substance_table_FL1{sub_id_FL1, 'tlfb_cal_scr_tot_su_ud'}],2);
%0 in FL1
drugs_use_FL2=any([substance_table_FL2{sub_id_FL2, 'tlfb_cal_scr_alc_td_cum'},substance_table_FL2{sub_id_FL2, 'su_tlfb_cal_scr_nic_days_cum'},...
    substance_table_FL2{sub_id_FL2, 'su_tlfb_cal_scr_mj_days_cum'},substance_table_FL2{sub_id_FL2, 'tlfb_cal_scr_tot_su_ud'}],2);
%49 in FL2
drugs_use_FL3=any([substance_table_FL3{sub_id_FL3, 'tlfb_cal_scr_alc_td_cum'},substance_table_FL3{sub_id_FL3, 'su_tlfb_cal_scr_nic_days_cum'},...
    substance_table_FL3{sub_id_FL3, 'su_tlfb_cal_scr_mj_days_cum'},substance_table_FL3{sub_id_FL3, 'tlfb_cal_scr_tot_su_ud'}],2);
%74 in FL3 



sub_id_FL23=sub_id_FL3(ismember(sub_id_FL3,sub_id_FL2)); % 3430 children that have both FL2 and FL3 substance use data


drugs_use_FL23=any([drugs_use_FL2(ismember(sub_id_FL2,sub_id_FL3)),drugs_use_FL3(ismember(sub_id_FL3,sub_id_FL2))],2);
%92

    [~,chi2,temp] = crosstab(drugs_use_FL23, psm_table{sub_id_FL23,'TS'});% chi2=15.9 p=0.0012

    



%% combat 

sorted_sites=dem_table{sub_id_BL,'abcd_site'};
sorted_sites(sorted_sites==20)=1;
sorted_sites(sorted_sites==21)=17;
sorted_sites(sorted_sites==22)=19;

net_BL_combat=combat([network_table_BL{sub_id_BL,23:191} sub_network_table_BL{sub_id_BL,23:269}]',sorted_sites',...
[dem_table{sub_id_BL,'sex'},dem_table{sub_id_BL,'race4'},network_table_BL{sub_id_BL,'interview_age'} ],1);% running combat to control the effects of those variables 

sub_network_table_BL_com=sub_network_table_BL;
sub_network_table_BL_com{sub_id_BL,23:269}=net_BL_combat(170:end,:)';

DTI_BL_combat=combat([DTI_table_BL{sub_id_BL,10:93}]',sorted_sites',...
[dem_table{sub_id_BL,'sex'},dem_table{sub_id_BL,'race4'},network_table_BL{sub_id_BL,'interview_age'} ],1);% running combat to control the effects of those variables 

DTI_table_BL_com=DTI_table_BL;
DTI_table_BL_com{sub_id_BL,10:93}=DTI_BL_combat';


temp_FL2=sum(isnan([[network_table_FL2{sub_id_FL2img,23:191} sub_network_table_FL2{sub_id_FL2img,23:269}],dem_table{sub_id_FL2img,'abcd_site'},...
[dem_table{sub_id_FL2img,'sex'},dem_table{sub_id_FL2img,'race4'},network_table_FL2{sub_id_FL2img,'interview_age'} ]]),2);

sub_id_FL2imgnan=sub_id_FL2img(temp_FL2==0);


sorted_FL2=dem_table{sub_id_FL2imgnan,'abcd_site'};
sorted_FL2(sorted_FL2==20)=1;
sorted_FL2(sorted_FL2==21)=17;
sorted_FL2(sorted_FL2==22)=19;



net_FL_combat=combat([network_table_FL2{sub_id_FL2imgnan,23:191} sub_network_table_FL2{sub_id_FL2imgnan,23:269}]',sorted_FL2',...
[dem_table{sub_id_FL2imgnan,[16 22]},network_table_FL2{sub_id_FL2imgnan,'interview_age'} ],1);% running combat to control the effects of those variables 

sub_network_table_FL_com=sub_network_table_FL2;
sub_network_table_FL_com{sub_id_FL2imgnan,23:269}=net_FL_combat(170:end,:)';


DTI_FL_combat=combat([DTI_table_FL2{sub_id_FL2imgnan,10:93}]',sorted_FL2',...
[dem_table{sub_id_FL2imgnan,'sex'},dem_table{sub_id_FL2imgnan,'race4'},network_table_FL2{sub_id_FL2imgnan,'interview_age'} ],1);% running combat to control the effects of those variables 

DTI_table_FL_com=DTI_table_FL2;
DTI_table_FL_com{sub_id_FL2imgnan,10:93}=DTI_FL_combat';




pubertal_f_FL2=sleep_parent_SDSC_FL2{sub_id_FL2imgnan,'pds_p_ss_female_category'};
pubertal_m_FL2=sleep_parent_SDSC_FL2{sub_id_FL2imgnan,'pds_p_ss_male_category'};
pubertal_FL2=sum([pubertal_f_FL2,pubertal_m_FL2],2,'omitnan');
pubertal_FL2(pubertal_FL2==0)=nan;
    

%% test DTI results dmdtifp1_28 dmdtifp1_29



 data_test=table;
    data_test.TS=psm_table{sub_id_BL,'TS'};
      data_test.TSD=log10(sleep_parent_SDSC_BL{sub_id_BL,'sds_p_ss_total'});
    data_test.IC=NIH_TB_BL{sub_id_BL,'nihtbx_flanker_uncorrected'};% inhibition control flaker task
    data_test.CS=NIH_TB_BL{sub_id_BL,'nihtbx_cardsort_uncorrected'};% cognitive flexibility card sorting task
    data_test.LS=NIH_TB_BL{sub_id_BL,'nihtbx_list_uncorrected'};%list sorting
    data_test.FA_sl_striatum=DTI_table_BL_com{sub_id_BL,'dmdtifp1_27'};
    data_test.FA_slp_striatum=DTI_table_BL_com{sub_id_BL,'dmdtifp1_31'};
    
    
    
    
    
    
    data_test.HI=psm_table{sub_id_BL,'HI'};
    data_test.PU=psm_table{sub_id_BL,'puberty'};
    data_test.BAS_drive=mental_health_BL{sub_id_BL,'bis_y_ss_bas_drive'};
    data_test.BAS_rr=mental_health_BL{sub_id_BL,'bis_y_ss_bas_rr'};% reward responsiness
    data_test.upps_nu=mental_health_BL{sub_id_BL,'upps_y_ss_negative_urgency'};
    data_test.upps_lplan=mental_health_BL{sub_id_BL,'upps_y_ss_lack_of_planning'};
    data_test.upps_ss=mental_health_BL{sub_id_BL,'upps_y_ss_sensation_seeking'};
    data_test.upps_pu=mental_health_BL{sub_id_BL,'upps_y_ss_positive_urgency'};
    data_test.upps_lper=mental_health_BL{sub_id_BL,'upps_y_ss_lack_of_perseverance'};

    DTI_motion=nan(length(sub_id_BL),1);
    for i=1:length(sub_id_BL)
        
    try
        DTI_motion(i,1)=DTI_roi_table_BL{sub_id_BL(i),'dmdtifp1_1183'};
    catch
        DTI_motion(i,1)=nan;
        continue
    end
    end
    data_test.DTI_motion=DTI_motion;
    data_test=[data_test,network_table_BL(sub_id_BL,'interview_age'),...
        dem_table(sub_id_BL,[16 23]),dem_table(sub_id_BL,[9 15])...
        ,dem_table(sub_id_BL,'edu'),network_table_BL(sub_id_BL,'rsfmri_c_ngd_ntpoints'),...
        network_table_BL(sub_id_BL,'rsfmri_c_ngd_meanmotion')...
        BMI_BL(sub_id_BL,'anthro_bmi_calc')];

    
    
    
    data_test.overall_SU_Parent=overall_SU;
    data_test.drug_use_pregenant=all_drug_pre;
    data_test.cerc_cdelh=sub_network_table_BL_com{sub_id_BL,'rsfmri_cor_ngd_cerc_scs_cdelh'};
    data_test.copa_plrh=sub_network_table_BL_com{sub_id_BL,'rsfmri_cor_ngd_copa_scs_plrh'};
    data_test.insomia=log10(sleep_parent_SDSC_BL{sub_id_BL,'sds_p_ss_dims'});
    data_test.sbd=log10(sleep_parent_SDSC_BL{sub_id_BL,'sds_p_ss_sbd'});
    
    
    
    data_test_FL2=table;
    data_test_FL2.TS=total_sleep_duration_FL2(ismember(sub_id_FL2img,sub_id_FL2imgnan));
  
    data_test_FL2.TSD=log10(sleep_parent_SDSC_FL2{sub_id_FL2imgnan,'sds_p_ss_total'});
    
%     data_test_FL2.HI=psm_table{sub_id_BL,'HI'};
    data_test_FL2.PU=pubertal_FL2;
    data_test_FL2.BAS_drive=mental_health_FL2{sub_id_FL2imgnan,'bis_y_ss_bas_drive'};
    data_test_FL2.BAS_rr=mental_health_FL2{sub_id_FL2imgnan,'bis_y_ss_bas_rr'};% reward responsiness
    data_test_FL2.upps_nu=mental_health_FL2{sub_id_FL2imgnan,'upps_y_ss_negative_urgency'};
    data_test_FL2.upps_lplan=mental_health_FL2{sub_id_FL2imgnan,'upps_y_ss_lack_of_planning'};
    data_test_FL2.upps_ss=mental_health_FL2{sub_id_FL2imgnan,'upps_y_ss_sensation_seeking'};
    data_test_FL2.upps_pu=mental_health_FL2{sub_id_FL2imgnan,'upps_y_ss_positive_urgency'};
    data_test_FL2.upps_lper=mental_health_FL2{sub_id_FL2imgnan,'upps_y_ss_lack_of_perseverance'};
    data_test_FL2.FA_sl_striatum=DTI_table_FL_com{sub_id_FL2imgnan,'dmdtifp1_27'};
    data_test_FL2.FA_slp_striatum=DTI_table_FL_com{sub_id_FL2imgnan,'dmdtifp1_31'};
    
    
    
    data_test_FL2=[data_test_FL2,network_table_FL2(sub_id_FL2imgnan,'interview_age'),...
        dem_table(sub_id_FL2imgnan,[16 23]),dem_table(sub_id_FL2imgnan,[9 15])...
        ,dem_table(sub_id_FL2imgnan,'edu'),network_table_FL2(sub_id_FL2imgnan,'rsfmri_c_ngd_ntpoints'),...
        network_table_FL2(sub_id_FL2imgnan,'rsfmri_c_ngd_meanmotion'),...
        data_test(sub_id_FL2imgnan,'overall_SU_Parent'),data_test(sub_id_FL2imgnan,'drug_use_pregenant')];
        
    %NDAR_INVTHXM3AU1
%     data_test_FL2.BMI=BMI_FL2(sub_id_FL2imgnan,'anthro_bmi_calc');
    
    data_test_FL2.cerc_cdelh=sub_network_table_FL_com{sub_id_FL2imgnan,'rsfmri_cor_ngd_cerc_scs_cdelh'};
    data_test_FL2.copa_plrh=sub_network_table_FL_com{sub_id_FL2imgnan,'rsfmri_cor_ngd_copa_scs_plrh'};
    data_test_FL2.PU_BL=psm_table{sub_id_FL2imgnan,'puberty'};
    data_test_FL2.BMI_BL=BMI_BL{sub_id_FL2imgnan,'anthro_bmi_calc'};
    data_test_FL2.HI=psm_table{sub_id_FL2imgnan,'HI'};
    
%%    

X_var=[data_test.TS,data_test.TSD];
mediators=[data_test.FA_sl_striatum,data_test.FA_slp_striatum,data_test.cerc_cdelh, data_test.copa_plrh];
Y_var=[data_test.BAS_rr,data_test.BAS_drive,...
    data_test.upps_nu,data_test.upps_pu,data_test.upps_ss,data_test.upps_lplan,data_test.upps_lper];
covs_fullvar=[data_test.interview_age,data_test.sex,data_test.race4,data_test.edu,data_test.HI,...
    data_test.PU,data_test.rsfmri_c_ngd_ntpoints,data_test.rsfmri_c_ngd_meanmotion,...
    data_test.abcd_site,data_test.rel_family_id,data_test.overall_SU_Parent,data_test.drug_use_pregenant,data_test.anthro_bmi_calc];
covs_var=[data_test.interview_age,data_test.sex,data_test.race4,data_test.edu,data_test.HI,...
    data_test.PU,...
    data_test.abcd_site,data_test.rel_family_id,data_test.overall_SU_Parent,data_test.drug_use_pregenant];
covs_fullvar_FL=[data_test_FL2{sub_id_FL2imgnan,'interview_age'},data_test{sub_id_FL2imgnan,'sex'},data_test{sub_id_FL2imgnan,'race4'},...
   data_test{sub_id_FL2imgnan,'HI'},data_test{sub_id_FL2imgnan,'edu'},...
    data_test{sub_id_FL2imgnan,'PU'},network_table_FL2{sub_id_FL2imgnan,'rsfmri_c_ngd_ntpoints'},...
        network_table_FL2{sub_id_FL2imgnan,'rsfmri_c_ngd_meanmotion'},...
    data_test{sub_id_FL2imgnan,'abcd_site'},data_test{sub_id_FL2imgnan,'rel_family_id'},...
    data_test{sub_id_FL2imgnan,'overall_SU_Parent'},data_test{sub_id_FL2imgnan,'drug_use_pregenant'},data_test_FL2.BMI_BL];


mediators_FL=[data_test_FL2.cerc_cdelh, data_test_FL2.copa_plrh];
Y_var_FL2=[data_test_FL2.BAS_rr,data_test_FL2.BAS_drive,...
    data_test_FL2.upps_nu,data_test_FL2.upps_pu,data_test_FL2.upps_ss,data_test_FL2.upps_lplan,data_test_FL2.upps_lper];



p_lmm=nan(2,14);
t_lmm=nan(2,14);
X_var_name={'TS','TSD'};
Y_var_name=data_test.Properties.VariableNames([6     7     29    30      3     5     4    11    10    12    15    15    13    16 ]);
for i=1:2
    for j=1:14
    lme=fitlme(data_test,[Y_var_name{j},'~',X_var_name{i} '+interview_age+sex+race4+edu+HI+PU+anthro_bmi_calc+rsfmri_c_ngd_ntpoints+rsfmri_c_ngd_meanmotion+overall_SU_Parent+drug_use_pregenant+(1|abcd_site)+(1|abcd_site:rel_family_id)']);
%     lme=fitlme(data_test,[Y_var_name{j},'~',X_var_name{i} '+interview_age+sex+race4+edu+HI+PU+anthro_bmi_calc+(1|abcd_site)+(1|abcd_site:rel_family_id)']);
    [~,~,stats]=fixedEffects(lme);
    test_table=dataset2table(stats);
    test_table.Properties.RowNames=test_table.Name;
    p_lmm(i,j)=test_table{X_var_name{i},'pValue'};
    
    t_lmm(i,j)=test_table{X_var_name{i},'tStat'};
    end
end

p_lmm_fc=nan(3,10);
for i=1:5
    for j=1:10
    lme=fitlme(data_test,[Y_var_name{j+9},'~',Y_var_name{i+4} '+interview_age+sex+race4+edu+HI+PU+anthro_bmi_calc+rsfmri_c_ngd_ntpoints+rsfmri_c_ngd_meanmotion+overall_SU_Parent+drug_use_pregenant+(1|abcd_site)+(1|abcd_site:rel_family_id)']);
    [~,~,stats]=fixedEffects(lme);
    test_table=dataset2table(stats);
    test_table.Properties.RowNames=test_table.Name;
    p_lmm_fc(i,j)=test_table{Y_var_name{i+4},'pValue'};
    
    end
end




p_lmm_FL=nan(2,9);% follow-up
t_lmm_FL=nan(2,9);
X_var_name={'TS','TSD'};
Y_var_name_FL=data_test_FL2.Properties.VariableNames([23 24 5 4 6 9 8 7 10  ]);
for i=1
    for j=1:9
    lme=fitlme(data_test_FL2,[Y_var_name_FL{j},'~',X_var_name{i} '+interview_age+sex+race4+edu+HI+PU_BL+BMI_BL+rsfmri_c_ngd_ntpoints+rsfmri_c_ngd_meanmotion+overall_SU_Parent+drug_use_pregenant+(1|abcd_site)+(1|abcd_site:rel_family_id)']);
%     lme=fitlme(data_test,[Y_var_name{j},'~',X_var_name{i} '+interview_age+sex+race4+edu+HI+PU+anthro_bmi_calc+(1|abcd_site)+(1|abcd_site:rel_family_id)']);
    [~,~,stats]=fixedEffects(lme);
    test_table=dataset2table(stats);
    test_table.Properties.RowNames=test_table.Name;
    p_lmm_FL(i,j)=test_table{X_var_name{i},'pValue'};
    
    t_lmm_FL(i,j)=test_table{X_var_name{i},'tStat'};
    end
end







p_med=nan(5,10,2);
beta_med=nan(5,10,2);
for i=1:2
    for j=1:4
        for k=1:7
                
[paths, stats_med]=mediation(X_var(:,i),...
Y_var(:,k),mediators(:,j)...
,'cov', [covs_fullvar],'verbose','boot','bootsamples',10000);
p_med(j,k,i)=stats_med.p(5);
beta_med(j,k,i)=paths(5);
        end
    end
end




p_med_FL2=nan(2,7);
beta_med_FL2=nan(2,7);
for i=1:2
    for j=1:7
        stats_med=[];paths=[];
                [paths, stats_med]=mediation(data_test_FL2.TS,...
Y_var_FL2(:,j),mediators_FL(:,i)...
,'cov', [covs_fullvar_FL],'verbose','boot','bootsamples',10000);
p_med_FL2(i,j)=stats_med.p(5);
beta_med_FL2(i,j)=paths(5);            
       
    end
end



p_med_FL=[];
for i=1
    for j=1:5
        
                [paths, stats_med]=mediation(X_var(ismember(sub_id_BL,sub_id_FL2imgnan),1),...
mental_health_FL2{sub_id_FL2imgnan,'upps_y_ss_positive_urgency'},data_test{sub_id_FL2imgnan,'copa_plrh'}...
,'cov', [covs_fullvar(ismember(sub_id_BL,sub_id_FL2imgnan),:),mental_health_BL{sub_id_FL2imgnan,'upps_y_ss_positive_urgency'}],'verbose','boot','bootsamples',10000);

p_med_FL(1,2)=stats_med.p(5);

beta_med_FL(1,2)=paths(5);

                [paths, stats_med]=mediation(X_var(ismember(sub_id_BL,sub_id_FL2imgnan),1),...
mental_health_FL2{sub_id_FL2imgnan,'upps_y_ss_positive_urgency'},data_test{sub_id_FL2imgnan,'cerc_cdelh'}...
,'cov', [covs_fullvar(ismember(sub_id_BL,sub_id_FL2imgnan),:),mental_health_BL{sub_id_FL2imgnan,'upps_y_ss_positive_urgency'}],'verbose','boot','bootsamples',10000);
p_med_FL(1,1)=stats_med.p(5);
beta_med_FL(1,1)=paths(5);           
       

 [paths, stats_med]=mediation(mental_health_BL{sub_id_FL2imgnan,'upps_y_ss_positive_urgency'},...
    data_test_FL2.TS,data_test{sub_id_FL2imgnan,'copa_plrh'}...
,'cov', [covs_fullvar(ismember(sub_id_BL,sub_id_FL2imgnan),:),...
data_test{ismember(sub_id_BL,sub_id_FL2imgnan),'TS'}],'verbose','boot','bootsamples',10000)
p_med_FL(2,2)=stats_med.p(5);
beta_med_FL(2,2)=paths(5);  

 [paths, stats_med]=mediation(mental_health_BL{sub_id_FL2imgnan,'upps_y_ss_positive_urgency'},...
    data_test_FL2.TS,data_test{sub_id_FL2imgnan,'cerc_cdelh'}...
,'cov', [covs_fullvar(ismember(sub_id_BL,sub_id_FL2imgnan),:),...
data_test{ismember(sub_id_BL,sub_id_FL2imgnan),'TS'}],'verbose','boot','bootsamples',10000)

p_med_FL(2,1)=stats_med.p(5);
beta_med_FL(2,1)=paths(5);       



%negative urgency
[paths, stats_med]=mediation(X_var(ismember(sub_id_BL,sub_id_FL2imgnan),1),...
mental_health_FL2{sub_id_FL2imgnan,'upps_y_ss_negative_urgency'},data_test{sub_id_FL2imgnan,'copa_plrh'}...
,'cov', [covs_fullvar(ismember(sub_id_BL,sub_id_FL2imgnan),:),mental_health_BL{sub_id_FL2imgnan,'upps_y_ss_negative_urgency'}],'verbose','boot','bootsamples',10000);

p_med_FL(3,2)=stats_med.p(5);

beta_med_FL(3,2)=paths(5);

                [paths, stats_med]=mediation(X_var(ismember(sub_id_BL,sub_id_FL2imgnan),1),...
mental_health_FL2{sub_id_FL2imgnan,'upps_y_ss_negative_urgency'},data_test{sub_id_FL2imgnan,'cerc_cdelh'}...
,'cov', [covs_fullvar(ismember(sub_id_BL,sub_id_FL2imgnan),:),mental_health_BL{sub_id_FL2imgnan,'upps_y_ss_negative_urgency'}],'verbose','boot','bootsamples',10000);
p_med_FL(3,1)=stats_med.p(5);
beta_med_FL(3,1)=paths(5);           
       

 [paths, stats_med]=mediation(mental_health_BL{sub_id_FL2imgnan,'upps_y_ss_negative_urgency'},...
    data_test_FL2.TS,data_test{sub_id_FL2imgnan,'copa_plrh'}...
,'cov', [covs_fullvar(ismember(sub_id_BL,sub_id_FL2imgnan),:),...
data_test{ismember(sub_id_BL,sub_id_FL2imgnan),'TS'}],'verbose','boot','bootsamples',10000)
p_med_FL(4,2)=stats_med.p(5);
beta_med_FL(4,2)=paths(5);  

 [paths, stats_med]=mediation(mental_health_BL{sub_id_FL2imgnan,'upps_y_ss_negative_urgency'},...
    data_test_FL2.TS,data_test{sub_id_FL2imgnan,'cerc_cdelh'}...
,'cov', [covs_fullvar(ismember(sub_id_BL,sub_id_FL2imgnan),:),...
data_test{ismember(sub_id_BL,sub_id_FL2imgnan),'TS'}],'verbose','boot','bootsamples',10000)

p_med_FL(4,1)=stats_med.p(5);
beta_med_FL(4,1)=paths(5);      








    end
end




p_med_fa=nan(4,10,2);
for i=1:2
    for j=1:4
        for k=1:10
[paths, stats_med]=mediation(X_var(:,i),...
Y_var(:,k),mediators(:,j)...
,'cov', [covs_var],'verbose','boot','bootsamples',10000);%DTI without motion
p_med_fa(j,k,i)=stats_med.p(5);

        end
    end
end



p_value_shaped=[reshape(p_value(1:169),13,13);reshape(p_value(170:end),9,13);];
p_value_shaped=[tril(p_value_shaped(1:13,1:13));p_value_shaped(14:end,:)];
sum(fdr_bh(p_value_shaped(p_value_shaped>0),0.05))
t_stats_reshaped=[reshape(t_stats(1:169),13,13);reshape(t_stats(170:end),9,13);];
t_stats_reshaped=[tril(t_stats_reshaped(1:13,1:13));t_stats_reshaped(14:end,:)];

t_stats_reshaped=t_stats_reshaped([1:6 8:end],[1:6 8:end]);     




%% follow-up
 


data_CLPM=table;

data_CLPM.TS_1=data_test{sub_id_FL2imgnan,'TS'};
data_CLPM.TS_2=data_test_FL2.TS;
data_CLPM.cdelh_1=data_test{sub_id_FL2imgnan,'cerc_cdelh'};
data_CLPM.cdelh_2=data_test_FL2.cerc_cdelh;
data_CLPM.plrh_1=data_test{sub_id_FL2imgnan,'copa_plrh'};
data_CLPM.plrh_2=data_test_FL2.copa_plrh;
data_CLPM.upps_pu_1=data_test{sub_id_FL2imgnan,'upps_pu'};
data_CLPM.upps_pu_2=data_test_FL2.upps_pu;
data_CLPM.upps_nu_1=data_test{sub_id_FL2imgnan,'upps_nu'};
data_CLPM.upps_nu_2=data_test_FL2.upps_nu;
data_CLPM.upps_lper_1=data_test{sub_id_FL2imgnan,'upps_lper'};
data_CLPM.upps_lper_2=data_test_FL2.upps_lper;
data_CLPM.BAS_drive_1=data_test{sub_id_FL2imgnan,'BAS_drive'};
data_CLPM.BAS_drive_2=data_test_FL2.BAS_drive;



covs_fullvar=[data_test.interview_age,data_test.sex,data_test.race4,data_test.HI,data_test.edu,...
    data_test.PU,data_test.rsfmri_c_ngd_ntpoints,data_test.rsfmri_c_ngd_meanmotion,...
    data_test.abcd_site,data_test.rel_family_id,data_test.overall_SU_Parent,data_test.drug_use_pregenant];


covs_fullvar_FL2=[data_test{sub_id_FL2imgnan,'interview_age'},data_test{sub_id_FL2imgnan,'sex'},data_test{sub_id_FL2imgnan,'race4'},...
   data_test{sub_id_FL2imgnan,'HI'},data_test{sub_id_FL2imgnan,'edu'},...
    data_test{sub_id_FL2imgnan,'PU'},network_table_FL2{sub_id_FL2imgnan,'rsfmri_c_ngd_ntpoints'},...
        network_table_FL2{sub_id_FL2imgnan,'rsfmri_c_ngd_meanmotion'},...
    data_test{sub_id_FL2imgnan,'abcd_site'},data_test{sub_id_FL2imgnan,'rel_family_id'},...
    data_test{sub_id_FL2imgnan,'overall_SU_Parent'},data_test{sub_id_FL2imgnan,'drug_use_pregenant'}];

data_CLPM_cov=data_CLPM;
covs_fullvar_BL=covs_fullvar(ismember(sub_id_BL,sub_id_FL2imgnan),:);
tempv=[];
for i=1:2:14
    tempv=data_CLPM{:,i};
ind_nn=~isnan(tempv);

b=regress(tempv(ind_nn),[ones(length(tempv(ind_nn)),1) covs_fullvar_BL(ind_nn,:)]);
Yhat=[ones(length(tempv(ind_nn)),1) covs_fullvar_BL(ind_nn,:)]*b;
resid(ind_nn)=tempv(ind_nn)-Yhat+mean(tempv(ind_nn));
resid(~ind_nn)=nan;
data_CLPM_cov{:,i}=resid';
end

tempv=[];
for i=2:2:14
    tempv=data_CLPM{:,i};
ind_nn=~isnan(tempv);

b=regress(tempv(ind_nn),[ones(length(tempv(ind_nn)),1) covs_fullvar_FL2(ind_nn,:)]);
Yhat=[ones(length(tempv(ind_nn)),1) covs_fullvar_FL2(ind_nn,:)]*b;
resid(ind_nn)=tempv(ind_nn)-Yhat+mean(tempv(ind_nn));
resid(~ind_nn)=nan;
data_CLPM_cov{:,i}=resid';
end



data_CLPM_cov=data_CLPM_cov(~isnan(sum(data_CLPM_cov{:,1:14},2)),:);

% data_CLPM_cov.upps_1=(data_CLPM_cov.upps_nu_1+data_CLPM_cov.upps_pu_1)/2;
% data_CLPM_cov.upps_2=(data_CLPM_cov.upps_nu_2+data_CLPM_cov.upps_pu_2)/2;



writetable(data_CLPM,'CLPM.csv')

writetable(data_CLPM_cov,'CLPM_cov.csv')


%% tstat behavioral plot
figure;
b=bar(t_lmm(1,[1 2 8 9 ]),'FaceColor','flat');

for i=1:2
b.CData(i,:) = [0 0 1];
end


for i=3:4
b.CData(i,:) = [0 1 0];
end

hold on
% yline([0.15],'--')
set(gca, 'FontName', 'Helvetica')
xticks(1:42);
xticklabels(strrep(Y_var_name([1 2 8 9 ]),'_','\_'))
ylabel({'t stats'},'FontName', 'Helvetica')
title({'Brain Connectivity'})
set(gca,'TickLength',[0 0])
set(gca,'LineWidth',2)
set(gca,'fontsize',18)
set(gcf,'units','inches','position',[0, 0, 9, 6])
box off
printeps(1,'brain');

figure;
b=bar(t_lmm(1,end-4:end),'FaceColor','flat');

hold on
% yline([0.15],'--')
set(gca, 'FontName', 'Helvetica')
xticks(1:42);
xticklabels(strrep(Y_var_name(end-4:end),'_','\_'))
ylabel({'t stats'},'FontName', 'Helvetica')
title({'Impulsivity'})
set(gca,'TickLength',[0 0])
set(gca,'LineWidth',2)
set(gca,'fontsize',18)
set(gcf,'units','inches','position',[0, 0, 9, 6])
box off
printeps(1,'Impulsivity');

figure;
b=bar(t_lmm_FL(1,[1 2 5 6 ]),'FaceColor','flat');

for i=1:2
b.CData(i,:) = [0 0 1];
end


for i=3:4
b.CData(i,:) = [0 1 0];
end

hold on
% yline([0.15],'--')
set(gca, 'FontName', 'Helvetica')
xticks(1:42);
xticklabels(strrep(Y_var_name_FL([1 2 5 6]),'_','\_'))
ylabel({'t stats'},'FontName', 'Helvetica')
title({'Brain Connectivity'})
set(gca,'TickLength',[0 0])
set(gca,'LineWidth',2)
set(gca,'fontsize',18)
set(gcf,'units','inches','position',[0, 0, 9, 6])
box off
printeps(1,'brain_FL');




figure;
y=[nanmean(data_test{sub_id_FL2imgnan,[14 17 16 15 18]});...
    nanmean(data_test_FL2{sub_id_FL2imgnan,[6 9 8 7 10]})]';
y_std=[nanstd(data_test{sub_id_FL2imgnan,[14 17 16 15 18]})./sqrt(length(sub_id_FL2imgnan));...
    nanstd(data_test_FL2{sub_id_FL2imgnan,[6 9 8 7 10]})./sqrt(length(sub_id_FL2imgnan))]';
barwitherr(y_std,y)
set(gca, 'FontName', 'Helvetica')
% xticks(1:42);
ylim([0 12])
xticklabels(strrep(Y_var_name_FL([5:9]),'_','\_'))
% ylabel({'t stats'},'FontName', 'Helvetica')
% title({'Brain Connectivity'})
set(gca,'TickLength',[0 0])
set(gca,'LineWidth',2)
set(gca,'fontsize',18)
set(gcf,'units','inches','position',[0, 0, 9, 6])
legend('Baseline','FL2')
box off
printeps(1,'Impulsivity_BL_FL');

    


figure;
bar(beta_med(3,[6:7],1),'w','LineWidth',2)
xticklabels({'upps\_nu','upps\_pu' });
ylabel({'beta of mediation effect (a*b)'},'fontsize',24)
ylim([-0.04 0])
set(gca,'TickLength',[0 0])
set(gca, 'XAxisLocation', 'top')
set(gca,'LineWidth',2)
set(gca,'fontsize',18)
set(gcf,'units','inches','position',[0, 0, 9, 6])
% title({'CON\_CAU'})
box off
printeps(1,'mediation_bl_cau')

figure;
bar(beta_med(4,[6:7],1),'w','LineWidth',2)
xticklabels({'upps\_nu','upps\_pu' });
ylabel({'beta of mediation effect (a*b)'},'fontsize',24)
ylim([-0.04 0])
set(gca,'TickLength',[0 0])
set(gca, 'XAxisLocation', 'top')
set(gca,'LineWidth',2)
set(gca,'fontsize',18)
set(gcf,'units','inches','position',[0, 0, 9, 6])
% title({'CON\_CAU'})

box off
printeps(1,'mediation_bl_gp')


figure;
bar(beta_med_FL(1,:),'w','LineWidth',2)
xticklabels({'cerc\_cdelh','copa\_plrh' });
ylabel({'beta of mediation effect (a*b)'},'fontsize',24)
ylim([-0.04 0])
set(gca,'TickLength',[0 0])
set(gca, 'XAxisLocation', 'top')
set(gca,'LineWidth',2)
set(gca,'fontsize',18)
set(gcf,'units','inches','position',[0, 0, 9, 6])
% title({'CON\_CAU'})
box off
printeps(1,'mediation_fl_ts_pu')

figure;
bar(beta_med_FL(2,:),'w','LineWidth',2)
xticklabels({'cerc\_cdelh','copa\_plrh' });
ylabel({'beta of mediation effect (a*b)'},'fontsize',24)
ylim([-0.003 0])
set(gca,'TickLength',[0 0])
set(gca, 'XAxisLocation', 'top')
set(gca,'LineWidth',2)
set(gca,'fontsize',18)
set(gcf,'units','inches','position',[0, 0, 9, 6])
% title({'CON\_CAU'})
box off
printeps(1,'mediation_fl_pu_ts')



figure;
bar(beta_med_FL2(1,3:4),'w','LineWidth',2)
xticklabels({'upps\_nu','upps\_pu' });
ylabel({'beta of mediation effect (a*b)'},'fontsize',24)
ylim([-0.04 0])
set(gca,'TickLength',[0 0])
set(gca, 'XAxisLocation', 'top')
set(gca,'LineWidth',2)
set(gca,'fontsize',18)
set(gcf,'units','inches','position',[0, 0, 9, 6])
% title({'CON\_CAU'})
box off
printeps(1,'mediation_fl_cau')

figure;
bar(beta_med_FL2(2,[3:4]),'w','LineWidth',2)
xticklabels({'upps\_nu','upps\_pu' });
ylabel({'beta of mediation effect (a*b)'},'fontsize',24)
ylim([-0.04 0])
set(gca,'TickLength',[0 0])
set(gca, 'XAxisLocation', 'top')
set(gca,'LineWidth',2)
set(gca,'fontsize',18)
set(gcf,'units','inches','position',[0, 0, 9, 6])
% title({'CON\_CAU'})

box off
printeps(1,'mediation_fl_gp')









figure1=figure;
axes1 = axes('Parent',figure1);
hold(axes1,'all');
 plot1=plot(cohend(1:20)',t_lmm(1:4)','Parent',axes1,'LineStyle','none','Marker','o','Markersize',9,'MarkerFaceColor','b','DisplayName','Structrual connec')
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







for i = 1: 4
    y(i,:)=mean(data_test{data_test.TS==i,[12 15]},'omitnan');
    y_std(i,:)=std(data_test{data_test.TS==i,[12 15]},'omitnan')./sqrt(sum(~isnan(data_test{data_test.TS==i,[12 15]})));
end
figure;
barwitherr(y_std',y')
set(gca, 'FontName', 'Helvetica')
xticks(1:2);
ylim([5 10])
xticklabels(strrep({'upps_nu','upps_pu'},'_','\_'))
% ylabel({'t stats'},'FontName', 'Helvetica')
% title({'Brain Connectivity'})
set(gca,'TickLength',[0 0])
set(gca,'LineWidth',2)
set(gca,'fontsize',18)
set(gcf,'units','inches','position',[0, 0, 9, 6])
legend('<7','7-8','8-9','>9')
box off
printeps(1,'PU_NU_Sleep_Duration');



for i = 1: 5
    n_pu(i,1)=sum(data_test{:,'PU'}==i);
    n_pu(i,2)=sum(data_test_FL2{:,'PU'}==i);
end

[t,p]=corr(data_test.cerc_cdelh, data_test.upps_pu,'row','pairwise')
[t,p]=corr(data_test.copa_plrh, data_test.upps_nu,'row','pairwise')




%% demographic table

data_table=table('Size',[9 4],'VariableTypes',["string" "string"  "string" "string"],...
    'VariableNames',{'Characteristic','Mean (SD)','Number of missing values','n'});

summary_table_SDS=data_test(:,[20 21 11  32 33 22]);
summary_table_SDS.overall_SU_Parent=summary_table_SDS.overall_SU_Parent>0;

for i = 1:5 
    
    
    
    data_table{i,2}={sprintf('%3.2f (%3.2f)',nanmean(summary_table_SDS{sub_id_BL,i}),nanstd(summary_table_SDS{sub_id_BL,i}))};
    data_table{i,3}={sprintf('%3.0f',sum(isnan(summary_table_SDS{sub_id_BL,i})))};
    data_table{i,4}={sprintf('%3.2f (%3.3f)',sum(summary_table_SDS{sub_id_BL,i}),sum(summary_table_SDS{sub_id_BL,i})/height(summary_table_SDS))};
   
end

 
% i=10;
%  data_table{i,2}={sprintf('%3.0f (%3.0f)',sum(summary_table_SDS{sub_id_BL,2}==0),sum(summary_table_SDS{sub_id_BL,2}==1))};


 data_table{6,2}={sprintf('%3.0f ',sum(summary_table_SDS{sub_id_BL,6}==1))};
 data_table{7,2}={sprintf('%3.0f ',sum(summary_table_SDS{sub_id_BL,6}==2))};
  data_table{8,2}={sprintf('%3.0f ',sum(summary_table_SDS{sub_id_BL,6}==3))};
  data_table{9,2}={sprintf('%3.0f ',sum(summary_table_SDS{sub_id_BL,6}==4))};

data_table{1:6,1}=summary_table_SDS.Properties.VariableNames';





%% 2 year follow-up



data_table=table('Size',[9 4],'VariableTypes',["string" "string"  "string" "string"],...
    'VariableNames',{'Characteristic','Mean (SD)','Number of missing values','n'});


summary_table_FL2=data_test_FL2(:,[11 12 3 19 20 13]);
summary_table_FL2.overall_SU_Parent=summary_table_FL2.overall_SU_Parent>0;



for i = 1:5
    
    
    data_table{i,2}={sprintf('%3.2f (%3.2f)',nanmean(summary_table_FL2{sub_id_FL2imgnan,i}),nanstd(summary_table_FL2{sub_id_FL2imgnan,i}))};
    data_table{i,3}={sprintf('%3.0f',sum(isnan(summary_table_FL2{sub_id_FL2imgnan,i})))};
     data_table{i,4}={sprintf('%3.2f (%3.3f)',sum(summary_table_FL2{sub_id_FL2imgnan,i}),sum(summary_table_FL2{sub_id_FL2imgnan,i})/height(summary_table_FL2))};
   
    

    
   
end

 
% i=10;
%  data_table{i,2}={sprintf('%3.0f (%3.0f)',sum(summary_table_FL2{sub_id_FL2imgnan,2}==0),sum(summary_table_FL2{sub_id_FL2imgnan,2}==1))};
% 

 data_table{6,2}={sprintf('%3.0f ',sum(summary_table_FL2{sub_id_FL2imgnan,6}==1))};
 data_table{7,2}={sprintf('%3.0f ',sum(summary_table_FL2{sub_id_FL2imgnan,6}==2))};
  data_table{8,2}={sprintf('%3.0f ',sum(summary_table_FL2{sub_id_FL2imgnan,6}==3))};
data_table{9,2}={sprintf('%3.0f ',sum(summary_table_FL2{sub_id_FL2imgnan,6}==4))}

data_table{1:6,1}=summary_table_FL2.Properties.VariableNames';