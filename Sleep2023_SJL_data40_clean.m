%% initial parameters

father_dir='/home/yangf7/Documents/Nils/ABCD';%working dir
cd(father_dir)
beh_dir=fullfile(father_dir,'Beh_tabulated_data_40/');%data folder

network_table=readtable(fullfile(beh_dir,'abcd_betnet02.txt'));%load netowrk data
network_table_FL2=network_table(strcmp(network_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
network_table_BL=network_table(strcmp(network_table{:,'eventname'},'baseline_year_1_arm_1'),:);
network_table_BL.Properties.RowNames=network_table_BL.src_subject_id;
network_table_FL2.Properties.RowNames=network_table_FL2.src_subject_id;
qc_inclusion=readtable(fullfile(beh_dir, 'abcd_imgincl01.txt'));
qc_inclusion_BL=qc_inclusion(strcmp(qc_inclusion{:,'eventname'},'baseline_year_1_arm_1'),:);
qc_inclusion_FL2=qc_inclusion(strcmp(qc_inclusion{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
qc_inclusion_BL.Properties.RowNames=qc_inclusion_BL.src_subject_id;
qc_inclusion_FL2.Properties.RowNames=qc_inclusion_FL2.src_subject_id;

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
DTI_table=readtable(fullfile(beh_dir,'abcd_dmdtifp101.txt'));
DTI_table_BL=DTI_table(strcmp(DTI_table{:,'eventname'},'baseline_year_1_arm_1'),:);
DTI_table_BL.Properties.RowNames=DTI_table_BL.src_subject_id;
DTI_table_FL2=DTI_table(strcmp(DTI_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
DTI_table_FL2.Properties.RowNames=DTI_table_FL2.src_subject_id;

NIH_TB=readtable(fullfile(beh_dir,'abcd_tbss01.txt'));% total score (sds_p_ss_total) for NIH toolbox 
NIH_TB_BL=NIH_TB(strcmp(NIH_TB{:,'eventname'},'baseline_year_1_arm_1'),:);
NIH_TB_BL.Properties.RowNames=NIH_TB_BL.src_subject_id;
NIH_TB_FL2=NIH_TB(strcmp(NIH_TB{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
NIH_TB_FL2.Properties.RowNames=NIH_TB_FL2.src_subject_id;
NIH_TB_FL1=NIH_TB(strcmp(NIH_TB{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
NIH_TB_FL1.Properties.RowNames=NIH_TB_FL1.src_subject_id;
NIH_TB_FL3=NIH_TB(strcmp(NIH_TB{:,'eventname'},'3_year_follow_up_y_arm_1'),:);
NIH_TB_FL3.Properties.RowNames=NIH_TB_FL3.src_subject_id;




sub_id_BL=NIH_TB_BL.Properties.RowNames;% subject ID for baseline



sub_network_table=readtable(fullfile(beh_dir,'mrirscor02.txt'));% network FC for subcortical-cortical connections
sub_network_table_FL2=sub_network_table(strcmp(sub_network_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
sub_network_table_BL=sub_network_table(strcmp(sub_network_table{:,'eventname'},'baseline_year_1_arm_1'),:);
sub_network_table_BL.Properties.RowNames=sub_network_table_BL.src_subject_id;
sub_network_table_FL2.Properties.RowNames=sub_network_table_FL2.src_subject_id;

cbcl_sum=readtable(fullfile(beh_dir,'abcd_cbcls01.txt'));% children behavior checklists
cbcl_sum_FL2=cbcl_sum(strcmp(cbcl_sum{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
cbcl_sum_BL=cbcl_sum(strcmp(cbcl_sum{:,'eventname'},'baseline_year_1_arm_1'),:);
cbcl_sum_BL.Properties.RowNames=cbcl_sum_BL.src_subject_id;
cbcl_sum_FL2.Properties.RowNames=cbcl_sum_FL2.src_subject_id;
cbcl_sum_FL1=cbcl_sum(strcmp(cbcl_sum{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
cbcl_sum_FL1.Properties.RowNames=cbcl_sum_FL1.src_subject_id;
cbcl_sum_FL3=cbcl_sum(strcmp(cbcl_sum{:,'eventname'},'3_year_follow_up_y_arm_1'),:);
cbcl_sum_FL3.Properties.RowNames=cbcl_sum_FL3.src_subject_id;



sleep_parent_SDSC=readtable(fullfile(beh_dir,'abcd_ssphp01.txt'));% total score (sds_p_ss_total) for SDSC 
sleep_parent_SDSC_BL=sleep_parent_SDSC(strcmp(sleep_parent_SDSC{:,'eventname'},'baseline_year_1_arm_1'),:);
sleep_parent_SDSC_BL.Properties.RowNames=sleep_parent_SDSC_BL.src_subject_id;
sleep_parent_SDSC_FL2=sleep_parent_SDSC(strcmp(sleep_parent_SDSC{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
sleep_parent_SDSC_FL2.Properties.RowNames=sleep_parent_SDSC_FL2.src_subject_id;
sleep_parent_SDSC_FL1=sleep_parent_SDSC(strcmp(sleep_parent_SDSC{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
sleep_parent_SDSC_FL1.Properties.RowNames=sleep_parent_SDSC_FL1.src_subject_id;
sleep_parent_SDSC_FL3=sleep_parent_SDSC(strcmp(sleep_parent_SDSC{:,'eventname'},'3_year_follow_up_y_arm_1'),:);
sleep_parent_SDSC_FL3.Properties.RowNames=sleep_parent_SDSC_FL3.src_subject_id;





sleep_parent_SDS=readtable(fullfile(beh_dir,'abcd_sds01.txt'));% total hours of sleep (sleep_1_p) for SDS 
%(Parent Sleep Disturbance Scale for Children)
sleep_parent_SDS_BL=sleep_parent_SDS(strcmp(sleep_parent_SDS{:,'eventname'},'baseline_year_1_arm_1'),:);
sleep_parent_SDS_BL.Properties.RowNames=sleep_parent_SDS_BL.src_subject_id;
sleep_parent_SDS_FL2=sleep_parent_SDS(strcmp(sleep_parent_SDS{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
sleep_parent_SDS_FL2.Properties.RowNames=sleep_parent_SDS_FL2.src_subject_id;
sleep_parent_SDS_FL1=sleep_parent_SDS(strcmp(sleep_parent_SDS{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
sleep_parent_SDS_FL1.Properties.RowNames=sleep_parent_SDS_FL1.src_subject_id;





dem_table=readtable(fullfile(beh_dir,'ABCD_sites_effect.xlsx'));% load demgraphic information
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




BMI=readtable(fullfile(beh_dir,'DEAP-data-download.csv'));% BMI
BMI_BL=BMI(strcmp(BMI{:,'event_name'},'baseline_year_1_arm_1'),:);
BMI_BL.Properties.RowNames=BMI_BL.src_subject_id;
BMI_FL2=BMI(strcmp(BMI{:,'event_name'},'2_year_follow_up_y_arm_1'),:);
BMI_FL2.Properties.RowNames=BMI_FL2.src_subject_id;
BMI_FL3=BMI(strcmp(BMI{:,'event_name'},'3_year_follow_up_y_arm_1'),:);
BMI_FL3.Properties.RowNames=BMI_FL3.src_subject_id;
BMI_BL.anthro_bmi_calc=str2double(BMI_BL.anthro_bmi_calc);
BMI_FL2.anthro_bmi_calc=str2double(BMI_FL2.anthro_bmi_calc);




mental_health=readtable(fullfile(beh_dir,'abcd_mhy02.txt'));% mental health variables
mental_health_BL=mental_health(strcmp(mental_health{:,'eventname'},'baseline_year_1_arm_1'),:);
mental_health_BL.Properties.RowNames=mental_health_BL.src_subject_id;
mental_health_FL1=mental_health(strcmp(mental_health{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
mental_health_FL1.Properties.RowNames=mental_health_FL1.src_subject_id;
mental_health_FL2=mental_health(strcmp(mental_health{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
mental_health_FL2.Properties.RowNames=mental_health_FL2.src_subject_id;
mental_health_FL3=mental_health(strcmp(mental_health{:,'eventname'},'3_year_follow_up_y_arm_1'),:);
mental_health_FL3.Properties.RowNames=mental_health_FL3.src_subject_id;




grade_table=readtable(fullfile(beh_dir,'abcd_ysaag01.txt'));% grades in school
grade_table_BL=grade_table(strcmp(grade_table{:,'eventname'},'baseline_year_1_arm_1'),:);
grade_table_BL.Properties.RowNames=grade_table_BL.src_subject_id;
grade_table_FL1=grade_table(strcmp(grade_table{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
grade_table_FL1.Properties.RowNames=grade_table_FL1.src_subject_id;
grade_table_FL2=grade_table(strcmp(grade_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
grade_table_FL2.Properties.RowNames=grade_table_FL2.src_subject_id;
grade_table_FL3=grade_table(strcmp(grade_table{:,'eventname'},'3_year_follow_up_y_arm_1'),:);
grade_table_FL3.Properties.RowNames=grade_table_FL3.src_subject_id;

%fitbit data processing

fitbit_PA=readtable(fullfile(beh_dir,'abcd_fbdpas01.txt'));
fitbit_sleep=readtable(fullfile(beh_dir,'abcd_fbdss01.txt'));
fitbit_sleep_BL=fitbit_sleep(strcmp(fitbit_sleep{:,'eventname'},'baseline_year_1_arm_1'),:);
% fitbit_sleep_week_BL.Properties.RowNames=fitbit_sleep_week_BL.src_subject_id;
fitbit_sleep_FL2=fitbit_sleep(strcmp(fitbit_sleep{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
fitbit_sleep_FL3=fitbit_sleep(strcmp(fitbit_sleep{:,'eventname'},'3_year_follow_up_y_arm_1'),:);

fitbit_sleep_week=readtable(fullfile(beh_dir,'abcd_fbwss01.txt'));
fitbit_sleep_week_BL=fitbit_sleep_week(strcmp(fitbit_sleep_week{:,'eventname'},'baseline_year_1_arm_1'),:);
% fitbit_sleep_week_BL.Properties.RowNames=fitbit_sleep_week_BL.src_subject_id;
fitbit_sleep_week_FL2=fitbit_sleep_week(strcmp(fitbit_sleep_week{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
% fitbit_sleep_week_FL2.Properties.RowNames=fitbit_sleep_week_FL2.src_subject_id;
sub_id_FB_bl=unique(fitbit_sleep_BL.src_subject_id);%79 have both timepoint
sub_id_FB_FL2=unique(fitbit_sleep_FL2.src_subject_id);
% 

fitbit_daily_BL=fitbit_sleep_BL;% baseline data only have 100-ish children

fitbit_daily_BL=fitbit_daily_BL(fitbit_daily_BL.fit_ss_protocol_wear==1,:);
sub_id_BL7=unique(fitbit_daily_BL.src_subject_id);
sub_id_7_BL=[];
for i=1:length(sub_id_BL7)
   if sum(ismember(fitbit_daily_BL.src_subject_id,sub_id_BL7(i)))>6
    ind_sid=ismember(fitbit_daily_BL.src_subject_id,sub_id_BL7(i));
    ind_sid(ind_sid==1)= hours(fitbit_daily_BL{ind_sid,'fit_ss_wakeup_minutes'} - fitbit_daily_BL{ind_sid,'fit_ss_first_sleep_minutes'})<24;
    
%     
%    ind_sid(ind_sid==1)= hours(fitbit_daily_test{ind_sid,'fit_ss_wakeup_minutes'} - fitbit_daily_test{ind_sid,'fit_ss_first_sleep_minutes'})<24&...
%         hours(fitbit_daily_test{ind_sid,'fit_ss_wakeup_minutes'} - fitbit_daily_test{ind_sid,'fit_ss_first_sleep_minutes'})>2;
    weekend_ind=ind_sid&fitbit_daily_BL.fit_ss_weekend_ind==1;
    weekday_ind=ind_sid&fitbit_daily_BL.fit_ss_weekend_ind==0;
    
    
   if sum(ind_sid)>6&sum(weekend_ind)>1&sum(weekday_ind)>1
       sub_id_7_BL=[sub_id_7_BL;sub_id_BL7(i)];
   end
       
       
     
   end
end
% sub_id_both=sub_id_7_BL(ismember(sub_id_7_BL,sub_id_7));%49

fitbit_daily_sum_BL=table;
fitbit_daily_sum_BL.sub_id=sub_id_7_BL;
fitbit_daily_sum_BL.Properties.RowNames=sub_id_7_BL;
fitbit_daily_sum_BL.weekday_sleep=nan(length(sub_id_7_BL),1);
fitbit_daily_sum_BL.weekend_sleep=nan(length(sub_id_7_BL),1);
for i=1:length(sub_id_7_BL)
    weekday_ind=ismember(fitbit_daily_BL.src_subject_id,sub_id_7_BL(i))&fitbit_daily_BL.fit_ss_weekend_ind==0;
    weekend_ind=ismember(fitbit_daily_BL.src_subject_id,sub_id_7_BL(i))&fitbit_daily_BL.fit_ss_weekend_ind==1;

    weekday_ind(weekday_ind==1)=(hours(fitbit_daily_BL{weekday_ind,'fit_ss_wakeup_minutes'} - fitbit_daily_BL{weekday_ind,'fit_ss_first_sleep_minutes'})<24);%make sure sleep duration less than 24 hours and longer than 2 hours. 
    weekend_ind(weekend_ind==1)=(hours(fitbit_daily_BL{weekend_ind,'fit_ss_wakeup_minutes'} - fitbit_daily_BL{weekend_ind,'fit_ss_first_sleep_minutes'})<24);%make sure sleep duration less than 24 hours and longer than 2 hours. 

    
    
    fitbit_daily_sum_BL{sub_id_7_BL(i),'weekday_sleep'}=mean(fitbit_daily_BL{weekday_ind,"fit_ss_sleepperiod_minutes"});
    fitbit_daily_sum_BL{sub_id_7_BL(i),'weekday_sleep_num'}=length(fitbit_daily_BL{weekday_ind,"fit_ss_sleepperiod_minutes"});
    fitbit_daily_sum_BL{sub_id_7_BL(i),'weekend_sleep'}=mean(fitbit_daily_BL{weekend_ind,"fit_ss_sleepperiod_minutes"});
    fitbit_daily_sum_BL{sub_id_7_BL(i),'weekend_sleep_num'}=length(fitbit_daily_BL{weekend_ind,"fit_ss_sleepperiod_minutes"});
    fitbit_daily_sum_BL{sub_id_7_BL(i),'avg_sleep'}=mean(fitbit_daily_BL{ismember(fitbit_daily_BL.src_subject_id,sub_id_7_BL(i)),"fit_ss_sleepperiod_minutes"});
    fitbit_daily_sum_BL{sub_id_7_BL(i),'avg_sleep_num'}=length(fitbit_daily_BL{ismember(fitbit_daily_BL.src_subject_id,sub_id_7_BL(i)),"fit_ss_sleepperiod_minutes"});
    fitbit_daily_sum_BL{sub_id_7_BL(i),'sd_sleep'}=std(fitbit_daily_BL{ismember(fitbit_daily_BL.src_subject_id,sub_id_7_BL(i)),"fit_ss_sleepperiod_minutes"});
    fitbit_daily_sum_BL{sub_id_7_BL(i),'weekday_wake'}=mean(fitbit_daily_BL{weekday_ind,"fit_ss_wake_minutes"});
    fitbit_daily_sum_BL{sub_id_7_BL(i),'weekday_wake_count'}=mean(fitbit_daily_BL{weekday_ind,"fit_ss_wake_count"});
    fitbit_daily_sum_BL{sub_id_7_BL(i),'weekend_wake'}=mean(fitbit_daily_BL{weekend_ind,"fit_ss_wake_minutes"});
    fitbit_daily_sum_BL{sub_id_7_BL(i),'weekend_wake_count'}=mean(fitbit_daily_BL{weekend_ind,"fit_ss_wake_count"});
    fitbit_daily_sum_BL{sub_id_7_BL(i),'avg_wake'}=mean(fitbit_daily_BL{ismember(fitbit_daily_BL.src_subject_id,sub_id_7_BL(i)),"fit_ss_wake_minutes"});
    fitbit_daily_sum_BL{sub_id_7_BL(i),'avg_wake_num'}=length(fitbit_daily_BL{ismember(fitbit_daily_BL.src_subject_id,sub_id_7_BL(i)),"fit_ss_wake_count"});
 
    fitbit_daily_sum_BL{sub_id_7_BL(i),'weekday_midpoint'}=mean(timeofday(fitbit_daily_BL{weekday_ind,'fit_ss_first_sleep_minutes'}+...
        (fitbit_daily_BL{weekday_ind,'fit_ss_wakeup_minutes'} - fitbit_daily_BL{weekday_ind,'fit_ss_first_sleep_minutes'})/2));
    fitbit_daily_sum_BL{sub_id_7_BL(i),'weekend_midpoint'}=mean(timeofday(fitbit_daily_BL{weekend_ind,'fit_ss_first_sleep_minutes'}+...
        (fitbit_daily_BL{weekend_ind,'fit_ss_wakeup_minutes'} - fitbit_daily_BL{weekend_ind,'fit_ss_first_sleep_minutes'})/2));

    
    
    
    
end

fitbit_daily_sum_BL.avg_timeinbed=fitbit_daily_sum_BL.avg_sleep+fitbit_daily_sum_BL.avg_wake;
fitbit_daily_sum_BL.sleep_diff=fitbit_daily_sum_BL.weekend_sleep-fitbit_daily_sum_BL.weekday_sleep;
fitbit_daily_sum_BL.midpoint_diff=hours(fitbit_daily_sum_BL.weekend_midpoint-fitbit_daily_sum_BL.weekday_midpoint);
% 
% % 
% % 
% 
% 
% 
% 


fitbit_daily_test=fitbit_sleep_FL2;% two-year follow-up data 

 fitbit_daily_test=fitbit_daily_test(fitbit_daily_test.fit_ss_wkno<4,:);
sub_id=unique(fitbit_daily_test.src_subject_id);
sub_id_7=[];
for i=1:length(sub_id)
    
   ind_sid=ismember(fitbit_daily_test.src_subject_id,sub_id(i));
   ind_sid(ind_sid==1)= hours(fitbit_daily_test{ind_sid,'fit_ss_wakeup_minutes'} - fitbit_daily_test{ind_sid,'fit_ss_first_sleep_minutes'})<24;
    

    weekend_ind=ind_sid&fitbit_daily_test.fit_ss_weekend_ind==1;
    weekday_ind=ind_sid&fitbit_daily_test.fit_ss_weekend_ind==0;
    
    
   if sum(ind_sid)>6&sum(weekend_ind)>2&sum(weekday_ind)>2
       sub_id_7=[sub_id_7;sub_id(i)];
   end
end

fitbit_daily_sum=table;
fitbit_daily_sum.sub_id=sub_id_7;
fitbit_daily_sum.Properties.RowNames=sub_id_7;
fitbit_daily_sum.weekday_sleep=nan(length(sub_id_7),1);
fitbit_daily_sum.weekend_sleep=nan(length(sub_id_7),1);
for i=1:length(sub_id_7)
    weekday_ind=ismember(fitbit_daily_test.src_subject_id,sub_id_7(i))&fitbit_daily_test.fit_ss_weekend_ind==0;
    weekend_ind=ismember(fitbit_daily_test.src_subject_id,sub_id_7(i))&fitbit_daily_test.fit_ss_weekend_ind==1;
    weekday_ind(weekday_ind==1)=(hours(fitbit_daily_test{weekday_ind,'fit_ss_wakeup_minutes'} - fitbit_daily_test{weekday_ind,'fit_ss_first_sleep_minutes'})<15&...
        hours(fitbit_daily_test{weekday_ind,'fit_ss_wakeup_minutes'} - fitbit_daily_test{weekday_ind,'fit_ss_first_sleep_minutes'})>3);%make sure sleep duration less than 15 hours and longer than 3 hours. 
    weekend_ind(weekend_ind==1)=(hours(fitbit_daily_test{weekend_ind,'fit_ss_wakeup_minutes'} - fitbit_daily_test{weekend_ind,'fit_ss_first_sleep_minutes'})<15&...
        hours(fitbit_daily_test{weekend_ind,'fit_ss_wakeup_minutes'} - fitbit_daily_test{weekend_ind,'fit_ss_first_sleep_minutes'})>3);%make sure sleep duration less than 15 hours and longer than 3 hours. 

    fitbit_daily_sum{sub_id_7(i),'weekday_sleep'}=mean(fitbit_daily_test{weekday_ind,"fit_ss_sleepperiod_minutes"});
    fitbit_daily_sum{sub_id_7(i),'weekday_sleep_num'}=length(fitbit_daily_test{weekday_ind,"fit_ss_sleepperiod_minutes"});
    fitbit_daily_sum{sub_id_7(i),'weekend_sleep'}=mean(fitbit_daily_test{weekend_ind,"fit_ss_sleepperiod_minutes"});
    fitbit_daily_sum{sub_id_7(i),'weekend_sleep_num'}=length(fitbit_daily_test{weekend_ind,"fit_ss_sleepperiod_minutes"});
    fitbit_daily_sum{sub_id_7(i),'avg_sleep'}=mean(fitbit_daily_test{ismember(fitbit_daily_test.src_subject_id,sub_id_7(i)),"fit_ss_sleepperiod_minutes"});
    fitbit_daily_sum{sub_id_7(i),'avg_sleep_num'}=length(fitbit_daily_test{ismember(fitbit_daily_test.src_subject_id,sub_id_7(i)),"fit_ss_sleepperiod_minutes"});
    fitbit_daily_sum{sub_id_7(i),'sd_sleep'}=std(fitbit_daily_test{ismember(fitbit_daily_test.src_subject_id,sub_id_7(i)),"fit_ss_sleepperiod_minutes"});
    fitbit_daily_sum{sub_id_7(i),'weekday_wake'}=mean(fitbit_daily_test{weekday_ind,"fit_ss_wake_minutes"});
    fitbit_daily_sum{sub_id_7(i),'weekday_wake_count'}=mean(fitbit_daily_test{weekday_ind,"fit_ss_wake_count"});
    fitbit_daily_sum{sub_id_7(i),'weekend_wake'}=mean(fitbit_daily_test{weekend_ind,"fit_ss_wake_minutes"});
    fitbit_daily_sum{sub_id_7(i),'weekend_wake_count'}=mean(fitbit_daily_test{weekend_ind,"fit_ss_wake_count"});
    fitbit_daily_sum{sub_id_7(i),'avg_wake'}=mean(fitbit_daily_test{ismember(fitbit_daily_test.src_subject_id,sub_id_7(i)),"fit_ss_wake_minutes"});
    fitbit_daily_sum{sub_id_7(i),'avg_wake_num'}=length(fitbit_daily_test{ismember(fitbit_daily_test.src_subject_id,sub_id_7(i)),"fit_ss_wake_count"});
    inbedtime=[];
    inbedtime=fitbit_daily_test{ismember(fitbit_daily_test.src_subject_id,sub_id_7(i)),"fit_ss_first_inbed_minutes"};
    inbedtime=inbedtime(hour(inbedtime)>6&hour(inbedtime)<24);
    fitbit_daily_sum{sub_id_7(i),'sd_inbedtime'}=std(hour(inbedtime)*60+minute(inbedtime));
    fitbit_daily_sum{sub_id_7(i),'weekday_midpoint'}=mean(timeofday(fitbit_daily_test{weekday_ind,'fit_ss_first_sleep_minutes'}+...
        (fitbit_daily_test{weekday_ind,'fit_ss_wakeup_minutes'} - fitbit_daily_test{weekday_ind,'fit_ss_first_sleep_minutes'})/2));
    fitbit_daily_sum{sub_id_7(i),'weekend_midpoint'}=mean(timeofday(fitbit_daily_test{weekend_ind,'fit_ss_first_sleep_minutes'}+...
        (fitbit_daily_test{weekend_ind,'fit_ss_wakeup_minutes'} - fitbit_daily_test{weekend_ind,'fit_ss_first_sleep_minutes'})/2));
    fitbit_daily_sum{sub_id_7(i),'sjlsc'}=hours(mean(timeofday(fitbit_daily_test{weekend_ind,'fit_ss_first_sleep_minutes'}))...
        -mean(timeofday(fitbit_daily_test{weekday_ind,'fit_ss_first_sleep_minutes'})));

%     fitbit_daily_sum{sub_id_7(i),'midpoint_diff'}=fitbit_daily_sum{sub_id_7(i),'weekend_midpoint'}-fitbit_daily_sum{sub_id_7(i),'weekday_midpoint'};
    
    
end

fitbit_daily_sum.avg_timeinbed=fitbit_daily_sum.avg_sleep+fitbit_daily_sum.avg_wake;
fitbit_daily_sum.sleep_diff=fitbit_daily_sum.weekend_sleep-fitbit_daily_sum.weekday_sleep;
fitbit_daily_sum.sjl=hours(fitbit_daily_sum.weekend_midpoint-fitbit_daily_sum.weekday_midpoint);
fitbit_daily_sum.sjl=abs(fitbit_daily_sum.sjl);
fitbit_daily_sum.sjl(fitbit_daily_sum.sjl>12)=fitbit_daily_sum.sjl(fitbit_daily_sum.sjl>12)-12;
fitbit_daily_sum.sjlsc=abs(fitbit_daily_sum.sjlsc);
fitbit_daily_sum.sjlsc(fitbit_daily_sum.sjlsc>12)=fitbit_daily_sum.sjlsc(fitbit_daily_sum.sjlsc>12)-12;






parent_dem_table=readtable(fullfile(beh_dir,'abcd_lpds01.txt'));% demographic for parents
parent_dem_table_BL=parent_dem_table(strcmp(parent_dem_table{:,'eventname'},'baseline_year_1_arm_1'),:);
parent_dem_table_BL.Properties.RowNames=parent_dem_table_BL.src_subject_id;
parent_dem_table_FL2=parent_dem_table(strcmp(parent_dem_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
parent_dem_table_FL2.Properties.RowNames=parent_dem_table_FL2.src_subject_id;

% network_conn_table=([network_table_BL(sub_id_BL,23:191) sub_network_table_BL(sub_id_BL,23:269)]);
% network_conn_table_FL2=([network_table_FL2(sub_id_FL2,23:191) sub_network_table_FL2(sub_id_FL2,23:269)]);

T1_table=readtable(fullfile(beh_dir,'abcd_mrisdp10201.txt'));% brain structural varaibles
T1_table_BL=T1_table(strcmp(T1_table{:,'eventname'},'baseline_year_1_arm_1'),:);
T1_table_BL.Properties.RowNames=T1_table_BL.src_subject_id;
T1_table_FL2=T1_table(strcmp(T1_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
T1_table_FL2.Properties.RowNames=T1_table_FL2.src_subject_id;


T1_sub=readtable(fullfile(beh_dir,'abcd_smrip201.txt'));%330:375 brain structural varaibles (subcortical)
T1_sub_BL=T1_sub(strcmp(T1_sub{:,'eventname'},'baseline_year_1_arm_1'),:);
T1_sub_BL.Properties.RowNames=T1_sub_BL.src_subject_id;
T1_sub_FL2=T1_sub(strcmp(T1_sub{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
T1_sub_FL2.Properties.RowNames=T1_sub_FL2.src_subject_id;


T1_qc_table=readtable(fullfile(beh_dir,'mriqcrp10301.txt'));
T1_qc_table_BL=T1_qc_table(strcmp(T1_qc_table{:,'eventname'},'baseline_year_1_arm_1'),:);
T1_qc_table_BL.Properties.RowNames=T1_qc_table_BL.src_subject_id;
T1_qc_table_FL2=T1_qc_table(strcmp(T1_qc_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
T1_qc_table_FL2.Properties.RowNames=T1_qc_table_FL2.src_subject_id;





ps_table=readtable(fullfile(beh_dir,'abcd_ps01.txt'));
ps_table_BL=ps_table(strcmp(ps_table{:,'eventname'},'baseline_year_1_arm_1'),:);
ps_table_BL.Properties.RowNames=ps_table_BL.src_subject_id;
ps_table_FL1=ps_table(strcmp(ps_table{:,'eventname'},'1_year_follow_up_y_arm_1'),:);
ps_table_FL1.Properties.RowNames=ps_table_FL1.src_subject_id;
ps_table_FL2=ps_table(strcmp(ps_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
ps_table_FL2.Properties.RowNames=ps_table_FL2.src_subject_id;


mcqc_table=readtable(fullfile(beh_dir,'abcd_mcqc01.txt'));% munich choronotype 
mcqc_table_FL3=mcqc_table(strcmp(mcqc_table{:,'eventname'},'3_year_follow_up_y_arm_1'),:);
mcqc_table_FL3.Properties.RowNames=mcqc_table_FL3.src_subject_id;
mcqc_table_FL2=mcqc_table(strcmp(mcqc_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
mcqc_table_FL2.Properties.RowNames=mcqc_table_FL2.src_subject_id;




sub_id_FL1=sub_id_BL(ismember(sub_id_BL,sleep_parent_SDS_FL1.Properties.RowNames));

sub_id_FL2=sub_id_BL(ismember(sub_id_BL,mcqc_table_FL2.Properties.RowNames));
sub_id_FL3=sub_id_FL2(ismember(sub_id_FL2,mcqc_table_FL3.Properties.RowNames));





pubertal_f_FL2=sleep_parent_SDSC_FL2{sub_id_FL2,'pds_p_ss_female_category_2'};% change to categor_2
pubertal_m_FL2=sleep_parent_SDSC_FL2{sub_id_FL2,'pds_p_ss_male_category_2'};%change to categor_2
pubertal_FL2=sum([pubertal_f_FL2,pubertal_m_FL2],2,'omitnan');
pubertal_FL2(pubertal_FL2==0)=nan;

sjlsc1=mcqc_table_FL2{sub_id_FL2,'mctq_sow_calc'}+0.5*mcqc_table_FL2{sub_id_FL2,'mctq_sdweek_calc'};
sjlsc1(sjlsc1>=24)=sjlsc1(sjlsc1>=24)-24;
sjlsc2=mcqc_table_FL2{sub_id_FL2,'mctq_sof_calc'}+0.5*mcqc_table_FL2{sub_id_FL2,'mctq_sdweek_calc'};
sjlsc2(sjlsc2>=24)=sjlsc2(sjlsc2>=24)-24;
sjlsc=abs(sjlsc2-sjlsc1);
sjlsc(sjlsc>12)=24-abs(sjlsc(sjlsc>12));





psm_table=[NIH_TB_BL(sub_id_BL,'interview_age'),...
        dem_table(sub_id_BL,'race4'),dem_table(sub_id_BL,'sex'),dem_table(sub_id_BL,'edu')...
        ,dem_table(sub_id_BL,'abcd_site'),dem_table(sub_id_BL,'rel_family_id'),...
       ...
        BMI_BL(sub_id_BL,'anthro_bmi_calc')];% data for propensity score matching

psm_table.Properties.VariableNames={'age','race','sex','prt_edu','sites','family','BMI'};
psm_table.age_sex=psm_table.age.*psm_table.sex;
% psm_table.puberty=pubertal;
psm_table.HI=House_in;
TS=double(total_sleep_duration);
psm_table.TS=TS;


sub_Fl2_bad_SD=sub_id_FL2(mcqc_table_FL2{sub_id_FL2,'mctq_sdw_calc'}>=15|mcqc_table_FL2{sub_id_FL2,'mctq_sdf_calc'}>=15|mcqc_table_FL2{sub_id_FL2,'mctq_sdw_calc'}<=3|mcqc_table_FL2{sub_id_FL2,'mctq_sdf_calc'}<=3);
% exclude children sleep less than 3 hour and longer than 15 hours


 [Handedness_BL,~,~,~,~]=abcd_load_var2table(fullfile(beh_dir,'abcd_ehis01.txt'));% handedness data



psm_table_FL2=[NIH_TB_FL2(sub_id_FL2,'interview_age'),...
        dem_table(sub_id_FL2,'race4'),dem_table(sub_id_FL2,'sex'),parent_dem_table_FL2(sub_id_FL2,'demo_prnt_ed_v2_2yr_l')...
        ,dem_table(sub_id_FL2,'abcd_site'),dem_table(sub_id_FL2,'rel_family_id'),...
       ...
        BMI_FL2(sub_id_FL2,'anthro_bmi_calc')];% data for propensity score matching two-year follow-up

psm_table_FL2.Properties.VariableNames={'age','race','sex','prt_edu','sites','family','BMI_BL'};
psm_table_FL2.age_sex=psm_table_FL2.age.*psm_table_FL2.sex;
psm_table_FL2.puberty=pubertal_FL2;% using baseline pubetal status as proximy
psm_table_FL2.HI=parent_dem_table_FL2{sub_id_FL2,'demo_comb_income_v2_l'};
psm_table_FL2.TS=mcqc_table_FL2{sub_id_FL2,'mctq_sdweek_calc'};
psm_table_FL2.SJLsc=sjlsc;
psm_table_FL2.SJL=mcqc_table_FL2{sub_id_FL2,'mctq_sjl_calc'};
psm_table_FL2.SDis=sleep_parent_SDSC_FL2{sub_id_FL2,'sds_p_ss_total'};

psm_table_FL2.HI(psm_table_FL2.HI==777|psm_table_FL2.HI==999)=nan;%999, Don't know ; 777, Refuse to answer 
psm_table_FL2.prt_edu(psm_table_FL2.prt_edu==777|psm_table_FL2.prt_edu==999)=nan;%999, Don't know ; 777, Refuse to answer 


%% runing psm in R
psm_table_FL2_nonan=psm_table_FL2(~isnan(sum(psm_table_FL2{:,1:14},2)),:);% exclued subs with missing values; n=8966 
%6565
psm_table_FL2_nonan(ismember(psm_table_FL2_nonan.Properties.RowNames,sub_Fl2_bad_SD),:)=[];%n=8675
%6335

sum(isnan(psm_table_FL2{:,1:14}))
psm_table_FL2_nonan.num=(1:height(psm_table_FL2_nonan))';
writetable(psm_table_FL2_nonan,'psm_table_FL2_nonan_40_nooutlier_r_bmi.csv')

psm_table_FL2_nonan.SJLsc_b=double(psm_table_FL2_nonan.SJLsc>1.000000001);

mdl = fitglm(psm_table_FL2_nonan,'SJLsc_b~age+sex+race+prt_edu+sites+HI+BMI_BL+puberty+age_sex+TS+SDis','Distribution','binomial');%adjR2 0.0903
%% follow up 2 behavior


psm_FL2_matched=readtable('/home/yangf7/Documents/Nils/ABCD/sjl/psm_table_FL2_nn_matched40_nooutlier_r_bmi.csv');
psm_FL2_matched.Properties.RowNames=psm_table_FL2_nonan.Properties.RowNames(psm_FL2_matched.num);
psm_FL2_matched.SJLsc_b=double(psm_FL2_matched.SJLsc>1.000000001);
psm_FL2_matched_table=sortrows(psm_FL2_matched,[17 20]);%group_diff and subclass
psm_sjl0=psm_FL2_matched_table.Properties.RowNames(1:height(psm_FL2_matched_table)/2);%sjlsc <=1
psm_sjl1=psm_FL2_matched_table.Properties.RowNames(height(psm_FL2_matched_table)/2+1:end);%sjlsc >1
sjl_psm=[psm_sjl0;psm_sjl1];

ind_fitbit=ismember(psm_sjl0,sub_id_7)&ismember(psm_sjl1,sub_id_7);%matched subject that both have fitbit data
sjl_psm_fitbit=[psm_sjl0(ind_fitbit);psm_sjl1(ind_fitbit)];
psm_sjl0_fitbit=psm_sjl0(ind_fitbit);
psm_sjl1_fitbit=psm_sjl1(ind_fitbit);


% [t,p]=corr(fitbit_daily_sum{sub_id_7,'sjlsc'},psm_table_FL2{sub_id_7,'SJLsc'},'rows','pairwise')

ind_FL3=ismember(psm_sjl0,sub_id_FL3)&ismember(psm_sjl1,sub_id_FL3);%matched subject that both have FL3 data
sjl_psm_FL3=[psm_sjl0(ind_FL3);psm_sjl1(ind_FL3)];
psm_sjl0_FL3=psm_sjl0(ind_FL3);
psm_sjl1_FL3=psm_sjl1(ind_FL3);



cbcl_raw_FL2=cbcl_sum_FL2(sjl_psm,10:4:86);
cbcl_raw_FL2=cbcl_raw_FL2(:,sort(cbcl_raw_FL2.Properties.VariableNames));

cbcl_raw_FL2=cbcl_raw_FL2(:,[1:18 20 19]);% put total pro to the last one


cognition_row_FL2=NIH_TB_FL2(sjl_psm,[12:5:42 46:4:54]);
mhy_row_FL2=mental_health_FL2(sjl_psm,[28 31 37 41 44 47 50 53 56:3:65]);
% screentime_row_FL2=screen_time_FL2(sleep_psc_FL2,[10:21]);
% 
% 


sjl_grade=grade_table_FL2{sjl_psm,'sag_grades_last_yr'};
sjl_grade(sjl_grade==-1|sjl_grade==777)=nan;
psm_FL2_matched_table.grades=13-sjl_grade;% reverse coding, so higher value means better grade
% 
cbcl_raw_FL2=[cbcl_raw_FL2,cognition_row_FL2,mhy_row_FL2,psm_FL2_matched_table(sjl_psm,'grades')];



cohend_FL2=nan(1,width(cbcl_raw_FL2));  
p_value=nan(1,width(cbcl_raw_FL2));
for i=1:width(cbcl_raw_FL2)
    if mean(isnan(cbcl_raw_FL2{psm_sjl0,i}))>0.9
          cohend_FL2(i)=0;
        p_value(i)=1;
    else
        
%         effect=meanEffectSize(cbcl_raw_FL2{fitbit_rs,i},cbcl_raw_FL2{fitbit_is,i},Effect="cohen")
    cohend_FL2(i)=computeCohen_d(cbcl_raw_FL2{psm_sjl0,i},cbcl_raw_FL2{psm_sjl1,i},'independent');
     [h,p_value(i)]=ttest2(cbcl_raw_FL2{psm_sjl0,i},cbcl_raw_FL2{psm_sjl1,i});
    
      
    end
end
% p_value=-log10(p_value);

[~,~,~,beh_fdr_p_FL2]=fdr_bh(p_value);
beh_fdr_p_FL2(beh_fdr_p_FL2<0.0001)=0.0001;


cohend_combined=cohend_FL2;
figure;
b=bar(cohend_combined,'FaceColor','flat');
for i=21:21+width(cognition_row_FL2)-1
b.CData(i,:) = [0 1 0];
end

for i=21+width(cognition_row_FL2):21+width(cognition_row_FL2)-1+width(mhy_row_FL2)
b.CData(i,:) = [1 1 0];
end
for i=43
b.CData(i,:) = [1 0 0];
end
hold on
yline([-0.15 0.15],'--')

% xticks(1:width(cbcl_raw_BL));
set(gca, 'FontName', 'Helvetica','FontSize',12)
xticks(1:width(cohend_combined));
xticklabels(strrep(cbcl_raw_FL2.Properties.VariableNames(1:43),'_','\_'))
ylabel({'Cohen''s d (Low - High SJLsc)'},'FontName', 'Helvetica','FontSize',24)


title({'Behavioral variables at FL2'},'FontSize',24)
set(gcf,'units','inches','position',[0,0,15,6])
set(gca,'TickLength',[0 0])
xtickangle(45)
box off

%% mediation analysis
mediator_FC=[network_conn_table_FL2(sjl_psm_FL2img,'rsfmri_cor_ngd_cerc_scs_hprh'),network_conn_table_FL2(sjl_psm_FL2img,'rsfmri_cor_ngd_cerc_scs_agrh'),...
    network_conn_table_FL2(sjl_psm_FL2img,'rsfmri_cor_ngd_copa_scs_vtdclh'),network_conn_table_FL2(sjl_psm_FL2img,'rsfmri_cor_ngd_smh_scs_cderh')...
    network_conn_table_FL2(sjl_psm_FL2img,'rsfmri_cor_ngd_smh_scs_ptrh')];
cov_med=[psm_FL2_matched_table(sjl_psm_FL2img,[2:6 8:12 15 ]),network_table_FL2(sjl_psm_FL2img,'rsfmri_c_ngd_ntpoints')...
        ,network_table_FL2(sjl_psm_FL2img,'rsfmri_c_ngd_meanmotion')];
dep_med=[cbcl_raw_FL2(sjl_psm_FL2img,abs(cohend_FL2)>0.15)];


p_med=nan(size(mediator_FC,2),size(dep_med,2));
beta_med=nan(size(mediator_FC,2),size(dep_med,2));
for i=1:size(mediator_FC,2)
    for j=1:size(dep_med,2)
    [paths, stats_med]=mediation(psm_FL2_matched_table{sjl_psm_FL2img,'SJLsc_b'},dep_med{sjl_psm_FL2img,j},...
    mediator_FC{sjl_psm_FL2img,i}...
    ,'cov', cov_med{sjl_psm_FL2img,:},'verbose','boot','bootsamples',10000);
    p_med(i,j)=stats_med.p(5);
    beta_med(i,j)=stats_med.mean(5);
    end
end

p_med_fitbit=nan(size(mediator_FC,2),size(dep_med,2));
beta_med_fitbit=nan(size(mediator_FC,2),size(dep_med,2));
for i=1:size(mediator_FC,2)
    for j=1:size(dep_med,2)
    [paths, stats_med]=mediation(fitbit_daily_sum{sub_id_7(ismember(sub_id_7,sjl_psm_FL2img)),'sjlsc'},dep_med{sub_id_7(ismember(sub_id_7,sjl_psm_FL2img)),j},...
    mediator_FC{sub_id_7(ismember(sub_id_7,sjl_psm_FL2img)),i}...
    ,'cov', cov_med{sub_id_7(ismember(sub_id_7,sjl_psm_FL2img)),:},'verbose','boot','bootsamples',10000);
    p_med_fitbit(i,j)=stats_med.p(5);
    beta_med_fitbit(i,j)=stats_med.mean(5);
    end
end

% mediation figure
p_med(p_med>=0.05)=1;
fig1=figure;
imagesc(-log10(p_med))
title('Mediation Effect (a*b)','FontSize',18);
set(gcf,'Position',[100 100 800 800])
set(gca, 'XTick', 1:size(dep_med,2)); % center x-axis ticks on bins
set(gca, 'YTick', 1:size(mediator_FC,2)); % center y-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel',{'picture vocabulary','reading','CI','PPS number','PU','grades'}); % set x-axis labels
set(gca, 'YTickLabel', {'cerc-hprh','cerc-agrh','copa-vtdclh','smh-cderh'...
    ,'smh-ptrh'},'FontSize',18); % set y-axis labels
set(gcf,'units','inches','position',[0,0,9,9])
c = colorbar;
c.Limits= [0 5];
c.TickLength = [0];
c.Label.String = '-log10(p) of mediation effect';
c.FontSize=16;
colormap(bluewhitered)
% figname='/home/yangf7/Documents/Nils/ABCD/sjl/mediation.jpg';
% exportfig(fig1,figname,'Format','png','Color','cmyk',...
%     'Resolution',300,'Renderer','opengl');

% mediation figure
figure;
imagesc(-log10(p_med_fitbit))
title('Follow-up 2 Fitbit','FontSize',18);
set(gcf,'Position',[100 100 800 800])
set(gca, 'XTick', 1:size(dep_med,2)); % center x-axis ticks on bins
set(gca, 'YTick', 1:size(mediator_FC,2)); % center y-axis ticks on bins
set(gca, 'TickLength',[0 0])
set(gca, 'XTickLabel',{'picture vocabulary','picture','reading','CI','PU','grades'}); % set x-axis labels
set(gca, 'YTickLabel', {'cerc-hprh','cerc-agrh','copa-vtdclh','smh-cderh'...
    ,'smh-ptrh'}); % set y-axis labels
set(gcf,'units','inches','position',[0,0,9,9])
c = colorbar;
c.Limits= [0 5];
c.TickLength = [0];
c.Label.String = '-log10(p) of mediation effect';
c.FontSize=16;
colormap(bluewhitered)









ind_FL3img=ismember(psm_sjl0_FL2img,sjl_psm_FL3)&ismember(psm_sjl1_FL2img,sjl_psm_FL3);%matched img subject that both have FL3 data
sjl_psm_FL3img=[psm_sjl0_FL2img(ind_FL3img);psm_sjl1_FL2img(ind_FL3img)];
% psm_sjl0_FL3=psm_sjl0(ind_FL3);
% psm_sjl1_FL3=psm_sjl1(ind_FL3);

%   length(sjl_psm_FL3(ismember(sjl_psm_FL3,sjl_psm_FL2img)))
  [paths, stats_med]=mediation(psm_FL2_matched_table{sjl_psm_FL3img,'SJLsc_b'},cbcl_raw_FL3{sjl_psm_FL3img,'grades'},...
    network_conn_table_FL2{sjl_psm_FL3img,'rsfmri_cor_ngd_cerc_scs_hprh'}...
    ,'cov', [psm_FL2_matched_table{sjl_psm_FL3img,[2:6 8:12 15 ]},cbcl_raw_FL2{sjl_psm_FL3img,'grades'}],'verbose','boot','bootsamples',10000);

p_med_FL3=nan(size(mediator_FC,2),1);
beta_med_FL3=nan(size(mediator_FC,2),1);
for i=1:size(mediator_FC,2)
    
    [paths, stats_med]=mediation(psm_FL2_matched_table{sjl_psm_FL3img,'SJLsc_b'},cbcl_raw_FL3{sjl_psm_FL3img,'grades'},...
    mediator_FC{sjl_psm_FL3img,i}...
    ,'cov', [cov_med{sjl_psm_FL3img,:},cbcl_raw_FL2{sjl_psm_FL3img,'grades'}.^4],'verbose','boot','bootsamples',10000);

    p_med_FL3(i)=stats_med.p(5);
    beta_med_FL3(i)=stats_med.mean(5);
    
end












%% FL3 behavior
cbcl_raw_FL3=cbcl_sum_FL3(sjl_psm_FL3,10:4:86);
cbcl_raw_FL3=cbcl_raw_FL3(:,sort(cbcl_raw_FL3.Properties.VariableNames));

cbcl_raw_FL3=cbcl_raw_FL3(:,[1:18 20 19]);% put total pro to the last one





sjl_grade_FL3=grade_table_FL3{sjl_psm_FL3,'sag_grades_last_yr'};
sjl_grade_FL3(sjl_grade_FL3==-1|sjl_grade_FL3==777)=nan;


TS_FL3=mcqc_table_FL3{sjl_psm_FL3,'mctq_sdweek_calc'};

SDis_FL3=sleep_parent_SDSC_FL3{sjl_psm_FL3,'sds_p_ss_total'};% psm_FL2_matched_table.grade=sjl_grade;


sjlsc1_FL3=mcqc_table_FL3{sjl_psm_FL3,'mctq_sow_calc'}+0.5*mcqc_table_FL3{sjl_psm_FL3,'mctq_sdweek_calc'};
sjlsc1_FL3(sjlsc1_FL3>=24)=sjlsc1_FL3(sjlsc1_FL3>=24)-24;
sjlsc2_FL3=mcqc_table_FL3{sjl_psm_FL3,'mctq_sof_calc'}+0.5*mcqc_table_FL3{sjl_psm_FL3,'mctq_sdweek_calc'};
sjlsc2_FL3(sjlsc2_FL3>=24)=sjlsc2_FL3(sjlsc2_FL3>=24)-24;
sjlsc_FL3=abs(sjlsc2_FL3-sjlsc1_FL3);
sjlsc_FL3(sjlsc_FL3>12)=24-abs(sjlsc_FL3(sjlsc_FL3>12));


cbcl_raw_FL3=[cbcl_raw_FL3];
cbcl_raw_FL3.grades=13-sjl_grade_FL3;%reverse coding, higher value means better grade
cbcl_raw_FL3.TS=TS_FL3;
cbcl_raw_FL3.SDis=SDis_FL3;
cbcl_raw_FL3.SJLsc=sjlsc_FL3;

cohend_FL3=nan(1,width(cbcl_raw_FL3));  
p_value=nan(1,width(cbcl_raw_FL3));
for i=1:width(cbcl_raw_FL3)
    if mean(isnan(cbcl_raw_FL3{psm_sjl0_FL3,i}))>0.9
          cohend_FL3(i)=0;
        p_value(i)=1;
    else
        
%         effect=meanEffectSize(cbcl_raw_FL2{fitbit_rs,i},cbcl_raw_FL2{fitbit_is,i},Effect="cohen")
    cohend_FL3(i)=computeCohen_d(cbcl_raw_FL3{psm_sjl0_FL3,i},cbcl_raw_FL3{psm_sjl1_FL3,i},'independent');
     [h,p_value(i)]=ttest2(cbcl_raw_FL3{psm_sjl0_FL3,i},cbcl_raw_FL3{psm_sjl1_FL3,i});
    
      
    end
end
% p_value=-log10(p_value);

[~,~,~,beh_fdr_p_FL3]=fdr_bh(p_value(1:21));
beh_fdr_p_FL3(beh_fdr_p_FL3<0.0001)=0.0001;


figure;
b=bar(cohend_FL3(1:21),'FaceColor','flat');
for i=21
b.CData(i,:) = [1 0 0];
end

set(gca, 'FontName', 'Helvetica','FontSize',12)
xticks(1:21);
xticklabels(strrep(cbcl_raw_FL3.Properties.VariableNames,'_','\_'))
ylabel({'Cohen''s d (Low - High SJLsc)'},'FontName', 'Helvetica','FontSize',24)


title({'Behavioral variables at FL3'},'FontSize',24)
set(gcf,'units','inches','position',[0,0,10,6])
set(gca,'TickLength',[0 0])
xtickangle(45)
box off


%% imaging FL2
sub_id_FL2img=sub_id_FL2(ismember(sub_id_FL2,network_table_FL2.Properties.RowNames));
ind_FL2=ismember(psm_sjl0,sub_id_FL2img)&ismember(psm_sjl1,sub_id_FL2img);%matched subject that both have FL3 data
sjl_psm_FL2img=[psm_sjl0(ind_FL2);psm_sjl1(ind_FL2)];
psm_sjl0_FL2img=psm_sjl0(ind_FL2);
psm_sjl1_FL2img=psm_sjl1(ind_FL2);

network_conn_table_FL2=([network_table_FL2(sjl_psm_FL2img,23:191) sub_network_table_FL2(sjl_psm_FL2img,23:269)]);


p_stats=nan(1,size(network_conn_table_FL2,2));
cohen_d_stats=nan(1,size(network_conn_table_FL2,2));
p_value=nan(1,size(network_conn_table_FL2,2));



for i=1:size(network_conn_table_FL2,2)
    tempv=network_conn_table_FL2{sjl_psm_FL2img,i};
    ind_nn=~isnan(tempv);
    b=regress(tempv(ind_nn),[ones(length(tempv(ind_nn)),1) network_table_FL2{sjl_psm_FL2img(ind_nn),'rsfmri_c_ngd_ntpoints'}...
        ,network_table_FL2{sjl_psm_FL2img(ind_nn),'rsfmri_c_ngd_meanmotion'}]);
    Yhat=[ones(length(tempv(ind_nn)),1) network_table_FL2{sjl_psm_FL2img(ind_nn),'rsfmri_c_ngd_ntpoints'}...
        ,network_table_FL2{sjl_psm_FL2img(ind_nn),'rsfmri_c_ngd_meanmotion'}]*b;
    resid(ind_nn)=tempv(ind_nn)-Yhat+mean(tempv(ind_nn));
    resid(~ind_nn)=nan;



    [h,p]=ttest2(resid(ismember(sjl_psm_FL2img,psm_sjl0_FL2img)),resid(ismember(sjl_psm_FL2img,psm_sjl1_FL2img)));
    d=computeCohen_d(resid(ismember(sjl_psm_FL2img,psm_sjl0_FL2img)),resid(ismember(sjl_psm_FL2img,psm_sjl1_FL2img)),'independent');
    
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
p_test=p_value;
p_test(p_value>0)=fdr_bh(p_value(p_value>0),0.05);
% cohen_d_shaped= cohen_d_shaped.*p_test;
cohen_d_shaped(abs(cohen_d_shaped)<0.15)=0;
cohen_d_shaped=[tril(cohen_d_shaped(1:12,1:12));cohen_d_shaped(13:end,:)];

net_names={'auditory network', 'cingulo-opercular network', 'cingulo-parietal network', 'default network'...
    ,'dorsal attention network', 'fronto-parietal network','retrosplenial temporal network',...
    'sensorimotor hand network','sensorimotor mouth network', 'salience network',...
    'ventral attention network','visual network','left-cerebellum-cortex','left-thalamus-proper',...
    'left-caudate',' left-putamen','left-pallidum','brain-stem','left-hippocampus',...
    'left-amygdala','left-accumbens-area','left-ventraldc','right-cerebellum-cortex',...
    'right-thalamus-proper','right-caudate','right-putamen','right-pallidum', 'right-hippocampus'...
    'right-amygdala','right-accumbens-area','right-ventraldc'};

% create_matrix_ABCD(p_stats_shaped,'Total Sleep Duration FL2',net_names,0)

% create_matrix_ABCD_cohensD(cohen_d_shaped,'Total Sleep Duration Cohen''s D FL2',net_names,0) 
FL_net.p_value=p_value;
FL_net.cohen_d_stats=cohen_d_stats;
FL_net.cohen_d_threhold=cohen_d_shaped;

create_matrix_ABCD_cohensD(FL_net.cohen_d_threhold,'SJLsc',net_names,0) 




%% ind t test on greymatter volume 463:613 with TIV regressed out FL2

p_stats=nan(1,151);
cohen_d_stats=nan(1,151);
p_values=nan(1,151);


% sub_id_FL2img=sub_id_FL2(ismember(sub_id_FL2,T1_table_FL2.Properties.RowNames));
ind=ismember(psm_sjl0,T1_table_FL2.Properties.RowNames)&ismember(psm_sjl1,T1_table_FL2.Properties.RowNames);

% T1_sub_FL2
sleep_temp_FL2=sjl_psm([ind;ind]);
% ind_FL2_gmv=ismember(sleep_NE,sleep_temp_FL2)&ismember(sleep_FS,sleep_temp_FL2);%matched subject that both have FL2 data
% sleep_temp_FL2=[sleep_FS(ind_FL2_gmv);sleep_NE(ind_FL2_gmv)];
% ind=qc_table_FL2{sleep_FS(ind_FL2_gmv),'imgincl_t1w_include'}==1&qc_table_FL2{sleep_NE(ind_FL2_gmv),'imgincl_t1w_include'}==1;%matched subject that both passed qc

for i=1:151
    tempv=T1_table_FL2{sleep_temp_FL2,i+462};
    ind_nn=~isnan(tempv);
    b=regress(tempv(ind_nn),[ones(length(tempv(ind_nn)),1) T1_table_FL2{sleep_temp_FL2(ind_nn),'mrisdp_604'}]);
    Yhat=[ones(length(tempv(ind_nn)),1) T1_table_FL2{sleep_temp_FL2(ind_nn),'mrisdp_604'}]*b;
    resid=nan(length(tempv),1);
    resid(ind_nn)=tempv(ind_nn)-Yhat+mean(tempv(ind_nn));
    
    
%     resid(~ind_nn)=nan;
    [h,p]=ttest2(resid(ismember(sleep_temp_FL2,psm_sjl0)),resid(ismember(sleep_temp_FL2,psm_sjl1)));
    d=computeCohen_d(resid(ismember(sleep_temp_FL2,psm_sjl0)),resid(ismember(sleep_temp_FL2,psm_sjl1)),'independent');

    p_stats(i)=-log10(p)*sign(d);
    p_values(i)=p;
     cohen_d_stats(i)=d;
end


table_GMV_TIV_FL2=T1_table_BL(1,463:613);
table_GMV_TIV_FL2{1,:}=p_stats;
table_GMV_TIV_FL2{2,:}=cohen_d_stats;
table_GMV_TIV_FL2{3,:}=p_values;
table_GMV_TIV_FL2(:,table_GMV_TIV_FL2{3,:}<0.05) % mrisdp_570 right temporal pole cohensd = 0.0936 did not survive fdr correction

%% subcortical volume

p_stats=nan(1,46);
cohen_d_stats=nan(1,46);
p_values=nan(1,46);
% sleep_temp_FL2=sleep_psc_FL2(ismember(sleep_psc_FL2,T1_table_FL2.Properties.RowNames));
% ind_FL2_gmv=ismember(sleep_NE,sleep_temp_FL2)&ismember(sleep_FS,sleep_temp_FL2);%matched subject that both have FL2 data
% sub_id_FL2img=sub_id_FL2(ismember(sub_id_FL2,T1_table_FL2.Properties.RowNames));
ind=ismember(psm_sjl0,T1_sub_FL2.Properties.RowNames)&ismember(psm_sjl1,T1_sub_FL2.Properties.RowNames);

% T1_sub_FL2
sleep_temp_sub_FL2=sjl_psm([ind;ind]);
% ind=qc_table_FL2{sleep_FS(ind_FL2_gmv),'imgincl_t1w_include'}==1&qc_table_FL2{sleep_NE(ind_FL2_gmv),'imgincl_t1w_include'}==1;%matched subject that both passed qc

for i=1:46
       try
    tempv=T1_sub_FL2{sleep_temp_sub_FL2,i+329};
    ind_nn=~isnan(tempv);
    b=regress(tempv(ind_nn),[ones(length(tempv(ind_nn)),1) T1_sub_FL2{sleep_temp_sub_FL2(ind_nn),'smri_vol_scs_intracranialv'}]);
    Yhat=[ones(length(tempv(ind_nn)),1) T1_sub_FL2{sleep_temp_sub_FL2(ind_nn),'smri_vol_scs_intracranialv'}]*b;
    resid(ind_nn)=tempv(ind_nn)-Yhat+mean(tempv(ind_nn));
    resid(~ind_nn)=nan;
    [h,p]=ttest2(resid(ismember(sleep_temp_sub_FL2,psm_sjl0)),resid(ismember(sleep_temp_sub_FL2,psm_sjl1)));
    d=computeCohen_d(resid(ismember(sleep_temp_sub_FL2,psm_sjl0)),resid(ismember(sleep_temp_sub_FL2,psm_sjl1)),'independent');

    
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
table_subvol_TIV_FL2(:,table_subvol_TIV_FL2{3,:}<0.01) 
% table_thk_TIV(2,table_thk_TIV{2,1:148}>0.09)




%% fitbit
sub_id_7_sjl=sub_id_7(ismember(sub_id_7,sjl_psm));

ind_fitbit=ismember(psm_sjl0,sub_id_7_sjl)&ismember(psm_sjl1,sub_id_7_sjl);
cohen_d_est=nan(width(cbcl_raw_FL2),1);
cov=[psm_FL2_matched_table{sub_id_7_sjl,[2:12 15]}    ];
p=nan(width(cbcl_raw_FL2),1);r=nan(width(cbcl_raw_FL2),1);
for i = 1:width(cbcl_raw_FL2)
    ind_nan=~isnan(cbcl_raw_FL2{sub_id_7_sjl,i});
    if sum(ind_nan)<1000
        r(i)=0;
        p(i)=1;
        continue
    end
    [r(i),p(i)]=partialcorr(fitbit_daily_sum{sub_id_7_sjl(ind_nan),'sjlsc'},cbcl_raw_FL2{sub_id_7_sjl(ind_nan),i},cov(ind_nan,:));
    [b]=regress(zscore(cbcl_raw_FL2{sub_id_7_sjl(ind_nan),i}),[ones(length(fitbit_daily_sum{sub_id_7_sjl(ind_nan),'sjlsc'}),1) fitbit_daily_sum{sub_id_7_sjl(ind_nan),'sjlsc'}]);
%     cohen_d_est(i)=r(i)*b(2)/(std(fitbit_daily_sum{sub_id_7_sjl(ind_nan),'sjlsc'})*sqrt(1-r(i)^2));
end


% cohend_combined=cohend_FL2;
figure;
% b=bar(log10(p).*sign(r),'FaceColor','flat');
b=bar(-r,'FaceColor','flat');

for i=21:21+width(cognition_row_FL2)-1
b.CData(i,:) = [0 1 0];
end

for i=21+width(cognition_row_FL2):21+width(cognition_row_FL2)-1+width(mhy_row_FL2)
b.CData(i,:) = [1 1 0];
end
for i=43
b.CData(i,:) = [1 0 0];
end
hold on
% yline([-0.15 0.15],'--')

% xticks(1:width(cbcl_raw_BL));
set(gca, 'FontName', 'Helvetica','FontSize',12)
xticks(1:width(cbcl_raw_FL2));
xticklabels(strrep(cbcl_raw_FL2.Properties.VariableNames(1:43),'_','\_'))
ylabel({'negative partial correlation coeffecient'},'FontName', 'Helvetica','FontSize',18)


title({'Behavioral variables at FL2'},'FontSize',24)
set(gcf,'units','inches','position',[0,0,15,6])
set(gca,'TickLength',[0 0])
xtickangle(45)
box off


[r1,p1]=corr(fitbit_daily_sum{sub_id_7_sjl,'sjlsc'},psm_FL2_matched_table{sub_id_7_sjl,'SJLsc'})
% [r1,p1]=corr(fitbit_daily_sum{sub_id_7_sjl,'sjlsc'},psm_table_FL2_nonan{sub_id_7_sjl,'SJLsc'})
% [r1,p1]=corr(fitbit_daily_sum{sub_id_7_sjl,'avg_sleep'},psm_table_FL2_nonan{sub_id_7_sjl,'TS'})


[r1,p1]=corr(fitbit_daily_sum{sub_id_7_sjl,'avg_sleep'},psm_FL2_matched_table{sub_id_7_sjl,'TS'})

[r1,p1]=corr(fitbit_daily_sum{[psm_sjl0(ind_fitbit);psm_sjl1(ind_fitbit)],'sjlsc'},psm_FL2_matched_table{[psm_sjl0(ind_fitbit);psm_sjl1(ind_fitbit)],'SJLsc'})

[r, LB, UB, F, df1, df2, p] = ICC([fitbit_daily_sum{sub_id_7_sjl,'sjlsc'},psm_FL2_matched_table{sub_id_7_sjl,'SJLsc'}], 'C-1', 0.05, 0)
[r, LB, UB, F, df1, df2, p] = ICC([fitbit_daily_sum{sub_id_7_sjl,'avg_sleep'}/60,psm_FL2_matched_table{sub_id_7_sjl,'TS'}], 'C-1', 0.05, 0)


(sum(fitbit_daily_sum{sub_id_7_sjl,'sjlsc'}<=1&psm_FL2_matched_table{sub_id_7_sjl,'SJLsc'}<=1)+sum(fitbit_daily_sum{sub_id_7_sjl,'sjlsc'}>1&psm_FL2_matched_table{sub_id_7_sjl,'SJLsc'}>1))/length(sub_id_7_sjl)
%58.25% in the similar categories. 

x1=[sum(fitbit_daily_sum{sub_id_7_sjl,'sjlsc'}<=1),sum(fitbit_daily_sum{sub_id_7_sjl,'sjlsc'}>1)];
x2=[sum(psm_FL2_matched_table{sub_id_7_sjl,'SJLsc'}<=1),sum(psm_FL2_matched_table{sub_id_7_sjl,'SJLsc'}>1)];
[t,chi2,p2]=crosstab(fitbit_daily_sum{sub_id_7_sjl,'sjlsc'}>1,psm_FL2_matched_table{sub_id_7_sjl,'SJLsc'}>1)

 CFig_ABCD_sleep(fitbit_daily_sum{sub_id_7_sjl,'sjlsc'},psm_FL2_matched_table{sub_id_7_sjl,'SJLsc'},"Objective SJLsc","Subjecitve SJLsc")
  CFig_ABCD_sleep(fitbit_daily_sum{[psm_sjl0(ind_fitbit);psm_sjl1(ind_fitbit)],'sjlsc'},psm_FL2_matched_table{[psm_sjl0(ind_fitbit);psm_sjl1(ind_fitbit)],'SJLsc'},"Objective SJLsc","Subjecitve SJLsc")
  CFig_ABCD_sleep(psm_FL2_matched_table{sub_id_7_sjl,'SJLsc'},cbcl_raw_FL2{sub_id_7_sjl,'nihtbx_cryst_uncorrected'},"Objective SJLsc","Subjecitve SJLsc")




%% imaging FL2 fitbit


sub_id_7_sjl_img=sub_id_7(ismember(sub_id_7,sjl_psm_FL2img));



network_conn_table_FL2=([network_table_FL2(sub_id_7_sjl_img,23:191) sub_network_table_FL2(sub_id_7_sjl_img,23:269)]);


p_stats=nan(1,size(network_conn_table_FL2,2));
cohen_d_stats=nan(1,size(network_conn_table_FL2,2));
p_value=nan(1,size(network_conn_table_FL2,2));
r_value=nan(1,size(network_conn_table_FL2,2));
cov=[psm_FL2_matched_table{sub_id_7_sjl_img,[2:6 8:12 15]}    ];

for i=1:size(network_conn_table_FL2,2)
    [r_value(i),p_value(i)]=partialcorr(fitbit_daily_sum{sub_id_7_sjl_img,'sjlsc'},network_conn_table_FL2{sub_id_7_sjl_img,i},[cov,network_table_FL2{sub_id_7_sjl_img,'rsfmri_c_ngd_ntpoints'}...
        ,network_table_FL2{sub_id_7_sjl_img,'rsfmri_c_ngd_meanmotion'}], 'rows','pairwise' );




    p_stats(i)=log10(p_value(i))*sign(r_value(i));
%     p_value(i)=p;
%     cohen_d_stats(i)=d;
end
p_stats_shaped=[reshape(p_stats(1:169),13,13);reshape(p_stats(170:end),19,13);];
r_value=[reshape(r_value(1:169),13,13);reshape(r_value(170:end),19,13);];
r_value=r_value([1:6 8:end],[1:6 8:end]);
p_value=[reshape(p_value(1:169),13,13);reshape(p_value(170:end),19,13);];
p_value=p_value([1:6 8:end],[1:6 8:end]);
p_value=[tril(p_value(1:12,1:12));p_value(13:end,:)];
p_test=p_value;
p_test(p_value>0)=fdr_bh(p_value(p_value>0),0.1);
p_stats_shaped=p_stats_shaped([1:6 8:end],[1:6 8:end]);
p_stats_shaped=[tril(p_stats_shaped(1:12,1:12));p_stats_shaped(13:end,:)];
% p_stats_shaped= p_stats_shaped.*p_test;
p_stats_shaped(FL_net.cohen_d_threhold==0)=0;
r_value(FL_net.cohen_d_threhold==0)=0;

net_names={'auditory network', 'cingulo-opercular network', 'cingulo-parietal network', 'default network'...
    ,'dorsal attention network', 'fronto-parietal network','retrosplenial temporal network',...
    'sensorimotor hand network','sensorimotor mouth network', 'salience network',...
    'ventral attention network','visual network','left-cerebellum-cortex','left-thalamus-proper',...
    'left-caudate',' left-putamen','left-pallidum','brain-stem','left-hippocampus',...
    'left-amygdala','left-accumbens-area','left-ventraldc','right-cerebellum-cortex',...
    'right-thalamus-proper','right-caudate','right-putamen','right-pallidum', 'right-hippocampus'...
    'right-amygdala','right-accumbens-area','right-ventraldc'};

create_matrix_ABCD(p_stats_shaped,'SJLsc Fitbit',net_names,0)

create_matrix_ABCD_cohensD(-r_value,'Objecitvie SJLsc',net_names,0) 


% 
% FL_net.p_value=p_value;
% FL_net.cohen_d_stats=cohen_d_stats;
% FL_net.cohen_d_threhold=cohen_d_shaped;



