%%  initial parameters
clear;clc;
father_dir='/home/yangf7/Documents/Nils/ABCD';
cd(father_dir)
beh_dir=fullfile(father_dir,'abcd-data-release-5.1/core/');
addpath(genpath('/misc/imeel/yangf7/matlab/matlab_code_Nils'));

%% loading  var

%imaging

[network_table]=abcd_load_var2cell(fullfile(beh_dir,'imaging/mri_y_rsfmr_cor_gp_gp.csv'));%functional connectivities
[sub_network_table]=abcd_load_var2cell(fullfile(beh_dir,'imaging/mri_y_rsfmr_cor_gp_aseg.csv'));%subcortical-cortical FCs
% [qc_inclusion]=abcd_load_var2cell(fullfile(beh_dir,'imaging/mri_y_qc_incl.csv'));%


qc_inclusion=readtable(fullfile(beh_dir, 'imaging/mri_y_qc_incl.csv'));
qc_inclusion_table{1}=qc_inclusion(strcmp(qc_inclusion{:,'eventname'},'baseline_year_1_arm_1'),:);
qc_inclusion_table{3}=qc_inclusion(strcmp(qc_inclusion{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
qc_inclusion_table{5}=qc_inclusion(strcmp(qc_inclusion{:,'eventname'},'4_year_follow_up_y_arm_1'),:);

length(unique(qc_inclusion_table{1}.src_subject_id));
indx=strcmp(qc_inclusion_table{1}.src_subject_id,'NDAR_INV2F729N9A');
qc_inclusion_table{1}(find(indx==1,1),:)=[];
qc_inclusion_table{1}.Properties.RowNames=qc_inclusion_table{1}.src_subject_id;
[v, w] = unique( qc_inclusion_table{3}.src_subject_id, 'stable' );
duplicate_indices = setdiff( 1:numel(qc_inclusion_table{3}.src_subject_id), w );
qc_inclusion_table{3}=qc_inclusion_table{3}(w,:);% delete duplicates (2nd one is deleted)
qc_inclusion_table{3}.Properties.RowNames=qc_inclusion_table{3}.src_subject_id;
qc_inclusion_table{5}.Properties.RowNames=qc_inclusion_table{5}.src_subject_id;



[qc_table]=abcd_load_var2cell(fullfile(beh_dir,'imaging/mri_y_qc_man_post_fmr.csv'));%not useful
[qc_motion]=abcd_load_var2cell(fullfile(beh_dir,'imaging/mri_y_qc_motion.csv'));%motion

[T1_table]=abcd_load_var2cell(fullfile(beh_dir,'imaging/mri_y_smr_vol_dst.csv'));%gmv

[T1_sub_table]=abcd_load_var2cell(fullfile(beh_dir,'imaging/mri_y_smr_vol_aseg.csv'));%gmv subcortical
[T1_qc]=abcd_load_var2cell(fullfile(beh_dir,'imaging/mri_y_qc_raw_smr_t1.csv'));%t1 qc
[rfmri_qc]=abcd_load_var2cell(fullfile(beh_dir,'imaging/mri_y_qc_raw_rsfmr.csv'));%rfmri qc

[NIH_TB]=abcd_load_var2cell(fullfile(beh_dir,'imaging/mri_y_qc_raw_rsfmr.csv'));%NIH toolbox 

% mental health

[bisbas]=abcd_load_var2cell(fullfile(beh_dir,'mental-health/mh_y_bisbas.csv'));%bisbas
[pps]=abcd_load_var2cell(fullfile(beh_dir,'mental-health/mh_y_pps.csv'));%mh_y_pps
[upps]=abcd_load_var2cell(fullfile(beh_dir,'mental-health/mh_y_upps.csv'));%mh_y_upps



[cbcl_sum]=abcd_load_var2cell(fullfile(beh_dir,'mental-health/mh_p_cbcl.csv'));%cbcl
[sleep_parent]=abcd_load_var2cell(fullfile(beh_dir,'physical-health/ph_p_sds.csv'));%total score (sds_p_ss_total) for SDSC 



%% general demographic
 [pdem_table]=abcd_load_var2cell(fullfile(beh_dir,'abcd-general/abcd_p_demo.csv'));

 [pdem_table_BL,pdem_table_FL1,pdem_table_FL2,pdem_table_FL3,pdem_table_FL4]=abcd_load_var2table(fullfile(beh_dir,'abcd-general/abcd_p_demo.csv'));
 [youth_table_BL,youth_table_FL1,youth_table_FL2,youth_table_FL3,youth_table_FL4]=abcd_load_var2table(fullfile(beh_dir,'abcd-general/abcd_y_lt.csv'));

sub_id_BL=pdem_table_BL.Properties.RowNames;%subject to change
dem_table=table;
dem_table.src_subject_id=sub_id_BL;
dem_table.Properties.RowNames=dem_table.src_subject_id;
dem_table=[dem_table,youth_table_BL(sub_id_BL,3:10)];
dem_table{sub_id_BL,'sex'}= pdem_table_BL{sub_id_BL,'demo_sex_v2'}; 
dem_table{sub_id_BL,'gender'}= pdem_table_BL{sub_id_BL,'demo_gender_id_v2'};       
dem_table{sub_id_BL,'race'}= pdem_table_BL{sub_id_BL,'race_ethnicity'};   %1 = White; 2 = Black; 3 = Hispanic; 4 = Asian; 5 = Other
dem_table{sub_id_BL,'HI'}= pdem_table_BL{sub_id_BL,'demo_comb_income_v2'}; %  1= Less than $5,000; 2=$5,000 through $11,999; 3=$12,000 through $15,999; 4=$16,000 through $24,999; 5=$25,000 through $34,999; 6=$35,000 through $49,999; 7=$50,000 through $74,999; 8= $75,000 through $99,999; 9=$100,000 through $199,999; 10=$200,000 and greater. 999 = Don't know No lo s√É¬© ; 777 = Refuse to answer No deseo responde
dem_table{dem_table{sub_id_BL,'HI'}>100,'HI'}=nan;
for i=1:size(dem_table,1)
    dem_table{i,'abcd_site'}=str2double([dem_table{i,'site_id_l'}{1}(5:end)]);

end

pdem_table=readtable(fullfile(beh_dir,'abcd-general/abcd_p_demo.csv'));

prnt_edu=pdem_table(~isnan(pdem_table{:,'demo_prnt_ed_v2_l'}),[1,170]);
prnt_edu{prnt_edu{:,2}>100,2}=nan;
[v, w] = unique( prnt_edu.src_subject_id, 'stable' );
prnt_edu=prnt_edu(w,:);% delete duplicates (2nd one is deleted)

prnt_edu.Properties.RowNames=prnt_edu.src_subject_id;
dem_table{prnt_edu.Properties.RowNames,'edu'}= prnt_edu{prnt_edu.Properties.RowNames,'demo_prnt_ed_v2_l'};   %0 = Never attended/Kindergarten only Nunca asist?/Kinder solamente; 1 = 1st grade 1.er grado; 2 = 2nd grade 2.? grado; 3 = 3rd grade 3.er grado;4 = 4th grade 4.? grado; 5 = 5th grade 5.? grado; 6 = 6th grade 6.? grado; 7 = 7th grade 7.? grado; 8 = 8th grade 8.? grado; 9 = 9th grade 9.? grado; 10 = 10th grade 10.? grado; 11 = 11th grade 11.? grado; 12 = 12th grade, no diploma 12.? grado, sin certificado/diploma; 13 = High school graduate Preparatoria terminada; 14 = GED or equivalent Diploma General de Equivalencia (GED) o equivalente; 22 = Less than 1 year of college credit/post-secondary education (or less than 10 classes) Menos de 1 a?o de cr?dito universitario/educaci?n post-secundaria (o menos de 10 clases); 23 = One year or more of college credit, no degree Un a?o o mas de cr?dito universitario, sin t?tulo/diploma; 16 = Associate degree: Occupational, Technical, or Vocational T?tulo de asociado: programa ocupacional, t?cnico o vocacional; 17 = Associate degree: Academic Program T?tulo de asociado: programa acad?mico; 18 = Bachelor's degree (ex. BA, AB, BS, BBS) Licenciatura (p. ej., BA, AB, BS, BBA); 19 = Master's degree (ex. MA, MS, MEng, MEd, MBA) Maestr?a (p. ej., MA, MS, MEng, MEd, MBA); 20 = Professional School degree (ex. MD, DDS, DVN, JD) T?tulo de escuela profesional (p. ej., MD, DDS, DVM, JD); 21 = Doctoral degree (ex. PhD, EdD) Doctorado (p. ej., PhD, EdD); 777 = Refused to answer Prefiero no responder; 999 = Don't Know No lo s?
dem_table{dem_table{sub_id_BL,'edu'}>100,'edu'}=nan;

[anthro_table_BL,anthro_table_FL1,anthro_table_FL2,anthro_table_FL3,anthro_table_FL4]=abcd_load_var2table(fullfile(beh_dir,'physical-health/ph_y_anthro.csv'));
anthro_table_BL.BMI=703*anthro_table_BL.anthroweightcalc./(anthro_table_BL.anthroheightcalc).^2;
anthro_table_FL1.BMI=703*anthro_table_FL1.anthroweightcalc./(anthro_table_FL1.anthroheightcalc).^2;
anthro_table_FL2.BMI=703*anthro_table_FL2.anthroweightcalc./(anthro_table_FL2.anthroheightcalc).^2;
anthro_table_FL3.BMI=703*anthro_table_FL3.anthroweightcalc./(anthro_table_FL3.anthroheightcalc).^2;
anthro_table_FL4.BMI=703*anthro_table_FL4.anthroweightcalc./(anthro_table_FL4.anthroheightcalc).^2;



dem_table.BMI=nan(height(dem_table),1);
dem_table{anthro_table_BL.Properties.RowNames,'BMI'}= anthro_table_BL{anthro_table_BL.Properties.RowNames,'BMI'};  %There are 20 children have BMI lower than 10, consider to change it to nan.

[pds_BL,pds_FL1,pds_FL2,pds_FL3,pds_FL4]=abcd_load_var2table(fullfile(beh_dir,'physical-health/ph_y_pds.csv'));% youth report puberty
[ppds_BL,ppds_FL1,ppds_FL2,ppds_FL3,ppds_FL4]=abcd_load_var2table(fullfile(beh_dir,'physical-health/ph_p_pds.csv'));%parent report puberty

[ppds]=abcd_load_var2cell(fullfile(beh_dir,'physical-health/ph_p_pds.csv'));%parent report puberty



dem_table.pubertal_status=sum([ppds_BL{sub_id_BL,'pds_p_ss_female_category_2'},ppds_BL{sub_id_BL,'pds_p_ss_male_category_2'}],2,'omitnan');%parent report puberty
dem_table.pubertal_status(dem_table.pubertal_status==0)=nan;








%% 

mental_health_BL=[bisbas{1}(sub_id_BL,[24:3:42]),pps{1}(sub_id_BL,66:3:78),upps{1}(sub_id_BL,23:3:35)];

% linked-external-data
[rhds]=abcd_load_var2cell(fullfile(beh_dir,'linked-external-data/led_l_densbld.csv'));%resident density

[adi]=abcd_load_var2cell(fullfile(beh_dir,'linked-external-data/led_l_adi.csv'));%resident density
dem_table{adi{1}.Properties.RowNames,'ADI'}= adi{1}{adi{1}.Properties.RowNames,'reshist_addr1_adi_perc'};  % ADI for address 1 






% culture enviornment
[grade]=abcd_load_var2cell(fullfile(beh_dir,'culture-environment/ce_y_sag.csv'));%ce_y_sag






%% novel-tech

% fitbit_PA=readtable(fullfile(beh_dir,'abcd_fbdpas01.txt'));
fitbit_sleep=readtable(fullfile(beh_dir,'novel-technologies/nt_y_fitb_slp_d.csv'));
fitbit_sleep_BL=fitbit_sleep(strcmp(fitbit_sleep{:,'eventname'},'baseline_year_1_arm_1'),:);
% fitbit_sleep_week_BL.Properties.RowNames=fitbit_sleep_week_BL.src_subject_id;
fitbit_sleep_FL2=fitbit_sleep(strcmp(fitbit_sleep{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
fitbit_sleep_FL4=fitbit_sleep(strcmp(fitbit_sleep{:,'eventname'},'4_year_follow_up_y_arm_1'),:);

fitbit_sleep_week=readtable(fullfile(beh_dir,'novel-technologies/nt_y_fitb_slp_w.csv'));
fitbit_sleep_week_BL=fitbit_sleep_week(strcmp(fitbit_sleep_week{:,'eventname'},'baseline_year_1_arm_1'),:);
% fitbit_sleep_week_BL.Properties.RowNames=fitbit_sleep_week_BL.src_subject_id;
fitbit_sleep_week_FL2=fitbit_sleep_week(strcmp(fitbit_sleep_week{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
% fitbit_sleep_week_FL2.Properties.RowNames=fitbit_sleep_week_FL2.src_subject_id;
fitbit_sleep_week_FL4=fitbit_sleep_week(strcmp(fitbit_sleep_week{:,'eventname'},'4_year_follow_up_y_arm_1'),:);
sub_id_FB_BL=unique(fitbit_sleep_BL.src_subject_id);%79 have both timepoint
sub_id_FB_FL2=unique(fitbit_sleep_FL2.src_subject_id);
sub_id_FB_FL4=unique(fitbit_sleep_FL4.src_subject_id);

% 

fitbit_daily_BL=fitbit_sleep_BL;

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


fitbit_daily_test=fitbit_sleep_FL2;

% fitbit_daily_test=fitbit_daily_test(fitbit_daily_test.fit_ss_protocol_wear==1,:);
 fitbit_daily_test=fitbit_daily_test(fitbit_daily_test.fit_ss_wkno<4,:);
sub_id=unique(fitbit_daily_test.src_subject_id);
sub_id_7=[];
for i=1:length(sub_id)
    
   ind_sid=ismember(fitbit_daily_test.src_subject_id,sub_id(i));
   ind_sid(ind_sid==1)= hours(fitbit_daily_test{ind_sid,'fit_ss_wakeup_minutes'} - fitbit_daily_test{ind_sid,'fit_ss_first_sleep_minutes'})<24;
    
%     
%    ind_sid(ind_sid==1)= hours(fitbit_daily_test{ind_sid,'fit_ss_wakeup_minutes'} - fitbit_daily_test{ind_sid,'fit_ss_first_sleep_minutes'})<24&...
%         hours(fitbit_daily_test{ind_sid,'fit_ss_wakeup_minutes'} - fitbit_daily_test{ind_sid,'fit_ss_first_sleep_minutes'})>2;
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
% 
%     weekday_ind(weekday_ind==1)=(hours(fitbit_daily_test{weekday_ind,'fit_ss_wakeup_minutes'} - fitbit_daily_test{weekday_ind,'fit_ss_first_sleep_minutes'})<24);%make sure sleep duration less than 24 hours and longer than 2 hours. 
%     weekend_ind(weekend_ind==1)=(hours(fitbit_daily_test{weekend_ind,'fit_ss_wakeup_minutes'} - fitbit_daily_test{weekend_ind,'fit_ss_first_sleep_minutes'})<24);%make sure sleep duration less than 24 hours and longer than 2 hours. 

    weekday_ind(weekday_ind==1)=(hours(fitbit_daily_test{weekday_ind,'fit_ss_wakeup_minutes'} - fitbit_daily_test{weekday_ind,'fit_ss_first_sleep_minutes'})<15&...
        hours(fitbit_daily_test{weekday_ind,'fit_ss_wakeup_minutes'} - fitbit_daily_test{weekday_ind,'fit_ss_first_sleep_minutes'})>3);%make sure sleep duration less than 24 hours and longer than 2 hours. 
    weekend_ind(weekend_ind==1)=(hours(fitbit_daily_test{weekend_ind,'fit_ss_wakeup_minutes'} - fitbit_daily_test{weekend_ind,'fit_ss_first_sleep_minutes'})<15&...
        hours(fitbit_daily_test{weekend_ind,'fit_ss_wakeup_minutes'} - fitbit_daily_test{weekend_ind,'fit_ss_first_sleep_minutes'})>3);%make sure sleep duration less than 24 hours and longer than 2 hours. 

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





%%



% subtanse use
[sui]=abcd_load_var2cell(fullfile(beh_dir,'substance-use/su_y_sui.csv'));%substance use interview
[tlfb]=abcd_load_var2cell(fullfile(beh_dir,'substance-use/su_y_tlfb.csv'));%time follow back interview


% physical health

[mctq]=abcd_load_var2cell(fullfile(beh_dir,'physical-health/ph_y_mctq.csv'));% mctq


[pdem_table]=abcd_load_var2cell(fullfile(beh_dir,'abcd-general/abcd_p_demo.csv'));

% 
sub_id{1}=pdem_table{1}.Properties.RowNames;
for ii=2:length(sleep_parent)

    sub_id{ii}=sub_id{1}(ismember(sub_id{1},pdem_table{ii}.Properties.RowNames));
end


%% loading vars

var_table=readtable(fullfile(father_dir,'ABCD5.0_race_differences.xlsx'));
var_table.var_count=ones(height(var_table),1);
for n_var=1:height(var_table)
    if ~isempty(var_table{n_var,'VariableName5'}{:})
        var_table.var_count(n_var)=5;
        continue;
    end
    if ~isempty(var_table{n_var,'VariableName4'}{:})
        var_table.var_count(n_var)=4;
        continue;
    end
    if ~isempty(var_table{n_var,'VariableName3'}{:})
        var_table.var_count(n_var)=3;
        continue;
    end
    if ~isempty(var_table{n_var,'VariableName2'}{:})
        var_table.var_count(n_var)=2;
        continue;
    end

end

var_race=cell(height(var_table),5);
for n_var=[1:37 41:height(var_table)] %38-40 fitbit with duplicate names
    [var_race{n_var,1},var_race{n_var,2},var_race{n_var,3},var_race{n_var,4},var_race{n_var,5}]=abcd_load_var2table(...
        fullfile(beh_dir,var_table{n_var,'ABCD_folder'}{:},[var_table{n_var,'table_name'}{:},'.csv']),',');%
end

sub_var=cellfun(@height,var_race);
time_point=[];
for ii=[1:37 41:height(var_table)]
    try
        time_point(ii,1)=find(sub_var(ii,:)>1000,1);% more than 1000 entry
    catch
        time_point(ii,1)=0;
    end

end

time_point(21)=2;% data not available at BL
time_point(23)=2;% data not available at BL

var_time_point=sub_var>=1000;
var_time_point([21,23],1)=0;% data not available at BL

var_race{18,4}{:,'erq_ss_reappraisal_pr'}=sum(var_race{18,4}{:,[5:2:9]},2);
var_race{18,4}{:,'erq_ss_suppress_pr'}=sum(var_race{18,4}{:,[4:2:8]},2);
var_race{18,4}{var_race{18,4}{:,'erq_ss_reappraisal_pr'}>100,'erq_ss_reappraisal_pr'}=nan;
var_race{18,4}{var_race{18,4}{:,'erq_ss_suppress_pr'}>100,'erq_ss_suppress_pr'}=nan;


var_race{19,4}=var_race{18,4};

var_table_cell={};         temp_var=nan(size(sub_id_BL));

for ii= 1:size(var_time_point,2)
    var_table_cell{ii}=dem_table(:,1);
    for n_var=[1:45 47:height(var_table)] 
         if var_time_point(n_var,ii)>0
            var_table_cell{ii}=addvars(var_table_cell{ii},temp_var,'NewVariableNames',var_table{n_var,'VariableName1'}{:});
            var_table_cell{ii}{var_race{n_var,ii}.Properties.RowNames,var_table{n_var,'VariableName1'}{:}}=...
                var_race{n_var,ii}{:,var_table{n_var,'VariableName1'}{:}};
            if ~isempty(var_table{n_var,'VariableName2'}{:})
                var_table_cell{ii}=addvars(var_table_cell{ii},temp_var,'NewVariableNames',var_table{n_var,'VariableName2'}{:});
                var_table_cell{ii}{var_race{n_var,ii}.Properties.RowNames,var_table{n_var,'VariableName2'}{:}}=...
                var_race{n_var,ii}{:,var_table{n_var,'VariableName2'}{:}};
            end
            
            if ~isempty(var_table{n_var,'VariableName3'}{:})
                var_table_cell{ii}=addvars(var_table_cell{ii},temp_var,'NewVariableNames',var_table{n_var,'VariableName3'}{:});
                var_table_cell{ii}{var_race{n_var,ii}.Properties.RowNames,var_table{n_var,'VariableName3'}{:}}=...
                var_race{n_var,ii}{:,var_table{n_var,'VariableName3'}{:}};
            end
    
            if ~isempty(var_table{n_var,'VariableName4'}{:})
                var_table_cell{ii}=addvars(var_table_cell{ii},temp_var,'NewVariableNames',var_table{n_var,'VariableName4'}{:});
                var_table_cell{ii}{var_race{n_var,ii}.Properties.RowNames,var_table{n_var,'VariableName4'}{:}}=...
                var_race{n_var,ii}{:,var_table{n_var,'VariableName4'}{:}};
            end
    
            if ~isempty(var_table{n_var,'VariableName5'}{:})
                var_table_cell{ii}=addvars(var_table_cell{ii},temp_var,'NewVariableNames',var_table{n_var,'VariableName5'}{:});
                var_table_cell{ii}{var_race{n_var,ii}.Properties.RowNames,var_table{n_var,'VariableName5'}{:}}=...
                var_race{n_var,ii}{:,var_table{n_var,'VariableName5'}{:}};
            end







         end
    end
end


    var_table_cell{1}{anthro_table_BL.Properties.RowNames,'BMI'}=anthro_table_BL{:,'BMI'}; 
   
    var_table_cell{2}{anthro_table_FL1.Properties.RowNames,'BMI'}=anthro_table_FL1{:,'BMI'}; 
    var_table_cell{3}{anthro_table_FL2.Properties.RowNames,'BMI'}=anthro_table_FL2{:,'BMI'}; 
    var_table_cell{4}{anthro_table_FL3.Properties.RowNames,'BMI'}=anthro_table_FL3{:,'BMI'}; 
    var_table_cell{5}{anthro_table_FL4.Properties.RowNames,'BMI'}=anthro_table_FL4{:,'BMI'}; 

for ii=1:length(var_table_cell)
    var_table_cell{ii}{ppds{ii}.Properties.RowNames,'pubertal_status'}=sum([ppds{ii}{:,'pds_p_ss_female_category_2'},ppds{ii}{:,'pds_p_ss_male_category_2'}],2,'omitnan');
   
end
    for ii=1:length(var_table_cell)
            var_table_cell{ii}{var_table_cell{ii}{:,'BMI'}<5,'BMI'}=nan;
            var_table_cell{ii}{var_table_cell{ii}{:,'pubertal_status'}==0,'pubertal_status'}=nan;
            var_table_cell{ii}{var_table_cell{ii}{:,"BMI"}>=60&var_table_cell{ii}{:,"anthro_waist_cm"}<60,"BMI"}=nan;
            var_table_cell{ii}{var_table_cell{ii}{:,"BMI"}<=10&var_table_cell{ii}{:,"anthro_waist_cm"}>20,"BMI"}=nan;
            var_table_cell{ii}{var_table_cell{ii}{:,"anthro_waist_cm"}>=60&var_table_cell{ii}{:,"BMI"}<60,"anthro_waist_cm"}=nan;
            var_table_cell{ii}{var_table_cell{ii}{:,"anthro_waist_cm"}<=10,"anthro_waist_cm"}=nan;
            try
            var_table_cell{ii}{var_table_cell{ii}{:,'mctq_msfsc_calc'}>=24,'mctq_msfsc_calc'}=var_table_cell{ii}{var_table_cell{ii}{:,'mctq_msfsc_calc'}>=24,'mctq_msfsc_calc'}-24;
            catch
                continue
            end

    end






variable_table=dem_table(:,1);
for n_var=[1:45 47:height(var_table)] 
% for n_var=[18 19] %38-40 fitbit with duplicate names

    if time_point(n_var,1)>0
        temp_var=nan(size(sub_id_BL));
        variable_table=addvars(variable_table,temp_var,'NewVariableNames',var_table{n_var,'VariableName1'}{:});
        variable_table{var_race{n_var,time_point(n_var)}.Properties.RowNames,var_table{n_var,'VariableName1'}{:}}=...
            var_race{n_var,time_point(n_var)}{var_race{n_var,time_point(n_var)}.Properties.RowNames,var_table{n_var,'VariableName1'}{:}};
    
        if ~isempty(var_table{n_var,'VariableName2'}{:})
            variable_table=addvars(variable_table,temp_var,'NewVariableNames',var_table{n_var,'VariableName2'}{:});
            variable_table{var_race{n_var,time_point(n_var)}.Properties.RowNames,var_table{n_var,'VariableName2'}{:}}=...
                var_race{n_var,time_point(n_var)}{var_race{n_var,time_point(n_var)}.Properties.RowNames,var_table{n_var,'VariableName2'}{:}};
        end
        
        if ~isempty(var_table{n_var,'VariableName3'}{:})
            variable_table=addvars(variable_table,temp_var,'NewVariableNames',var_table{n_var,'VariableName3'}{:});
            variable_table{var_race{n_var,time_point(n_var)}.Properties.RowNames,var_table{n_var,'VariableName3'}{:}}=...
                var_race{n_var,time_point(n_var)}{var_race{n_var,time_point(n_var)}.Properties.RowNames,var_table{n_var,'VariableName3'}{:}};
        end

        if ~isempty(var_table{n_var,'VariableName4'}{:})
            variable_table=addvars(variable_table,temp_var,'NewVariableNames',var_table{n_var,'VariableName4'}{:});
            variable_table{var_race{n_var,time_point(n_var)}.Properties.RowNames,var_table{n_var,'VariableName4'}{:}}=...
                var_race{n_var,time_point(n_var)}{var_race{n_var,time_point(n_var)}.Properties.RowNames,var_table{n_var,'VariableName4'}{:}};
        end

        if ~isempty(var_table{n_var,'VariableName5'}{:})
            variable_table=addvars(variable_table,temp_var,'NewVariableNames',var_table{n_var,'VariableName5'}{:});
            variable_table{var_race{n_var,time_point(n_var)}.Properties.RowNames,var_table{n_var,'VariableName5'}{:}}=...
                var_race{n_var,time_point(n_var)}{var_race{n_var,time_point(n_var)}.Properties.RowNames,var_table{n_var,'VariableName5'}{:}};
        end
    end
end

variable_table.sleep_fitbit=nan(size(sub_id_BL));
variable_table{sub_id_BL(ismember(sub_id_BL,fitbit_daily_sum.sub_id)),'sleep_fitbit'}=fitbit_daily_sum{sub_id_BL(ismember(sub_id_BL,fitbit_daily_sum.sub_id)),'avg_sleep'};

variable_table.timeinbed_fitbit=nan(size(sub_id_BL));
variable_table{sub_id_BL(ismember(sub_id_BL,fitbit_daily_sum.sub_id)),'timeinbed_fitbit'}=fitbit_daily_sum{sub_id_BL(ismember(sub_id_BL,fitbit_daily_sum.sub_id)),'avg_timeinbed'};





psm_table=dem_table(:,[3,8,10,12:18]);

% psm_table.Properties.VariableNames={'age','race','sex','prt_edu','sites','family','BMI'};
psm_table.age_sex=psm_table.interview_age.*psm_table.sex;
psm_table{psm_table.Properties.RowNames,'BMI'}=var_table_cell{1}{psm_table.Properties.RowNames,'BMI'};
psm_table_race=psm_table(psm_table.race==1|psm_table.race==2,:);%n=7957
% % psm_table.puberty=pubertal;
% psm_table.HI=House_in;
% TS=double(total_sleep_duration);
% psm_table.TS=TS;



%abcd_mid02
%tfmri_mid_all_beh_srw_mrt	Float		Recommended	Average reaction time for small reward trials for run 1 and run 2
%tfmri_mid_all_beh_srwpfb_mrt               Average reaction time for small reward trials with positive feedback for run 1 and run 2 
% tfmri_mid_all_beh_lrw_mrt	Float		Recommended	Average reaction time for large reward trials for run 1 and run 2
%tfmri_mid_all_beh_lrwpfb_mrt
% tfmri_mid_all_beh_hrw_mrt	Float		Recommended	Average reaction time for large and small reward trials for run 1 and run 2
% tfmri_mid_all_beh_sl_mrt	Float		Recommended	Average reaction time for small loss trials for run 1 and run 2
%tfmri_mid_all_beh_ll_mrt
%tfmri_mid_all_beh_hlpfb_mrt
% tfmri_mid_all_beh_hl_mrt	Float		Recommended	Average reaction time for large and small loss trials for run 2
%tfmri_mid_all_beh_nt_mrt	Float		Recommended	Average reaction time for neutral trials for run 1 and run 2
%tfmri_mid_all_beh_ntpfb_mrt

%% runing psm
psm_table_nonan=psm_table_race(~isnan(sum(psm_table_race{:,2:11},2)),:);% exclued subs with missing values; n=6988

psm_table_nonan(psm_table_nonan.sex==3,:)=[];%n=6986
race_org=psm_table_nonan.Properties.RowNames;
race_org1=race_org(psm_table_nonan{race_org,'race'}==2);
race_org0=race_org(psm_table_nonan{race_org,'race'}==1);

% sum(isnan(psm_table_FL2{:,1:14}))
psm_table_nonan.num=(1:height(psm_table_nonan))';
writetable(psm_table_nonan,'psm_table_race_50_nooutlier.csv')


% mdl = fitglm(psm_table_FL2_nonan,'SJLsc_b~age+sex+race+prt_edu+sites+HI+BMI_BL+puberty_BL+age_sex+TS+SDis','Distribution','binomial');
%% load matched table


% psm_matched=readtable('/home/yangf7/Documents/Nils/ABCD/psm_table_race_matched50_nooutlier.csv');
psm_matched=readtable('/home/yangf7/Documents/Nils/ABCD/race/psm_table_race_matched50_nooutlier_nopuberty.csv');

psm_matched.Properties.RowNames=psm_table_nonan.Properties.RowNames(psm_matched.num);
psm_matched_table=sortrows(psm_matched,[5 16]);%group_diff and subclass
psm_race0=psm_matched_table.Properties.RowNames(1:height(psm_matched_table)/2);%white
psm_race1=psm_matched_table.Properties.RowNames(height(psm_matched_table)/2+1:end);%black
race_psm=[psm_race0;psm_race1];

psm_matched_withpub=readtable('/home/yangf7/Documents/Nils/ABCD/race/psm_table_race_matched50_nooutlier_withpuberty.csv');
psm_matched_withpub.Properties.RowNames=psm_table_nonan.Properties.RowNames(psm_matched_withpub.num);
psm_matched_withpub=sortrows(psm_matched_withpub,[5 16]);%group_diff and subclass
psmwp_race0=psm_matched_withpub.Properties.RowNames(1:height(psm_matched_withpub)/2);%white
psmwp_race1=psm_matched_withpub.Properties.RowNames(height(psm_matched_withpub)/2+1:end);%black
race_psmwp=[psmwp_race0;psmwp_race1];





%% behavioral plots
options=statset("bootci");options.UseParallel=1;

cohend_d=nan(2,width(var_table_cell{1})-1,length(var_table_cell));  cohend_d_wp=nan(2,width(var_table_cell{1})-1,length(var_table_cell)); 
cohend_d_CI=nan(2,width(var_table_cell{1})-1,length(var_table_cell));cohend_d_wp_CI=nan(2,width(var_table_cell{1})-1,length(var_table_cell));
p_value=nan(2,width(var_table_cell{1})-1,length(var_table_cell));p_value_wp=nan(2,width(var_table_cell{1})-1,length(var_table_cell));

numofsubs=nan(1,width(var_table_cell{1})-1,length(var_table_cell));numofsubs_wp=nan(1,width(var_table_cell{1})-1,length(var_table_cell));

cohend_d_CI_psm=nan(2,width(var_table_cell{1})-1,length(var_table_cell));cohend_d_wp_CI_psm=nan(2,width(var_table_cell{1})-1,length(var_table_cell));

for ii=1:length(var_table_cell)
    for n_var=1:width(var_table_cell{ii})-1
        if mean(isnan(var_table_cell{ii}{psm_race0,n_var+1}))>0.9
              cohend_d(1,n_var,ii)=0;
              p_value(1,n_var,ii)=1;
              
        else
            ind_nn=~isnan(var_table_cell{ii}{psm_race0,n_var+1})&~isnan(var_table_cell{ii}{psm_race1,n_var+1});% only both PSM pairs have valid values
            %         effect=meanEffectSize(beh_BL{fitbit_rs,i},beh_BL{fitbit_is,i},Effect="cohen")
            cohen_temp=[];
            cohen_temp=meanEffectSize(var_table_cell{ii}{psm_race1(ind_nn),n_var+1},var_table_cell{ii}{psm_race0(ind_nn),n_var+1},...
            BootstrapOptions=options, ConfidenceIntervalType="exact" ,Effect="cohen");
            cohend_d(1,n_var,ii)=cohen_temp.Effect;
            cohend_d_CI_psm(:,n_var,ii)=cohen_temp.ConfidenceIntervals;
            
            [h,p_value(1,n_var,ii)]=ttest2(var_table_cell{ii}{psm_race1(ind_nn),n_var+1},var_table_cell{ii}{psm_race0(ind_nn),n_var+1});
            numofsubs(1,n_var,ii)=sum(ind_nn);
            tempv=var_table_cell{ii}{race_org,n_var+1};cohen_temp=[];
            ind_nn=~isnan(tempv);
            b=regress(tempv(ind_nn),[ones(length(tempv(ind_nn)),1) ...
            psm_table_nonan{race_org(ind_nn),[2,3,5:7 10:11]}-mean(psm_table_nonan{race_org(ind_nn),[2,3,5:7 10:11]},'omitnan')]);
            Yhat=[ones(length(tempv(ind_nn)),1)  psm_table_nonan{race_org(ind_nn),[2,3,5:7 10:11]}-mean(psm_table_nonan{race_org(ind_nn),[2,3,5:7 10:11]},'omitnan')]*b;
            resid(ind_nn)=tempv(ind_nn)-Yhat+mean(tempv(ind_nn));
            resid(~ind_nn)=nan;
            
            
            cohen_temp=meanEffectSize(resid(ismember(race_org,race_org1)),resid(ismember(race_org,race_org0)),...
            BootstrapOptions=options, ConfidenceIntervalType="exact" ,Effect="cohen");
            cohend_d(2,n_var,ii)=cohen_temp.Effect;
            cohend_d_CI(:,n_var,ii)=cohen_temp.ConfidenceIntervals;
            [h,p_value(2,n_var,ii)]=ttest2(resid(ismember(race_org,race_org1)),resid(ismember(race_org,race_org0)));



            ind_nn_wp=~isnan(var_table_cell{ii}{psmwp_race0,n_var+1})&~isnan(var_table_cell{ii}{psmwp_race1,n_var+1});% only both PSM pairs have valid values
            %         effect=meanEffectSize(beh_BL{fitbit_rs,i},beh_BL{fitbit_is,i},Effect="cohen")
            cohen_temp=[];
            cohen_temp=meanEffectSize(var_table_cell{ii}{psmwp_race1(ind_nn_wp),n_var+1},var_table_cell{ii}{psmwp_race0(ind_nn_wp),n_var+1},...
            BootstrapOptions=options, ConfidenceIntervalType="exact" ,Effect="cohen");
            cohend_d_wp(1,n_var,ii)=cohen_temp.Effect;
            cohend_d_wp_CI_psm(:,n_var,ii)=cohen_temp.ConfidenceIntervals;
            
            [h,p_value_wp(1,n_var,ii)]=ttest2(var_table_cell{ii}{psmwp_race1(ind_nn_wp),n_var+1},var_table_cell{ii}{psmwp_race0(ind_nn_wp),n_var+1});
            numofsubs_wp(1,n_var,ii)=sum(ind_nn_wp);
            tempv=var_table_cell{ii}{race_org,n_var+1};cohen_temp=[];
            ind_nn_wp=~isnan(tempv);
            b=regress(tempv(ind_nn_wp),[ones(length(tempv(ind_nn_wp)),1) ...
            psm_table_nonan{race_org(ind_nn_wp),[2,3,5:7 9:11]}-mean(psm_table_nonan{race_org(ind_nn_wp),[2,3,5:7 9:11]},'omitnan')]);
            Yhat=[ones(length(tempv(ind_nn_wp)),1)  psm_table_nonan{race_org(ind_nn_wp),[2,3,5:7 9:11]}-mean(psm_table_nonan{race_org(ind_nn_wp),[2,3,5:7 9:11]},'omitnan')]*b;
            resid(ind_nn_wp)=tempv(ind_nn_wp)-Yhat+mean(tempv(ind_nn_wp));
            resid(~ind_nn_wp)=nan;
            
            
            cohen_temp=meanEffectSize(resid(ismember(race_org,race_org1)),resid(ismember(race_org,race_org0)),...
            BootstrapOptions=options, ConfidenceIntervalType="exact" ,Effect="cohen");
            cohend_d_wp(2,n_var,ii)=cohen_temp.Effect;
            cohend_d_wp_CI(:,n_var,ii)=cohen_temp.ConfidenceIntervals;
            [h,p_value_wp(2,n_var,ii)]=ttest2(resid(ismember(race_org,race_org1)),resid(ismember(race_org,race_org0)));

            
            % cohend_d(2,n_var,ii)=computeCohen_d(resid(ismember(race_org,race_org1)),resid(ismember(race_org,race_org0)),'independent');
            

    % cohend_BL(3,i)=computeCohen_d(beh_BL_mean{race_org1,i},beh_BL_mean{race_org0,i},'independent');

        end
    end

end
% p_value=-log10(p_value);
% 
% [~,~,~,beh_fdr_p_FL2]=fdr_bh(p_value);
% beh_fdr_p_FL2(beh_fdr_p_FL2<0.0001)=0.0001;

data_points={'Baseline','One-year Follow-up','Two-year Follow-up','Three-year Follow-up','Four-year Follow-up'};

for ii=1:length(var_table_cell)
    ind_nz=cohend_d(1,:,ii)~=0&numofsubs(1,:,ii)>=10;% num of psm paris higher than 10
    indx=1:sum(ind_nz);
    indx=[indx;cohend_d(:,ind_nz,ii)];
    [B,I2]=sort(indx(3,:));
% [~,I2]=sort(indx(4,:));
    temp_cohen=cohend_d(:,ind_nz,ii);
    temp_cohen_CI=cohend_d_CI(:,ind_nz,ii);
    temp_cohen_CI_psm=cohend_d_CI_psm(:,ind_nz,ii);

    figure;
    % errorbar(temp_cohen(1,I2)',indx(1,:),temp_cohen_CI_psm(1,:),temp_cohen_CI_psm(2,:),"horizontal",'*','color','r','DisplayName','PSM-matched')

    plot(temp_cohen(1,I2)',indx(1,:),'.',"MarkerSize",32,'color','r','DisplayName','PSM')
    xlim([-1,1])
    hold on 
    errorbar(temp_cohen(2,I2)',indx(1,:),abs(temp_cohen_CI(1,:)-temp_cohen_CI(2,:))/2,abs(temp_cohen_CI(1,:)-temp_cohen_CI(2,:))/2,"horizontal",'.',...
        'color','b','LineWidth',2,"MarkerSize",18,'capsize',0,'DisplayName','Regression w/ 95%CI')
    % plot(temp_cohen(2,I2)',indx(1,:),'o','color','b','DisplayName','Covariates-corrected')
    hold on 
    % plot(cohend_BL(3,I2)',1:size(cohend_BL,2),'o','color','magenta','DisplayName','No-control')
    set(gca, 'YTick', indx(1,:)); % center y-axis ticks on bins
    set(gca, 'TickLength',[0 0])
    var_names=var_table_cell{ii}.Properties.VariableNames(logical([0,ind_nz(1:width(var_table_cell{ii})-1)]));
    set(gca, 'YTickLabel',strrep( var_names(I2),'_','\_')              ,'FontSize',18); % set y-axis labels
    xline([ 0 ],'--')
    xline([ 0.15],'-','0.15')
        xline([ -0.15],'-','-0.15')

    yline([indx(1,:)],'-.','Alpha',0.2)
    f=get(gca,'Children');
    legend([f(end-1),f(end)])
    ylim([0,indx(1,end)+1])
    title(data_points{ii},'FontSize',24)
    
    legend('location','best')
    set(gcf,'Position',[2500 400 1200 1200])
    print(['/home/yangf7/Documents/Nils/ABCD/race/',data_points{ii}],'-depsc')
    % printeps(1,strrep(['/home/yangf7/Documents/Nils/ABCD/race/',data_points{ii}],' ','_'))
    sum(abs(temp_cohen(1,:))>abs((temp_cohen(2,:))))./sum(temp_cohen(2,:)~=0)
    sum((abs(temp_cohen(1,:)-temp_cohen(2,:)))>abs(temp_cohen_CI(1,:)-temp_cohen_CI(2,:))/2)./sum(temp_cohen(2,:)~=0)

    close;
end



for ii=1:length(var_table_cell)
    ind_nz=cohend_d_wp(1,:,ii)~=0&numofsubs_wp(1,:,ii)>=10;% num of psm paris higher than 10
    ind_nz(find(ind_nz,1,'last'))=0;
    indx=1:sum(ind_nz);
    indx=[indx;cohend_d_wp(:,ind_nz,ii)];
    [B,I2]=sort(indx(3,:));
% [~,I2]=sort(indx(4,:));
    temp_cohen=cohend_d_wp(:,ind_nz,ii);
    temp_cohen_CI=cohend_d_wp_CI(:,ind_nz,ii);
    temp_cohen_CI_psm=cohend_d_wp_CI_psm(:,ind_nz,ii);

    figure;
    % errorbar(temp_cohen(1,I2)',indx(1,:),temp_cohen_CI_psm(1,:),temp_cohen_CI_psm(2,:),"horizontal",'*','color','r','DisplayName','PSM-matched')

    plot(temp_cohen(1,I2)',indx(1,:),'.',"MarkerSize",24,'color','r','DisplayName','PSM')
    xlim([-1,1])
    hold on 
    errorbar(temp_cohen(2,I2)',indx(1,:),abs(temp_cohen_CI(1,:)-temp_cohen_CI(2,:))/2,abs(temp_cohen_CI(1,:)-temp_cohen_CI(2,:))/2,"horizontal",'.',...
        'color','b','LineWidth',2,"MarkerSize",18,'capsize',0,'DisplayName','Regression w/ 95%CI')
    % plot(temp_cohen(2,I2)',indx(1,:),'o','color','b','DisplayName','Covariates-corrected')
    hold on 
    % plot(cohend_BL(3,I2)',1:size(cohend_BL,2),'o','color','magenta','DisplayName','No-control')
    set(gca, 'YTick', indx(1,:)); % center y-axis ticks on bins
    set(gca, 'TickLength',[0 0])
    var_names=var_table_cell{ii}.Properties.VariableNames(logical([0,ind_nz(1:width(var_table_cell{ii})-1)]));
    set(gca, 'YTickLabel',strrep( var_names(I2),'_','\_')              ,'FontSize',18); % set y-axis labels
    xline([ 0 ],'--')
    xline([ 0.15],'-','0.15')
        xline([ -0.15],'-','-0.15')

    yline([indx(1,:)],'-.','Alpha',0.2)
    f=get(gca,'Children');
    legend([f(end-1),f(end)])
    ylim([0,indx(1,end)])
    title([data_points{ii},'\_with\_puberty'],'FontSize',24)
    
    legend('location','best')
    set(gcf,'Position',[2500 400 1600 1200])
    print(['/home/yangf7/Documents/Nils/ABCD/race/','with_puberty_',data_points{ii}],'-depsc')
    a(ii)=sum(abs(temp_cohen(1,:))>abs((temp_cohen(2,:))))./sum(temp_cohen(2,:)~=0)
    b(ii)=sum((abs(temp_cohen(1,:)-temp_cohen(2,:)))>abs(temp_cohen_CI(1,:)-temp_cohen_CI(2,:))/2)./sum(temp_cohen(2,:)~=0)
    close;
end




% dep_med=[beh_BL(race_psm,:)];

cov_med=psm_table_nonan{race_psm,[2,3,5:7 10:11]}-mean(psm_table_nonan{race_psm,[2,3,5:7 10:11]},'omitnan');

p_med=cell(5,1);
beta_med=cell(5,1);
for ii=1:length(var_table_cell)
        ind_nz=cohend_d(1,:,ii)~=0&numofsubs(1,:,ii)>=20;
        indx=1:sum(ind_nz);
        mediator=var_table_cell{ii}{race_psm,'pubertal_status'};
        dep_med=var_table_cell{ii}{race_psm,logical([0,ind_nz(1:width(var_table_cell{ii})-2)])};
    for jj=1:sum(ind_nz)-1
        [paths, stats_med]=mediation(psm_table_nonan{race_psm,'race'},dep_med(:,jj),...
        mediator...
        ,'cov', cov_med,'verbose','boot','bootsamples',10000);
        p_med{ii}(jj,:)=stats_med.p;
        beta_med{ii}(jj,:)=stats_med.mean;
    end




end

cov_med=psm_table_nonan{race_org,[2,3,5:7 10:11]}-mean(psm_table_nonan{race_org,[2,3,5:7 10:11]},'omitnan');

p_med=cell(5,1);
beta_med=cell(5,1);
for ii=1:length(var_table_cell)
        ind_nz=cohend_d(1,:,ii)~=0&numofsubs(1,:,ii)>=20;
        indx=1:sum(ind_nz);
        % [coeff,score,latent,tsquared,explained,mu]=pca([var_table_cell{ii}{race_psm,'pubertal_status'},var_table_cell{ii}{race_psm,'hormone_scr_dhea_mean'},...
        %     var_table_cell{ii}{race_psm,'hormone_scr_ert_mean'}],'rows','pairwise');
        % mediator=score(:,1);
        mediator=var_table_cell{ii}{race_org,'pubertal_status'};
        dep_med=var_table_cell{ii}{race_org,logical([0,ind_nz(1:width(var_table_cell{ii})-2)])};
    for jj=1:sum(ind_nz)-1
        [paths, stats_med]=mediation(psm_table_nonan{race_org,'race'},dep_med(:,jj),...
        mediator...
        ,'cov', cov_med,'verbose','boot','bootsamples',10000);
        p_med{ii}(jj,:)=stats_med.p;
        beta_med{ii}(jj,:)=stats_med.mean;
    end




end
% 
% 
% 





for ii=1:length(var_table_cell)

    percentofmed(ii)=sum(p_med{ii}(:,4)<=0.05&p_med{ii}(:,5)<=0.05)./sum(p_med{ii}(:,4)<=0.05);
end


cov_med=psm_table_nonan{race_psm,[2,3,5:7 10:11]}-mean(psm_table_nonan{race_psm,[2,3,5:7 10:11]},'omitnan');

p_med_psm=cell(5,1);
beta_med=cell(5,1);
for ii=1:length(var_table_cell)
        ind_nz=cohend_d(1,:,ii)~=0&numofsubs(1,:,ii)>=20;
        indx=1:sum(ind_nz);
        % [coeff,score,latent,tsquared,explained,mu]=pca([var_table_cell{ii}{race_psm,'pubertal_status'},var_table_cell{ii}{race_psm,'hormone_scr_dhea_mean'},...
        %     var_table_cell{ii}{race_psm,'hormone_scr_ert_mean'}],'rows','pairwise');
        % mediator=score(:,1);
        mediator=var_table_cell{ii}{race_psm,'pubertal_status'};
        dep_med_psm=var_table_cell{ii}{race_psm,logical([0,ind_nz(1:width(var_table_cell{ii})-2)])};
    for jj=1:sum(ind_nz)-1
        [paths, stats_med]=mediation(psm_table_nonan{race_psm,'race'},dep_med_psm(:,jj),...
        mediator...
        ,'cov', cov_med,'verbose','boot','bootsamples',10000);
        p_med_psm{ii}(jj,:)=stats_med.p;
        beta_med{ii}(jj,:)=stats_med.mean;
    end




end
% 




for ii=1:length(var_table_cell)

    percentofmed_psm(ii)=sum(p_med_psm{ii}(:,4)<=0.05&p_med_psm{ii}(:,5)<=0.05)./sum(p_med_psm{ii}(:,4)<=0.05);
end


%% errorbar with shade

ind_pub=[43,26,53,39,52]-1;
ind_BMI=ind_pub-1;
ind_Waist=ind_pub-2;
cohen_pub=[cohend_d(:,ind_pub(1),1),cohend_d(:,ind_pub(2),2),cohend_d(:,ind_pub(3),3),cohend_d(:,ind_pub(4),4),cohend_d(:,ind_pub(5),5)];
cohen_pub_CI=[cohend_d_CI(:,ind_pub(1),1),cohend_d_CI(:,ind_pub(2),2),cohend_d_CI(:,ind_pub(3),3),cohend_d_CI(:,ind_pub(4),4),cohend_d_CI(:,ind_pub(5),5)];
cohen_pub_CI=abs(cohen_pub_CI(1,:)-cohen_pub_CI(2,:))/2;

cohen_BMI=[cohend_d(:,ind_BMI(1),1),cohend_d(:,ind_BMI(2),2),cohend_d(:,ind_BMI(3),3),cohend_d(:,ind_BMI(4),4),cohend_d(:,ind_BMI(5),5)];
cohen_BMI_CI=[cohend_d_CI(:,ind_BMI(1),1),cohend_d_CI(:,ind_BMI(2),2),cohend_d_CI(:,ind_BMI(3),3),cohend_d_CI(:,ind_BMI(4),4),cohend_d_CI(:,ind_BMI(5),5)];
cohen_BMI_CI=abs(cohen_BMI_CI(1,:)-cohen_BMI_CI(2,:))/2;

cohen_BMI_wp=[cohend_d_wp(:,ind_BMI(1),1),cohend_d_wp(:,ind_BMI(2),2),cohend_d_wp(:,ind_BMI(3),3),cohend_d_wp(:,ind_BMI(4),4),cohend_d_wp(:,ind_BMI(5),5)];
cohen_BMI_wp_CI=[cohend_d_wp_CI(:,ind_BMI(1),1),cohend_d_wp_CI(:,ind_BMI(2),2),cohend_d_wp_CI(:,ind_BMI(3),3),cohend_d_wp_CI(:,ind_BMI(4),4),cohend_d_wp_CI(:,ind_BMI(5),5)];
cohen_BMI_wp_CI=abs(cohen_BMI_wp_CI(1,:)-cohen_BMI_wp_CI(2,:))/2;

cohen_Waist=[cohend_d(:,ind_Waist(1),1),cohend_d(:,ind_Waist(2),2),cohend_d(:,ind_Waist(3),3),cohend_d(:,ind_Waist(4),4),cohend_d(:,ind_Waist(5),5)];
cohen_Waist_CI=[cohend_d_CI(:,ind_Waist(1),1),cohend_d_CI(:,ind_Waist(2),2),cohend_d_CI(:,ind_Waist(3),3),cohend_d_CI(:,ind_Waist(4),4),cohend_d_CI(:,ind_Waist(5),5)];
cohen_Waist_CI=abs(cohen_Waist_CI(1,:)-cohen_Waist_CI(2,:))/2;

cohen_Waist_wp=[cohend_d_wp(:,ind_Waist(1),1),cohend_d_wp(:,ind_Waist(2),2),cohend_d_wp(:,ind_Waist(3),3),cohend_d_wp(:,ind_Waist(4),4),cohend_d_wp(:,ind_Waist(5),5)];
cohen_Waist_wp_CI=[cohend_d_wp_CI(:,ind_Waist(1),1),cohend_d_wp_CI(:,ind_Waist(2),2),cohend_d_wp_CI(:,ind_Waist(3),3),cohend_d_wp_CI(:,ind_Waist(4),4),cohend_d_wp_CI(:,ind_Waist(5),5)];
cohen_Waist_wp_CI=abs(cohen_Waist_wp_CI(1,:)-cohen_Waist_wp_CI(2,:))/2;




figure;
h=shadedErrorBar(1:5,cohen_pub(2,:),cohen_pub_CI,'lineProps',{'.b-','MarkerSize',32},'patchSaturation',0.05);
hold on
plot(1:5,cohen_pub(1,:),'r.-',"MarkerSize",32)
set(gca, 'TickLength',[0 0])
% yline([ 0.15],'-','0.15')
ylim([0 0.7])
xlim([0.5 5.5])
set(gca, 'XTick', 1:5); % center x-axis ticks on bins
% set(gca, 'YTick', 1:5,'FontWeight','Normal'); % center x-axis ticks on bins

set(gca, 'XTickLabel',data_points  ,'FontSize',24,'FontWeight','Bold'); % set y-axis labels
ylabel('Cohen''s d'  ,'FontSize',28,'FontWeight','Bold')
legend('Regression w/ 95%CI','PSM','FontSize',24,'FontWeight','Normal')
title('Pubertal Status','FontSize',24)

legend('location','best')
set(gcf,'Position',[2500 400 2000 1200])
print(['/home/yangf7/Documents/Nils/ABCD/race/','pubertalstatus'],'-depsc')
% printeps(1,['/home/yangf7/Documents/Nils/ABCD/race/','pubertalstatus'])

close


figure;
subplot(1,2,1)
h=shadedErrorBar(1:5,cohen_BMI(2,:),cohen_BMI_CI,'lineProps',{'.b-','MarkerSize',32},'patchSaturation',0.05);
hold on
plot(1:5,cohen_BMI(1,:),'r.-',"MarkerSize",32)
set(gca, 'TickLength',[0 0])
% yline([ 0.15],'-','0.15')
ylim([0 0.7])
xlim([0.5 5.5])
set(gca, 'XTick', 1:5); % center x-axis ticks on bins
set(gca, 'XTickLabel',data_points  ,'FontSize',24,'FontWeight','Bold'); % set y-axis labels
ylabel('Cohen''s d'  ,'FontSize',28,'FontWeight','Bold')
legend('Regression w/ 95%CI','PSM','FontSize',24,'FontWeight','Normal')
title('BMI','FontSize',24)

legend('location','best')
% set(gcf,'Position',[2500 400 1600 1200])


subplot(1,2,2)

h=shadedErrorBar(1:5,cohen_BMI_wp(2,:),cohen_BMI_wp_CI,'lineProps',{'.b-','MarkerSize',32},'patchSaturation',0.05);
hold on
plot(1:5,cohen_BMI_wp(1,:),'r.-',"MarkerSize",32)
set(gca, 'TickLength',[0 0])
% yline([ 0.15],'-','0.15')
ylim([0 0.7])
xlim([0.5 5.5])
set(gca, 'XTick', 1:5); % center x-axis ticks on bins
set(gca, 'XTickLabel',data_points  ,'FontSize',24,'FontWeight','Bold'); % set y-axis labels
ylabel('Cohen''s d'  ,'FontSize',28,'FontWeight','Bold')
legend('Regression w/ 95%CI','PSM','FontSize',24,'FontWeight','Normal')
title('BMI after controlling for pubertal status','FontSize',24)

legend('location','best')
set(gcf,'Position',[2500 400 2800 1200])

print(['/home/yangf7/Documents/Nils/ABCD/race/','BMI_wPub'],'-depsc')
% ['/home/yangf7/Documents/Nils/ABCD/race/','BMI_wPub']
close



figure;
subplot(1,2,1)
h=shadedErrorBar(1:5,cohen_Waist(2,:),cohen_Waist_CI,'lineProps','ob-');
hold on
plot(1:5,cohen_Waist(1,:),'rs-',"MarkerSize",10)
set(gca, 'TickLength',[0 0])
% yline([ 0.15],'-','0.15')
ylim([-0.16 0.7])
xlim([0.5 5.5])
set(gca, 'XTick', 1:5); % center x-axis ticks on bins
set(gca, 'XTickLabel',data_points  ,'FontSize',18); % set y-axis labels

legend('Regression','PSM','FontSize',24)
title('Waist Size','FontSize',24)

legend('location','best')
% set(gcf,'Position',[2500 400 1600 1200])


subplot(1,2,2)

h=shadedErrorBar(1:5,cohen_Waist_wp(2,:),cohen_Waist_wp_CI,'lineProps','ob-');
hold on
plot(1:5,cohen_Waist_wp(1,:),'rs-',"MarkerSize",10)
set(gca, 'TickLength',[0 0])
% yline([ 0],'-','0')
ylim([-0.16 0.7])
xlim([0.5 5.5])
set(gca, 'XTick', 1:5); % center x-axis ticks on bins
set(gca, 'XTickLabel',data_points  ,'FontSize',18); % set y-axis labels

legend('Regression','PSM','FontSize',24)
title('Waist Size after controlling for pubertal status','FontSize',24)

legend('location','best')
set(gcf,'Position',[2500 400 2800 1200])
print(['/home/yangf7/Documents/Nils/ABCD/race/','Waist_wPub'],'-djpeg')
close




%% imaging BL
sub_id_img=cell(5);ind_img=cell(5);race_img=cell(5);network_conn_table=cell(5);
for ii=1:2:5
    sub_id_img{ii}=sub_id{ii}(ismember(sub_id{ii},network_table{ii}.Properties.RowNames));
    sub_id_img{ii}=sub_id_img{ii}(ismember(sub_id_img{ii},psm_table_nonan.Properties.RowNames));
    sub_id_img{ii}=sub_id_img{ii}(ismember(sub_id_img{ii},qc_inclusion_table{ii}.Properties.RowNames));
    sub_id_img{ii}=sub_id_img{ii}(ismember(sub_id_img{ii},qc_motion{ii}.Properties.RowNames));
    sub_id_img{ii}=sub_id_img{ii}(qc_inclusion_table{ii}{sub_id_img{ii},'imgincl_rsfmri_include'}==1);%passed qc
    ind_img{ii}=ismember(psm_race0,sub_id_img{ii})&ismember(psm_race1,sub_id_img{ii});%matched subject that both have img data
    race_img{ii}=[psm_race0(ind_img{ii});psm_race1(ind_img{ii})];
    network_conn_table{ii}=[network_table{ii}(sub_id_img{ii},3:171) sub_network_table{ii}(sub_id_img{ii},3:249)];
    % network_conn_table{ii}=[network_table{ii}(sub_id_img{ii},3:171)];


end


% network_conn_table=([network_table_BL(sub_id_BLimg,3:171) sub_network_table_BL(sub_id_BLimg,3:249)]);



p_stats=nan(2,size(network_table_BL,2));
cohend_img{ii}=nan(2,size(network_table_BL,2));
p_value=nan(2,size(network_table_BL,2));



for ii=1:length(var_table_cell)
cohend_img{ii}=nan(2,size(network_table{ii},2));


    for jj=1:size(network_conn_table{ii},2)
        tempv=[];resid=[];p=[];d=[];
        tempv=network_conn_table{ii}{race_img{ii},jj};
        ind_nn=~isnan(tempv);
        b=regress(tempv(ind_nn),[ones(length(tempv(ind_nn)),1) qc_motion{ii}{race_img{ii}(ind_nn),'rsfmri_ntpoints'}...
            ,qc_motion{ii}{race_img{ii}(ind_nn),'rsfmri_meanmotion'}]);
        Yhat=[ones(length(tempv(ind_nn)),1) qc_motion{ii}{race_img{ii}(ind_nn),'rsfmri_ntpoints'}...
            ,qc_motion{ii}{race_img{ii}(ind_nn),'rsfmri_meanmotion'}]*b;
        resid(ind_nn)=tempv(ind_nn)-Yhat+mean(tempv(ind_nn));
        resid(~ind_nn)=nan;
        
    
    
    
    
    
    
    
    
    
        % [h,p]=ttest2(resid(ismember(race_BLimg,psm_race1_BLimg)),resid(ismember(race_BLimg,psm_race0_BLimg)));
        d=computeCohen_d(resid(ismember(race_img{ii},psm_race1)),resid(ismember(race_img{ii},psm_race0)),'independent');
    
    
    
    
        
        % p_stats(1,i)=-log10(p)*sign(d);
        % p_value(1,i)=p;
        cohend_img{ii}(1,jj)=d;
    
        tempv=[];resid=[];p=[];d=[];
        tempv=network_conn_table{ii}{sub_id_img{ii},jj};
        ind_nn=~isnan(tempv);
        cov_temp=psm_table_nonan{sub_id_img{ii}(ind_nn),[2,3,5:7 10:11]}-mean(psm_table_nonan{sub_id_img{ii}(ind_nn),[2,3,5:7 10:11]},'omitnan');
    
        b=regress(tempv(ind_nn),[ones(length(tempv(ind_nn)),1) qc_motion{ii}{sub_id_img{ii}(ind_nn),'rsfmri_ntpoints'}...
            ,qc_motion{ii}{sub_id_img{ii}(ind_nn),'rsfmri_meanmotion'},cov_temp]);
        Yhat=[ones(length(tempv(ind_nn)),1) qc_motion{ii}{sub_id_img{ii}(ind_nn),'rsfmri_ntpoints'}...
            ,qc_motion{ii}{sub_id_img{ii}(ind_nn),'rsfmri_meanmotion'},cov_temp]*b;
        resid(ind_nn)=tempv(ind_nn)-Yhat+mean(tempv(ind_nn));
        resid(~ind_nn)=nan;
    
    
        % [h,p]=ttest2(resid(psm_table_nonan{sub_id_BLimg,'race'}==2),resid(psm_table_nonan{sub_id_BLimg,'race'}==1));
        d=computeCohen_d(resid(ismember(sub_id_img{ii},race_org1)),resid(ismember(sub_id_img{ii},race_org0)),'independent');
    
    
            
        % p_stats(2,i)=-log10(p)*sign(d);
        % p_value(2,i)=p;
         cohend_img{ii}(2,jj)=d;
    
    
    
    
    
    
    
    end
end
% 
% include_indx=abs(cohend_BL(2,:))>=0.1711|abs(cohend_BL(1,:))>=0.1711;% 0.8 statistical power

% indx=1:size(cohend_BL,2);
% indx=[indx;cohend_BL];
% [B,I2]=sort(indx(3,:));
% [~,I2]=sort(indx(4,:));


% figure;
% plot(cohend_BL(1,I2)',1:size(cohend_BL,2),'o','color','r','DisplayName','PSM-matched')
% xlim([-1,1])
% hold on 
% plot(cohend_BL(2,I2)',1:size(cohend_BL,2),'o','color','b','DisplayName','Covariates-corrected')
% hold on 
% % plot(cohend_BL(3,I2)',1:size(cohend_BL,2),'o','color','magenta','DisplayName','No-control')
% set(gca, 'YTick', 1:size(cohend_BL,2)); % center y-axis ticks on bins
% set(gca, 'TickLength',[0 0])
% set(gca, 'YTickLabel',strrep( beh_BL.Properties.VariableNames(I2),'_','\_')              ,'FontSize',18); % set y-axis labels
% xline([ 0 ],'--')
% xline([-0.1303  0.1303],'-')
% 
% yline([1:size(cohend_BL,2)],'-.','Alpha',0.2)
% f=get(gca,'Children');
% legend([f(end-1),f(end)])
% ylim([0,size(cohend_BL,2)])
% 
% legend('location','best')
% set(gcf,'Position',[1500 100 1000 800])




% 
% cov_med_img=[qc_motion_BL{race_BLimg,'rsfmri_ntpoints'}...
%         ,qc_motion_BL{race_BLimg,'rsfmri_meanmotion'},psm_table_nonan{race_BLimg,[2,3,5:7 10:11]}]...
%         -mean([qc_motion_BL{race_BLimg,'rsfmri_ntpoints'}...
%         ,qc_motion_BL{race_BLimg,'rsfmri_meanmotion'},psm_table_nonan{race_BLimg,[2,3,5:7 10:11]}],'omitnan');

p_med_img=cell(5,1);
beta_med_img=cell(5,1);
% stand_ab=cell(5,1);
for ii=1:2:length(var_table_cell)

        mediator=var_table_cell{ii}{sub_id_img{ii},'pubertal_status'};
        % mediator=var_table_cell{ii}{race_BLimg,'hormone_scr_dhea_mean'};

        cov_med_img=[qc_motion{ii}{sub_id_img{ii},'rsfmri_ntpoints'}...
        ,qc_motion{ii}{sub_id_img{ii},'rsfmri_meanmotion'},psm_table_nonan{sub_id_img{ii},[2,3,5:7 10:11]}]...
        -mean([qc_motion{ii}{sub_id_img{ii},'rsfmri_ntpoints'}...
        ,qc_motion{ii}{sub_id_img{ii},'rsfmri_meanmotion'},psm_table_nonan{sub_id_img{ii},[2,3,5:7 10:11]}],'omitnan');
        dep_med=network_conn_table{ii}{sub_id_img{ii},:};
    for jj=1:size(dep_med,2)
        [paths, stats_med]=mediation(psm_table_nonan{sub_id_img{ii},'race'},dep_med(:,jj),...
        mediator...
        ,'cov', cov_med_img,'verbose','boot','bootsamples',10000);
        p_med_img{ii}(jj,:)=stats_med.p;
        beta_med_img{ii}(jj,:)=stats_med.mean;
        % stand_ab{ii}(jj)=stats_med.mean(5).*(std(psm_table_nonan{race_BLimg,'race'},'omitnan'))./std(dep_med(:,jj),'omitnan');

    end




end



% percentofmed_img=[];
for ii=1:2:length(var_table_cell)

    percentofmed_img(ii)=sum(p_med_img{ii}(:,4)<=0.05&p_med_img{ii}(:,5)<=0.05)./sum(p_med_img{ii}(:,4)<=0.05);
end




p_med_img_psm=cell(5,1);
beta_med_img=cell(5,1);
% stand_ab=cell(5,1);
for ii=1:2:length(var_table_cell)

        mediator=var_table_cell{ii}{race_img{ii},'pubertal_status'};
        % mediator=var_table_cell{ii}{race_BLimg,'hormone_scr_dhea_mean'};

        cov_med_img=[qc_motion{ii}{race_img{ii},'rsfmri_ntpoints'}...
        ,qc_motion{ii}{race_img{ii},'rsfmri_meanmotion'},psm_table_nonan{race_img{ii},[2,3,5:7 10:11]}]...
        -mean([qc_motion{ii}{race_img{ii},'rsfmri_ntpoints'}...
        ,qc_motion{ii}{race_img{ii},'rsfmri_meanmotion'},psm_table_nonan{race_img{ii},[2,3,5:7 10:11]}],'omitnan');
        dep_med=network_conn_table{ii}{race_img{ii},:};
    for jj=1:size(dep_med,2)
        [paths, stats_med]=mediation(psm_table_nonan{race_img{ii},'race'},dep_med(:,jj),...
        mediator...
        ,'cov', cov_med_img,'verbose','boot','bootsamples',10000);
        p_med_img_psm{ii}(jj,:)=stats_med.p;
        beta_med_img{ii}(jj,:)=stats_med.mean;
        % stand_ab{ii}(jj)=stats_med.mean(5).*(std(psm_table_nonan{race_BLimg,'race'},'omitnan'))./std(dep_med(:,jj),'omitnan');

    end




end



% percentofmed_img=[];
for ii=1:2:length(var_table_cell)

    percentofmed_img_psm(ii)=sum(p_med_img_psm{ii}(:,4)<=0.05&p_med_img_psm{ii}(:,5)<=0.05)./sum(p_med_img_psm{ii}(:,4)<=0.05);
end



% p_stats_shaped=[reshape(p_stats(1:169),13,13);reshape(p_stats(170:end),19,13);];
% cohen_d_shaped=[reshape(cohen_d_stats(1:169),13,13);reshape(cohen_d_stats(170:end),19,13);];
% cohen_d_shaped=cohen_d_shaped([1:6 8:end],[1:6 8:end]);
% p_value=[reshape(p_value(1:169),13,13);reshape(p_value(170:end),19,13);];
% p_value=p_value([1:6 8:end],[1:6 8:end]);
% p_value=[tril(p_value(1:12,1:12));p_value(13:end,:)];
% p_test=p_value;
% p_test(p_value>0)=fdr_bh(p_value(p_value>0),0.05);
% % cohen_d_shaped= cohen_d_shaped.*p_test;
% cohen_d_shaped(abs(cohen_d_shaped)<0.15)=0;
% cohen_d_shaped=[tril(cohen_d_shaped(1:12,1:12));cohen_d_shaped(13:end,:)];
% 
% net_names={'auditory network', 'cingulo-opercular network', 'cingulo-parietal network', 'default network'...
%     ,'dorsal attention network', 'fronto-parietal network','retrosplenial temporal network',...
%     'sensorimotor hand network','sensorimotor mouth network', 'salience network',...
%     'ventral attention network','visual network','left-cerebellum-cortex','left-thalamus-proper',...
%     'left-caudate',' left-putamen','left-pallidum','brain-stem','left-hippocampus',...
%     'left-amygdala','left-accumbens-area','left-ventraldc','right-cerebellum-cortex',...
%     'right-thalamus-proper','right-caudate','right-putamen','right-pallidum', 'right-hippocampus'...
%     'right-amygdala','right-accumbens-area','right-ventraldc'};
% 
% % create_matrix_ABCD(p_stats_shaped,'Total Sleep Duration FL2',net_names,0)
% 
% create_matrix_ABCD_cohensD(cohen_d_shaped,'Total Sleep Duration Cohen''s D FL2',net_names,0) 
% 
% create_matrix_ABCD_cohensD(FL_net.cohen_d_threhold,'SJLsc',net_names,0) 
% 
% FL_net.p_value=p_value;
% FL_net.cohen_d_stats=cohen_d_stats;
% FL_net.cohen_d_threhold=cohen_d_shaped;


%% ind t test on greymatter volume with TIV regressed out

% p_stats=nan(1,151);
% cohen_d_stats=nan(1,151);
% p_values=nan(1,151);
sub_id_t1img=cell(5);race_t1img=cell(5);
for ii=1:2:5


    sub_id_t1img{ii}=sub_id_img{ii}(qc_inclusion_table{ii}{sub_id_img{ii},'imgincl_t1w_include'}==1);%passed t1 qc
    ind_BL=ismember(psm_race0,sub_id_t1img{ii})&ismember(psm_race1,sub_id_t1img{ii});%matched subject that both have BL img data
    
    race_t1img{ii}=[psm_race0(ind_BL);psm_race1(ind_BL)];

end
% T1_sub_FL2
% sleep_temp_FL2=sjl_psm([ind;ind]);
% ind_FL2_gmv=ismember(sleep_NE,sleep_temp_FL2)&ismember(sleep_FS,sleep_temp_FL2);%matched subject that both have FL2 data
% sleep_temp_FL2=[sleep_FS(ind_FL2_gmv);sleep_NE(ind_FL2_gmv)];
% ind=qc_table_FL2{sleep_FS(ind_FL2_gmv),'imgincl_t1w_include'}==1&qc_table_FL2{sleep_NE(ind_FL2_gmv),'imgincl_t1w_include'}==1;%matched subject that both passed qc

% for i=1:151
%     tempv=T1_table_FL2{sleep_temp_FL2,i+462};
%     ind_nn=~isnan(tempv);
%     b=regress(tempv(ind_nn),[ones(length(tempv(ind_nn)),1) T1_table_FL2{sleep_temp_FL2(ind_nn),'mrisdp_604'}]);
%     Yhat=[ones(length(tempv(ind_nn)),1) T1_table_FL2{sleep_temp_FL2(ind_nn),'mrisdp_604'}]*b;
%     resid=nan(length(tempv),1);
%     resid(ind_nn)=tempv(ind_nn)-Yhat+mean(tempv(ind_nn));
% 
% 
% %     resid(~ind_nn)=nan;
%     [h,p]=ttest2(resid(ismember(sleep_temp_FL2,psm_race0)),resid(ismember(sleep_temp_FL2,psm_race1)));
%     d=computeCohen_d(resid(ismember(sleep_temp_FL2,psm_race0)),resid(ismember(sleep_temp_FL2,psm_race1)),'independent');
% 
%     p_stats(i)=-log10(p)*sign(d);
%     p_values(i)=p;
%      cohen_d_stats(i)=d;
% end

p_med_t1=cell(5,1);
beta_med_t1=cell(5,1);
for ii=1:2:5

cov_med_img=[T1_table{ii}{sub_id_t1img{ii},'mrisdp_604'}, qc_motion{ii}{sub_id_t1img{ii},'rsfmri_ntpoints'}...
        ,qc_motion{ii}{sub_id_t1img{ii},'rsfmri_meanmotion'},psm_table_nonan{sub_id_t1img{ii},[2,3,5:7 10:11]}]...
        -mean([T1_table{ii}{sub_id_t1img{ii},'mrisdp_604'}, qc_motion{ii}{sub_id_t1img{ii},'rsfmri_ntpoints'}...
        ,qc_motion{ii}{sub_id_t1img{ii},'rsfmri_meanmotion'},psm_table_nonan{sub_id_t1img{ii},[2,3,5:7 10:11]}],'omitnan');

stand_ab=cell(5,1);


        mediator=var_table_cell{ii}{sub_id_t1img{ii},'pubertal_status'};
        % mediator=var_table_cell{ii}{sub_id_t1img{ii},'hormone_scr_dhea_mean'};

        dep_med=[T1_table{ii}{sub_id_t1img{ii},3:150},T1_sub_table{ii}{sub_id_t1img{ii},3:34}];
    for jj=1:size(dep_med,2)
        [paths, stats_med]=mediation(psm_table_nonan{sub_id_t1img{ii},'race'},dep_med(:,jj),...
        mediator...
        ,'cov', cov_med_img,'verbose','boot','bootsamples',10000);
        p_med_t1{ii}(jj,:)=stats_med.p;
        beta_med_t1{ii}(jj,:)=stats_med.mean;
        % stand_ab{ii}(jj)=stats_med.mean(5).*(std(psm_table_nonan{sub_id_BLimg_str,'race'},'omitnan'))./std(dep_med(:,jj),'omitnan');

    end




end


for ii=1:2:length(var_table_cell)

    percentofmed_t1(ii)=sum(p_med_t1{ii}(:,4)<=0.05&p_med_t1{ii}(:,5)<=0.05)./sum(p_med_t1{ii}(:,4)<=0.05);
end

p_med_t1_psm=cell(5,1);
beta_med_t1=cell(5,1);
for ii=1:2:5

cov_med_img=[T1_table{ii}{race_t1img{ii},'mrisdp_604'}, qc_motion{ii}{race_t1img{ii},'rsfmri_ntpoints'}...
        ,qc_motion{ii}{race_t1img{ii},'rsfmri_meanmotion'},psm_table_nonan{race_t1img{ii},[2,3,5:7 10:11]}]...
        -mean([T1_table{ii}{race_t1img{ii},'mrisdp_604'}, qc_motion{ii}{race_t1img{ii},'rsfmri_ntpoints'}...
        ,qc_motion{ii}{race_t1img{ii},'rsfmri_meanmotion'},psm_table_nonan{race_t1img{ii},[2,3,5:7 10:11]}],'omitnan');

stand_ab=cell(5,1);


        mediator=var_table_cell{ii}{race_t1img{ii},'pubertal_status'};
        % mediator=var_table_cell{ii}{race_t1img{ii},'hormone_scr_dhea_mean'};

        dep_med=[T1_table{ii}{race_t1img{ii},3:150},T1_sub_table{ii}{race_t1img{ii},3:34}];
    for jj=1:size(dep_med,2)
        [paths, stats_med]=mediation(psm_table_nonan{race_t1img{ii},'race'},dep_med(:,jj),...
        mediator...
        ,'cov', cov_med_img,'verbose','boot','bootsamples',10000);
        p_med_t1_psm{ii}(jj,:)=stats_med.p;
        beta_med_t1{ii}(jj,:)=stats_med.mean;
        % stand_ab{ii}(jj)=stats_med.mean(5).*(std(psm_table_nonan{sub_id_BLimg_str,'race'},'omitnan'))./std(dep_med(:,jj),'omitnan');

    end




end


for ii=1:2:length(var_table_cell)

    percentofmed_t1psm(ii)=sum(p_med_t1_psm{ii}(:,4)<=0.05&p_med_t1_psm{ii}(:,5)<=0.05)./sum(p_med_t1_psm{ii}(:,4)<=0.05);
end

%% pie chart


figure

pie_beh=[sum(p_med{1}(:,4)>=0.05),sum(p_med{1}(:,4)<0.05)-sum(p_med{1}(:,4)<0.05&p_med{1}(:,5)<0.05),...
    sum(p_med{1}(:,4)<0.05&p_med{1}(:,5)<0.05)];
pie_fc=[sum(p_med_img{1}(:,4)>=0.05),sum(p_med_img{1}(:,4)<0.05)-sum(p_med_img{1}(:,4)<0.05&p_med_img{1}(:,5)<0.05),...
    sum(p_med_img{1}(:,4)<0.05&p_med_img{1}(:,5)<0.05)];

pie_t1=[sum(p_med_t1{1}(:,4)>=0.05),sum(p_med_t1{1}(:,4)<0.05)-sum(p_med_t1{1}(:,4)<0.05&p_med_t1{1}(:,5)<0.05),...
    sum(p_med_t1{1}(:,4)<0.05&p_med_t1{1}(:,5)<0.05)];
labels={'c >= 0.05','c < 0.05 & a*b >= 0.05','c < 0.05 & a*b < 0.05'};




t = tiledlayout(1,3,'TileSpacing','compact');

% Create pie charts
ax1 = nexttile;
pie(ax1,pie_beh)
title('Behavior')

ax2 = nexttile;
pie(ax2,pie_fc)
title('Functional connectivity')

ax3 = nexttile;
pie(ax3,pie_t1)
title('GMV')

% Create legend
lgd = legend(labels,'FontSize',18);
lgd.Layout.Tile = 'east';
set(gcf,'Position',[1600 400 1600 1200])
