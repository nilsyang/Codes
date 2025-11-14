% this is code used in paper https://doi.org/10.1002/hbm.25772
%% initial parameters

father_dir='/home/zwang/Documents/Nils/ABCD';
cd(father_dir)
beh_dir='/home/zwang/Documents/Nils/ABCD/Beh_tabulated_data/';

network_table=readtable('/Beh_tabulated_data/abcd_betnet02.txt');% network connectivity 
network_table_FL2=network_table(strcmp(network_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
network_table_BL=network_table(strcmp(network_table{:,'eventname'},'baseline_year_1_arm_1'),:);
network_table_BL.Properties.RowNames=network_table_BL.src_subject_id;
network_table_FL2.Properties.RowNames=network_table_FL2.src_subject_id;
qc_table=readtable(fullfile(beh_dir, 'abcd_imgincl01.txt'));% quality control of mri data
qc_table_BL=qc_table(strcmp(qc_table{:,'eventname'},'baseline_year_1_arm_1'),:);
qc_table_FL2=qc_table(strcmp(qc_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
qc_table_BL.Properties.RowNames=qc_table_BL.src_subject_id;
qc_table_FL2.Properties.RowNames=qc_table_FL2.src_subject_id;
sub_id_BL=network_table_BL.Properties.RowNames;% subjects' id in baseline


sub_id_BL=sub_id_BL(ismember(sub_id_BL,qc_table_BL.Properties.RowNames));% all subs that have qc values

sub_id_BL=sub_id_BL(qc_table_BL{sub_id_BL,'imgincl_rsfmri_include'}==1);% all subs that passed rsfmri qc


sub_network_table=readtable('/Beh_tabulated_data/mrirscor02.txt');% sub cortical network connecitivity
sub_network_table_FL2=sub_network_table(strcmp(sub_network_table{:,'eventname'},'2_year_follow_up_y_arm_1'),:);
sub_network_table_BL=sub_network_table(strcmp(sub_network_table{:,'eventname'},'baseline_year_1_arm_1'),:);
sub_network_table_BL.Properties.RowNames=sub_network_table_BL.src_subject_id;
sub_network_table_FL2.Properties.RowNames=sub_network_table_FL2.src_subject_id;

cbcl_sum=readtable('/Beh_tabulated_data/abcd_cbcls01.txt');% child behavioral checklist
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

total_sleep_duration=sleep_parent_SDS_BL{sub_id_BL,'sleepdisturb1_p'};%1 = 9-11 hours/ 9 a 11 horas; 2 = 8-9 hours /8 a 9 horas; 3 = 7-8 hours /7 a 8 horas; 4 = 5-7 hours /5 a 7 horas; 5 = Less than 5 hours/ Menos de 5 horas// Consider each question pertaining to the PAST 6 MONTHS of the child's life 
total_sleep_duration(total_sleep_duration==5)=4; % combine 4 & 5 as only 26 subjects have a value of 5
total_sleep_duration=5-total_sleep_duration; %reverse sleep duration to 4 are recommed sleep time.



dem_table=readtable(fullfile(beh_dir,'ABCD_sites_effect.xlsx'));% this excel file is downloaded from deap website : https://deap.nimhda.org/
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
dem_table.race4=dem_table.white+2*dem_table.black+3*dem_table.asian+4*dem_table.mixed;% recoded race
dem_table.race4(dem_table.race4==0)=4;
pdem_table=readtable(fullfile(beh_dir,'pdem02.txt'));
pdem_table.Properties.RowNames=pdem_table.src_subject_id;

dem_table.edu1=double(strcmp(dem_table.high_educ,'< HS Diploma'));% recoded parent educational level
dem_table.edu2=double(strcmp(dem_table.high_educ,'HS Diploma/GED'));
dem_table.edu3=double(strcmp(dem_table.high_educ,'Some College'));
dem_table.edu4=double(strcmp(dem_table.high_educ,'Bachelor'));
dem_table.edu5=double(strcmp(dem_table.high_educ,'Post Graduate Degree'));
dem_table.edu=dem_table.edu1+2*dem_table.edu2+3*dem_table.edu3+4*dem_table.edu4+5*dem_table.edu5;


sub_id_BL=sub_id_BL(sum(isnan(network_table_BL{sub_id_BL,23:191}'))==0); %exculded 14 subs that has missing values on network connectivity
total_sleep_disturbance=sleep_parent_SDSC_BL{sub_id_BL,'sds_p_ss_total'};% total sleep disturbance score
tot_pro=cbcl_sum_BL{sub_id_BL,'cbcl_scr_syn_totprob_t'};% total behaviral problems from cbcl aka mental problems (mp) in the paper


sub_id_BL=sub_id_BL(~isnan(total_sleep_disturbance)&~isnan(tot_pro)); % excluded 23 subs that have missing values on sleep disturbance and tot_pro



% combat
net_BL_combat=combat([network_table_BL{sub_id_BL,23:191} sub_network_table_BL{sub_id_BL,23:269}]',dem_table{sub_id_BL,'abcd_site'}',...
[dem_table{sub_id_BL,[16 20]},network_table_BL{sub_id_BL,'interview_age'}],1);% running combat to control the effects of those variables 

pubertal_f=sleep_parent_SDSC_BL{sub_id_BL,'pds_p_ss_female_category'};% puberty statues of each child
pubertal_m=sleep_parent_SDSC_BL{sub_id_BL,'pds_p_ss_male_category'};
pubertal=sum([pubertal_f,pubertal_m],2,'omitnan');
pubertal(pubertal==0)=nan;
total_comp_fc=NIH_TB_BL{sub_id_BL,'nihtbx_totalcomp_fc'};% total cognition score
total_sleep_disturbance=sleep_parent_SDSC_BL{sub_id_BL,'sds_p_ss_total'};% total sleep disturbance score
tot_pro=cbcl_sum_BL{sub_id_BL,'cbcl_scr_syn_totprob_t'};% total behaviral problems from cbcl aka mental problems (mp) in the paper


summary_table_SDS=[network_table_BL(sub_id_BL,'interview_age'),...
        dem_table(sub_id_BL,[16 20]),dem_table(sub_id_BL,[9 15])];
summary_table_SDS.PubertalStatus=pubertal;
summary_table_SDS.TotalProblems=tot_pro;
summary_table_SDS.TotalSleepDisturbance=total_sleep_disturbance;
summary_table_SDS.TotalCognitiveCompositeScores=total_comp_fc;


%% netowork connecitiviy analyses (combined left and right subcortical regions into one region)
total_sleep_disturbance=sleep_parent_SDSC_BL{sub_id_BL,'sds_p_ss_total'};
total_sleep_duration=sleep_parent_SDS_BL{sub_id_BL,'sleepdisturb1_p'};%1 = 9-11 hours/ 9 a 11 horas; 2 = 8-9 hours /8 a 9 horas; 3 = 7-8 hours /7 a 8 horas; 4 = 5-7 hours /5 a 7 horas; 5 = Less than 5 hours/ Menos de 5 horas// Consider each question pertaining to the PAST 6 MONTHS of the child's life 
total_sleep_duration(total_sleep_duration==5)=4; % combine 4 & 5 as only 26 subjects have a value of 5
total_sleep_duration=5-total_sleep_duration; %reverse sleep duration to 4 are recommed sleep time.

House_in=pdem_table{sub_id_BL,'demo_comb_income_v2'};% household income data
House_in(House_in>100)=nan;






net_BL_combat_combined=net_BL_combat(1:169,:);%(combined left and right subcortical regions into one region)
net_BL_sub=net_BL_combat(170:end,:);
net_BL_sub=reshape(net_BL_sub,13,19,9350);
net_BL_sub=net_BL_sub(:,[1:5,7:19],:);
net_BL_sub1=(net_BL_sub(:,1:9,:)+net_BL_sub(:,10:18,:))/2;
net_BL_combat_combined=[net_BL_combat_combined;reshape(net_BL_sub1,13*9,9350)];

p_stats=nan(1,286);
p_value=nan(1,286);
t_stats=nan(1,286);
r2_stats=nan(1,286);
 data_test=table;
    data_test.TSD=log10(total_sleep_disturbance);
    data_test.FC=net_BL_combat_combined(1,:)';
    data_test.MP=tot_pro;
    data_test.HI=House_in;
    data_test.PU=pubertal;
    data_test=[data_test,network_table_BL(sub_id_BL,'interview_age'),...
        dem_table(sub_id_BL,[16 20]),dem_table(sub_id_BL,[9 15])...
        ,dem_table(sub_id_BL,'edu'),network_table_BL(sub_id_BL,'rsfmri_c_ngd_ntpoints'),...
        network_table_BL(sub_id_BL,'rsfmri_c_ngd_meanmotion')...
        BMI_BL(sub_id_BL,'anthro_bmi_calc')];
for i=1:size(net_BL_combat_combined,1)% run linear-mixed effect models
   data_test.FC=net_BL_combat_combined(i,:)';
    lme=fitlme(data_test,'FC~TSD+interview_age+sex+race+rsfmri_c_ngd_meanmotion+PU+(1|abcd_site)+(1|abcd_site:rel_family_id)');
    [~,~,stats]=fixedEffects(lme);
    r2=lme.Rsquared;
    p_stats(i)=-log10(stats.pValue(2))*sign(stats.tStat(2));
    p_value(i)=stats.pValue(2);
    t_stats(i)=stats.tStat(2);
    r2_stats(i)=r2.Adjusted;
end

p_value_shaped=[reshape(p_value(1:169),13,13);reshape(p_value(170:end),9,13);];
p_value_shaped=[tril(p_value_shaped(1:13,1:13));p_value_shaped(14:end,:)];
sum(fdr_bh(p_value_shaped(p_value_shaped>0),0.05))
t_stats_reshaped=[reshape(t_stats(1:169),13,13);reshape(t_stats(170:end),9,13);];
t_stats_reshaped=[tril(t_stats_reshaped(1:13,1:13));t_stats_reshaped(14:end,:)];

t_stats_reshaped=t_stats_reshaped([1:6 8:end],[1:6 8:end]);


net_names={'auditory network', 'cingulo-opercular network', 'cingulo-parietal network', 'default network'...
    ,'dorsal attention network', 'fronto-parietal network','retrosplenial temporal network',...
    'sensorimotor hand network','sensorimotor mouth network', 'salience network',...
    'ventral attention network','visual network','cerebellum-cortex','thalamus-proper',...
    'caudate',' putamen','pallidum','hippocampus',...
    'amygdala','accumbens-area','ventraldc'};

create_matrix_ABCD_tstats(t_stats_reshaped,'Total Sleep Disturbance',net_names,0)


% MP
p_stats=nan(1,286);
p_value=nan(1,286);
t_stats=nan(1,286);
r2_stats=nan(1,286);
 data_test=table;
    data_test.TSD=log10(total_sleep_disturbance);
    data_test.FC=net_BL_combat_combined(1,:)';
    data_test.MP=tot_pro;
    data_test.HI=House_in;
    data_test.PU=pubertal;
    data_test=[data_test,network_table_BL(sub_id_BL,'interview_age'),...
        dem_table(sub_id_BL,[16 20]),dem_table(sub_id_BL,[9 15])...
        ,dem_table(sub_id_BL,'edu'),network_table_BL(sub_id_BL,'rsfmri_c_ngd_ntpoints'),...
        network_table_BL(sub_id_BL,'rsfmri_c_ngd_meanmotion')...
        BMI_BL(sub_id_BL,'anthro_bmi_calc')];
for i=1:size(net_BL_combat_combined,1)
   data_test.FC=net_BL_combat_combined(i,:)';
    lme=fitlme(data_test,'FC~MP+interview_age+sex+race+rsfmri_c_ngd_meanmotion+PU+(1|abcd_site)+(1|abcd_site:rel_family_id)');
    [~,~,stats]=fixedEffects(lme);
    r2=lme.Rsquared;
    p_stats(i)=-log10(stats.pValue(2))*sign(stats.tStat(2));
    p_value(i)=stats.pValue(2);
    t_stats(i)=stats.tStat(2);
    r2_stats(i)=r2.Adjusted;
end

p_value_shaped=[reshape(p_value(1:169),13,13);reshape(p_value(170:end),9,13);];
p_value_shaped=[tril(p_value_shaped(1:13,1:13));p_value_shaped(14:end,:)];
sum(fdr_bh(p_value_shaped(p_value_shaped>0),0.05))
t_stats_reshaped=[reshape(t_stats(1:169),13,13);reshape(t_stats(170:end),9,13);];
t_stats_reshaped=[tril(t_stats_reshaped(1:13,1:13));t_stats_reshaped(14:end,:)];

t_stats_reshaped=t_stats_reshaped([1:6 8:end],[1:6 8:end]);
p_test=p_value_shaped;
p_test(p_value_shaped>0)=fdr_bh(p_value_shaped(p_value_shaped>0),0.05);
t_stats_threshold= t_stats_reshaped.*p_test;
p_value_shaped=p_value_shaped([1:6 8:end],[1:6 8:end]);
t_stats_threshold=t_stats_threshold([1:6 8:end],[1:6 8:end]);
net_names={'auditory network', 'cingulo-opercular network', 'cingulo-parietal network', 'default network'...
    ,'dorsal attention network', 'fronto-parietal network','retrosplenial temporal network',...
    'sensorimotor hand network','sensorimotor mouth network', 'salience network',...
    'ventral attention network','visual network','cerebellum-cortex','thalamus-proper',...
    'caudate',' putamen','pallidum','hippocampus',...
    'amygdala','accumbens-area','ventraldc'};

create_matrix_ABCD_tstats(t_stats_reshaped,'Mental Problems',net_names,0)
create_matrix_ABCD_tstats(t_stats_threshold,'Mental Problems',net_names,0)


% full covariate
p_stats=nan(1,286);
p_value=nan(1,286);
t_stats=nan(1,286);
r2_stats=nan(1,286);
 data_test=table;
 
    data_test.TSD=log10(total_sleep_disturbance);
    data_test.FC=net_BL_combat_combined(1,:)';
    data_test.HI=House_in;
    data_test.PU=pubertal;
    data_test=[data_test,network_table_BL(sub_id_BL,'interview_age'),...
        dem_table(sub_id_BL,[16 20]),dem_table(sub_id_BL,[9 15])...
        ,dem_table(sub_id_BL,'edu'),network_table_BL(sub_id_BL,'rsfmri_c_ngd_ntpoints'),...
        network_table_BL(sub_id_BL,'rsfmri_c_ngd_meanmotion')...
        BMI_BL(sub_id_BL,'anthro_bmi_calc')];
for i=1:size(net_BL_combat_combined,1)
   data_test.FC=net_BL_combat_combined(i,:)';
    lme=fitlme(data_test,'FC~TSD+interview_age+sex+race+edu+HI+PU+rsfmri_c_ngd_ntpoints+rsfmri_c_ngd_meanmotion+anthro_bmi_calc+(1|abcd_site)+(1|abcd_site:rel_family_id)');
    [~,~,stats]=fixedEffects(lme);
    r2=lme.Rsquared;
    p_stats(i)=-log10(stats.pValue(2))*sign(stats.tStat(2));
    p_value(i)=stats.pValue(2);
    t_stats(i)=stats.tStat(2);
    r2_stats(i)=r2.Adjusted;
end

p_value_shaped=[reshape(p_value(1:169),13,13);reshape(p_value(170:end),9,13);];
p_value_shaped=[tril(p_value_shaped(1:13,1:13));p_value_shaped(14:end,:)];
sum(fdr_bh(p_value_shaped(p_value_shaped>0),0.05))
t_stats_reshaped=[reshape(t_stats(1:169),13,13);reshape(t_stats(170:end),9,13);];
t_stats_reshaped=[tril(t_stats_reshaped(1:13,1:13));t_stats_reshaped(14:end,:)];

t_stats_reshaped=t_stats_reshaped([1:6 8:end],[1:6 8:end]);


net_names={'auditory network', 'cingulo-opercular network', 'cingulo-parietal network', 'default network'...
    ,'dorsal attention network', 'fronto-parietal network','retrosplenial temporal network',...
    'sensorimotor hand network','sensorimotor mouth network', 'salience network',...
    'ventral attention network','visual network','cerebellum-cortex','thalamus-proper',...
    'caudate',' putamen','pallidum','hippocampus',...
    'amygdala','accumbens-area','ventraldc'};

create_matrix_ABCD_tstats(t_stats_reshaped,'Total Sleep Disturbance full covariate',net_names,1)

p_stats=nan(1,286);
p_value=nan(1,286);
t_stats=nan(1,286);
r2_stats=nan(1,286);
 data_test=table;
    data_test.TSD=log10(total_sleep_disturbance);
    data_test.MP=tot_pro;
    data_test.FC=net_BL_combat_combined(1,:)';
    data_test.HI=House_in;
    data_test.PU=pubertal;
    data_test=[data_test,network_table_BL(sub_id_BL,'interview_age'),...
        dem_table(sub_id_BL,[16 20]),dem_table(sub_id_BL,[9 15])...
        ,dem_table(sub_id_BL,'edu'),network_table_BL(sub_id_BL,'rsfmri_c_ngd_ntpoints'),...
        network_table_BL(sub_id_BL,'rsfmri_c_ngd_meanmotion')...
        BMI_BL(sub_id_BL,'anthro_bmi_calc')];
for i=1:size(net_BL_combat_combined,1)
   data_test.FC=net_BL_combat_combined(i,:)';
    lme=fitlme(data_test,'FC~MP+interview_age+sex+race+edu+HI+PU+rsfmri_c_ngd_ntpoints+rsfmri_c_ngd_meanmotion+anthro_bmi_calc+(1|abcd_site)+(1|abcd_site:rel_family_id)');
    [~,~,stats]=fixedEffects(lme);
    r2=lme.Rsquared;
    p_stats(i)=-log10(stats.pValue(2))*sign(stats.tStat(2));
    p_value(i)=stats.pValue(2);
    t_stats(i)=stats.tStat(2);
    r2_stats(i)=r2.Adjusted;
end
p_value_shaped=[reshape(p_value(1:169),13,13);reshape(p_value(170:end),9,13);];
p_value_shaped=[tril(p_value_shaped(1:13,1:13));p_value_shaped(14:end,:)];
sum(fdr_bh(p_value_shaped(p_value_shaped>0),0.05))
t_stats_reshaped=[reshape(t_stats(1:169),13,13);reshape(t_stats(170:end),9,13);];
t_stats_reshaped=[tril(t_stats_reshaped(1:13,1:13));t_stats_reshaped(14:end,:)];

t_stats_reshaped=t_stats_reshaped([1:6 8:end],[1:6 8:end]);



net_names={'auditory network', 'cingulo-opercular network', 'cingulo-parietal network', 'default network'...
    ,'dorsal attention network', 'fronto-parietal network','retrosplenial temporal network',...
    'sensorimotor hand network','sensorimotor mouth network', 'salience network',...
    'ventral attention network','visual network','cerebellum-cortex','thalamus-proper',...
    'caudate',' putamen','pallidum','hippocampus',...
    'amygdala','accumbens-area','ventraldc'};
create_matrix_ABCD_tstats(t_stats_reshaped,'Mental problems full covariate',net_names,1)






%%  1 year follow up data

sub_id_FL1=sub_id_BL(ismember(sub_id_BL,sleep_parent_SDSC_FL1.Properties.RowNames));

total_sleep_disturbance_FL1=sleep_parent_SDSC_FL1{sub_id_FL1,'sds_p_ss_total'};
tot_pro_FL1=cbcl_sum_FL1{sub_id_FL1,'cbcl_scr_syn_totprob_t'};


covs_FL1=[covs{:,1:5},pubertal];
covs_FL1=covs_FL1(ismember(sub_id_BL,sub_id_FL1));

pubertal_f_FL1=sleep_parent_SDSC_FL1{sub_id_FL1,'pds_p_ss_female_category'};
pubertal_m_FL1=sleep_parent_SDSC_FL1{sub_id_FL1,'pds_p_ss_male_category'};
pubertal_FL1=sum([pubertal_f_FL1,pubertal_m_FL1],2,'omitnan');


% total_comp_fc_FL1=NIH_TB_FL1{sub_id_FL1,'nihtbx_totalcomp_fc'};


summary_table_FL1=[sleep_parent_SDSC_FL1(sub_id_FL1,'interview_age'),...
        dem_table(sub_id_FL1,[16 20]),dem_table(sub_id_FL1,[9 15])];
summary_table_FL1.PubertalStatus=pubertal_FL1;
summary_table_FL1.TotalProblems=tot_pro_FL1;
summary_table_FL1.TotalSleepDisturbance=total_sleep_disturbance_FL1;
summary_table_FL1.TotalCognitiveCompositeScores=total_comp_fc(ismember(sub_id_BL,sub_id_FL1));



%% mediation with mean FD 
% need mediation toolbox
paths=nan(6,5);
[paths(1,:), stats_med]=mediation(total_sleep_disturbance,tot_pro,net_BL_combat_table.rsfmri_c_ngd_dt_ngd_dla...
    ,'cov',[data_test{:,[5:10 13]}],'doCIs','boot','bootsamples',10000) % p=9.5424e-05

[paths(2,:), stats_med]=mediation(total_sleep_disturbance,tot_pro,net_BL_combat_table.rsfmri_c_ngd_dt_ngd_dt...
    ,'cov',[data_test{:,[5:10 13]}],'doCIs','boot','bootsamples',10000)% 0.0051

[paths(3,:), stats_med]=mediation(total_sleep_disturbance,tot_pro,net_BL_combat_table.rsfmri_c_ngd_dla_ngd_dla...
    ,'cov',[data_test{:,[5:10 13]}],'doCIs','boot','bootsamples',10000) %1.9689e-04

[paths(4,:), stats_med]=mediation(tot_pro,total_sleep_disturbance,net_BL_combat_table.rsfmri_c_ngd_dt_ngd_dla...
    ,'cov',[data_test{:,[5:10 13]}],'doCIs','boot','bootsamples',10000) % p=0.1019

[paths(5,:), stats_med]=mediation(tot_pro,total_sleep_disturbance,net_BL_combat_table.rsfmri_c_ngd_dt_ngd_dt...
    ,'cov',[data_test{:,[5:10 13]}],'doCIs','boot','bootsamples',10000)% 0.0131

[paths(6,:), stats_med]=mediation(tot_pro,total_sleep_disturbance,net_BL_combat_table.rsfmri_c_ngd_dla_ngd_dla...
    ,'cov',[data_test{:,[5:10 13]}],'doCIs','boot','bootsamples',10000) %0.1424


paths_FL=nan(6,5);

[paths_FL(1,:), stats_med]=mediation(total_sleep_disturbance(ismember(sub_id_BL,sub_id_FL1))...
    ,tot_pro_FL1,net_BL_combat_table.rsfmri_c_ngd_dt_ngd_dla(ismember(sub_id_BL,sub_id_FL1))...
    ,'cov',[data_test{ismember(sub_id_BL,sub_id_FL1),[5:10 13]},tot_pro(ismember(sub_id_BL,sub_id_FL1))],'doCIs','boot','bootsamples',10000)
%0.0022
[paths_FL(2,:), stats_med]=mediation(total_sleep_disturbance(ismember(sub_id_BL,sub_id_FL1))...
    ,tot_pro_FL1,net_BL_combat_table.rsfmri_c_ngd_dt_ngd_dt(ismember(sub_id_BL,sub_id_FL1))...
    ,'cov',[data_test{ismember(sub_id_BL,sub_id_FL1),[5:10 13]},tot_pro(ismember(sub_id_BL,sub_id_FL1))],'doCIs','boot','bootsamples',10000)
% 0.0150
[paths_FL(3,:), stats_med]=mediation(total_sleep_disturbance(ismember(sub_id_BL,sub_id_FL1))...
    ,tot_pro_FL1,net_BL_combat_table.rsfmri_c_ngd_dla_ngd_dla(ismember(sub_id_BL,sub_id_FL1))...
    ,'cov',[data_test{ismember(sub_id_BL,sub_id_FL1),[5:10 13]},tot_pro(ismember(sub_id_BL,sub_id_FL1))],'doCIs','boot','bootsamples',10000)

%0.0138

[paths_FL(4,:), stats_med]=mediation(tot_pro(ismember(sub_id_BL,sub_id_FL1)),...
    total_sleep_disturbance_FL1...
    ,net_BL_combat_table.rsfmri_c_ngd_dt_ngd_dla(ismember(sub_id_BL,sub_id_FL1))...
    ,'cov',[data_test{ismember(sub_id_BL,sub_id_FL1),[5:10 13]},total_sleep_disturbance(ismember(sub_id_BL,sub_id_FL1))],'doCIs','boot','bootsamples',10000)
%0.0056
[paths_FL(5,:), stats_med]=mediation(tot_pro(ismember(sub_id_BL,sub_id_FL1)),...
    total_sleep_disturbance_FL1...
    ,net_BL_combat_table.rsfmri_c_ngd_dt_ngd_dt(ismember(sub_id_BL,sub_id_FL1))...
    ,'cov',[data_test{ismember(sub_id_BL,sub_id_FL1),[5:10 13]},total_sleep_disturbance(ismember(sub_id_BL,sub_id_FL1))],'doCIs','boot','bootsamples',10000)
%0.0427
[paths_FL(6,:), stats_med]=mediation(tot_pro(ismember(sub_id_BL,sub_id_FL1)),...
    total_sleep_disturbance_FL1...
    ,net_BL_combat_table.rsfmri_c_ngd_dla_ngd_dla(ismember(sub_id_BL,sub_id_FL1))...
    ,'cov',[data_test{ismember(sub_id_BL,sub_id_FL1),[5:10 13]},total_sleep_disturbance(ismember(sub_id_BL,sub_id_FL1))],'doCIs','boot','bootsamples',10000)

%0.0160



% plot 
figure;
bar(paths(1:3,5),'w')
xticklabels({'DMN-DAN','DMN-DMN','DAN-DAN' });
ylabel({'beta of mediation effect (a*b)'},'fontsize',18)
ylim([0 0.004])
set(gca,'TickLength',[0 0])
box off


figure;
bar(paths(4:6,5),'w')
xticklabels({'DMN-DAN','DMN-DMN','DAN-DAN' });
ylabel({'beta of mediation effect (a*b)'},'fontsize',18)
ylim([0 0.002])
set(gca,'TickLength',[0 0])
box off

figure;
bar(paths_FL(1:3,5),'w')
xticklabels({'DMN-DAN','DMN-DMN','DAN-DAN' });
ylabel({'beta of mediation effect (a*b)'},'fontsize',18)
ylim([0 0.004])
set(gca,'TickLength',[0 0])
box off


figure;
bar(paths_FL(4:6,5),'w')
xticklabels({'DMN-DAN','DMN-DMN','DAN-DAN' });
ylabel({'beta of mediation effect (a*b)'},'fontsize',18)
ylim([0 0.002])
set(gca,'TickLength',[0 0])
box off





