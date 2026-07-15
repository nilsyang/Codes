addpath(genpath('/misc/imeel/yangf7/matlab/'))

rmpath(genpath('/misc/imeel/yangf7/matlab/BCT'))

load('/raid/common/sleep1/derivatives/nils/ROI_list.mat');%ROI list
load('/misc/imeel/yangf7/matlab/SW_detection/label_table_0705.mat')% 



ROI_list.Network_Label_v4=nan(300,1);
for ii=1:height(ROI_list)
    if ROI_list.Network_Label(ii)<18
        ROI_list.Network_Label_v4(ii)=ROI_list.Network_Label(ii);
    elseif ROI_list.Network_Label_v2(ii)==18 && ROI_list.Network_Label(ii)==37
        ROI_list.Network_Label_v4(ii)=18;
    elseif ROI_list.Network_Label_v2(ii)==18 && ROI_list.Anatomical_lable(ii)==3
        ROI_list.Network_Label_v4(ii)=19;
    elseif ROI_list.Network_Label_v2(ii)==18 && (ROI_list.Anatomical_lable(ii)==5||ROI_list.Anatomical_lable(ii)==4)
        ROI_list.Network_Label_v4(ii)=20;
    elseif ROI_list.Network_Label_v2(ii)==18 && ROI_list.Anatomical_lable(ii)==6
        ROI_list.Network_Label_v4(ii)=21;
    elseif ROI_list.Network_Label(ii)>=39 &&ROI_list.Network_Label(ii)<60
        ROI_list.Network_Label_v4(ii)=22;
    end

end


ROI_list.Network_Label_v5=ROI_list.Network_Label_v4;
ROI_list.Network_Label_v5(244:245)=27;


Network_Label=[ROI_list.Network_Label_v5;(23:26)';(23:26)'];
net_n_sep_hip={'Un','DMN','VIS','FPN','REW','DAN','VAN','SAL','CON','dSMN','lSMN','AUD','PMN','MTL','aHIP','pHIP','BG','THAL','CB','aHIPi','pHIPi','ECi','paraHIPi'};

% 
% prepro={
%           'procbasic_ric_rvt_ppga' % with RETROICOR 
%         };

prepro={

        'procbasic_ric_gsrwb' % with RETROICOR and whole-brain GS
 };



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

        % if label_table{i,'diff_count'}==1|label_table{i,'EEG_missing'}==1
        if label_table{i,'diff_count'}==1

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
         test_mat{i,2}=readtable(['/raid/common/sleep1/derivatives/nils/dfc_results/HP_individulized_left_', label_table{i,'label_fmri'},'_',char(prepro{j}),'_000.netts'],'FileType','text');
         test_mat{i,2}=table2array(test_mat{i,2}([1 2 5 8],censor_ind));%with censor
         test_mat{i,3}=readtable(['/raid/common/sleep1/derivatives/nils/dfc_results/HP_individulized_right_', label_table{i,'label_fmri'},'_',char(prepro{j}),'_000.netts'],'FileType','text');
         test_mat{i,3}=table2array(test_mat{i,3}([1 2 5 8],censor_ind));%with censor
         test_mat{i,1}=[test_mat{i,1};test_mat{i,2};test_mat{i,3}];



         sleep_score_censor{i}=sleep_score{i}(censor_ind+3);% first 3 volumes were removed


    end
end



test_mat=test_mat(:,1);

sleep_score_censor_subs=[];sleep_score_censor_subs_n1=[];
rawFC_mat=[];rawFC_mat_n1=[];
sleep_score_table={};sleep_score_table_n1={};
for ii=1:length(subs)
    sleep_score_censor_subs{ii}=cat(2,sleep_score_censor{label_table{:,'diff_count'}~=1&label_table{:,'subs'}==ii&label_table{:,'sess'}==2});
% figure;plot( sleep_score_censor_n2_subs{ii})

rawFC_mat{ii}=cat(2,test_mat{label_table{:,'diff_count'}~=1&label_table{:,'subs'}==ii&label_table{:,'sess'}==2});

        i=[find(diff(sleep_score_censor_subs{ii}))];% find when the number (sleep stages) changes
        n=[i numel(sleep_score_censor_subs{ii})]-[0 i];% how long the sleep stages last
        i=[i numel(sleep_score_censor_subs{ii})];
%         startpoint=(i(n>=10)-n(n>=10)+1)';
%         endpoint=i(n>=10)';
        startpoint=(i-n+1)';
        endpoint=i';
sleep_score_table{ii}=table(startpoint,endpoint);
sleep_score_table{ii}.Score=sleep_score_censor_subs{ii}(startpoint)';% stages 5=rem 4=N1;
sleep_score_table{ii}.num_volume=n';
sleep_score_table{ii}=sleep_score_table{ii}(sleep_score_table{ii}.num_volume>=30,:);
sleep_score_table{ii}.subs=ii.*ones(height(sleep_score_table{ii}),1);

end


for ii=1:length(subs)
    sleep_score_censor_subs_n1{ii}=cat(2,sleep_score_censor{label_table{:,'diff_count'}~=1&label_table{:,'subs'}==ii&label_table{:,'sess'}==1});
% figure;plot( sleep_score_censor_n2_subs{ii})

rawFC_mat_n1{ii}=cat(2,test_mat{label_table{:,'diff_count'}~=1&label_table{:,'subs'}==ii&label_table{:,'sess'}==1});

        i=[find(diff(sleep_score_censor_subs_n1{ii}))];% find when the number (sleep stages) changes
        n=[i numel(sleep_score_censor_subs_n1{ii})]-[0 i];% how long the sleep stages last
        i=[i numel(sleep_score_censor_subs_n1{ii})];
%         startpoint=(i(n>=10)-n(n>=10)+1)';
%         endpoint=i(n>=10)';
        startpoint=(i-n+1)';
        endpoint=i';
sleep_score_table_n1{ii}=table(startpoint,endpoint);
sleep_score_table_n1{ii}.Score=sleep_score_censor_subs_n1{ii}(startpoint)';% stages 5=rem 4=N1;
sleep_score_table_n1{ii}.num_volume=n';
sleep_score_table_n1{ii}=sleep_score_table_n1{ii}(sleep_score_table_n1{ii}.num_volume>=30,:);
sleep_score_table_n1{ii}.subs=ii.*ones(height(sleep_score_table_n1{ii}),1);

end





win_len=30;%30TR
step_len=10;
sleep_score_table_dfc={};
for ii=1:length(subs)
sleep_score_table_dfc{ii}=sleep_score_table{ii}(1,:);
n=1;
    for jj=1:height(sleep_score_table{ii})

        FC_steps=1+floor((sleep_score_table{ii}.num_volume(jj)-win_len)./step_len);
        for kk=1:FC_steps
            sleep_score_table_dfc{ii}.startpoint(n)=sleep_score_table{ii}.startpoint(jj)+(kk-1)*step_len;
            sleep_score_table_dfc{ii}.endpoint(n)=sleep_score_table_dfc{ii}.startpoint(n)+win_len-1;
            sleep_score_table_dfc{ii}.num_volume(n)=win_len;
            sleep_score_table_dfc{ii}.Score(n)=sleep_score_table{ii}.Score(jj);
            sleep_score_table_dfc{ii}.subs(n)=sleep_score_table{ii}.subs(jj);
            sleep_score_table_dfc{ii}.nums(n)=n;

            n=n+1;
        end    
    end

end

for  ii=1:length(subs)
    sleep_score_table_dfc{ii}.unique_FC(1)=1;
    temp=sleep_score_table_dfc{ii}.endpoint(1);
      for jj=2:height(sleep_score_table_dfc{ii})
          if sleep_score_table_dfc{ii}.startpoint(jj)<temp
               sleep_score_table_dfc{ii}.unique_FC(jj)=0;
          else
              sleep_score_table_dfc{ii}.unique_FC(jj)=1;
              temp=sleep_score_table_dfc{ii}.endpoint(jj);
          end
      end
end


sleep_score_table_dfc_n1={};
for ii=1:length(subs)
sleep_score_table_dfc_n1{ii}=sleep_score_table_n1{ii}(1,:);
n=1;
    for jj=1:height(sleep_score_table_n1{ii})

        FC_steps=1+floor((sleep_score_table_n1{ii}.num_volume(jj)-win_len)./step_len);
        for kk=1:FC_steps
            sleep_score_table_dfc_n1{ii}.startpoint(n)=sleep_score_table_n1{ii}.startpoint(jj)+(kk-1)*step_len;
            sleep_score_table_dfc_n1{ii}.endpoint(n)=sleep_score_table_dfc_n1{ii}.startpoint(n)+win_len-1;
            sleep_score_table_dfc_n1{ii}.num_volume(n)=win_len;
            sleep_score_table_dfc_n1{ii}.Score(n)=sleep_score_table_n1{ii}.Score(jj);
            sleep_score_table_dfc_n1{ii}.subs(n)=sleep_score_table_n1{ii}.subs(jj);
            sleep_score_table_dfc_n1{ii}.nums(n)=n;

            n=n+1;
        end    
    end

end






subjects = {'00003','00003','00006','00006','00014','00014','00029','00029','00047','00047','00056','00056','00060','00060'...
    ,'00063','00063','00065','00065','00078','00078','00086','00086','00105','00105'};
sessione = {'1','2','1','2','1','2','1','2','1','2','1','2','1','2','1','2','1','2','1','2','1','2','1','2'};




%% dfc
sleep_score_table_dfc_2_5=sleep_score_table_dfc;




dfc_table_256=[];
score_kmeans=[];
for ii=1:length(subs)

    score_kmeans=[score_kmeans;sleep_score_table_dfc_2_5{ii}.Score];
    dfc_table_256=[dfc_table_256;sleep_score_table_dfc_2_5{ii}];
end


%% 302*302
net_code=unique(Network_Label);
ROI_num=1:length(Network_Label);
corr_mat_sep={};corr_mat_sep_all=[];
indx_mat=tril(ones(302,302),-1);
% indx_mat(1:294,:)=0;%only hpc connections
for ii=1:length(sleep_score_table_dfc_2_5)
% corr_mat_sep{ii}=nan(302*301/2,height(sleep_score_table_dfc_2_5{ii}));
corr_mat_sep{ii}=nan(sum(indx_mat==1,'all'),height(sleep_score_table_dfc_2_5{ii}));
for jj=1:height(sleep_score_table_dfc_2_5{ii})
% temp_mat=atanh(partialcorr(rawFC_mat{ii}(:,sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj))',mean(rawFC_mat{ii}(262:273,sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj)))'));
temp_mat=atanh(corr(rawFC_mat{ii}([1:233 240:308],sleep_score_table_dfc_2_5{ii}.startpoint(jj):sleep_score_table_dfc_2_5{ii}.endpoint(jj))'));%remove 234-239 MTL + hippo
temp_mat(logical(eye(size(temp_mat,1))))=4;
corr_mat_sep{ii}(:,jj)=temp_mat(indx_mat>0);
end
corr_mat_sep_all=cat(2,corr_mat_sep_all,corr_mat_sep{ii});
end






%% 302*302 night1
net_code=unique(Network_Label);
ROI_num=1:length(Network_Label);
corr_mat_sep_n1={};corr_mat_sep_all_n1=[];
indx_mat=tril(ones(302,302),-1);
% indx_mat(1:294,:)=0;%only hpc connections
for ii=1:length(sleep_score_table_dfc_n1)
% corr_mat_sep{ii}=nan(302*301/2,height(sleep_score_table_dfc_2_5{ii}));
corr_mat_sep_n1{ii}=nan(sum(indx_mat==1,'all'),height(sleep_score_table_dfc_n1{ii}));
for jj=1:height(sleep_score_table_dfc_n1{ii})
% temp_mat=atanh(partialcorr(rawFC_mat{ii}(:,sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj))',mean(rawFC_mat{ii}(262:273,sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj)))'));
temp_mat=atanh(corr(rawFC_mat_n1{ii}([1:233 240:308],sleep_score_table_dfc_n1{ii}.startpoint(jj):sleep_score_table_dfc_n1{ii}.endpoint(jj))'));%remove 234-239 MTL + hippo
temp_mat(logical(eye(size(temp_mat,1))))=4;
corr_mat_sep_n1{ii}(:,jj)=temp_mat(indx_mat>0);
end
corr_mat_sep_all_n1=cat(2,corr_mat_sep_all,corr_mat_sep_n1{ii});
end



%% identity  matrix
stages=[6,4,3,2,5];
sleep_stages_names={'Wake','N1','N2','N3','REM'};
id_matrix=nan(12,12,5);
sub_id=[1:12];%
id_diag=eye(length(sub_id),length(sub_id));
id_nodiag=~id_diag;
id_self=nan(length(sub_id),length(stages),length(stages));
id_other=nan(length(sub_id),length(stages),length(stages));
id_diff_all=nan(length(sub_id),11,length(stages),length(stages));
id_matrix_stg=[];
for trget=1:5
    id_matrix=nan(12,12,5);

    for i_st=1:length(stages)
        for x_sub=1:length(sub_id)
            for y_sub=1:length(sub_id)
                if i_st ==trget
                    try
                    temp_fc=corr(zscore(corr_mat_sep{sub_id(x_sub)}(:,sleep_score_table_dfc_2_5{sub_id(x_sub)}.Score==stages(i_st)&...
                        sleep_score_table_dfc_2_5{sub_id(x_sub)}.unique_FC==1)),...
                        zscore(corr_mat_sep{sub_id(y_sub)}(:,sleep_score_table_dfc_2_5{sub_id(y_sub)}.Score==stages(trget)&...
                        sleep_score_table_dfc_2_5{sub_id(y_sub)}.unique_FC==1)));
                    id_matrix(x_sub,y_sub,i_st)=mean((temp_fc(temp_fc~=1)));
                                        catch me
                        disp(me)
                    end
                % end
                else
                    try
                        temp_fc=corr(zscore(corr_mat_sep{sub_id(x_sub)}(:,sleep_score_table_dfc_2_5{sub_id(x_sub)}.Score==stages(i_st))),...
                            zscore(corr_mat_sep{sub_id(y_sub)}(:,sleep_score_table_dfc_2_5{sub_id(y_sub)}.Score==stages(trget))));
                        id_matrix(x_sub,y_sub,i_st)=mean((temp_fc(temp_fc~=1)));
                    catch me
                        disp(me)
                    end
                end
            end
        end
        temp_var=id_matrix(:,:,i_st);
       
        id_self(:,trget,i_st)=temp_var(id_diag>0);


        id_other(:,trget,i_st)=mean(reshape(temp_var(id_nodiag>0),11,12),'omitnan');
        id_diff_all(:,:,trget,i_st)=id_self(:,trget,i_st)-reshape(temp_var(id_nodiag>0),11,12)';
    end
     id_matrix_stg{trget}=id_matrix;
    
end


id_diff=id_self-id_other;

 [t,p]=ttest2(reshape(id_self(:,:,1),[],1),reshape(id_self(:,:,3),[],1))
  [t,p]=ttest2(reshape(id_other(:,:,1),[],1),reshape(id_other(:,:,4),[],1))
  [t,p]=ttest2(reshape(id_diff(:,3,1),[],1),reshape(id_diff(:,3,4),[],1))

    [t,p]=ranksum(reshape(id_diff_all(:,:,2:5,1),[],1),reshape(id_diff_all(:,:,[1:3,5],4),[],1))

ranksum(reshape(id_diff_all(:,:,:,1),[],1),reshape(id_diff_all(:,:,:,4),[],1),'tail','left')

ranksum(reshape(id_other(:,:,1),[],1),reshape(id_other(:,:,4),[],1),'tail','left')

figure;
subplot(1,2,1)
h=plot_box_scatter([reshape(id_other(:,:,1),[],1);reshape(id_other(:,:,4),[],1)]...
    ,[ones(60,1);ones(60,1).*2],'color',{'g','b'},'symbol',{'o','s'});
set(gca, 'XTick', 1:2)
set(gca, 'XTickLabel', {'Wake','N3'}); % set x-axis labels
ylabel('Self-Other','FontSize',28,'FontWeight','Bold'); % set x-axis labels
set(gcf,'Position',[800,147,1600,1023])

% set(h, {'DisplayName'}, sleep_stages_names')
% legend()

subplot(1,2,2)
h=plot_box_scatter([reshape(id_diff_all(:,:,:,1),[],1);reshape(id_diff_all(:,:,:,4),[],1)]...
    ,[ones(660,1);ones(660,1).*2],'color',{'g','b'},'symbol',{'o','s'});
set(gca, 'XTick', 1:2)
set(gca, 'XTickLabel', {'Wake','N3'}); % set x-axis labels
ylabel('Self-Self minus Self-Other','FontSize',28,'FontWeight','Bold'); % set x-axis labels
set(gcf,'Position',[800,147,1600,1023])

% % set(h, {'DisplayName'}, sleep_stages_names')
% legend()





figure;
h=bar(1:5,squeeze(nanmean(id_diff,1)));

figure;
h=bar(1:5,squeeze(nanmean(id_self,1)));

set(gca, 'XTick', 1:5)
set(gca, 'XTickLabel', sleep_stages_names); % set x-axis labels
ylabel('self-self minus self-other','FontSize',28,'FontWeight','Bold'); % set x-axis labels
set(gcf,'Position',[800,147,1600,1023])

set(h, {'DisplayName'}, sleep_stages_names')
legend()

figure;
h=bar(1:5,squeeze(nanmean(id_other,1)));

set(gca, 'XTick', 1:5)
set(gca, 'XTickLabel', sleep_stages_names); % set x-axis labels
ylabel('self-self minus self-other','FontSize',28,'FontWeight','Bold'); % set x-axis labels
set(gcf,'Position',[800,147,1600,1023])

set(h, {'DisplayName'}, sleep_stages_names')
legend()




figure;
h=bar(1:5,squeeze(nanmean(id_self,1)));

set(gca, 'XTick', 1:5)
set(gca, 'XTickLabel', sleep_stages_names); % set x-axis labels
ylabel('self-self','FontSize',28,'FontWeight','Bold'); % set x-axis labels
set(gcf,'Position',[800,147,1600,1023])

set(h, {'DisplayName'}, sleep_stages_names')
legend()


figure;
h=bar(1:5,squeeze(nanmean(id_other,1)));

set(gca, 'XTick', 1:5)
set(gca, 'XTickLabel', sleep_stages_names); % set x-axis labels
ylabel('self-other','FontSize',28,'FontWeight','Bold'); % set x-axis labels
set(gcf,'Position',[800,147,1600,1023])

set(h, {'DisplayName'}, sleep_stages_names')
legend()






figure
for jj=[1 4]
% for jj=1:5
for ii=1:length(stages)

        if jj~=1
            subplot(5,5,ii+5*(jj-1))
        else
            subplot(5,5,ii)
        end
        imagesc(id_matrix_stg{jj}(:,:,ii));
        c = colorbar;
        
        c.FontSize=12;
         clim([0.1 0.7])
        colormap(bluewhitered)
        
        c.TickLength=0;
        c.Box='off';
        % xtickangle(90)
        % set(gca, 'YTick', 1:18)
        % set(gca, 'YTickLabel', net_n_sub); % set x-axis labels
        % set(gca, 'XTick', 1:18)
        xlabel(['Target (',sleep_stages_names{ii},')'],'FontSize',24,'FontWeight','Bold'); % set x-axis labels
        ylabel(['Database (',sleep_stages_names{jj},')'],'FontSize',24,'FontWeight','Bold'); % set x-axis labels
        % ylabel('REM','FontSize',28,'FontWeight','Bold'); % set x-axis labels


        set(gca,'box','off') 
        
        % title(sleep_stages_names(ii),'FontSize',28,'FontWeight','Bold')


end
end
set(gcf,'Position',[800,147,2400,1823])

id_rate=nan(5,5);
for itrgt=1:5
    for istg=1:5
        
        diag_greater = false(12,1); 
        for isub=1:12
            row_values = id_matrix_stg{itrgt}(isub, :,istg);  % Extract the row
            diag_value = id_matrix_stg{itrgt}(isub, isub,istg);  % Extract the diagonal element
            row_values(isub) = -Inf;  % Exclude diagonal element from comparison
            if diag_value > max(row_values)
                diag_greater(isub) = true;
            end
        end
        id_rate(itrgt,istg)=sum(diag_greater)/sum(~isnan(diag(id_matrix_stg{itrgt}(:, :,istg))));
    end
end




figure

imagesc(id_rate);
c = colorbar;

c.FontSize=12;
 clim([0.6 1])
colormap(jet)

c.TickLength=0;
c.Box='off';
% xtickangle(90)
set(gca, 'YTick', 1:5)
set(gca, 'YTickLabel', sleep_stages_names,'FontSize',18); % set x-axis labels
set(gca, 'XTick', 1:5)
set(gca, 'XTickLabel', sleep_stages_names,'FontSize',18); % set x-axis labels

xlabel('Testing','FontSize',28,'FontWeight','Bold'); % set x-axis labels
ylabel('Database','FontSize',28,'FontWeight','Bold'); % set x-axis labels
% ylabel('REM','FontSize',28,'FontWeight','Bold'); % set x-axis labels


set(gca,'box','off') 
set(gcf,'Position',[800,147,2000,1823])


figure

imagesc(id_rate([1,4],:));
c = colorbar;

c.FontSize=12;
 clim([0.6 1])
colormap(jet)

c.TickLength=0;
c.Box='off';
% xtickangle(90)
set(gca, 'YTick', 1:2)
set(gca, 'YTickLabel', sleep_stages_names(1,[1 4]),'FontSize',18); % set x-axis labels
set(gca, 'XTick', 1:5)
set(gca, 'XTickLabel', sleep_stages_names,'FontSize',18); % set x-axis labels

xlabel('Target','FontSize',28,'FontWeight','Bold'); % set x-axis labels
ylabel('Database','FontSize',28,'FontWeight','Bold'); % set x-axis labels
% ylabel('REM','FontSize',28,'FontWeight','Bold'); % set x-axis labels


set(gca,'box','off') 
set(gcf,'Position',[800,147,1000,423])




%% manhatton distance
sub_id=1:12;stages=[6,4,3,2,5];fc_stage=[];
for i_st=1:length(stages)
    for x_sub=1:length(sub_id)

        j_ind=sleep_score_table_dfc_2_5{sub_id(x_sub)}.Score==stages(i_st)&...
                        sleep_score_table_dfc_2_5{sub_id(x_sub)}.unique_FC==1;
        fc_stage=corr_mat_sep{sub_id(x_sub)}(:,j_ind)';
        % x=covarianceShrinkage(fc_stage);
        ManhDistance=pdist2(fc_stage,mean(fc_stage,1),"cityblock");
        try
        avg_man_dis(i_st,x_sub)=mean(ManhDistance(:,1));
        std_man_dis(i_st,x_sub)=std(ManhDistance(:,1));
        catch
            continue
        end
    end
end
avg_man_dis(avg_man_dis==0)=nan;
avg_man_dis(5,2)=nan;
std_man_dis(5,2)=nan;
ttest(avg_man_dis(:,1)',avg_man_dis(:,4)')
ttest(avg_man_dis(:,1)',avg_man_dis(:,3)')






%% night 1




id_diag=eye(length(subs),length(subs));
id_nodiag=~id_diag;

id_self_n1=nan(length(subs),length(stages),length(stages));
id_other_n1=nan(length(subs),length(stages),length(stages));


for trget=1:5
    id_matrix=nan(12,12,5);
    for i_st=1:length(stages)
        for x_sub=1:length(subs)
            for y_sub=1:length(subs)
                % if i_st ==trget
                %     temp_fc=corr(zscore(corr_mat_sep_n1{sub_id(x_sub)}(:,sleep_score_table_dfc_n1{sub_id(x_sub)}.Score==stages(i_st))),...
                %         zscore(corr_mat_sep{sub_id(y_sub)}(:,sleep_score_table_dfc_2_5{sub_id(y_sub)}.Score==stages(i_st))));
                %     id_matrix(x_sub,y_sub,i_st)=mean((temp_fc(temp_fc~=1)));
                % else
                    try
                        temp_fc=corr(zscore(corr_mat_sep_n1{x_sub}(:,sleep_score_table_dfc_n1{x_sub}.Score==stages(i_st))),...
                            zscore(corr_mat_sep_n1{y_sub}(:,sleep_score_table_dfc_n1{y_sub}.Score==stages(trget))));
                        id_matrix(x_sub,y_sub,i_st)=nanmean((temp_fc(temp_fc~=1)));
                    catch me
                        disp(me)
                    end
                % end
            end
        end
        temp_var=id_matrix(:,:,i_st);
        id_self_n1(:,trget,i_st)=temp_var(id_diag>0);


        id_other_n1(:,trget,i_st)=nanmean(reshape(temp_var(id_nodiag>0),length(subs)-1,length(subs)));
    end
 
    
end


id_diff_n1=id_self_n1-id_other_n1;

id_self_n1n2=nan(length(subs),length(stages),length(stages));
id_other_n1n2=nan(length(subs),length(stages),length(stages));


id_matrix_stg_n1n2=[];

for trget=1:5
    id_matrix=nan(12,12,5);
    for i_st=1:length(stages)
        for x_sub=1:length(subs)
            for y_sub=1:length(subs)
                % if i_st ==trget
                %     temp_fc=corr(zscore(corr_mat_sep_n1{sub_id(x_sub)}(:,sleep_score_table_dfc_n1{sub_id(x_sub)}.Score==stages(i_st))),...
                %         zscore(corr_mat_sep{sub_id(y_sub)}(:,sleep_score_table_dfc_2_5{sub_id(y_sub)}.Score==stages(i_st))));
                %     id_matrix(x_sub,y_sub,i_st)=mean((temp_fc(temp_fc~=1)));
                % else
                    try
                        temp_fc=corr(zscore(corr_mat_sep{x_sub}(:,sleep_score_table_dfc_2_5{x_sub}.Score==stages(i_st))),...
                            zscore(corr_mat_sep_n1{y_sub}(:,sleep_score_table_dfc_n1{y_sub}.Score==stages(trget))));
                        id_matrix(x_sub,y_sub,i_st)=nanmean((temp_fc(temp_fc~=1)));
                    catch me
                        disp(me)
                    end
                % end
            end
        end
        temp_var=id_matrix(:,:,i_st);
        id_self_n1n2(:,trget,i_st)=temp_var(id_diag>0);


        id_other_n1n2(:,trget,i_st)=nanmean(reshape(temp_var(id_nodiag>0),length(subs)-1,length(subs)));
    end
    id_matrix_stg_n1n2{trget}=id_matrix;

    
end


id_diff_n1n2=id_self_n1n2-id_other_n1n2;


for jj=1:5
figure
for ii=1:length(stages)


        subplot(2,3,ii)
        imagesc(id_matrix_stg_n1n2{jj}(:,:,ii));
        c = colorbar;
        
        c.FontSize=12;
         clim([0.1 0.7])
        colormap(bluewhitered)
        
        c.TickLength=0;
        c.Box='off';
        % xtickangle(90)
        % set(gca, 'YTick', 1:18)
        % set(gca, 'YTickLabel', net_n_sub); % set x-axis labels
        % set(gca, 'XTick', 1:18)
        xlabel(sleep_stages_names(ii),'FontSize',28,'FontWeight','Bold'); % set x-axis labels
        ylabel(sleep_stages_names(jj),'FontSize',28,'FontWeight','Bold'); % set x-axis labels
        % ylabel('REM','FontSize',28,'FontWeight','Bold'); % set x-axis labels


        set(gca,'box','off') 
        set(gcf,'Position',[800,147,3200,1823])
        % title(sleep_stages_names(ii),'FontSize',28,'FontWeight','Bold')


end
end



id_rate=nan(5,5);
for itrgt=1:5
    for istg=1:5
        
        diag_greater = false(12,1); 
        for isub=1:12
            row_values = id_matrix_stg_n1n2{itrgt}(isub, :,istg);  % Extract the row
            diag_value = id_matrix_stg_n1n2{itrgt}(isub, isub,istg);  % Extract the diagonal element
            row_values(isub) = -Inf;  % Exclude diagonal element from comparison
            if diag_value > max(row_values)
                diag_greater(isub) = true;
            end
        end
        id_rate(itrgt,istg)=sum(diag_greater)/sum(~isnan(diag(id_matrix_stg_n1n2{itrgt}(:, :,istg))));
    end
end




figure

imagesc(id_rate);
c = colorbar;

c.FontSize=12;
 clim([0.6 1])
colormap(jet)

c.TickLength=0;
c.Box='off';
% xtickangle(90)
set(gca, 'YTick', 1:5)
set(gca, 'YTickLabel', sleep_stages_names,'FontSize',18); % set x-axis labels
set(gca, 'XTick', 1:5)
set(gca, 'XTickLabel', sleep_stages_names,'FontSize',18); % set x-axis labels

xlabel('Testing (night 1)','FontSize',28,'FontWeight','Bold'); % set x-axis labels
ylabel('Database (night 2)','FontSize',28,'FontWeight','Bold'); % set x-axis labels
% ylabel('REM','FontSize',28,'FontWeight','Bold'); % set x-axis labels


set(gca,'box','off') 
set(gcf,'Position',[800,147,2000,1823])





figure;
bar(1:5,squeeze(nanmean(id_diff,1)))

set(gca, 'XTick', 1:5)
set(gca, 'XTickLabel', sleep_stages_names); % set x-axis labels
ylabel('self-slef minus self-other','FontSize',28,'FontWeight','Bold'); % set x-axis labels
        set(gcf,'Position',[800,147,3200,1823])


figure;
bar(1:5,squeeze(nanmean(id_diff_n1,1)))

figure;
bar(1:5,squeeze(nanmean(id_diff_n1n2,1)))
set(gca, 'XTick', 1:5)
set(gca, 'XTickLabel', sleep_stages_names); % set x-axis labels
ylabel('self-slef minus self-other','FontSize',28,'FontWeight','Bold'); % set x-axis labels
        set(gcf,'Position',[800,147,3200,1823])




figure
for ii=1:length(stages)


        subplot(2,3,ii)
        imagesc(id_matrix(:,:,ii));
        c = colorbar;
        
        c.FontSize=12;
         clim([0.1 0.7])
        colormap(bluewhitered)
        
        c.TickLength=0;
        c.Box='off';
        % xtickangle(90)
        % set(gca, 'YTick', 1:18)
        % set(gca, 'YTickLabel', net_n_sub); % set x-axis labels
        % set(gca, 'XTick', 1:18)
        xlabel(sleep_stages_names(ii),'FontSize',28,'FontWeight','Bold'); % set x-axis labels
        ylabel(sleep_stages_names(trget),'FontSize',28,'FontWeight','Bold'); % set x-axis labels
        % ylabel('REM','FontSize',28,'FontWeight','Bold'); % set x-axis labels


        set(gca,'box','off') 
        set(gcf,'Position',[800,147,3200,1823])
        % title(sleep_stages_names(ii),'FontSize',28,'FontWeight','Bold')


end

figure

imagesc(id_rate([1,4],:));
c = colorbar;

c.FontSize=12;
 clim([0.6 1])
colormap(jet)

c.TickLength=0;
c.Box='off';
% xtickangle(90)
set(gca, 'YTick', 1:2)
set(gca, 'YTickLabel', sleep_stages_names(1,[1 4]),'FontSize',18); % set x-axis labels
set(gca, 'XTick', 1:5)
set(gca, 'XTickLabel', sleep_stages_names,'FontSize',18); % set x-axis labels

xlabel('Target (Night 1)','FontSize',28,'FontWeight','Bold'); % set x-axis labels
ylabel('Database','FontSize',28,'FontWeight','Bold'); % set x-axis labels
% ylabel('REM','FontSize',28,'FontWeight','Bold'); % set x-axis labels


set(gca,'box','off') 
set(gcf,'Position',[800,147,1000,423])




%% Dp and theta

stages=[6,4,3,2,5];
sleep_stages_names={'Wake','N1','N2','N3','REM'};
dp_matrix=nan(size(corr_mat_sep{1},1),5);

theta_ii=[];
for istg=1:length(stages)
    for isub=1:length(subs)
        temp_fc=zscore(corr_mat_sep{isub}(:,sleep_score_table_dfc_2_5{isub}.Score==stages(istg)));
        for iedge=1:size(corr_mat_sep{1},1)
                thea_temp=temp_fc(iedge,:).*temp_fc(iedge,:)';
                thea_temp=tril(thea_temp,-1);
                theta_ii(iedge,isub,istg)=mean(thea_temp(thea_temp~=0));
        end
    end
end


theta_ij=nan(size(corr_mat_sep{1},1),12,12,5);
for istg=1:length(stages)
    for isub=1:length(subs)
        temp_fc=zscore(corr_mat_sep{isub}(:,sleep_score_table_dfc_2_5{isub}.Score==stages(istg)));
        for jsub=1:length(subs)
            if jsub~=isub
                temp_fc_j=zscore(corr_mat_sep{jsub}(:,sleep_score_table_dfc_2_5{jsub}.Score==stages(istg)));
            else
                continue
            end

            for iedge=1:size(corr_mat_sep{1},1)
                    theta_ij(iedge,isub,jsub,istg)=mean(temp_fc(iedge,:).*temp_fc_j(iedge,:)','all','omitnan');
             
            end


        end
    end
end
pij=nan(size(corr_mat_sep{1},1),12,5);

for istg=1:length(stages)
    for isub=1:length(subs) 
        for iedge=1:size(corr_mat_sep{1},1)

            temp_thea_ij=squeeze(theta_ij(iedge,isub,:,istg));
            if isnan(theta_ii(iedge,isub,istg))
                pij(iedge,isub,istg)=nan;
            else
                pij(iedge,isub,istg)=sum(theta_ii(iedge,isub,istg)<temp_thea_ij)/sum(~isnan(temp_thea_ij));
            end
        end
    end
end

pij(pij==0)=1/exp(4);

DP=squeeze(mean(-log(pij),2,'omitnan'));



Network_Label_v2=Network_Label([1:233 240:308]);
Network_Label_v2(295:end)=23;%combine hippocampal rois
Network_Label_v2(238:239)=24;
net_code=unique(Network_Label_v2);

ROI_num=1:length(Network_Label_v2);
net_num=1:18;
net_n_hip={'Un','DMN','VIS','FPN','REW','DAN','VAN','SAL','CON','dSMN','lSMN','AUD','PMN','BG','THAL','CB','HIPi','Amy'};
% Network_Label_v2=Network_Label([1:233 240:308]);
% net_code=unique(Network_Label_v2);
% net_num=1:21;
% ROI_num=1:length(Network_Label_v2);
% net_n_hip={'Un','DMN','VIS','FPN','REW','DAN','VAN','SAL','CON','dSMN','lSMN','AUD','PMN','BG','THAL','CB','aHIPi','pHIPi','ECi','paraHIPi','Amy'};

net_mat=tril(ones(18,18));
net_mat(net_mat>0)=1:sum(net_mat>0,'all');
net_mat=net_mat+tril(net_mat,-1)';
edge_net=tril(ones(302,302),-1);
for iedge=1:size(edge_net,1)
    for jedge=1:size(edge_net,2)
        if indx_mat(iedge,jedge)>0
            edge_net(iedge,jedge)=net_mat(net_num(net_code==Network_Label_v2(iedge)),net_num(net_code==Network_Label_v2(jedge)));
        end
    end

end

edge_net=edge_net(edge_net>0);



DP_net=nan(length(net_n_hip),5);
for istg=1:length(stages)
    for inet=1:length(unique(net_mat)) 
        DP_net(inet,istg)=sum(DP(edge_net==inet,istg)>prctile(DP(:,istg),95))./sum(edge_net==inet) ;

    end
end


DP_mat=nan(18,18,5);
for ii=1:5
    id_mat=tril(ones(18,18));
    id_mat(id_mat>0)=DP_net(:,ii);
    t_mat=id_mat;
      t_mat=t_mat+tril(t_mat,-1)';
            t_mat=tril(t_mat([1:2 4:9 3 10:18 ],[1:2 4:9 3 10:18 ]));
    DP_mat(:,:,ii)=t_mat;
end

% t_mat=t_mat(1:18,1:18);
%   t_mat=t_mat+tril(t_mat,-1)';
%             t_mat=tril(t_mat([1:2 4:9 3 10:18 ],[1:2 4:9 3 10:18 ]));




figure
for ii=1:length(stages)


        subplot(2,3,ii)
        imagesc(DP_mat(:,:,ii));
        c = colorbar;
        
        c.FontSize=12;
         clim([0.05 0.2])
        colormap(bluewhitered)
        
        c.TickLength=0;
        c.Box='off';
        xtickangle(90)
        set(gca, 'YTick', 1:21)
        set(gca, 'YTickLabel', net_n_hip([1:2 4:9 3 10:18 ])); % set x-axis labels
        set(gca, 'XTick', 1:21)
        set(gca, 'XTickLabel', net_n_hip([1:2 4:9 3 10:18 ])); % set x-axis labels

        % xlabel(sleep_stages_names(ii),'FontSize',28,'FontWeight','Bold'); % set x-axis labels
        % ylabel('Wake','FontSize',28,'FontWeight','Bold'); % set x-axis labels

        set(gca,'box','off') 
        set(gcf,'Position',[800,147,3200,1823])
        title(sleep_stages_names(ii),'FontSize',28,'FontWeight','Bold')
% 

end


phi_net=nan(length(net_n_hip),5);
thea=squeeze(mean(theta_ii,2,'omitnan'));
for istg=1:length(stages)
    for inet=1:length(unique(net_mat)) 
        phi_net(inet,istg)=sum(thea(edge_net==inet,istg)>prctile(thea(:,istg),95))./sum(edge_net==inet) ;

    end
end




phi_mat=nan(21,21,5);
for ii=1:5
    id_mat=tril(ones(21,21));
    id_mat(id_mat>0)=phi_net(:,ii);
    phi_mat(:,:,ii)=id_mat;
end


figure
for ii=1:length(stages)


        subplot(2,3,ii)
        imagesc(phi_mat(:,:,ii));
        c = colorbar;
        
        c.FontSize=12;
         clim([0.05 0.5])
        colormap(bluewhitered)
        
        c.TickLength=0;
        c.Box='off';
        xtickangle(90)
        set(gca, 'YTick', 1:21)
        set(gca, 'YTickLabel', net_n_hip); % set x-axis labels
        set(gca, 'XTick', 1:21)
        set(gca, 'XTickLabel', net_n_hip); % set x-axis labels

        % xlabel(sleep_stages_names(ii),'FontSize',28,'FontWeight','Bold'); % set x-axis labels
        % ylabel('Wake','FontSize',28,'FontWeight','Bold'); % set x-axis labels

        set(gca,'box','off') 
        set(gcf,'Position',[800,147,3200,1823])
        title(sleep_stages_names(ii),'FontSize',28,'FontWeight','Bold')
% 

end





