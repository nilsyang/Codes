%% setup
clearvars;clc
addpath(genpath('/misc/imeel/yangf7/matlab/'))
rmpath(genpath('/misc/imeel/yangf7/matlab/matlab_code_Nils/eeglab2023.0'))
rmpath(genpath('/misc/imeel/dezwart/matlab/from_others/eeglab2024.2'))
rmpath(genpath('/misc/imeel/yangf7/matlab/BCT'))
load('/raid/common/sleep1/derivatives/nils/ROI_list.mat');%ROI list
load('/misc/imeel/yangf7/matlab/SW_detection/label_table_0705.mat')%
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
prepro={
         'procbasic_ric_rvt_ppga' % with RETROICOR and ppga
       };
% prepro={
%
%         'procbasic_ric_gsrwb' % with RETROICOR and whole-brain GS
%  };
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
sleep_score_censor_subs=[];
rawFC_mat=[];
sleep_score_table={};
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
win_len=30;%windows length = 30TR
sleep_score_table_dfc={};
for ii=1:length(subs)
sleep_score_table_dfc{ii}=sleep_score_table{ii}(1,:);
n=1;
   for jj=1:height(sleep_score_table{ii})
       FC_steps=sleep_score_table{ii}.num_volume(jj)-win_len+1;
       for kk=1:FC_steps
           sleep_score_table_dfc{ii}.startpoint(n)=sleep_score_table{ii}.startpoint(jj)+kk-1;
           sleep_score_table_dfc{ii}.endpoint(n)=sleep_score_table_dfc{ii}.startpoint(n)+win_len-1;
           sleep_score_table_dfc{ii}.num_volume(n)=win_len;
           sleep_score_table_dfc{ii}.Score(n)=sleep_score_table{ii}.Score(jj);
           sleep_score_table_dfc{ii}.subs(n)=sleep_score_table{ii}.subs(jj);
           sleep_score_table_dfc{ii}.nums(n)=n;
           n=n+1;
       end
  
  
   end
end
%% dfc
dfc_table_256=[];
score_kmeans=[];
for ii=1:length(subs)
   score_kmeans=[score_kmeans;sleep_score_table_dfc{ii}.Score];
   dfc_table_256=[dfc_table_256;sleep_score_table_dfc{ii}];
end
%% 302*302
% net_code=unique(Network_Label);
% ROI_num=1:length(Network_Label);
% corr_mat_sep={};corr_mat_sep_all=[];
% indx_mat=tril(ones(302,302),-1);
% % indx_mat(1:294,:)=0;%only hpc connections
% for ii=1:length(sleep_score_table_dfc)
% % corr_mat_sep{ii}=nan(302*301/2,height(sleep_score_table_dfc{ii}));
% corr_mat_sep{ii}=nan(sum(indx_mat==1,'all'),height(sleep_score_table_dfc{ii}));
% for jj=1:height(sleep_score_table_dfc{ii})
% % temp_mat=atanh(partialcorr(rawFC_mat{ii}(:,sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj))',mean(rawFC_mat{ii}(262:273,sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj)))'));
% temp_mat=atanh(corr(rawFC_mat{ii}([1:233 240:308],sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj))'));%remove 234-239 MTL + hippo
% temp_mat(logical(eye(size(temp_mat,1))))=4;
% corr_mat_sep{ii}(:,jj)=temp_mat(indx_mat>0);
% end
% corr_mat_sep_all=cat(2,corr_mat_sep_all,corr_mat_sep{ii});
% end
%% subcortical
net_code=unique(Network_Label);
ROI_num=1:length(Network_Label);
corr_mat_sep={};corr_mat_sep_all=[];
indx_mat=tril(ones(38,38),-1);
% indx_mat(1:294,:)=0;%only hpc connections
for ii=1:length(sleep_score_table_dfc)
% corr_mat_sep{ii}=nan(302*301/2,height(sleep_score_table_dfc{ii}));
corr_mat_sep{ii}=nan(sum(indx_mat==1,'all'),height(sleep_score_table_dfc{ii}));
for jj=1:height(sleep_score_table_dfc{ii})
% temp_mat=atanh(partialcorr(rawFC_mat{ii}(:,sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj))',mean(rawFC_mat{ii}(262:273,sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj)))'));
temp_mat=atanh(corr(rawFC_mat{ii}([244:273 301:308],sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj))'));%remove 234-239 MTL + hippo
temp_mat(logical(eye(size(temp_mat,1))))=4;
corr_mat_sep{ii}(:,jj)=temp_mat(indx_mat>0);
end
corr_mat_sep_all=cat(2,corr_mat_sep_all,corr_mat_sep{ii});
end
t_mat=zeros(size(indx_mat));
sub_roi_names={'AMYl-REW','AMYr-REW','NAccl-REW','NAccr-REW','CAUheadl-SAL','CAUheadr-SAL','CAUtaill-FPN','CAUtailr-FPN','PUT1l-SMN','PUT1r-SMN'...
   'PUT2l-SMN','PUT2r-SMN','PUT3l-VAN','PUT3r-VAN','PUT4l-CON','PUT4r-CON','PALl-CON','PALr-CON',...
   'THAL1l-DMN','THAL1r-DMN','THAL2l-VIS','THAL2r-VIS','THAL3l-CON','THAL3r-CON','THAL4l-CON','THAL4r-CON','THAL5l-SMN','THAL5r-SMN','THAL6-SMN','THAL6r-SMN',...
   'aHIPl','aHIPr','pHIPl','pHIPr',...
   'ECl','ECr','paraHIPl','paraHIPr'};
sleep_stages=[2,3,4,5,6];temp_fc=zeros(703,5,12);
figure
for ii=1:length(sleep_stages)
   for jj=1:length(subs)
       temp_fc(:,ii,jj)=mean(corr_mat_sep_all(:,dfc_table_256.Score==sleep_stages(ii)&dfc_table_256.subs==jj),2);
      
   end
       t_mat(indx_mat>0)=nanmean(temp_fc(:,ii,:),3);
       subplot(2,3,ii)
       imagesc(t_mat);
       c = colorbar;
      
       c.FontSize=12;
        clim([-0.5 1])
       colormap(bluewhitered)
      
       c.TickLength=0;
       c.Box='off';
       xtickangle(45)
       set(gca, 'YTick', 1:38)
       set(gca, 'YTickLabel', sub_roi_names); % set x-axis labels
       set(gca, 'XTick', 1:38)
       set(gca, 'XTickLabel', sub_roi_names); % set x-axis labels
       set(gca,'box','off')
       set(gcf,'Position',[8000 8000 2000 2000])
       title(['Sleep score ', num2str(sleep_stages(ii))])
end
Network_Label_sub=[1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,12,12,13,13,14,14,15,15,16,16,17,17,18,18]';
net_code_sub=unique(Network_Label_sub);
net_n_sub={'AMY_R_E_W','NAcc_R_E_W','CAUhead_S_A_L','CAUtail_F_P_N','PUT_d_S_M_N','PUT_l_S_M_N',...
   'PUT_V_A_N','PUT_C_O_N','PAL_C_O_N',...
   'THAL_D_M_N','THAL_V_I_S','THAL_C_O_N','THAL_l_S_M_N','THAL_d_S_M_N',...
   'aHIP','pHIP',...
   'EC','paraHIP'};
ROI_num=1:length(Network_Label_sub);
ind_net=[tril(ones(18,18))];
corr_mat_net={};corr_mat_net_all=[];
for ii=1:length(sleep_score_table_dfc)
corr_mat_net{ii}=nan(sum(ind_net==1,'all'),height(sleep_score_table_dfc{ii}));
   for jj=1:height(sleep_score_table_dfc{ii})
       % temp_mat=atanh(partialcorr(rawFC_mat{ii}(:,sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj))',mean(rawFC_mat{ii}(262:273,sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj)))'));
       temp_mat=atanh(corr(rawFC_mat{ii}([244:273 301:308],sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj))'));%remove 234-239 MTL + hippo
       temp_mat(logical(eye(size(temp_mat,1))))=4;
      
       w_wch=zeros(length(net_code_sub),1);%
       % b_wch=zeros(length(net_code_sub),1);
       bw_wch=zeros(length(net_code_sub)*(length(net_code_sub)-1)/2,1);
      
      
       net_mat=ones(length(net_code_sub),length(net_code_sub));
       net_mat_tril=tril(net_mat,-1);
       % net_mat_tril=[net_mat_tril,zeros(length(net_code),1)];
       net_mat_diag=eye(length(net_code_sub),length(net_code_sub));
       % net_mat_diag=[net_mat_diag,zeros(length(net_code),1)];
      
     
       % for k=1:size(temp_mat,3)
      
      
           for n_net = 1:length(net_code_sub)
               net_nodes=ROI_num(Network_Label_sub==net_code_sub(n_net));
                   m_net=temp_mat(net_nodes,net_nodes);
                   m_net=tril(m_net,-1);
           %         m_net_s(:,n_sub)=m_net(m_net>0);
                   ind_mat_temp=ones(numel(net_nodes),numel(net_nodes));
                   ind_mat_temp=tril(ind_mat_temp,-1);
                   ind_mat_temp=ind_mat_temp.*~isnan(m_net);
                   w_wch(n_net)=nansum(nansum(m_net))/sum(sum(ind_mat_temp));
              
           end
      
      
    
      
      
      
      
         n=1;
           for n_net = 1:length(net_code_sub)-1
               for n_net2= n_net+1:length(net_code_sub)
                   net_nodes=ROI_num(Network_Label_sub==net_code_sub(n_net));
                   net_nodes_2=ROI_num(Network_Label_sub==net_code_sub(n_net2));
                       m_net=temp_mat(net_nodes,net_nodes_2);
                       ind_mat_temp=ones(numel(net_nodes),numel(net_nodes_2));
                       ind_mat_temp=ind_mat_temp.*~isnan(m_net);
      
       %                 ind_mat_temp=tril(ind_mat,-1);
          
                       bw_wch(n)=nansum(nansum(m_net))/sum(ind_mat_temp(:));
                    
                  
                n=n+1;
               end
              
           end
      
      
      
      
              mat_test=zeros(length(net_code_sub),length(net_code_sub));
              mat_test(net_mat_diag>0)=w_wch;
              mat_test(net_mat_tril>0)=bw_wch;
              % mat_test(:,length(net_code)+1)=b_wch;
             
              corr_mat_net{ii}(:,jj)=mat_test(ind_net>0);
      
   end
   corr_mat_net_all=cat(2,corr_mat_net_all,corr_mat_net{ii});
end
t_mat=zeros(size(ind_net));
t_mat_cell=[];
sleep_stages=[6,4,3,2,5];temp_fc=zeros(size(corr_mat_net_all,1),5,12);
sleep_stages_names={'Wake','N1','N2','N3','REM','REM-Wake'};
figure
for ii=1:length(sleep_stages)+1
   if ii==length(sleep_stages)+1
       t_mat=t_mat_cell{5}-t_mat_cell{1};
   else
   for jj=1:length(subs)
       temp_fc(:,ii,jj)=mean(corr_mat_net_all(:,dfc_table_256.Score==sleep_stages(ii)&dfc_table_256.subs==jj),2);
      
   end
       t_mat(ind_net>0)=nanmean(temp_fc(:,ii,:),3);
       t_mat_cell{ii}=t_mat;
   end
       subplot(2,3,ii)
       imagesc(t_mat);
       c = colorbar;
      
       c.FontSize=12;
        clim([-0.5 1])
       colormap(bluewhitered)
      
       c.TickLength=0;
       c.Box='off';
       xtickangle(90)
       set(gca, 'YTick', 1:18)
       set(gca, 'YTickLabel', net_n_sub); % set x-axis labels
       set(gca, 'XTick', 1:18)
       set(gca, 'XTickLabel', net_n_sub); % set x-axis labels
       set(gca,'box','off')
       set(gcf,'Position',[800,147,3200,1823])
       title(sleep_stages_names(ii),'FontSize',28,'FontWeight','Bold')
end
for ii=1:length(sleep_stages)
   for jj=1:length(subs)
           temp_fc(:,ii,jj)=mean(corr_mat_net_all(:,dfc_table_256.Score==sleep_stages(ii)&dfc_table_256.subs==jj),2);
          
  
   end
   t_mat(ind_net>0)=nanmean(temp_fc(:,ii,:),3);
           t_mat_cell{ii}=t_mat;
end
for ii=1:4
   t_mat_cell{ii}=tril(t_mat_cell{5}-t_mat_cell{ii})';
end
figure
for ii=1:length(sleep_stages)
      
       subplot(2,3,ii)
       imagesc(t_mat_cell{ii});
       c = colorbar;
      
       c.FontSize=14;
        clim([-0.5 1.5])
       colormap(bluewhitered)
      
       c.TickLength=0;
       c.Box='off';
       xtickangle(90)
       set(gca, 'YTick', 1:19)
       set(gca, 'YTickLabel', net_n_sub); % set x-axis labels
       set(gca, 'XTick', 1:19)
       set(gca, 'XTickLabel', net_n_sub); % set x-axis labels
       set(gca,'box','off')
       set(gcf,'Position',[800,147,3200,1850])
       title(sleep_stages_names(ii),'FontSize',28,'FontWeight','Bold')
end
%% thal and cortical network
Network_Label_v2=Network_Label([1:233 240:308]);
net_code=unique(Network_Label_v2);
ROI_num=1:length(Network_Label_v2);
net_n_hip={'Un','DMN','VIS','FPN','REW','DAN','VAN','SAL','CON','dSMN','lSMN','AUD','PMN','BG','THAL','CB','aHIPi','pHIPi','ECi','paraHIPi','Amy'};
ind_net=[tril(ones(21,21)), ones(21,1)];
corr_mat_net={};corr_mat_net_all=[];
for ii=1:length(sleep_score_table_dfc)
corr_mat_net{ii}=nan(sum(ind_net==1,'all'),height(sleep_score_table_dfc{ii}));
   for jj=1:height(sleep_score_table_dfc{ii})
       % temp_mat=atanh(partialcorr(rawFC_mat{ii}(:,sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj))',mean(rawFC_mat{ii}(262:273,sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj)))'));
       temp_mat=atanh(corr(rawFC_mat{ii}([1:233 240:308],sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj))'));%remove 234-239 MTL + hippo
       temp_mat(logical(eye(size(temp_mat,1))))=4;
      
       w_wch=zeros(length(net_code),1);%
       b_wch=zeros(length(net_code),1);
       bw_wch=zeros(length(net_code)*(length(net_code)-1)/2,1);
      
      
       net_mat=ones(length(net_code),length(net_code));
       net_mat_tril=tril(net_mat,-1);
       % net_mat_tril=[net_mat_tril,zeros(length(net_code),1)];
       net_mat_diag=eye(length(net_code),length(net_code));
       % net_mat_diag=[net_mat_diag,zeros(length(net_code),1)];
      
     
       % for k=1:size(temp_mat,3)
      
      
           for n_net = 1:length(net_code)
               net_nodes=ROI_num(Network_Label_v2==net_code(n_net));
                   m_net=temp_mat(net_nodes,net_nodes);
                   m_net=tril(m_net,-1);
           %         m_net_s(:,n_sub)=m_net(m_net>0);
                   ind_mat_temp=ones(numel(net_nodes),numel(net_nodes));
                   ind_mat_temp=tril(ind_mat_temp,-1);
                   ind_mat_temp=ind_mat_temp.*~isnan(m_net);
                   w_wch(n_net)=nansum(nansum(m_net))/sum(sum(ind_mat_temp));
              
           end
      
      
           for n_net = 1:length(net_code)
               net_nodes=ROI_num(Network_Label_v2==net_code(n_net));
               net_nodes_2=ROI_num(Network_Label_v2<17&~(Network_Label_v2==net_code(n_net)));%all the cortical ROIs excluding the current one.
                   m_net=temp_mat(net_nodes,net_nodes_2);
                   ind_mat_temp=ones(numel(net_nodes),numel(net_nodes_2));
                   ind_mat_temp=ind_mat_temp.*~isnan(m_net);
       %             ind_mat_temp=tril(ind_mat,-1);
                   b_wch(n_net)=nansum(nansum(m_net))/sum(ind_mat_temp(:));
              
           end
      
      
      
      
         n=1;
           for n_net = 1:length(net_code)-1
               for n_net2= n_net+1:length(net_code)
                   net_nodes=ROI_num(Network_Label_v2==net_code(n_net));
                   net_nodes_2=ROI_num(Network_Label_v2==net_code(n_net2));
                       m_net=temp_mat(net_nodes,net_nodes_2);
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
             
              corr_mat_net{ii}(:,jj)=mat_test(ind_net>0);
      
   end
   corr_mat_net_all=cat(2,corr_mat_net_all,corr_mat_net{ii});
end
Network_Label_sub=[ROI_list.netWorkbenchLabel(1:243);[20,20,21,21,22,22,22,22,23,23,24,24]'];
net_code_sub=unique(Network_Label_sub);
net_n_sub={'Un','DMN','VIS','FPN','REW','DAN','VAN','SAL','CON','dSMN','lSMN','AUD','PMN','MTL'...
   'THAL_D_M_N','THAL_V_I_S','THAL_C_O_N','THAL_l_S_M_N','THAL_d_S_M_N'};
ROI_num=1:length(Network_Label_sub);
ind_net=[tril(ones(19,19))];
corr_mat_net={};corr_mat_net_all=[];
for ii=1:length(sleep_score_table_dfc)
corr_mat_net{ii}=nan(sum(ind_net==1,'all'),height(sleep_score_table_dfc{ii}));
   for jj=1:height(sleep_score_table_dfc{ii})
       % temp_mat=atanh(partialcorr(rawFC_mat{ii}(:,sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj))',mean(rawFC_mat{ii}(262:273,sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj)))'));
       temp_mat=atanh(corr(rawFC_mat{ii}([1:243 262:273],sleep_score_table_dfc{ii}.startpoint(jj):sleep_score_table_dfc{ii}.endpoint(jj))'));%remove 234-239 MTL + hippo
       temp_mat(logical(eye(size(temp_mat,1))))=4;
      
       w_wch=zeros(length(net_code_sub),1);%
       % b_wch=zeros(length(net_code_sub),1);
       bw_wch=zeros(length(net_code_sub)*(length(net_code_sub)-1)/2,1);
      
      
       net_mat=ones(length(net_code_sub),length(net_code_sub));
       net_mat_tril=tril(net_mat,-1);
       % net_mat_tril=[net_mat_tril,zeros(length(net_code),1)];
       net_mat_diag=eye(length(net_code_sub),length(net_code_sub));
       % net_mat_diag=[net_mat_diag,zeros(length(net_code),1)];
      
     
       % for k=1:size(temp_mat,3)
      
      
           for n_net = 1:length(net_code_sub)
               net_nodes=ROI_num(Network_Label_sub==net_code_sub(n_net));
                   m_net=temp_mat(net_nodes,net_nodes);
                   m_net=tril(m_net,-1);
           %         m_net_s(:,n_sub)=m_net(m_net>0);
                   ind_mat_temp=ones(numel(net_nodes),numel(net_nodes));
                   ind_mat_temp=tril(ind_mat_temp,-1);
                   ind_mat_temp=ind_mat_temp.*~isnan(m_net);
                   w_wch(n_net)=nansum(nansum(m_net))/sum(sum(ind_mat_temp));
              
           end
      
      
    
      
      
      
      
         n=1;
           for n_net = 1:length(net_code_sub)-1
               for n_net2= n_net+1:length(net_code_sub)
                   net_nodes=ROI_num(Network_Label_sub==net_code_sub(n_net));
                   net_nodes_2=ROI_num(Network_Label_sub==net_code_sub(n_net2));
                       m_net=temp_mat(net_nodes,net_nodes_2);
                       ind_mat_temp=ones(numel(net_nodes),numel(net_nodes_2));
                       ind_mat_temp=ind_mat_temp.*~isnan(m_net);
      
       %                 ind_mat_temp=tril(ind_mat,-1);
          
                       bw_wch(n)=nansum(nansum(m_net))/sum(ind_mat_temp(:));
                    
                  
                n=n+1;
               end
              
           end
      
      
      
      
              mat_test=zeros(length(net_code_sub),length(net_code_sub));
              mat_test(net_mat_diag>0)=w_wch;
              mat_test(net_mat_tril>0)=bw_wch;
              % mat_test(:,length(net_code)+1)=b_wch;
             
              corr_mat_net{ii}(:,jj)=mat_test(ind_net>0);
      
   end
   corr_mat_net_all=cat(2,corr_mat_net_all,corr_mat_net{ii});
end
t_mat=zeros(size(ind_net));t_mat_cell=[];
sleep_stages=[6,4,3,2,5];temp_fc=zeros(size(corr_mat_net_all,1),5,12);
sleep_stages_names={'Wake','N1','N2','N3','REM','REM-Wake'};
for ii=1:length(sleep_stages)
   for jj=1:length(subs)
           temp_fc(:,ii,jj)=mean(corr_mat_net_all(:,dfc_table_256.Score==sleep_stages(ii)&dfc_table_256.subs==jj),2);
          
  
   end
   t_mat(ind_net>0)=nanmean(temp_fc(:,ii,:),3);
       t_mat=t_mat+tril(t_mat,-1)';
           t_mat_cell{ii}=tril(t_mat([1:2 4:9 3 10:19 ],[1:2 4:9 3 10:19 ]));
           t_mat=zeros(size(ind_net));
end
%
for ii=1:4
   t_mat_cell{ii}=tril(t_mat_cell{5}-t_mat_cell{ii})';
end
net_n_sub_reorder={'Un','DMN','FPN','REW','DAN','VAN','SAL','CON','VIS','dSMN','lSMN','AUD','PMN','MTL'...
   'THAL_D_M_N','THAL_V_I_S','THAL_C_O_N','THAL_l_S_M_N','THAL_d_S_M_N'};
figure
for ii=1:length(sleep_stages)
      
       subplot(2,3,ii)
       imagesc(t_mat_cell{ii});
       c = colorbar;
      
       c.FontSize=14;
        clim([-0.5 1.5])
       colormap(bluewhitered)
      
       c.TickLength=0;
       c.Box='off';
       xtickangle(90)
       set(gca, 'YTick', 1:19)
       set(gca, 'YTickLabel', net_n_sub_reorder); % set x-axis labels
       set(gca, 'XTick', 1:19)
       set(gca, 'XTickLabel', net_n_sub_reorder); % set x-axis labels
       set(gca,'box','off')
       set(gcf,'Position',[800,147,3200,1850])
       title(sleep_stages_names(ii),'FontSize',28,'FontWeight','Bold')
end
dfc_table_256.stages=categorical(dfc_table_256.Score);
% statistical test
p_mat=ones(size(ind_net));
p_mat_test=1:190;
for ii=1:size(corr_mat_net_all,1)
   dfc_table_256.fc=corr_mat_net_all(ii,:)';
   dfc_table_56=dfc_table_256(dfc_table_256.Score>4,:);%only wake or rem
   lme=fitlme(dfc_table_256,'fc~stages+(1|subs)','DummyVarCoding','effects');
       % lme=fitlme(dfc_table_56,'fc~stages+(1|subs)');
   % anova_table=anova(lme,'DFmethod','Satterthwaite');
   [beta,~,stats]=fixedEffects(lme);% fixed effects
  test_table=dataset2table(stats);
  test_table.Properties.RowNames=test_table.Name;
  p_mat_test(ii)=test_table{'stages_5','pValue'};
  % t_value(jj)=test_table{'Score','tStat'};
                  % anova_table2=dataset2table(anova_table);
                  % anova_table2.Properties.RowNames=anova_table2.Term;
                  %
                  % p_mat_test(ii)=anova_table2{'stages','pValue'};
end
p_mat(ind_net>0)=p_mat_test<0.0001;
%% 20*21
dfc_table_256=[];
score_kmeans=[];
for ii=1:length(subs)
   score_kmeans=[score_kmeans;sleep_score_table_dfc{ii}.Score];
   dfc_table_256=[dfc_table_256;sleep_score_table_dfc{ii}];
end
t_mat=zeros(size(ind_net));
temp_fc=zeros(252,5,12);
sleep_stages=[6,4,3,2,5];
sleep_stages_names={'Wake','N1','N2','N3','REM'};
figure
for ii=1:length(sleep_stages)
   for jj=1:length(subs)
       temp_fc(:,ii,jj)=mean(corr_mat_net_all(:,dfc_table_256.Score==sleep_stages(ii)&dfc_table_256.subs==jj),2);
      
   end
       t_mat(ind_net>0)=nanmean(temp_fc(:,ii,:),3);
       subplot(2,3,ii)
       imagesc(t_mat);
       c = colorbar;
      
       c.FontSize=12;
        clim([-0.3 0.6])
       colormap(bluewhitered)
      
       c.TickLength=0;
       c.Box='off';
       xtickangle(45)
       set(gca, 'YTick', 1:21)
       set(gca, 'YTickLabel', net_n_hip); % set x-axis labels
       set(gca, 'XTick', 1:22)
       set(gca, 'XTickLabel', [net_n_hip,'one\_vs\_cor']); % set x-axis labels
       set(gca,'box','off')
       set(gcf,'Position',[8000 8000 2000 2000])
       title(sleep_stages_names(ii))
end
for ii=[2,3,4,5,6]
   median_fc(ii)=median(corr_mat_net_all(204,dfc_table_256.Score==ii));
end
FC_table=dfc_table_256(ismember(dfc_table_256.Score,[2,3,4,5]),:);
FC_table.swd=sw_freq_table_all(ismember(dfc_table_256.Score,[2,3,4,5]));
corr_mat_net_n2=corr_mat_net_all(:,ismember(dfc_table_256.Score,[2,3,4,5]));
p_mat=nan(size(ind_net));t_mat=nan(size(ind_net));t_value=[];p_value=[];
for jj=1:size(corr_mat_net_all,1)
 
      
       FC_table.FC=corr_mat_net_n2(jj,:)';
% FC_table.Score=categorical( FC_table.Score);
  lme=fitlme(FC_table,'FC~Score+(1|subs)');% linear mixed effect model
  [beta,~,stats]=fixedEffects(lme);% fixed effects
  test_table=dataset2table(stats);
  test_table.Properties.RowNames=test_table.Name;
  p_value(jj)=test_table{'Score','pValue'};
  t_value(jj)=test_table{'Score','tStat'};
end
t_mat(ind_net>0)=t_value;
t_mat(ind_net==0)=0;
figure
imagesc(t_mat);
c = colorbar;
c.FontSize=12;
colormap(bluewhitered)
% clim([-1 1])
c.TickLength=0;
c.Box='off';
xtickangle(45)
set(gca, 'YTick', 1:20)
set(gca, 'YTickLabel', net_n_hip); % set x-axis labels
set(gca, 'XTick', 1:20)
set(gca, 'XTickLabel', [net_n_hip,'one_vs_cor']); % set x-axis labels
set(gca,'box','off')
set(gcf,'Position',[8000 8000 2000 2000])
FC_table=dfc_table_256(dfc_table_256.Score==2|dfc_table_256.Score==6,:);
FC_table.swd=sw_freq_table_all(dfc_table_256.Score==2|dfc_table_256.Score==6);
corr_mat_net_n2=corr_mat_net_all(:,dfc_table_256.Score==2|dfc_table_256.Score==6);
p_mat=nan(size(ind_net));t_mat=nan(size(ind_net));t_value=[];p_value=[];
for jj=1:size(corr_mat_net_all,1)
 
      
       FC_table.FC=corr_mat_net_n2(jj,:)';
  lme=fitlme(FC_table,'FC~Score+(1|subs)');% linear mixed effect model
  [beta,~,stats]=fixedEffects(lme);% fixed effects
  test_table=dataset2table(stats);
  test_table.Properties.RowNames=test_table.Name;
  p_value(jj)=test_table{'Score','pValue'};
  t_value(jj)=test_table{'Score','tStat'};
end
t_mat(ind_net>0)=t_value;
t_mat(ind_net==0)=0;
figure
imagesc(t_mat);
c = colorbar;
c.FontSize=12;
colormap(bluewhitered)
% clim([-1 1])
c.TickLength=0;
c.Box='off';
xtickangle(45)
set(gca, 'YTick', 1:21)
set(gca, 'YTickLabel', net_n_hip); % set x-axis labels
set(gca, 'XTick', 1:22
set(gca, 'XTickLabel', [net_n_hip,'one_vs_cor']); % set x-axis labels
set(gca,'box','off')
set(gcf,'Position',[8000 8000 2000 2000])
temp=dfc_table_256.Score;
temp(temp==6)=1;
temp(temp==2)=99;
temp(temp==4)=2;
temp(temp==99)=4;
x=unique(temp(temp>0));
for ii=1:5
   data_temp(ii)=mean(corr_mat_net_all(204,temp==x(ii)));
   err_temp(ii)=std(corr_mat_net_all(204,temp==x(ii)))./sqrt(length(corr_mat_net_all(204,temp==x(ii))));
       % err_temp(ii)=std(corr_mat_net_all(204,temp==x(ii)));
end
figure
bar(x,data_temp)               
hold on
er = errorbar(x,data_temp,err_temp,err_temp);   
er.Color = [0 0 0];                           
er.LineStyle = 'none';
figure
for ii=1:12
subplot(2,6,ii)
try
boxplot(corr_mat_net_all(204,temp>0&dfc_table_256.subs==ii),temp(temp>0&dfc_table_256.subs==ii),...
    'Labels',{'Wake','N1','N2','N3','REM'})
title(['Participant ',num2str(ii,'%02d')])
catch
boxplot([corr_mat_net_all(204,temp>0&dfc_table_256.subs==ii),0],[temp(temp>0&dfc_table_256.subs==ii);5],...
    'Labels',{'Wake','N1','N2','N3','REM'})
title(['Participant ',num2str(ii,'%02d')])
end
% set(gca,'FontSize',18)
end
set(gcf,'Position',[8000 8000 2000 2000])
figure
for ii=1:12
subplot(2,6,ii)
boxplot(sw_freq_table_all(dfc_table_256.Score>0&dfc_table_256.subs==ii),dfc_table_256.Score(dfc_table_256.Score>0&dfc_table_256.subs==ii))
end
set(gcf,'Position',[8000 8000 2000 2000])
[t,p]=corr(corr_mat_net_all(204,ismember(dfc_table_256.Score,[2,3,4]))',sw_freq_table_all(ismember(dfc_table_256.Score,[2,3,4])),'rows','pairwise')
[t,p]=corr(corr_mat_net_all(129,ismember(dfc_table_256.Score,[2,3,4]))',sw_freq_table_all(ismember(dfc_table_256.Score,[2,3,4])),'rows','pairwise')
[t,p]=corr(corr_mat_net_all(204,dfc_table_256.Score==2)',sw_amp_table_all(dfc_table_256.Score==2),'rows','pairwise')
[t,p]=corr(corr_mat_net_all(129,dfc_table_256.Score==2)',sw_amp_table_all(dfc_table_256.Score==2),'rows','pairwise')
% [t,p]=corr(sw_freq_table_all(dfc_table_256.Score==2),sw_amp_table_all(dfc_table_256.Score==2),'rows','pairwise')
figure;
for ii=1:12
subplot(2,6,ii)
% boxplot(corr_mat_net_all(190,dfc_table_256.Score>0&dfc_table_256.subs==ii),dfc_table_256.Score(dfc_table_256.Score>0&dfc_table_256.subs==ii))
scatter(sw_freq_table_all(ismember(dfc_table_256.Score,[2])&dfc_table_256.subs==ii),corr_mat_net_all(204,ismember(dfc_table_256.Score,[2])&dfc_table_256.subs==ii))
[t(ii),p]=corr(corr_mat_net_all(204,ismember(dfc_table_256.Score,[2])&dfc_table_256.subs==ii)',sw_freq_table_all(ismember(dfc_table_256.Score,[2])&dfc_table_256.subs==ii),'rows','pairwise')
mean_bg_fc(ii)=median(corr_mat_net_all(204,ismember(dfc_table_256.Score,[2])&dfc_table_256.subs==ii)');
mean_freq(ii)=median(sw_freq_table_all(ismember(dfc_table_256.Score,[2])&dfc_table_256.subs==ii));
end
set(gcf,'Position',[8000 8000 2000 2000])
figure;
scatter(sw_freq_table_all(dfc_table_256.Score==2),corr_mat_net_all(204,dfc_table_256.Score==2))
set(gcf,'Position',[8000 8000 2000 2000])
%% shaded error bar
median_FC=nan(12,6);iqr_FC=nan(12,6);
for ii=1:12
   for jj=[2:6]
       var1=corr_mat_net_all(204,dfc_table_256.Score==jj&dfc_table_256.subs==ii);
       median_FC(ii,jj)=median(var1);
       iqr_FC(ii,jj)=1.57*iqr(var1)/sqrt(length(var1));
   end
end
data_points={'Wake','N1','N2','N3','REM'};
figure;
col_bar=jet(12);
for ii=1:12
   plot(1:5,median_FC(ii,[6,4:-1:2,5]),'Color',col_bar(ii,:),'LineWidth',2);
   hold on
set(gca, 'TickLength',[0 0])
% yline([ 0.15],'-','0.15')
% ylim([0 0.7])
xlim([0.5 5.5])
set(gca, 'XTick', 1:5); % center x-axis ticks on bins
set(gca, 'XTickLabel',data_points  ,'FontSize',24,'FontWeight','Bold');
end
ylabel('Intra-thalamic FC'  ,'FontSize',28,'FontWeight','Bold')
set(gcf,'Position',[2500 400 1600 1200])

