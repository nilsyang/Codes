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

