%clear清楚普通变量，不会清楚全局变量
clear;
%设置工作path
path(path,'../3_1_code(unsupervise)/CMNMF/old_result'); 
path(path,'../2_useful_data/human/');
%加载.mat文件中的变量
load('CMNMF_hsa_result_pathway_20161029T094839.mat','CMNMF_result_cell');
%result_cell赋值，{}用于cell类型数组赋值
result_cell = CMNMF_result_cell{2,1};
W_out = result_cell{2,1};    
load('g_p_network_2016_2.mat','gene_id');
load('hsa_pathway_data_2016_2.mat','pathway_old');
Z_threshold = 3;
[W_Z_filtered,predicted_pathway_gene] = predicted_pathway(W_out,Z_threshold,gene_id);
pathway_predicted_matrix = zeros(25,size(W_Z_filtered,2));

%pathway 210 214 215 230 231 232 233 234  264:285
%select 30 old pathway data as example,each row of it represents genes in this pathway
W_Z_selected = zeros(25,size(pathway_old,2));
W_Z_selected(1,:) = pathway_old(210,:);
W_Z_selected(2,:) = pathway_old(214,:);
W_Z_selected(3,:) = pathway_old(215,:);
for i = 4:8
    W_Z_selected(i,:) = pathway_old(226+i,:);
end
for i = 9:30
    W_Z_selected(i,:) = pathway_old(255+i,:);
end
for i = 1:30
    for j = 1:size(W_out,2)
        genes = W_Z_selected(i,W_Z_selected(i,:)>0);
        genes2 = gene_id(W_Z_filtered(:,j)>0);
        [C,~,~] = intersect(genes,genes2);
        pathway_predicted_matrix(i,j) = length(C);
    end
end
% frequency_array = zeros(30,1);
% for i = 1:size(W_out,2)
%    array = pathway_predicted_matrix(:,i);
%    [B,I] = sort(array,'descend'); 
%     for j = 1:6
%        if(B(j)>0)
%            frequency_array(I(j),1) = frequency_array(I(j),1)+1;
%        end
%     end
% end


%233 230 231 274 281 279 pathway
%140 173 2 239 127 143 cluster

%array2---
array = [140,173,2,239,127,143];
array2 = [233,230,231,274,281,279];
array3 = [7,4,5,19,26,24];
result = cell(6,4);
load('hsa_pathway_data_2016_2&9.mat');
for i = 1:6
   result{i,1} = array(i);
   result{i,2} = array2(i);
   genes = W_Z_selected(array3(i),W_Z_selected(array3(i),:)>0);
   genes2 = gene_id(W_Z_filtered(:,array(i))>0);
   [C,~,~] = intersect(genes,genes2);
   symbol_array = cell(length(C),1);
   for j = 1:length(C)
      [r,c] = find(hsa_pathway_ncbi_id_2016_2==C(j)) ;
      if(isempty(r)~=1)
          symbol_array{j,1} = hsa_pathway_2016_2_ncbi_symbol{r(1),c(1)};
      end
   end
   result{i,3} = C';
   result{i,4} = symbol_array';
end

rmpath('../3_1_code(unsupervise)/CMNMF'); 
rmpath('../2_useful_data/human/');
