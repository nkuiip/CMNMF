clear;
path(path,'../3_1_code(unsupervise)/CMNMF'); 
path(path,'../2_useful_data/mouse/');
load('CMNMF_mmu_result_pathway_20161109T094403.mat','CMNMF_result_cell');
result_cell = CMNMF_result_cell{2,1};
W_out = result_cell{2,1};    
load('g_p_network_2016_2.mat','gene_id');
load('mmu_pathway_data_2016_2.mat','pathway_old');
Z_threshold = 3;
[W_Z_filtered,predicted_pathway_gene] = predicted_pathway(W_out,Z_threshold,gene_id);
pathway_predicted_matrix = zeros(size(W_out,2),size(W_out,2));
for i = 1:size(W_out,2)
    for j = 1:size(W_out,2)
        genes = pathway_old(i,pathway_old(i,:)>0);
        genes2 = gene_id(W_Z_filtered(:,j)>0);
        [C,~,~] = intersect(genes,genes2);
        pathway_predicted_matrix(i,j) = length(C);
    end
end

%接下来根据pathway_predicted_matrix来进行分析
%先是看交集最大的
intersect_num = max(max(pathway_predicted_matrix));
[old_cluster,predict_cluster] = find(pathway_predicted_matrix == intersect_num);
%发现是259簇，pathways in cancer这个pathway，没什么分析意义
precision_vector = zeros(1,size(W_out,2));
for i = 1 : size(W_out,2)
   [num,row] = max(pathway_predicted_matrix(i,:)); 
   [X,~,~] = intersect(unique(pathway_old(i,:)),gene_id);
   old_num = length(X);
   precision_vector(i) = num/old_num;
end
%这个precision_vector希望能作一张图或者别的更清楚的看。


[~,useful_predict] = find(precision_vector>=0.15); %old pathway重合高
useful_predict_genes = zeros(length(useful_predict),1);
for i = 1 : length(useful_predict)
    vector = pathway_old(useful_predict(i),:);
    useful_predict_genes(i) = sum(vector~=0);
end
%275pathway 慢性骨髓白血病符合度高，0.15，可用进一步看


for i = 1 : size(W_out,2)
    [X,~,~] = intersect(unique(pathway_old(i,:)),gene_id);
    sum_vector(i) = length(X);
end
