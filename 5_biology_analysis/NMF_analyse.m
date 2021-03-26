clear;
path(path,'../3_1_code(unsupervise)/baselines/NMF/old_result'); 
path(path,'../2_useful_data/human/');
path(path,'../3_1_code(unsupervise)/common_tool_functions');
load('NMF_hsa_result_pathway_20161103T105114.mat','NMF_result_cell');
result_cell = NMF_result_cell{2,1};
W_out = result_cell{2,1};    
load('g_p_network_2016_2.mat','gene_id');
load('hsa_pathway_data_2016_9.mat','pathway_new');
load('hsa_pathway_data_2016_2&9.mat','hsa_pathway_name');
load('g_p_network_2016_2.mat','g_p_network_first','g_p_network_second');
load('hsa_pathway_2016_0209_difference.mat','hsa_pathway_2016_0209_difference');
Z_threshold = 2.5;
[W_Z_filtered,predicted_pathway_gene] = predicted_pathway(W_out,Z_threshold,gene_id);
cluster_num=size(W_out,2);
initial_pathway_new=pathway_new;

initial_gene_num=length(gene_id);
all_genes=unique(pathway_new(:));
delete_index=[];
for i=initial_gene_num:-1:1
    if ~ismember(gene_id(i),all_genes)
        predicted_pathway_gene(predicted_pathway_gene==gene_id(i))=0;
        delete_index=[delete_index,i];
    end
end
gene_id(delete_index)=[];
W_Z_filtered(delete_index,:)=[];
g_p_network_first(delete_index,:)=[];
g_p_network_second(delete_index,:)=[];
for i=1:cluster_num
    arranged_genes=unique(predicted_pathway_gene(:,i));
    arranged_genes(1)=[];
    predicted_pathway_gene(1:length(arranged_genes),i)=arranged_genes;
    predicted_pathway_gene(length(arranged_genes)+1:end,i)=0;
end
predicted_pathway_gene(sum(predicted_pathway_gene,2)==0,:)=[];


for i=1:cluster_num
    common_genes=intersect(pathway_new(i,:),gene_id);
    gene_num=length(common_genes);
    pathway_new(i,1:gene_num)=common_genes;
    pathway_new(i,gene_num+1:end)=0;
end

%common_gene_num_matrix---the num of common genes of each true pathway and predicted pathway
%common_gene_rate_matrix---the rate of common genes of each true pathway and predicted pathway
common_gene_num_matrix = zeros(cluster_num,cluster_num);
jaccard_matrix=zeros(cluster_num,cluster_num);
for i = 1:cluster_num
    for j = 1:cluster_num
        new_genes = pathway_new(i,pathway_new(i,:)>0);
        predicted_genes = gene_id(W_Z_filtered(:,j)>0);
        common_genes = intersect(new_genes,predicted_genes);
        all_genes=union(new_genes,predicted_genes);
        common_gene_num_matrix(i,j) = length(common_genes);
        jaccard_matrix(i,j)=length(common_genes)/length(all_genes);
    end
end
jaccard_matrix(isnan(jaccard_matrix))=0;

%pathway_predicted---predicted pathway
%cluster_distribution---corresponding cluster index of each pathway
pathway_predicted=zeros(cluster_num,max(sum(W_Z_filtered)));
cluster_distribution=zeros(1,cluster_num);
matrix_copy=jaccard_matrix;
exist_index=1:cluster_num;
for i=1:cluster_num
    [~,index]=max(matrix_copy(i,:));
    cluster_index=exist_index(index);
    predicted_genes=predicted_pathway_gene(:,cluster_index);
    pathway_predicted(i,1:length(predicted_genes))=predicted_genes;
    cluster_distribution(i)=cluster_index;
    matrix_copy(:,index)=[];
    exist_index(index)=[];
end

%接下来根据pathway_predicted_matrix来进行分析
%先是看交集最大的
intersect_num = max(max(common_gene_num_matrix));
[old_cluster,predict_cluster] = find(common_gene_num_matrix == intersect_num);
%发现是263簇，pathways in cancer这个pathway，没什么分析意义

%num of genes of each new pathway
gene_num_new_pathway=zeros(1,cluster_num);
for i = 1 : cluster_num
    gene_num_new_pathway(i)=nnz(pathway_new(i,:));
end
avg_gene_num_new_pathway=mean(gene_num_new_pathway);
std_gene_num_new_pathway=std(gene_num_new_pathway);

%num of genes of each predicted pathway
gene_num_predicted_pathway=zeros(1,cluster_num);
for i = 1 : cluster_num
    gene_num_predicted_pathway(i)=nnz(pathway_predicted(i,:));
end
avg_gene_num_predicted_pathway=mean(gene_num_predicted_pathway);
std_gene_num_predicted_pathway=std(gene_num_predicted_pathway);

%prediction recall of each new pathway
recall_vector = zeros(1,cluster_num);
precision_vector=zeros(1,cluster_num);
f1_vector=zeros(1,cluster_num);
predicted_gene_num_sum=0;
common_gene_num_sum=0;
new_gene_num_sum=0;
for i = 1 : cluster_num
   predicted_gene_num=nnz(pathway_predicted(i,:));
   new_gene_num=nnz(pathway_new(i,:));
   common_gene_num=nnz(intersect(pathway_predicted(i,:),pathway_new(i,:)));
   recall_vector(i)=common_gene_num/new_gene_num;
   precision_vector(i)=common_gene_num/predicted_gene_num;
   f1_vector(i)=2*recall_vector(i)*precision_vector(i)/(recall_vector(i)+precision_vector(i));
   predicted_gene_num_sum=predicted_gene_num_sum+predicted_gene_num;
   common_gene_num_sum=common_gene_num_sum+common_gene_num;
   new_gene_num_sum=new_gene_num_sum+new_gene_num;
end
recall_vector(isnan(recall_vector))=0;
precision_vector(isnan(precision_vector))=0;
f1_vector(isnan(f1_vector))=0;
avg_recall=mean(recall_vector);
avg_precision=mean(precision_vector);
avg_f1=mean(f1_vector);
std_recall=std(recall_vector);
std_precision=std(precision_vector);
std_f1=std(f1_vector);
global_avg_recall=common_gene_num_sum/new_gene_num_sum;
global_avg_precision=common_gene_num_sum/predicted_gene_num_sum;
global_avg_f1=2*global_avg_recall*global_avg_precision/(global_avg_recall+global_avg_precision);

%useful_predict---pathway whose prediction precision is more than 0.3
%useful_predict_genes---num of truely predicted genes of each useful pathway
[~,useful_predict] = find(recall_vector>=0.3); %old pathway重合高
useful_predict_genes = zeros(length(useful_predict),1);
for i = 1 : length(useful_predict)
    vector = pathway_new(useful_predict(i),:);
    useful_predict_genes(i) = sum(vector~=0);
end
%275pathway 甲状腺病符合度高，0.33，可用进一步看

%通过前后数据差距，可以发现134,137簇多了很多新数据
%recall_new_gene_pathway134---prediction recall of new genes of 134th pathway
%recall_new_gene_pathway137---prediction recall of new genes of 137th pathway
added_genes = unique(initial_pathway_new(134,:).*hsa_pathway_2016_0209_difference(134,:));
added_genes(:,1) = [];
predicted_genes=pathway_predicted(134,:);
common_genes=intersect(added_genes,predicted_genes);
recall_new_gene_pathway134=nnz(common_genes)/nnz(added_genes);

added_genes = unique(initial_pathway_new(137,:).*hsa_pathway_2016_0209_difference(137,:));
added_genes(:,1) = [];
predicted_genes=pathway_predicted(137,:);
common_genes=intersect(added_genes,predicted_genes);
recall_new_gene_pathway137=nnz(common_genes)/nnz(added_genes);

%recall of new genes of each new pathway
new_gene_recall_vector = zeros(cluster_num,1);
for i = 1 : cluster_num
    added_genes = unique(initial_pathway_new(i,:).*hsa_pathway_2016_0209_difference(i,:));
    added_genes(:,1) = [];
    predicted_genes = pathway_predicted(i,:);
    common_genes = intersect(added_genes,predicted_genes);
    new_gene_recall_vector(i) = nnz(common_genes)/nnz(added_genes);
end
new_gene_recall_vector(isnan(new_gene_recall_vector))=0;
avg_new_gene_recall=mean(new_gene_recall_vector);
std_new_gene_recall=std(new_gene_recall_vector);

%接下来看跟父层子层HPO之间的情况
%father_level_intersect---num of father level phenetype related to each new pathway
%son_level_intersect---num of son level phenetype related to each new pathway
pathway_new_father_level_intersect = zeros(cluster_num,1);
pathway_new_son_level_intersect = zeros(cluster_num,1);
for i = 1 : cluster_num
    new_genes = pathway_new(i,:);
    [~,~,gene_index] = intersect(new_genes,gene_id);
    pathway_new_father_level_intersect(i) = nnz(sum(g_p_network_first(gene_index,:)));
    pathway_new_son_level_intersect(i) = nnz(sum(g_p_network_second(gene_index,:)));
end
avg_pathway_new_father_level_intersect=mean(pathway_new_father_level_intersect);
avg_pathway_new_son_level_intersect=mean(pathway_new_son_level_intersect);

%CMNMF_father_level_intersect---num of father level phenetype related to each predicted pathway
%CMNMF_son_level_intersect---num of son level phenetype related to each predicted pathway
pathway_predicted_father_level_intersect = zeros(cluster_num,1);
pathway_predicted_son_level_intersect = zeros(cluster_num,1);
for i = 1 : cluster_num
    predicted_genes = pathway_predicted(i,:);
    [~,~,gene_index] = intersect(predicted_genes,gene_id);
    pathway_predicted_father_level_intersect(i) = nnz(sum(g_p_network_first(gene_index,:)));
    pathway_predicted_son_level_intersect(i) = nnz(sum(g_p_network_second(gene_index,:)));
end
avg_pathway_predicted_father_level_intersect=mean(pathway_predicted_father_level_intersect);
avg_pathway_predicted_son_level_intersect=mean(pathway_predicted_son_level_intersect);

%num of covered pathway of all genes of each pathway
intersect_cluster_vector = zeros(cluster_num,1);
for i = 1 : cluster_num
   vector = pathway_new(i,:);
   vector = vector(vector~=0);
   clusters = [];
   for j = 1 : length(vector)
       gene = vector(j);
       [row,~] = find(pathway_new == gene);
       clusters = unique([clusters;row]);
   end
    intersect_cluster_vector(i) = length(clusters);
end

evaluation_result={
    'Z_threshold',Z_threshold;
    'pathway_new',pathway_new;
    'pathway_predicted',pathway_predicted;
    'cluster_distribution',cluster_distribution;
    'common_gene_num_matrix',common_gene_num_matrix;
    'jaccard_matrix',jaccard_matrix;
    'gene_num_new_pathway',gene_num_new_pathway;
    'avg_gene_num_new_pathway',avg_gene_num_new_pathway;
    'std_gene_num_new_pathway',std_gene_num_new_pathway;
    'gene_num_predicted_pathway',gene_num_predicted_pathway;
    'avg_gene_num_predicted_pathway',avg_gene_num_predicted_pathway;
    'std_gene_num_predicted_pathway',std_gene_num_predicted_pathway;
    'f1_vector',f1_vector;
    'recall_vector',recall_vector;
    'precision_vecctor',precision_vector;
    'avg_f1',avg_f1;
    'avg_recall',avg_recall;
    'avg_precision',avg_precision;
    'std_f1',std_f1;
    'std_recall',std_recall;
    'std_precision',std_precision;
    'global_avg_f1',global_avg_f1;
    'global_avg_recall',global_avg_recall;
    'global_avg_precision',global_avg_precision;
    'useful_predict',useful_predict;
    'useful_predict_genes',useful_predict_genes;
    'recall_new_gene_pathway134',recall_new_gene_pathway134;
    'recall_new_gene_pathway137',recall_new_gene_pathway137;
    'new_gene_recall_vector',new_gene_recall_vector;
    'avg_new_gene_recall',avg_new_gene_recall;
    'std_new_gene_recall',std_new_gene_recall;
    'pathway_new_father_level_intersect',pathway_new_father_level_intersect;
    'pathway_new_son_level_intersect',pathway_new_son_level_intersect;
    'avg_pathway_new_father_level_intersect',avg_pathway_new_father_level_intersect;
    'avg_pathway_new_son_level_intersect',avg_pathway_new_son_level_intersect;
    'pathway_predicted_father_level_intersect',pathway_predicted_father_level_intersect;
    'pathway_predicted_son_level_intersect',pathway_predicted_son_level_intersect;
    'avg_pathway_predicted_father_level_intersect',avg_pathway_predicted_father_level_intersect;
    'avg_pathway_predicted_son_level_intersect',avg_pathway_predicted_son_level_intersect;
    };

result_file_name=['NMF.mat' ];
save(result_file_name,'evaluation_result');