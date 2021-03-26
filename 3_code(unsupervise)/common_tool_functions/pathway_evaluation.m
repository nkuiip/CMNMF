function [ micro_f,macro_f,micro_precision,macro_precision,micro_recall,macro_recall ] =...,
    pathway_evaluation(pathway,gene_id,W_Z_filtered,predicted_pathway_gene )

cluster_num=size(W_Z_filtered,2);

initial_gene_num=length(gene_id);
all_genes=unique(pathway(:));
delete_index=[];
for i=initial_gene_num:-1:1
    if ~ismember(gene_id(i),all_genes)
        predicted_pathway_gene(predicted_pathway_gene==gene_id(i))=0;
        delete_index=[delete_index,i];
    end
end
gene_id(delete_index)=[];
W_Z_filtered(delete_index,:)=[];
for i=1:cluster_num
    arranged_genes=unique(predicted_pathway_gene(:,i));
    arranged_genes(1)=[];
    predicted_pathway_gene(1:length(arranged_genes),i)=arranged_genes;
    predicted_pathway_gene(length(arranged_genes)+1:end,i)=0;
end
predicted_pathway_gene(sum(predicted_pathway_gene,2)==0,:)=[];

for i=1:cluster_num
    common_genes=intersect(pathway(i,:),gene_id);
    gene_num=length(common_genes);
    pathway(i,1:gene_num)=common_genes;
    pathway(i,gene_num+1:end)=0;
end

common_gene_num_matrix = zeros(cluster_num,cluster_num);
jaccard_matrix=zeros(cluster_num,cluster_num);
for i = 1:cluster_num
    for j = 1:cluster_num
        new_genes = pathway(i,pathway(i,:)>0);
        predicted_genes = gene_id(W_Z_filtered(:,j)>0);
        common_genes = intersect(new_genes,predicted_genes);
        all_genes=union(new_genes,predicted_genes);
        common_gene_num_matrix(i,j) = length(common_genes);
        jaccard_matrix(i,j)=length(common_genes)/length(all_genes);
    end
end
jaccard_matrix(isnan(jaccard_matrix))=0;

pathway_predicted=zeros(cluster_num,max(sum(W_Z_filtered)));
matrix_copy=jaccard_matrix;
exist_index=1:cluster_num;
for i=1:cluster_num
    [~,index]=max(matrix_copy(i,:));
    cluster_index=exist_index(index);
    predicted_genes=predicted_pathway_gene(:,cluster_index);
    pathway_predicted(i,1:length(predicted_genes))=predicted_genes;
    matrix_copy(:,index)=[];
    exist_index(index)=[];
end

recall_vector = zeros(1,cluster_num);
precision_vector=zeros(1,cluster_num);
f1_vector=zeros(1,cluster_num);
predicted_gene_num_sum=0;
common_gene_num_sum=0;
new_gene_num_sum=0;
for i = 1 : cluster_num
   predicted_gene_num=nnz(pathway_predicted(i,:));
   new_gene_num=nnz(pathway(i,:));
   common_gene_num=nnz(intersect(pathway_predicted(i,:),pathway(i,:)));
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
macro_recall=mean(recall_vector);
macro_precision=mean(precision_vector);
macro_f=mean(f1_vector);
micro_recall=common_gene_num_sum/new_gene_num_sum;
micro_precision=common_gene_num_sum/predicted_gene_num_sum;
micro_f=2*micro_recall*micro_precision/(micro_recall+micro_precision);

end

