clear;
path(path,'../3_1_code(unsupervise)/CMNMF/old_result'); 
path(path,'../2_useful_data/human/');
path(path,'../3_1_code(unsupervise)/common_tool_functions');
load('CMNMF_hsa_result_pathway_20161029T094839.mat','CMNMF_result_cell');
result_cell = CMNMF_result_cell{2,1};
W_out = result_cell{2,1};    
load('g_p_network_2016_2.mat','gene_id','M','g_p_network_first','g_p_network_second');
load('hsa_pathway_data_2016_9.mat','pathway_new');

vector = sum(M);
[~,more_father] = find(vector >= 2);
g_p_analyse = g_p_network_second(:,more_father);
vector = sum(g_p_analyse,2);
[row,~] = find(vector ~= 0);
%gene id of genes whose ralated phenetype has more than one father
more_father_genes = gene_id(row);

Z_threshold = 3;
[W_Z_filtered,predicted_pathway_gene] = predicted_pathway(W_out,Z_threshold,gene_id);
%the num of concordinate genes of each true pathway and predicted pathway
pathway_predicted_matrix = zeros(size(W_out,2),size(W_out,2));
for i = 1:size(W_out,2)
    for j = 1:size(W_out,2)
        genes = pathway_new(i,pathway_new(i,:)>0);
        genes2 = gene_id(W_Z_filtered(:,j)>0);
        [C,~,~] = intersect(genes,genes2);
        pathway_predicted_matrix(i,j) = length(C);
    end
end

%the num of truely predicted more-father-genes of each pathway
more_father_HPO_vector = zeros(296,1);
for i = 1 : size(W_out,2)
   [num,I] = max(pathway_predicted_matrix(i,:)); 
   predict_genes = unique(predicted_pathway_gene(:,I));
   predict_genes(1,:) = [];
   [X,~,~] = intersect(unique(pathway_new(i,:)),predict_genes);
   old_num = length(X);
   [Y,~,~] = intersect(X,more_father_genes);
   more_father_HPO_vector(i) = length(Y);
end