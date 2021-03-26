clear;
path(path,'../3_1_code(unsupervise)/CMNMF/old_result'); 
path(path,'../2_useful_data/human/');
path(path,'../3_1_code(unsupervise)/common_tool_functions');
load('g_p_network_2016_2.mat','gene_id');
load('hsa_ppi_2016_9.mat','ppi_new','ppi_new_index');
load('CMNMF_hsa_result_ppi_20161027T090505.mat','CMNMF_result_cell');
result_cell = CMNMF_result_cell{2,1};
W_out = result_cell{2,1};   
Z_threshold = 11;
[W_Z_filtered,predicted_pathway_gene] = predicted_pathway(W_out,Z_threshold,gene_id);
predict_GG_matrix = W_Z_filtered*W_Z_filtered';
predict_GG_matrix(predict_GG_matrix>0) = 1;
[~,ia,ib] = intersect(gene_id,ppi_new_index);
g_g_predict = predict_GG_matrix(ia,ia);
for i = 1 : size(g_g_predict,1)
   g_g_predict(i,i) = 0; 
end
g_g_standard = ppi_new(ib,ib);
g_g_standard = full(g_g_standard);
TP = nnz(g_g_predict+g_g_standard==2);
TN = nnz(g_g_predict+g_g_standard==0);
FP = nnz(g_g_predict-g_g_standard==1);
FN = nnz(g_g_predict-g_g_standard==-1);

RD = (TP+TN)/(TP+TN+FP+FN);

Precision = TP/(TP+FP);

Recall = TP/(TP+FN);

F=2*Recall*Precision/(Recall+Precision);

jaccard = TP/nnz(g_g_predict+g_g_standard);
%如果把score值升到11，最后预测的和标准的总量接近，最后除了召回率都提高了许多


%现在是把gene分为296类，如果对于每一类gene，在与g_p_network相交时对应的phenotype
%越少，那么分类效果更好。
for i = 1 : 296
    vector = predicted_pathway_gene(:,i);
    vector = vector(vector~=0);
    
end

