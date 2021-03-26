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
pathway_predicted_matrix = zeros(25,size(W_Z_filtered,2));

%pathway 210 214 215 230 231 232 233 234 260:281
W_Z_selected = zeros(30,size(pathway_old,2));
W_Z_selected(1,:) = pathway_old(210,:);
W_Z_selected(2,:) = pathway_old(214,:);
W_Z_selected(3,:) = pathway_old(215,:);
for i = 4:8
    W_Z_selected(i,:) = pathway_old(226+i,:);
end
for i = 9:30
    W_Z_selected(i,:) = pathway_old(251+i,:);
end
%W_Z_selected = pathway_old;
for i = 1:25
    for j = 1:size(W_out,2)
        genes = W_Z_selected(i,W_Z_selected(i,:)>0);
        genes2 = gene_id(W_Z_filtered(:,j)>0);
        [C,~,~] = intersect(genes,genes2);
        pathway_predicted_matrix(i,j) = length(C);
    end
end
rmpath('../3_1_code(unsupervise)/CMNMF'); 
rmpath('../2_useful_data/mouse/');

