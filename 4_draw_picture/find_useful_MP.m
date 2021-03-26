clear;
path(path,'../2_useful_data/mouse');
load('mmu_pathway_data_2016_2&9.mat');
load('g_p_network_2016_2.mat');
[r1,~] = find(mmu_pathway_id==4930);
[r2,~] = find(mmu_pathway_id==5215);
[r3,~] = find(mmu_pathway_id==5205);
pathway1_ids = mmu_pathway_mgi_id_2016_2(r1,mmu_pathway_mgi_id_2016_2(r1,:)>0);
pathway2_ids = mmu_pathway_mgi_id_2016_2(r2,mmu_pathway_mgi_id_2016_2(r2,:)>0);
pathway3_ids = mmu_pathway_mgi_id_2016_2(r3,mmu_pathway_mgi_id_2016_2(r3,:)>0);
[C1,~,i1] = intersect(pathway1_ids,gene_id);
[C2,~,i2] = intersect(pathway2_ids,gene_id);
[C3,~,i3] = intersect(pathway3_ids,gene_id);
%C为对应pathway的gene_id和我们选定的gene_id的公共部分
gene_MP_matrix1 = zeros(length(C1),60);
gene_MP_matrix2 = zeros(length(C2),120);
gene_MP_matrix3 = zeros(length(C3),65);
for i = 1:length(C1)
    row = i1(i);
    array = first_and_second_idx(g_p_network(row,:)==1,1);
    gene_MP_matrix1(i,1:length(array)) = array;

end
for i = 1:length(C2)
    row = i2(i);
    array = first_and_second_idx(g_p_network(row,:)==1,1);
    gene_MP_matrix2(i,1:length(array)) = array;

end
for i = 1:length(C3)
    row = i3(i);
    array = first_and_second_idx(g_p_network(row,:)==1,1);
    gene_MP_matrix3(i,1:length(array)) = array;

end
%以上把每个gene的所有相关MP找到

gene_MP_matrix1_temp = gene_MP_matrix1;
gene_MP_matrix1_temp(gene_MP_matrix1_temp>0) = 1;
array1 = sum(gene_MP_matrix1_temp,2);
[~,I] = sort(array1,'descend');
selected_MP_matrix1 = zeros(3,60);
selected_genes1 = [C1(15);C1(3);C1(14)];
selected_MP_matrix1 = [gene_MP_matrix1(15,:);gene_MP_matrix1(3,:);gene_MP_matrix1(14,:)];

gene_MP_matrix2_temp = gene_MP_matrix2;
gene_MP_matrix2_temp(gene_MP_matrix2_temp>0) = 1;
array1 = sum(gene_MP_matrix2_temp,2);
[~,I] = sort(array1,'descend');
selected_MP_matrix2 = zeros(3,120);
selected_genes2 = [C2(23);C2(19);C2(4)];
selected_MP_matrix2 = [gene_MP_matrix2(23,:);gene_MP_matrix2(19,:);gene_MP_matrix2(4,:)];

gene_MP_matrix3_temp = gene_MP_matrix3;
gene_MP_matrix3_temp(gene_MP_matrix3_temp>0) = 1;
array1 = sum(gene_MP_matrix3_temp,2);
[~,I] = sort(array1,'descend');
selected_MP_matrix3 = zeros(3,65);
selected_genes3 = [C3(14);C3(35);C3(24)];
selected_MP_matrix3 = [gene_MP_matrix3(14,:);gene_MP_matrix3(35,:);gene_MP_matrix3(2,:)];

%选出3组gene_id各3个，并列出所有相关MP。


MPO_id = unique(gene_MP_matrix1);
MPO_id(MPO_id==0) = [];
[father_MPO_id,~,i1] = intersect(MPO_id,first_level_id);
[son_MPO_id,~,i2] = intersect(MPO_id,second_level_id);
father_son_matrix = zeros(length(father_MPO_id),length(son_MPO_id));
for i = 1:length(father_MPO_id)
    for j = 1:length(son_MPO_id)
        father_son_matrix(i,j) = M(i1(i),i2(j));
    end
end
[r,c,~] = find(father_son_matrix==1);
selected_MP_father = unique(father_MPO_id(r));
selected_MP_son = unique(son_MPO_id(c));
%找到了有父子关系的MP
a  = zeros(0,0);
for i = 1:length(selected_MP_father)
   id =selected_MP_father(i);
   [r1,~] = find(selected_MP_matrix1==id);
   [r2,~] = find(selected_MP_matrix2==id);
   [r3,~] = find(selected_MP_matrix3==id);
   if(length(r1)+length(r2)+length(r3)==1)
       a(length(a)+1,1) = id;
   end
end
a(a==0) = [];


rmpath('../2_useful_data/mouse');