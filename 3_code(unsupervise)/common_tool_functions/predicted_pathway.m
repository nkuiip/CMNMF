function [Z_filtered,predicted_pathway_gene] = predicted_pathway(W,K,mgi_id)
%����Z-score
z_threshold =K;
Z=zscore(W');
%������С��threshold��ֵ��0
Z(Z<z_threshold)=0;
Z =Z';

%�۴�
[m,n]=size(Z);
%k(ÿ������ÿ��ɸѡ�ĸ�����
k = zeros(m,1);
%Z_filteredΪԤ�����gene-pathway��ϵ����Ϊ0,1����
Z_filtered=zeros(m,n);
Z(Z>0) = 1;
Z_filtered = Z;
    
%predicted_pathway_mgi_idÿһ��Ϊһ��pathway,���а�����pathway�л����mgi_id
%��Z_filter����ת��Ϊpredicted_pathway_mgi_id����

max_pathway_gene = max(sum(Z_filtered));
gene_location = zeros(max_pathway_gene,n);
predicted_pathway_gene = zeros(max_pathway_gene,n);
     
   for i = 1:n
         gene_location = find(Z_filtered(:,i)==1);
         gene_mgi_id = mgi_id(gene_location);
         predicted_pathway_gene(1:length(gene_mgi_id),i)=gene_mgi_id;
   end
   

     
