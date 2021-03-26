function [groundTruth_pathway, gene_idx_pathway] = GeneGene_GroundTruth_Pathway(pathway_hsa_ncbi)
%pathway_standard:��2016��2�·ݵ�pathway���ɻ���-�����Ƿ���ͬһpathway��
%���֣�����������Ϊ��֤����ѡ�������
    gene_idx_pathway = unique(pathway_hsa_ncbi);
    pathway_ncbi = zeros(length(gene_idx_pathway),size(pathway_hsa_ncbi,1));
    for i = 1:size(pathway_hsa_ncbi,1)
        [~,~,ib] = intersect(pathway_hsa_ncbi(i,:),gene_idx_pathway);
        pathway_ncbi(ib,i) = 1;
    end
    groundTruth_pathway = pathway_ncbi*pathway_ncbi';
end