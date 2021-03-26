function [ z_threshold ] = get_best_z( W,gene_id,pathway )
pair_num=0;
all_genes=unique(pathway(:));
delete_index=[];
for i=1:length(gene_id)
    if ~ismember(gene_id(i),all_genes)
        delete_index=[delete_index,i];
    end
end
gene_id(delete_index)=[];
W(delete_index,:)=[];
for i=1:size(pathway,1)
    pair_num=pair_num+length(intersect(pathway(i,:),gene_id));
end
W=zscore(W')';
left=0;
right=6;
for i=1:10
    mid=(left+right)/2;
    n=nnz(W>mid);
    if n>pair_num
        left=(left+right)/2;
    else
        right=(left+right)/2;
    end
end
z_threshold=mid;
end

