ground_truth_str = 'ppi';
speices = 'mouse';
if strcmp(ground_truth_str,'ppi') ==1
   if strcmp(speices,'human') == 1
       file_name  = 'CMNMF_hsa_result_ppi_20201230T225412.mat';
       %file_name = 'CMNMF_hsa_result_ppi_20161027T090505.mat';
   else
       file_name  = 'CMNMF_mmu_result_ppi_20201231T002102.mat';
       %file_name = 'CMNMF_mmu_result_ppi_20161027T173946.mat';
   end
else
   if strcmp(speices,'human') == 1
       file_name  = 'CMNMF_hsa_result_pathway_20201230T235133.mat';
       %file_name  = 'CMNMF_hsa_result_pathway_20161029T094839.mat'; 
   else
       file_name = 'CMNMF_mmu_result_pathway_20201230T105416.mat';
       %file_name = 'CMNMF_mmu_result_pathway_20161109T094403.mat';
   end
end
parameter_cell = {ground_truth_str;speices;file_name};
draw_colormap(parameter_cell);