function [CMNMF_result_cell] = CMNMF(parameter_cell)    
    %parameter_cell = {method_name};
    %%%%%%%%%%%%%%%parse parameter%%%%%%%%%%%%%%%%%%%%%%%%%%%
    method_dir = parameter_cell{1,1};
    species = parameter_cell{2,1};     
    alpha_set = parameter_cell{3,1};
    beta_set  = parameter_cell{4,1};    
    max_ites = parameter_cell{5,1};
    top_n_set = parameter_cell{6,1};
    cv_criteria = parameter_cell{7,1};
    validation_groundTruth_str = parameter_cell{8,1};
    file_num=parameter_cell{9,1};
    normalization=parameter_cell{10,1};
    test_groundTruth_str = validation_groundTruth_str;
    %%%%%%%%%%%%%%%parse parameter%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(species,'human')
        species_acronym = 'hsa';        
    elseif strcmp(species,'mouse')
        species_acronym = 'mmu';       
    end
    %%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    data_dir = ['../2_useful_data/' species '/'];
    path(path,data_dir); %invoke rmpath() to remove added path.
    load([data_dir 'g_p_network_2016_9.mat'],'g_p_network_first','g_p_network_second','M','gene_id');      
    load([data_dir [species_acronym '_pathway_data_2016_9_9.mat']], 'pathway_old');
    load([data_dir [species_acronym '_pathway_data_2016_9.mat']], 'pathway_new');
    load([data_dir [species_acronym '_ppi_2016_9_9.mat']], 'ppi_old','ppi_old_index'); 
    % ����ע��
    load([data_dir [species_acronym '_ppi_201609_202011_difference.mat']], 'ppi_difference','ppi_new_index');  
    ppi_new = ppi_difference;
    %load([data_dir [species_acronym '_ppi_2016_9.mat']], 'ppi_new','ppi_new_index');    
    %%%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V1 = g_p_network_first;
    V2 = g_p_network_second;
    [total_gene_num, total_p1_num] = size(V1);
    [~,total_p2_num] = size(V2);     
    [dimension_K_set,~] = size(pathway_old);
    %%%%%%%%%%%%%%%%%%%%%%%initialize matrix with random value%%%%%%%%%%%%%%%%%%%
    dir = [method_dir 'data/' species '/'];
    name_prefix = 'initial_matrix_fixed';
    get_files_parameter_cell = {file_num;dir;name_prefix;total_gene_num;...,
        dimension_K_set(1);total_p1_num;total_p2_num};
    file_name_cell = get_files(get_files_parameter_cell);
    
    %%%%%%%%%%%%%%%%%%%%%%%%initialize matrix with random value%%%%%%%%%%%%%%%%%%%
    gene_idx = gene_id;
    matrix_cell_train = {V1; V2; M};
    initial_matrixFileName_cell = file_name_cell;

    %%%%%%%%%%%%%%%%%%%%%%%%%validation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    [groundTruth_pathway_old, gene_idx_pathway] = GeneGene_GroundTruth_Pathway(pathway_old);
    matrix_cell_validation  = {gene_idx;groundTruth_pathway_old;gene_idx_pathway;ppi_old;ppi_old_index;validation_groundTruth_str;pathway_old};
    input_parameter_cell = {alpha_set; beta_set; max_ites; top_n_set; cv_criteria; method_dir; species;normalization};
    %��ѵ����fold���������㷨����������������ָ��c��w1��w2Ϊ�������飬cΪ��Ӧ����
    tic;
    [learned_matrix_cell,best_parameter_array,evaluation_result,loss] = learn(input_parameter_cell, matrix_cell_train, matrix_cell_validation, ...,
        initial_matrixFileName_cell);
    toc;
    %%%%%%%%%%%%%%%%%%%%%%%%%validation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %�����µ�pathway ��PPI����Ԥ��      
    %%%%%%%%%test%%%%%%%%%%%%%
    tic; 
    [groundTruth_pathway_new, gene_idx_pathway] = GeneGene_GroundTruth_Pathway(pathway_new);
    matrix_cell_test =  {gene_idx;groundTruth_pathway_new;gene_idx_pathway;ppi_new;ppi_new_index;test_groundTruth_str;pathway_new};
    test_parameter_cell = {best_parameter_array(1);best_parameter_array(2); dimension_K_set; max_ites ; top_n_set};
    evaluation_test = test(matrix_cell_test,learned_matrix_cell,test_parameter_cell);
    CMNMF_result_cell = {loss;evaluation_test;learned_matrix_cell;best_parameter_array;evaluation_result;max_ites;alpha_set;beta_set;normalization};
    toc;
    %%%%%%%%%test%%%%%%%%%%%%%%
    
    if strcmp(test_groundTruth_str,'ppi') == 1
        result_file_name = [method_dir 'old_result/CMNMF_' species_acronym '_result_ppi_difference_' datestr(now,30) '.mat' ];
        % ����ע��
        %result_file_name = [method_dir 'old_result/CMNMF_' species_acronym '_result_ppi_' datestr(now,30) '.mat' ];  
        save(result_file_name,'CMNMF_result_cell');
    elseif strcmp(test_groundTruth_str,'pathway') == 1
        result_file_name = [method_dir 'old_result/CMNMF_' species_acronym '_result_pathway_' datestr(now,30) '.mat' ];  
        save(result_file_name,'CMNMF_result_cell');
    end  
    
    rmpath(data_dir);
end

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

function [file_name_cell] = get_files(get_files_parameter_cell)
    %file_name_cell = get_files(10,'data','initial_matrix_fixed_2016');
    file_num_need = get_files_parameter_cell{1,1};
    dir_path = get_files_parameter_cell{2,1};
    name_prefix = get_files_parameter_cell{3,1};
    total_gene_num = get_files_parameter_cell{4,1};
    dimension_K_set(1) = get_files_parameter_cell{5,1};
    total_p1_num = get_files_parameter_cell{6,1};
    total_p2_num = get_files_parameter_cell{7,1};
    
    fileFolder=fullfile(dir_path);
    dirOutput=dir(fullfile(fileFolder,'*.mat'));
    fileNames={dirOutput.name}';
    count =0;
    for i =1:length(fileNames)
       if(strncmp(fileNames(i),name_prefix,7) == 1)
           count = count+1;
           file_name_cell{count,1} = fileNames{i,1};
       end
    end
    if count > file_num_need
       file_name_cell_temp = file_name_cell(end-file_num_need+1:end,:); 
       file_name_cell = file_name_cell_temp;
    elseif count < file_num_need
        count_temp = count;
        for i = 1:file_num_need - count
            file_name_new = ['initial_matrix_fixed_' datestr(now,30) '.mat' ];            
            W = rand(total_gene_num, dimension_K_set(1));
            H1 = rand(dimension_K_set(1), total_p1_num);
            H2 = rand(dimension_K_set(1), total_p2_num);    
            save(file_name_new,'W','H1','H2');
            pause(0.5);%pause 1 seconds
            count_temp = count_temp + 1;
            file_name_cell{count_temp,1} = file_name_new;
        end
    end
end