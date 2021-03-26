function [ CMNMF_result_cell ] = CMNMF_main()
%CMNMF_MAIN Summary of this function goes here
%   Detailed explanation goes here   
    method_dir = 'CMNMF/';
    species = 'mouse';    
    validation_groundTruth_str = 'pathway';%two choices: 'pathway' or 'ppi'
    %alpha_set = 10.^(-5:-5);           % pathway
    %beta_set  = 10.^(-1:-1);
    alpha_set = 0.0001;          % ppi[0.001, 0.01, 0.1, 1, 10, 100, 1000]
    beta_set  = 0.001;
    %dimension_K_set = [20];
    max_ites = 80;%80 is a proper number
    top_n_set = [200;600;1000];
    cv_criteria = 'micro_f';    
    file_num=1;
    normalization=7;     % 7 1 
    parameter_cell = {method_dir;species;alpha_set;
                        beta_set;max_ites;top_n_set;
                        cv_criteria;validation_groundTruth_str;file_num;normalization};
    path(path,method_dir);
    path(path,'common_tool_functions');    
    CMNMF_result_cell = CMNMF(parameter_cell);
    rmpath(method_dir);
    rmpath('common_tool_functions');
end

