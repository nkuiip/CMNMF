function [learned_matrix_cell,best_parameter_array,evaluation_result,loss] = learn(input_parameter_cell, matrix_cell_train,...,
    matrix_cell_validation, initial_matrixFileName_cell)
%LEARN Summary of this function goes here
%   Detailed explanation goes here
%input_parameter_cell = {alpha_set; beta_set; max_ites; top_n_set; cv_criteria};
alpha_set = input_parameter_cell{1,1};
beta_set = input_parameter_cell{2,1};
max_ites = input_parameter_cell{3,1};
top_n_set = input_parameter_cell{4,1};
cv_criteria = input_parameter_cell{5,1};
method_dir = input_parameter_cell{6,1};
species = input_parameter_cell{7,1};
normalization=input_parameter_cell{8,1};
top_n_cv = top_n_set(end);
i=1;
%each row:[RD,F,Precision,Recall,jaccard]
%evaluation_result:每组参数组合下，十次随机初始化学习到十次结果的平均值
if strcmp(species,'human')
    evaluation_num=7;
elseif strcmp(species,'mouse')
    evaluation_num=6;
end
evaluation_result = zeros(length(alpha_set)*length(beta_set),evaluation_num+2);
file_num=length(initial_matrixFileName_cell);
arg_num=length(alpha_set)*length(beta_set);
evaluation_result_all=zeros(length(alpha_set)*length(beta_set),file_num,evaluation_num+2);
loss=cell(arg_num,file_num+2);
learned_matrix_cells=cell(arg_num,file_num);
for alpha = alpha_set
    for beta = beta_set
        disp(['arg_index= ' num2str(i)]);
        cv_parameter_cell = {alpha; beta; max_ites; top_n_cv; cv_criteria; method_dir; species;normalization};
        [evaluation_result(i,1:end-2),evaluation_result_all(i,:,1:end-2),loss(i,:),learned_matrix_cells(i,:)] =...
            cv_train(cv_parameter_cell, matrix_cell_train, matrix_cell_validation, initial_matrixFileName_cell);
        evaluation_result(i,end-1:end) = [alpha,beta];
        i = i+1;
    end
end

[arg_index,best_parameter_array] = get_best_parameter(evaluation_result,cv_criteria);

%获取到最优参数后，我们选择十次随机初始化结果中，性能最好的一次初始化结果最为最后评价用的结果
%cv_parameter_cell = {best_parameter_array(1); best_parameter_array(2); max_ites; top_n_cv; cv_criteria; method_dir; species};
%[~,evaluation_result_best_parameter] = cv_train(cv_parameter_cell, matrix_cell_train, matrix_cell_validation,...,
%            initial_matrixFileName_cell);
%evaluation_result = [RD,F,Precision,Recall,jaccard],F1 value is in the
%second column

%     [M,I]  = max(evaluation_result_all(arg_index,:,1));
%     method_data_dir = [method_dir 'data/' species '/'];
%     file_name = [method_data_dir initial_matrixFileName_cell{I,1}];
%     load(file_name);
%     initial_matrix_cell = {W;H1;H2};
%     [learned_matrix_cell] =  CMNMF_Train(cv_parameter_cell, matrix_cell_train, initial_matrix_cell);
[~,file_index]=max(evaluation_result_all(arg_index,:,1));
learned_matrix_cell=learned_matrix_cells{arg_index,file_index};
end
%cv_train 返回在一组参数组合情况下，十次随机初始化学习到的十次结果
function  [evaluation_result_average,evaluation_result,loss,learned_matrix_cells] = ...,
    cv_train(cv_parameter_cell, matrix_cell_train, matrix_cell_validation, initial_matrixFileName_cell)
train_parameter_cell = cv_parameter_cell;
method_dir = cv_parameter_cell{6,1};
species = cv_parameter_cell{7,1};
method_data_dir = [method_dir 'data/' species '/'];
%given a fixed parameters, return the average result of ten times' different initial matrix values
%initial_matrix_cell = {file_name_cell;W; H1; H2};
if strcmp(species,'human')
    evaluation_num=7;
elseif strcmp(species,'mouse')
    evaluation_num=6;
end
file_num=length(initial_matrixFileName_cell);
evaluation_result = zeros(file_num,evaluation_num);
loss=cell(1,file_num+2);
learned_matrix_cells=cell(file_num,1);
for i=1:file_num
    disp(['file_index= ' num2str(i)]);
    file_name = [method_data_dir initial_matrixFileName_cell{i,1}];
    load(file_name);
    initial_matrix_cell = {W;H1;H2};
    [learned_matrix_cell] = CMNMF_Train(train_parameter_cell, matrix_cell_train, initial_matrix_cell);
    evaluation_result(i,:) = CMNMF_Evaluate(learned_matrix_cell, matrix_cell_validation,train_parameter_cell);
    loss{1,i}=learned_matrix_cell{1,1};
    learned_matrix_cells{i,1}=learned_matrix_cell;
end
loss{1,file_num+1}=cv_parameter_cell{1,1};
loss{1,file_num+2}=cv_parameter_cell{2,1};
evaluation_result_average = mean(evaluation_result,1);
end

function[arg_index,best_parameter_array] = get_best_parameter(evaluation_result,cv_criteria)

if strcmp(cv_criteria,'micro_f') == 1
    [M,I]  = max(evaluation_result(:,1));
    best_parameter_array = evaluation_result(I,end-1:end);
    arg_index=I;
elseif strcmp(cv_criteria,'macro_f') == 1
    [M,I]  = max(evaluation_result(:,2));
    best_parameter_array = evaluation_result(I,end-1:end);
    arg_index=I;
elseif strcmp(cv_criteria,'micro_precision') == 1
    [M,I]  = max(evaluation_result(:,3));
    arg_index=I;
    best_parameter_array = evaluation_result(I,end-1:end);
elseif strcmp(cv_criteria,'macro_precision') == 1
    [M,I]  = max(evaluation_result(:,4));
    best_parameter_array = evaluation_result(I,end-1:end);
    arg_index=I;
elseif strcmp(cv_criteria,'micro_recall') == 1
    [M,I]  = max(evaluation_result(:,5));
    best_parameter_array = evaluation_result(I,end-1:end);
    arg_index=I;
elseif strcmp(cv_criteria,'macro_recall') == 1
    [M,I]  = max(evaluation_result(:,6));
    best_parameter_array = evaluation_result(I,end-1:end);
    arg_index=I;
end

end