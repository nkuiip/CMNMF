function[] = draw_colormap(parameter_cell)

    ground_truth_str = parameter_cell{1,1};%ppi
    speices = parameter_cell{2,1}; 
    file_name  = parameter_cell{3,1};
    path(path,'../3_1_code(unsupervise)/CMNMF/');
    load(file_name,'CMNMF_result_cell');
    matrix = CMNMF_result_cell{5,1};         % CMNMF_result_cell{4,1}
    A=[matrix(1:7,1),matrix(8:14,1),matrix(15:21,1),matrix(22:28,1),matrix(29:35,1),...,
       matrix(36:42,1),matrix(43:49,1)];
    mat = A;           %# A n-by-n matrix of random values from 0 to 1
    save A.mat A
    imagesc(mat);            %# Create a colored plot of the matrix values
    colormap(flipud(gray));  %# Change the colormap to gray (so higher values are
                             %#   black and lower values are white)
    colorbar;
    xlabel('\alpha');
    ylabel('\beta');
    set(get(gca, 'YLabel'), 'Rotation', 0);
    set(gca,'XTick',1:7,...                         %# Change the axes tick marks
            'XTickLabel',{'0.001','0.01','0.1','1','10','100','1000'},...  %#   and tick labels
            'YTick',1:7,...
            'YTickLabel',{'0.001','0.01','0.1','1','10','100','1000'},...
            'TickLength',[0 0]);
    rmpath('../3_1_code(unsupervise)/CMNMF/');
    saveas(gca,[speices '_' ground_truth_str '_result'],'pdf');
end
