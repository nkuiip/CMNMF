setwd('E:\\zhangyaogong\\Matlab_Workspace\\CMNMF_2016_10_20\\4_draw_picture');
library('pheatmap')
library('plyr')


win_width = 100
win_height =200
row_name <- c('Mafa','Ins2','Abcc8','Ikbkg','Nfkbia','Ctnnb1','Tlr4','Vegfa','Tgfb2')
V3 = matrix(0, 9, 10)
V3[1, 1] = 1
V3[4, 2] = 1
V3[7, 3] = 1
V3[2, 4] = 1
V3[3, 5] = 1
V3[5, 6] = 1    
V3[6, 7] = 1
V3[8, 8] = 1
V3[8, 9] = 1
V3[9, 10] = 1
Mp_name <- c('MP5215\nabnormal pancreatic islet morphology','MP6042\nincreased apoptosis','MP2703\nabnormal renal tubule morphology','MP5216\nabnormal pancreatic alpha cell morphology','MP5217\nabnormal pancreatic beta cell morphology','MP6413\nincreased T cell apoptosis','MP11195\nincreased hair follicle apoptosis','MP4755\nabnormal loop of Henle morphology','MP4756\nabnormal proximal convoluted tubule morphology','MP9640\nabnormal renal tubule epithelium morphology')
RampPlette<-colorRampPalette(c('#ffffff','#9933ff'))(256)
color = RampPlette
cellheight = 15
cellwidth = 30
cluster_cols = 0
cluster_rows = 0
fontsize = 8
fontsize_row = 9
save_filename = 'v_a.pdf'
cluster_cols = 0
draw(V3, Mp_name, row_name, win_width, win_height, color, 
	cellheight, cellwidth, cluster_cols, cluster_rows, fontsize,
	 fontsize_row, save_filename)

