setwd('E:\\zhangyaogong\\Matlab_Workspace\\CMNMF_2016_10_20\\4_draw_picture');
library('plyr')
library('pheatmap')
source("draw.r")

win_width = 50
win_height =100
row_name <- c('Mafa','Ins2','Abcc8','Ikbkg','Nfkbia','Ctnnb1','Tlr4','Vegfa','Tgfb2')
Mp_name <- c('MP5215','MP6042','MP2703','MP5216','MP5217','MP6413','MP11195','MP4755','MP4756','MP9640')
col_name <- c('cluster e1','cluster e2','cluster e3')
col_name2 <- c('cluster d1','cluster d2','cluster d3')
col_name3 <- c('cluster c1','cluster c2','cluster c3')
V3=matrix(0,9,3)
V3[1,1] = 1
V3[2,1] = 1
V3[6,1] = 0.0057
V3[9,1] = 1
V3[4,2] = 1
V3[5,2] = 1
V3[6,3] = 0.9955
V3[7,3] = 1
V3[3,2] = 1
V3[8,3] = 1
t1=cbind(matrix(1,3,1),matrix(0,3,1),matrix(0,3,1))
t2=cbind(matrix(0,3,1),matrix(1,3,1),matrix(0,3,1))
t3=cbind(matrix(0,3,1),matrix(0,3,1),matrix(1,3,1))
V4 =rbind(t1,t2,t3)

V5=matrix(0,9,3)
V5[1,1] = 1
V5[2,1] = 1
V5[3,1] = 1
V5[6,1] = 0.9955
V5[4,2] = 1
V5[5,2] = 1
V5[6,3] = 0.0057
V5[7,3] = 1
V5[8,2] = 1
V5[9,3] = 1
RampPlette<-colorRampPalette(c('#ffffff','#9933ff'))(256)
color = RampPlette
cellheight = 15
cellwidth = 30
cluster_cols = 0
cluster_rows = 0
fontsize = 9
fontsize_row = 9
save_filename = 'v_d.pdf'
draw(V3, col_name2, row_name, win_width, win_height, color, 
	cellheight, cellwidth, cluster_cols, cluster_rows, fontsize,
	 fontsize_row, save_filename)
save_filename = 'v_e.pdf'
cluster_cols = 0
draw(V4, col_name, row_name, win_width, win_height, color, 
	cellheight, cellwidth, cluster_cols, cluster_rows, fontsize,
	 fontsize_row, save_filename)

save_filename = 'v_c.pdf'
draw(V5, col_name3, row_name, win_width, win_height, color, 
	cellheight, cellwidth, cluster_cols, cluster_rows, fontsize,
	 fontsize_row, save_filename)