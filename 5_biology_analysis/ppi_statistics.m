%����ppi������protein������protein�������Խ������ͳ��
clear;
path(path,'../2_useful_data/human/');
%load·���µ��������
load('hsa_ppi_2016_9.mat','ppi_new');
x=sum(ppi_new,2);
y=unique(x);
for i=1:length(y)
    a(i)=sum(x==y(i));
end
y1=full(y);
a1=full(a);
result(1,:)=y1(:,1);
result(2,:)=a1(1,:);
write_infile('a.txt',result);
save('ppi_sparsity.mat','result');

