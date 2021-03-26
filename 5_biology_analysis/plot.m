clear;
clc;
load('precision_vector.mat','CMNMF_precision','NMF_precision');

X = CMNMF_precision;
X(isnan(X)) = [];
Y = NMF_precision;
Y(isnan(Y)) = [];
scatter(X,Y,'.');
