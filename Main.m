clear all
close all
clc

%System
load('HE1_Dyn.mat');

%Controller Order
nc = 1;

%Options(1) Stop codition
%Options(2) Max number of iterations
options = [1e-2 500];

[Jopt,Jfeas,Zopt,Zfeas] = GKN_Method(A,B,E,C,D,G,H,nc,options);
