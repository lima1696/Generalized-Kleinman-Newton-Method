clear all
close all
clc

%System
load('HE1_dyn.mat');

%Controller Order
nc = 4;

%Options(1) Stop codition
%Options(2) Max number of iterations
%Options(3) Value of Tol (Z'Z < tol^2) 

options = [1e-2 5000 1e2];

[Jopt,Jfeas,Zopt,Zfeas] = GKN_Method(A,B,E,C,D,G,H,nc,options);