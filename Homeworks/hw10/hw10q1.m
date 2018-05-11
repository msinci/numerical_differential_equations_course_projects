% Mehmet Sinan INCI 
% MA 512 - HW 10 Q1
clc;
clear all;  % clear variables
close all;  % close old figures
format long; 

m=[32 64 128];it_method_demo(m,'Jacobi',@rhsfun);
m=[32 64 128];it_method_demo(m,'Gauss-Seidel',@rhsfun);
m=[32 64 128];it_method_demo(m,'SOR',@rhsfun);



