% Mehmet Sinan INCI 
% MA 512 - HW 10 Q1
clc;
clear all;  % clear variables
close all;  % close old figures
format long; 

m=[32 64 128];gmres_demo(m,[10 10 10],@rhsfun,'none');
m=[32 64 128];gmres_demo(m,[10 10 10],@rhsfun,'block Jacobi');
m=[32 64 128];gmres_demo(m,[10 10 10],@rhsfun,'ILU');
m=[32 64 128];gmres_demo(m,[10 10 10],@rhsfun,'ADI');


