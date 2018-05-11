% Mehmet Sinan INCI 
% MA 512 - HW 11 Q1
clc;
clear all;  % clear variables
close all;  % close old figures
format long; 


pcg_demo([32 64 128],'none','Poisson',@rhsfun);
pcg_demo([32 64 128],'block Jacobi','Poisson',@rhsfun);
pcg_demo([32 64 128],'SSOR','Poisson',@rhsfun);
pcg_demo([32 64 128],'incomplete Cholesky','Poisson',@rhsfun);
pcg_demo([32 64 128],'ADI','Poisson',@rhsfun);

