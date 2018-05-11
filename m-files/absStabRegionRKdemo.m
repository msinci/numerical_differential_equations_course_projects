%
% Plots absolute-stability regions for explicit RK methods of orders 1-8. 
% 
% Requires abs_stab_region.m. 
%

p = 6;

range = [-5 3 -5 5];

figure(1); clf; hold on; 

for k = 1:8
    switch k
        case 1
            range = [-4 2 -4 4];
            g = @(Z)1+Z;
            abs_stab_region(g,range); pause;
        case 2
            range = [-4 2 -4 4];
            g = @(Z)1+Z+2\Z.^2;
            abs_stab_region(g,range); pause;
        case 3 
            range = [-4 2 -4 4];
            g = @(Z)1+Z+2\Z.^2+6\Z.^3;
            abs_stab_region(g,range); pause;
        case 4
            range = [-4 2 -4 4];
            g = @(Z)1+Z+2\Z.^2+6\Z.^3+24\Z.^4;
            abs_stab_region(g,range); pause;
        case 5
            range = [-5 2 -5 5];
            g = @(Z)1+Z+2\Z.^2+6\Z.^3+24\Z.^4+120\Z.^5;
            clf; abs_stab_region(g,range); pause;
        case 6
            range = [-5 2 -5 5];
            g = @(Z)1+Z+2\Z.^2+6\Z.^3+24\Z.^4+120\Z.^5+720\Z.^6;
            clf; abs_stab_region(g,range); pause;
        case 7
            range = [-5 2 -5 5];
            g = @(Z)1+Z+2\Z.^2+6\Z.^3+24\Z.^4+120\Z.^5+720\Z.^6+5040\Z.^7;
            clf; abs_stab_region(g,range); pause;
        case 8
            range = [-5 3 -5 5];
            g = @(Z)1+Z+2\Z.^2+6\Z.^3+24\Z.^4+120\Z.^5+720\Z.^6+5040\Z.^7+40320\Z.^8;
            clf; abs_stab_region(g,range); 
    end
              
end

hold off;