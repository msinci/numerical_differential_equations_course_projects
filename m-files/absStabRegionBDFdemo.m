%% 
% Plots absolute-stability regions for BDF methods of orders 
% one through six. 
% 
%

p = 6;

xmin = -10; xmax = 30; 
ymin = -20; ymax = 20;

n = 500;
h = (n+1)\2*pi;
theta = 0:h:2*pi;


figure(2); clf; axis square; hold on;
plot([xmin xmax],[0 0],'-k','LineWidth',1);
plot([0 0],[ymin ymax],'-k','LineWidth',1);
axis([xmin xmax ymin ymax]);grid on;axis equal;axis square;

ASRT = 'Absolute Stability Regions for BDF Methods';
title(ASRT,'FontSize',16);

for k = 1:p
    switch k
        case 1
% Backward Euler
plot( exp(i*theta).\(exp(i*theta) - 1),'LineWidth',2);
title({ASRT;'Order 1 (Backward Euler, exterior of circle)'});
pause;
        case 2
% 2nd order
plot( ((3\2)*exp(i*(2*theta))).\(exp(i*(2*theta)) - (3\4)*exp(i*theta) ...
    + 3\1),'LineWidth',2);
title({ASRT;['Orders 1-2 (Exteriors of Closed Curves)']});
pause;
        case 3 
% 3rd order
plot( ((11\6)*exp(i*(3*theta))).\...
    (exp(i*(3*theta)) - (11\18)*exp(i*(2*theta)) + (11\9)*exp(i*theta) ... 
    - 11\2),'LineWidth',2); 
title({ASRT;'Orders 1-3 (Exteriors of Closed Curves)'});
pause;
        case 4
% 4th order
angles = [theta(86:417) theta(86)];
plot( ((25\12)*exp(i*(4*angles))).\...
    (exp(i*(4*angles)) - (25\48)*exp(i*(3*angles)) ...
    + (25\36)*exp(i*(2*angles)) - (25\16)*exp(i*angles) - 25\3), ...
    'LineWidth',2);
title({ASRT;'Orders 1-4 (Exteriors of Closed Curves)'});
pause;
        case 5
% 5th order
angles = theta;
plot( ((137\60)*exp(i*(5*angles))).\...
    (exp(i*(5*angles)) - (137\300)*exp(i*(4*angles)) ...
    + (137\300)*exp(i*(3*angles)) - (137\200)*exp(i*(2*angles)) ... 
    + (137\75)*exp(i*angles) - 137\12), ...
    'LineWidth',2);
title({ASRT;'Orders 1-5 (Exteriors of Closed Curves)'});
pause;
        case 6
%6th order
angles = theta;
plot( ((147\60)*exp(i*(6*angles))).\...
    (exp(i*(6*angles)) - (147\360)*exp(i*(5*angles)) ...
    + (147\450)*exp(i*(4*angles)) - (147\400)*exp(i*(3*angles)) ... 
    + (147\225)*exp(i*(2*angles)) - (147\72)*exp(i*angles) + 147\10), ...
    'LineWidth',2);
title({ASRT;'Orders 1-6 (Exteriors of Closed Curves)'});
%pause;

    end % End of switch. 
              
end % End of outer loop 

hold off;
