% This produces phase-plane quiver plots for several of the autonomous
% ODEs considered in the course. 

%% Predator-prey equation

% Set the alpha parameter. 
alpha = 1;

% Define the right-hand side functions. 
F1 = @(r,f)2*r-alpha*r.*f;
F2 = @(r,f)-f+alpha*r.*f;

% Create nxn arrays of r and f values. 
% Set n. 
n = 12; 
% Create the r values. 
rmin = 0; rmax = 2.5; hr = n\(rmax-rmin); rvals = rmin:hr:rmax;
% Create the f values. 
fmin = 0; fmax = 4; hf = n\(fmax-fmin); fvals = fmax:-hf:fmin;
% Create the arrays. 
[R,F] = meshgrid(rvals,fvals); 

% Evaluate the functions on the arrays. 
rprimevals = F1(R,F);
fprimevals = F2(R,F);

% Create the quiver plots. 
figure(1);
quiver(rvals,fvals,rprimevals,fprimevals,1.5);
axis([rmin-.1,rmax+.1,fmin-.1,fmax+.1]);axis('square');

%% Damped oscillator. 

% Set the D parameter. 
D = 1;

% Define the right-hand side functions. 
f1 = @(y1,y2) y2;
f2 = @(y1,y2) -y1 - D*y2;

% Create nxn arrays of y1 and y2 values. 
% Set n. 
n = 12;
% Create the y1 values. 
y1min = -2; y1max = 2; hy1 = n\(y1max-y1min); y1vals = y1min:hy1:y1max;
% Create the y2 values. 
y2min = -2; y2max = 2; hy2 = n\(y2max-y2min); y2vals = y2max:-hy2:y2min;
% Create the arrays. 
[Y1,Y2] = meshgrid(y1vals,y2vals); 

% Evaluate the functions on the arrays. 
y1primevals = f1(Y1,Y2);
y2primevals = f2(Y1,Y2);

% Create the quiver plots. 
figure(1);
quiver(y1vals,y2vals,y1primevals,y2primevals);
axis([y1min-.1,y1max+.1,y2min-.1,y2max+.1]);axis('square');

%% Clarinet reed. 

% Set the a, b, and k parameters. 
a = 5; b = 4; k = 5;

% Define the right-hand side functions. 
f1 = @(y1,y2) y2;
f2 = @(y1,y2) -k*y1 - b*y2.^3 + a*y2;

% Create nxn arrays of y1 and y2 values. 
% Set n. 
n = 12;
% Create the y1 values. 
y1min = -.7; y1max = .7; hy1 = n\(y1max-y1min); y1vals = y1min:hy1:y1max;
% Create the y2 values. 
y2min = -.7; y2max = .7; hy2 = n\(y2max-y2min); y2vals = y2max:-hy2:y2min;
% Create the arrays.
[Y1,Y2] = meshgrid(y1vals,y2vals); 

% Evaluate the functions on the arrays. 
y1primevals = f1(Y1,Y2);
y2primevals = f2(Y1,Y2);

% Create the quiver plots. 
figure(1);
quiver(y1vals,y2vals,y1primevals,y2primevals);
axis([y1min-.1,y1max+.1,y2min-.1,y2max+.1]);axis('square');


