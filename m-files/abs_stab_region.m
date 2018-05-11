function abs_stab_region(g,range)
%
% Draws the boundary of the region of absolute stability 
% of a method, given the function g that characterizes the 
% region by abs(g(z))<1.
%
% Inputs: 
%  g      = function handle (@g). The function g should be 
%           set up to return a matrix of values g(Z) when 
%           Z is a matrix argument. 
%  range  = 1x4 matrix [xmin xmax ymin ymax] that specifies 
%           the range to be considered. 
%
% Set up the mesh of z-values.
xmin = range(1,1); xmax = range(1,2); hx = 99\(xmax - xmin);
ymin = range(1,3); ymax = range(1,4); hy = 99\(ymax - ymin);
x = xmin:hx:xmax; y = ymin:hy:ymax; 
[X,Y]=meshgrid(x,y);
Z = complex(X,Y); 
% Evaluate g over the mesh and take absolute values.
ABSGVALS= abs(feval(g,Z));
% Draw the boundary of the region. 
contour(X,Y,ABSGVALS,[1 1],'-b');
axis([xmin xmax ymin ymax]);grid on;axis equal;
