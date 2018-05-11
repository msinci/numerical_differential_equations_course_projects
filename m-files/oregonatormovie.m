function ormv = oregonatormovie
%
% Creates and returns a movie file for the Oregonator problem. 
% The movie shows the solution point moving around a log/log/log 
% plot of the solution curve in R^3. 
%
% Must be run with files oregonatorfun.m and oregonatorjac.m 
% in the same directory. These evaluate, respectively, the 
% right-hand side function and its Jacobian with respect to y. 
%
% To play the movie, type "movie(ormv)" to play it once or 
% "movie(ormv,n)" to play it n times. 
%

% Set up for ode15s. 
options = odeset('Jacobian',@oregonatorjac);
y_0 = [1;1.8;1.8];             % Initial value. 
tlast = 290;		       % Final time.
tint = 2;		       % Time between frames. 

% Integrate and plot the reference curve. 
[T,Y]=ode15s(@oregonatorfun,[0 tlast],y_0,options);
plot3(log(Y(:,1)),log(Y(:,2)),log(Y(:,3)));
grid on; axis([-5 15 -10 10 0 12]); hold on;

% Integrate to create the movie points. 
[Tmov,Ymov]=ode15s(@oregonatorfun,[0:tint:tlast],y_0,options);

% Create the movie file. 
for i = 1:size(Ymov,1)-2, 
   plot3(log(Ymov(i,1)),log(Ymov(i,2)),log(Ymov(i,3)),'o',... 
         'MarkerSize',8,'MarkerFaceColor','r',...
	 'MarkerEdgeColor','k','LineWidth',1.8);
  ormv(i) = getframe;
  delete(findobj('Type','line','MarkerFaceColor','r'));
end
hold off; 

