close all; clear all; clc; 
n=10;
x_start=zeros(n,1);

% analytical p* computation
p_rb=0;
p_wc=-0.5*sum(linspace(1,n,n))+1; 
p_ic=sum(linspace(1,n*100,n))/200-sum(linspace(1/10,n*10,n))/10+1;

%% Gradient Descent
% remark: in orser to calculate error according to best numerical result
% instead of anlytic change ~ into p_rb,p_wc,p_ic  

[~,f_rb]=gradint_descent(@rosenbrock,x_start);
gd_rosenbrock_error=f_rb-0;

[~,f_wc]=gradint_descent(@well_cond,x_start);
gd_well_cond_error=f_wc-p_wc;

[~,f_ic]=gradint_descent(@ill_cond,x_start);
gd_ill_cond_error=f_ic-p_ic;



%% Newton Method

[~,f_rb]=newton_method(@rosenbrock, x_start);
nm_rosenbrock_error=f_rb-0;
turning_point_rb=find_max_diff(nm_rosenbrock_error);

[~,f_wc]=newton_method(@well_cond, x_start);
nm_well_cond_error=f_wc-p_wc;
% turning_point_wc=find_max_diff(nm_rosenbrock_error);

[~,f_ic]=newton_method(@ill_cond, x_start);
nm_ill_cond_error=f_ic-p_ic;
% turning_point_ic=find_max_diff(nm_rosenbrock_error);



%% plots

% Graient Descent
figure;
semilogy(1:length(gd_rosenbrock_error),gd_rosenbrock_error);
xlabel('Iteration Number');
ylabel('Error on f(x)');
title('Gradient Descent - Rosenbrock');


figure;
semilogy(1:1:length(gd_well_cond_error),gd_well_cond_error);
xlabel('Iteration Number');
ylabel('Error on f(x)');
title(['Gradient Descent - Well Conditioned Function CR:',num2str(1-1/n)]);


figure;
semilogy(1:1:length(gd_ill_cond_error),gd_ill_cond_error);
xlabel('Iteration Number');
ylabel('Error on f(x)');
title(['Gradient Descent - Ill Conditioned Function CR:',num2str(1-1/(n*10))]);

% Newton Method
figure;
semilogy(1:length(nm_rosenbrock_error),nm_rosenbrock_error);
xlabel('Iteration Number');
ylabel('Error on f(x)');
title('Newton Method - Rosenbrock');
text(turning_point_rb,(nm_rosenbrock_error(turning_point_rb)),'\fontsize{10}\leftarrow changes from linear to quadratic');

% graphs have been removed since in NM for quadratic functions it takes only one iteration

% figure;
% semilogy(1:1:length(nm_well_cond_error),nm_well_cond_error);
% xlabel('Iteration Number');
% ylabel('Error on f(x)');
% title('Gradient Descent - Well Conditioned Function');
% %text(turning_point_wc,(nm_rosenbrock_error(turning_point_wc)),'\fontsize{10}\leftarrow changes from linear to quadratic');
% 
% 
% figure;
% semilogy(1:1:length(nm_ill_cond_error),nm_ill_cond_error);
% xlabel('Iteration Number');
% ylabel('Error on f(x)');
% title('Gradient Descent - Ill Conditioned Function');
% %text(turning_point_ic,(nm_rosenbrock_error(turning_point_ic)),'\fontsize{10}\leftarrow changes from linear to quadratic');



% plot rosenbrock function convergence in 2D
banana = @(x,y) (1-x).^2 + 100*(y-x.^2).^2;

x = linspace(-1.5,1.5);
y = linspace(-1,1.5);
[xx,yy] = meshgrid(x,y); 
ff = banana(xx,yy);

levels =[1:1:10,10:10:100];
figure;
contour(x,y,ff,levels);
colorbar;
axis([-1.5 1.5 -1 1.5]); axis square; hold on;
title('Rosenbrock function 2D')

[x_gd,~]=gradint_descent(@rosenbrock,[-0.5;1.5]);
[x_nm,~]=newton_method(@rosenbrock,[-0.5;1.5]);

plot(x_gd(:,1),x_gd(:,2),'k','LineWidth',2);
plot(x_nm(:,1),x_nm(:,2),'r','LineWidth',2);

legend('Rosenbrock Function','Gradient Desecnt','Newton Method');