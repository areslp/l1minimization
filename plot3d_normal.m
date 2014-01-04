function [ output_args ] = plot3d_normal( input_args )

figure;
[X,Y] = meshgrid(-2:.5:2); 
R = sqrt(X.^2 + Y.^2) + eps;
% Z = sin(R)./R;
Z = exp( -(R^2/2) );
mesh(X,Y,Z,'EdgeColor','black')
surfnorm(X,Y,Z);
set(gcf, 'Color', 'w');
axis tight;



figure;
[X,Y] = meshgrid(-2:.5:2); 
R = sqrt(X.^2 + Y.^2) + eps;
% Z = sin(R)./R;
Z = exp( -(R^2/100000) );
mesh(X,Y,Z,'EdgeColor','black')
surfnorm(X,Y,Z);
set(gcf, 'Color', 'w');
axis tight;


% figure;
% X=linspace(0,1);
% Y=linspace(0,1);
% Z1=ones(100,100);
% Z2=zeros(100,100);
% mesh(X,Y,Z1);hold on;
% mesh (X,Y,Z2);hold on;
% mesh (Y,Z1,X);hold on;
% mesh (Y,Z2,X);hold on;
% mesh (Z1,Y,X);hold on;
% mesh (Z2,Y,X);
% % surfnorm(X,Y,Z);
% set(gcf, 'Color', 'w');

end

