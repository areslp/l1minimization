function [] = plot_1D_cloud(nin,nout)
% 可视化曲线实验里的法矢，曲线实验里点是按顺序来排列的
dim=size(nin,2);
n=size(nin,1);
for i=1:dim
    vin=nin(:,i);
    vout=nout(:,i);
    figure;
    hold on;
    plot(vin,'.r-');
    plot(vout,'.b-');
    hold off;
    axis([-1 n+1 -0.2 1.2]);
end
