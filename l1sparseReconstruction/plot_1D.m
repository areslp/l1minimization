function [ ] = plot_1D( I,OUT,row )
[m n k]=size(I);

if k==3
    rI=I(row,:,:);
    rI0=OUT(row,:,:);
    for i=1:k
        figure;
        hold on;
        plot(rI(:,:,i),'.r-');
        plot(rI0(:,:,i),'.b-');
        hold off;
        axis([-1 n+1 -0.2 1.2]);
    end
end
if k==1
    rI=I(row,:);
    rI0=OUT(row,:);
    figure;
    hold on;
    plot(rI,'.r-');
    plot(rI0,'.b-');
    hold off;
    axis([-1 n+1 -0.2 1.2]);
end


end

