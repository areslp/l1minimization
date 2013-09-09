function [] = image_as_cloud()
% 把图像作为cloud来处理，然后转回image看结果

I=imread('color.jpg');
[m,n,k]=size(I);
points_in=zeros(m*n,3);
normals_in=zeros(m*n,3);
for i=1:m
    for j=1:n
        points_in((i-1)*n+j,:)=[i,j,0];
        normals_in((i-1)*n+j,:)=I(i,j,:);
    end
end

[points_out,normals_out]=cloud_process(points_in,normals_in);

nI=I;
for i=1:m
    for j=1:n
        nI(i,j,:)=normals_out((i-1)*n+j,:);
    end
end


figure;
subplot(1,2,1);
imshow(I);
subplot(1,2,2);
imshow(nI);
