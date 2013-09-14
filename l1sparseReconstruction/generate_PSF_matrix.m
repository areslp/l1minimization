function [M] = generate_PSF_matrix(sigma,mn)
% 生成一个只对一个维度模糊的circle convolution matrix，只针对一个channel
% 对图像进行模糊时，先对I进行，再对I'进行，就相当于先对行进行模糊，再对列进行模糊
off=-ceil(sigma*3):ceil(sigma*3);  % psf pixel offsets
[a,b]=size(off);
off=reshape(off,a*b,1);
psf=exp( -off.^2/sigma^2 );        % psf values for offsets
psf=psf/sum(psf);                  % normalize to 1.0
 
% generate psf matrix by setting diagonals based on 1D psf above
tnum=mn;
M=sparse([],[],[],tnum,tnum,tnum*a*b);                   
for i=1:numel(psf),
    M = M + spdiags(psf(i)*ones(tnum,1),off(i),tnum,tnum);
end
