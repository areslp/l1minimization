function [A] = generate_PSF_matrix(n,m)
% 生成一个只对2维模糊的circle convolution matrix，没有考虑彩色图像的问题
mn=m*n;
off=9; % A矩阵的子块矩阵有off个需要计算
b0=fspecial('gaussian', [off off], 2);
A=sparse([],[],[],mn,mn,mn*off*off);
% 计算子矩阵
center_idx=floor(off/2)+1;
subA=cell(off,1);
for i=1:off
    M=sparse([],[],[],n,n,off*n);
    row=b0(i,:);
    for j=1:length(row) % length(row)==off
        M=M+spdiags(row(j)*ones(n,1),center_idx-j,n,n);
    end
    subA{i}=M;
end
% 将子矩阵叠到A中去
for i=1:off
    c=cell(m,1);
    znum=i-center_idx;
    zabs=abs(znum); % 需要填充的nxn 0矩阵的个数
    
    % Toeplitz
    % [c{1:m-zabs}]=deal(subA{i});
    % [c{m-zabs+1:end}]=deal(sparse(n,n,0));

    % Periodic end
    c(:)={subA{i}};

    D=blkdiag(c{:});
    D=circshift(D,[0,n*znum]);
    A=A+D;
end

function [r]=ni(i,n)
    r=mod(i-1,n);

function [r]=mi(i,n)
    r=floor((i-1)/n);

function [r]=b(aa,bb,b0)
    [m,n]=size(b0);
    if aa>=1 && aa<=m && bb>=1 && bb<=n
        r=b0(aa,bb);
    else
        r=0;
    end
