function [D,R]=floyd(W)
% @function: floyd算法计算两两顶点之间的最短路
% @illustrate: 输入边权图W，返回最短路径矩阵D，和路径回溯矩阵R
% @param W: 边权图 NxN
% @return D: 最短路径矩阵 NxN，D(i，j)表示由i到j的最短路径
% @return R: 最短路径的轨迹追随 NxN
%       若1-7的最短路径为1->2->3->7
%            R(1,7) = iter1 = 2
%            R(iter1,7) = iter2 = 3 
%            R(iter2,7) = iter3 = 7
%

    [m,n] = size(W);
    if m~=n
        error('W is not square');
    end
    
    D = W;
    R = repmat(1:n,n,1); % 将1：n复制成 n 行 1 列,记录回溯路径
    for iter = 1:n % 插点
        for i=1:n
            for j=1:n
                if D(i,j)>D(i,iter)+D(iter,j) % 两边之和小于第三边
                    D(i,j)=D(i,iter)+D(iter,j);
                    R(i,j)=iter;
                end
            end
        end
    end         
end