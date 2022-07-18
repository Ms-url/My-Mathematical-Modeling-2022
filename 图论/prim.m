function [edge,weight] = prim(W)
% @function: 普里母算法求最小生成树
% @illustrate: 
% @param W: 边权图 NxN
% @return edge: 最小生成树边集，2x(N-1)
%           每一列对应一条边的两个顶点
%               a1,a2,...,an
%               b1,b2,...,bn
% @return weight: 最小生成树的权之和
%
    [m,n] = size(W);
    if m~=n
        error('W is not square');
    end
    % 初始化
    U = 1;
    Ubar = 2:n ;
    [v1,pos1]=min(W(U,Ubar));
    
    weight=v1;
    edge=[1;Ubar(pos1)];
    
    U=[U,Ubar(pos1)];
    Ubar(pos1) = [];
    for i=1:n-2
        [v1,pos1]=min(W(U,Ubar));
        [v2,pos2]=min(v1);
        
        edge = [edge,[U(pos1(pos2));Ubar(pos2)]];
        
        weight = weight+v2; % 更新
        U = [U,Ubar(pos2)]; 
        Ubar(pos2)=[];
    end
end


