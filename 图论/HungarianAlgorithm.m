function count = HungarianAlgorithm(map,M,N)
%%%%%%%%%%%%%% 匈牙利算法 %%%%%%%%%%%%%%%%%
% 匈牙利算法主要用于解决一些与二分图匹配有关的问题用来寻找：
%       最大匹配数 或 最小覆盖点
%
% 最大匹配数
%       二分图左右两边的顶点1对1配对，找到最大配对数
%
% 最小覆盖点
%       另外一个关于二分图的问题是求最小点覆盖：我们想找到最少的一些点，使二分图所有的边都至少有一个端点在这些点之中。
%       倒过来说就是，删除包含这些点的边，可以删掉所有边。
%
% 关键点在于二分图的构建和转化，
%
% https://zhuanlan.zhihu.com/p/96229700 详情查看该链接 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @function ：寻找最大匹配数
% @param map: 边权图（邻接矩阵图），MxN
% @param M: 左侧元素数量
% @param N: 右侧元素数量
% @return: 最大匹配数
% 
    p = zeros([1,N]); % 记录当前右侧元素对应的左侧元素
    count = 0;
    for i=1:M
        vis = zeros([1,N]); % 记录右侧元素是否已被访问过
        [tag,p] = match(i,map,N,vis,p); % 寻找匹配点
        if tag==1
            count = count + 1;
        end
    end
end

function [tag,p] = match(i,map,N,vis,p)
%
% @function：递归寻找合适的匹配
% @return tag: 判断是否匹配成功，1为成功，0为未成功
% @return p：更新p的值
%
    tag = 0;
    for j=1:N
        if map(i,j)==1 && vis(j)==0 % 有边且未被访问
            vis(j)=1; % 标记为访问过
            [temp,~] = match(p(j),map,N,vis,p); % 判断该右侧点匹配的左侧点是否有其他匹配点
            if p(j)==0||temp == 1 
                if p(j)==1
                    [~,p] = match(p(j),map,N,vis,p);
                end
                p(j)=i; % 当前左侧元素成为当前右侧元素新匹配
                tag =1; % 匹配成功
                return
            end
        end
    end
end
