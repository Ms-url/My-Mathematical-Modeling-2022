function [edge,weight] = prim(W)
% @function: ����ĸ�㷨����С������
% @illustrate: 
% @param W: ��Ȩͼ NxN
% @return edge: ��С�������߼���2x(N-1)
%           ÿһ�ж�Ӧһ���ߵ���������
%               a1,a2,...,an
%               b1,b2,...,bn
% @return weight: ��С��������Ȩ֮��
%
    [m,n] = size(W);
    if m~=n
        error('W is not square');
    end
    % ��ʼ��
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
        
        weight = weight+v2; % ����
        U = [U,Ubar(pos2)]; 
        Ubar(pos2)=[];
    end
end


