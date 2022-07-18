function [D,R]=floyd(W)
% @function: floyd�㷨������������֮������·
% @illustrate: �����ȨͼW���������·������D����·�����ݾ���R
% @param W: ��Ȩͼ NxN
% @return D: ���·������ NxN��D(i��j)��ʾ��i��j�����·��
% @return R: ���·���Ĺ켣׷�� NxN
%       ��1-7�����·��Ϊ1->2->3->7
%            R(1,7) = iter1 = 2
%            R(iter1,7) = iter2 = 3 
%            R(iter2,7) = iter3 = 7
%

    [m,n] = size(W);
    if m~=n
        error('W is not square');
    end
    
    D = W;
    R = repmat(1:n,n,1); % ��1��n���Ƴ� n �� 1 ��,��¼����·��
    for iter = 1:n % ���
        for i=1:n
            for j=1:n
                if D(i,j)>D(i,iter)+D(iter,j) % ����֮��С�ڵ�����
                    D(i,j)=D(i,iter)+D(iter,j);
                    R(i,j)=iter;
                end
            end
        end
    end         
end