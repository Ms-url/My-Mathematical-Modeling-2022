%
% ţ�ٷ�����jacobi����(����)
%     �ŵ㣺�����������ٶȾ޿�(ƽ������)����У�����ᴫ�ݣ�
%     ȱ�㣺���Է��£�ֻҪ��ֵ���jacobi�������ϡ���ԡ������ԡ���̬�ȣ��͹��ˣ�
% 
% ��ţ�ٷ������߷�˼�룬�ý��ƾ�������jacobi����
%     �ŵ㣺jacobi��������������ﶼ�������⣡����ŵ㼫����߽ⷨ���ȶ��ԣ�����
%     ȱ�㣺�����ٶȽ���ƽ��������ֱ������֮�䣬����һ������
%
% ��������һ�����ٶ�ѡ����ţ�ٷ�����ţ�ٷ����ȶ���������һ��������
%
function x_result = SolveNonlinearEquations(F,x,x0,mode)
% @function:
%       ʹ�� ��ţ�ٷ� �� ţ�ٷ� �������Է���
% @illustrate��
%       �ú����е�ţ�ٷ�ʹ���˸�˹-���¶�������������pre_seidel����������
%
% @param mode:
%       1 -> ��ţ�ٷ�
%       2 -> ţ�ٷ�
% @example:
%       syms x1 x2; 
%       f1 = x1^2 - 10*x1 + x2^2 + 8;
%       f2 = x1*x2^2 + x1 - 10*x2 + 8;
%       F = [f1;f2]; % ������:������ 
%       x = [x1;x2]; % δ֪��:������
%       x0 = [2;3];  % ��ʼֵ:������
%       x_result = SolveNonlinearEquations(f,x,x0,1)
% 

    if mode==1
        disp("��ţ�ٷ�")
        x_result = ni_Newton(F,x,x0);
    elseif mode==2
        disp("ţ�ٷ�")
        x_result = Newton(F,x,x0);
    else
        disp("mode�����������")
    end
    
end

function x_result = Newton(F,x,x0)
% @@illustrate:
%       F = [f1,f2,...,fn]'
%       x(k+1) = x(k) - (F'^-1)F  
% @function��
%       ʹ��ţ�ٷ��������Է�����
%           f1(x1,x2,...,xn)=0
%           f2(x1,x2,...,xn)=0
%                   |
%           fn(x1,x2,...,xn)=0
%       �Ľ�,��̬����
%           1��dx(k)�����ľ���
%           2��F(x(k))�����ľ���
%           3����������
% @param F:
%       ���ű��������飬������
%       F = [f1 ; f2 ; ... ; fn ]
%
% @param x:
%       ������δ֪������������
%       x = [x1 ; x2 ; ... ; xn ]
%
% @param x0:
%       ������ʼֵ��������
%       x0 = [a1 ; a2 ; ... ; an]
%
% @return x_result: ������Ľ⣬������
%
% @example:
%       syms x1 x2; 
%       f1 = x1^2 - 10*x1 + x2^2 + 8;
%       f2 = x1*x2^2 + x1 - 10*x2 + 8;
%       F = [f1;f2]; % ������:������ 
%       x = [x1;x2]; % δ֪��:������
%       x0 = [2;3];  % ��ʼֵ:������
%       x_result = Newton(f,x,x0)
% 
    % x(k+1) = x(k) - dx(k)
    min_dxk = double( input('dx(k)�����ľ���:') ); % sprt(sum( ( dx(k) - 0 ).^2  )
    min_fkk = double( input('F(x(k))�����ľ���:') ); % sprt(sum( ( F(x(k)) - 0 ).^2  )
    num = input('��������:');
    
    % jacobi1 = [diff(f1,x1) diff(f1,x2);diff(f2,x1) diff(f2,x2)]
    % ֱ�����Դ��������ſ˱Ⱦ���:
    jacobi = jacobian(F',x');
    
    for k = 1:num
        dFk = double( subs(jacobi, x, x0) );
        Fk = double( subs(F, x, x0) );
        
        % ��˹-���¶�����������- ((F')^-1)*F
        % ת���ɽ����Է��̣�����������
        dxk = pre_seidel(dFk,-Fk, k);  
        
        x0 = x0 + dxk;
        fxk = double( subs(F, x, x0) );  % fk+1���������ж�
        if norm(dxk) < min_dxk || norm(fxk) < min_fkk
            break;
        end
    end
    
    if k < num
        x_result = x0;
        fprintf('�����Ѵ�Ҫ�󣬵�����ǰ����!\n');
        fprintf('%d�ε�����, ���ƽ�:\n',k);
    else
        x_result = x0;
        fprintf("���������Ѵ�����!\n");
        fprintf("���ƽ�Ϊ:\n",k);
    end
    
    fprintf('F(x)���Ϊ:%f\n',double( subs(F,x,x0) ));
    fprintf('dx(k)����:%f\n',norm(dxk));
    fprintf('F(x(k))����:%f\n',norm(fxk));
end

function  x_result = ni_Newton(F,x,x0)
% @illustrate:
%       B(0) = (F'(x(0)))^-1
%       s(k) = -B(k)*F(x(k))
%       x(k+1) = x(k) + s(k)
%       y(k) = F(x(k+1)) - F(x(k))
%       B(k+1) = B(k) + ( (s(k)-B(k)*y(k))*s(k)'*B(k) )/( s(k)'*B(k)*y(k) )
%   �����б���˵��
%   F(x(k)) -> Fk 
%   F(x(k+1)) -> Fkk
%
% @function��
%       ʹ����Broyden��1�ڶ������������Է�����
%           f1(x1,x2,...,xn)=0
%           f2(x1,x2,...,xn)=0
%                   |
%           fn(x1,x2,...,xn)=0
%       �Ľ�,��̬����
%           1��dx(k)�����ľ���
%           2��F(x(k))�����ľ���
%           3����������
%
% @param F:
%       ���ű��������飬������
%       F = [f1 ; f2 ; ... ; fn ]
%
% @param x:
%       ������δ֪������������
%       x = [x1 ; x2 ; ... ; xn ]
%
% @param x0:
%       ������ʼֵ��������
%       x0 = [a1 ; a2 ; ... ; an]

    min_sk = double( input('s��k�������ľ���:') );
    min_fkk = double( input('F��x(k)�������ľ���:') );
    num = input('��������:');
    
    % ֱ�����Դ��������ſ˱Ⱦ���:
    jacobi = jacobian(F,x);
    
    % ��ʼ���ǰ�ĳ�ֵ:
    Bk = inv( double( subs(jacobi,x,x0) ) );
    Fk = double( subs(F,x,x0) );
    % ѭ����ȫ����������
    for k = 1:num
        sk = -Bk*Fk;
        x0 = x0 + sk;
        Fkk = double( subs(F,x,x0) );
        if norm(sk) < min_sk || norm(Fkk) < min_fkk
            break;
        end
        yk = Fkk - Fk;
        Bkk = Bk + (sk-Bk*yk)*(sk')*Bk/( sk'*Bk*yk );  % ��ͬ������������ʽ����
        Fk = Fkk;
        Bk = Bkk;
    end
    
    if k < num
        x_result = x0;
        fprintf('�����Ѵ�Ҫ�󣬵�����ǰ����!\n');
        fprintf('%d�ε�����, ���ƽ�Ϊ:\n',k);
    else
        x_result = x0;
        fprintf('���������Ѵ�����!\n');
        fprintf("���ƽ�Ϊ:\n",k);
    end

fprintf('F(x)���Ϊ:%f\n',double( subs(F,x,x0) ));
fprintf('sk����:%f\n',norm(sk));
fprintf('fkk����:%f\n',norm(Fkk));

end

function x = pre_seidel(dF,F,n,varargin)
% @illustrate:
%       ���Է�����Ԥ������"��������"��˹-���¶�������
%       ���Ӹ������⣬ֻҪ�н⣬һ���ܽ����
%
% TDOT:ԭ��δ֪
%

    % Ԥ���� 
    F = dF'*F;
    dF = dF'*dF;
    % �����������ĸ�˹-���¶�����:
    D = diag(diag(dF));
    L = tril(dF,-1);     % ��������һ�����������;
    U = triu(dF,1);      % ��������һ�����������;
    B2 = -(D+L)\U;   % ���¶���������, ����������ж��������;
    g2 = (D+L)\F;
    
    radius = max(abs(eig(B2)));  % ����ֵ�п���Ϊ����, absȡ����ֵ + ȡģ
    fprintf('��%d��������Է�����, ���¶����������װ뾶Ϊ(ԽСԽ��): %.4f\n',n,radius);
    
    % �������㲿��:
    x = zeros(length(F),1);  % ��ʼ����4x1�ľ���(�о���): ��ʼֵȫΪ0
    
    if lenght(varargin)==0
        error = 0.0001; % ͳһ��һ�����ȼ���,���޸�
    else
        error = varargin(1);
    end
    
    count = 0;      % ����������
    while 1
        tmp = B2*x + g2;
        if max(abs(tmp - x)) < error
            break;
        end
        x = tmp;
        count = count + 1;
    end

end

