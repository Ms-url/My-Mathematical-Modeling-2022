%
% 牛顿法：用jacobi矩阵(导数)
%     优点：导数法收敛速度巨快(平方收敛)；自校正误差不会传递；
%     缺点：求导稍费事；只要赋值后的jacobi矩阵存在稀疏性、奇异性、病态等，就跪了；
% 
% 拟牛顿法：割线法思想，用近似矩阵趋近jacobi矩阵
%     优点：jacobi矩阵的问题在这里都不是问题！这个优点极大提高解法的稳定性！！！
%     缺点：收敛速度介于平方收敛和直线收敛之间，稍慢一丢丢。
%
% 建议牺牲一丢的速度选择拟牛顿法！拟牛顿法的稳定性真的提高一个量级。
%
function x_result = SolveNonlinearEquations(F,x,x0,mode)
% @function:
%       使用 拟牛顿法 或 牛顿法 求解非线性方程
% @illustrate：
%       该函数中的牛顿法使用了高斯-赛德尔迭代法（函数pre_seidel），无求逆
%
% @param mode:
%       1 -> 拟牛顿法
%       2 -> 牛顿法
% @example:
%       syms x1 x2; 
%       f1 = x1^2 - 10*x1 + x2^2 + 8;
%       f2 = x1*x2^2 + x1 - 10*x2 + 8;
%       F = [f1;f2]; % 方程组:列向量 
%       x = [x1;x2]; % 未知数:列向量
%       x0 = [2;3];  % 初始值:列向量
%       x_result = SolveNonlinearEquations(f,x,x0,1)
% 

    if mode==1
        disp("拟牛顿法")
        x_result = ni_Newton(F,x,x0);
    elseif mode==2
        disp("牛顿法")
        x_result = Newton(F,x,x0);
    else
        disp("mode参数输入错误")
    end
    
end

function x_result = Newton(F,x,x0)
% @@illustrate:
%       F = [f1,f2,...,fn]'
%       x(k+1) = x(k) - (F'^-1)F  
% @function：
%       使用牛顿法求解非线性方程组
%           f1(x1,x2,...,xn)=0
%           f2(x1,x2,...,xn)=0
%                   |
%           fn(x1,x2,...,xn)=0
%       的解,动态接收
%           1、dx(k)范数的精度
%           2、F(x(k))范数的精度
%           3、迭代次数
% @param F:
%       符号变量方程组，列向量
%       F = [f1 ; f2 ; ... ; fn ]
%
% @param x:
%       变量（未知数），列向量
%       x = [x1 ; x2 ; ... ; xn ]
%
% @param x0:
%       变量初始值，列向量
%       x0 = [a1 ; a2 ; ... ; an]
%
% @return x_result: 方程组的解，列向量
%
% @example:
%       syms x1 x2; 
%       f1 = x1^2 - 10*x1 + x2^2 + 8;
%       f2 = x1*x2^2 + x1 - 10*x2 + 8;
%       F = [f1;f2]; % 方程组:列向量 
%       x = [x1;x2]; % 未知数:列向量
%       x0 = [2;3];  % 初始值:列向量
%       x_result = Newton(f,x,x0)
% 
    % x(k+1) = x(k) - dx(k)
    min_dxk = double( input('dx(k)范数的精度:') ); % sprt(sum( ( dx(k) - 0 ).^2  )
    min_fkk = double( input('F(x(k))范数的精度:') ); % sprt(sum( ( F(x(k)) - 0 ).^2  )
    num = input('迭代次数:');
    
    % jacobi1 = [diff(f1,x1) diff(f1,x2);diff(f2,x1) diff(f2,x2)]
    % 直接用自带函数求雅克比矩阵:
    jacobi = jacobian(F',x');
    
    for k = 1:num
        dFk = double( subs(jacobi, x, x0) );
        Fk = double( subs(F, x, x0) );
        
        % 高斯-赛德尔迭代法计算- ((F')^-1)*F
        % 转换成解线性方程，避免多次求逆
        dxk = pre_seidel(dFk,-Fk, k);  
        
        x0 = x0 + dxk;
        fxk = double( subs(F, x, x0) );  % fk+1单纯用来判断
        if norm(dxk) < min_dxk || norm(fxk) < min_fkk
            break;
        end
    end
    
    if k < num
        x_result = x0;
        fprintf('精度已达要求，迭代提前结束!\n');
        fprintf('%d次迭代后, 近似解:\n',k);
    else
        x_result = x0;
        fprintf("迭代次数已达上限!\n");
        fprintf("近似解为:\n",k);
    end
    
    fprintf('F(x)结果为:%f\n',double( subs(F,x,x0) ));
    fprintf('dx(k)范数:%f\n',norm(dxk));
    fprintf('F(x(k))范数:%f\n',norm(fxk));
end

function  x_result = ni_Newton(F,x,x0)
% @illustrate:
%       B(0) = (F'(x(0)))^-1
%       s(k) = -B(k)*F(x(k))
%       x(k+1) = x(k) + s(k)
%       y(k) = F(x(k+1)) - F(x(k))
%       B(k+1) = B(k) + ( (s(k)-B(k)*y(k))*s(k)'*B(k) )/( s(k)'*B(k)*y(k) )
%   程序中变量说明
%   F(x(k)) -> Fk 
%   F(x(k+1)) -> Fkk
%
% @function：
%       使用逆Broyden秩1第二方法求解非线性方程组
%           f1(x1,x2,...,xn)=0
%           f2(x1,x2,...,xn)=0
%                   |
%           fn(x1,x2,...,xn)=0
%       的解,动态接收
%           1、dx(k)范数的精度
%           2、F(x(k))范数的精度
%           3、迭代次数
%
% @param F:
%       符号变量方程组，列向量
%       F = [f1 ; f2 ; ... ; fn ]
%
% @param x:
%       变量（未知数），列向量
%       x = [x1 ; x2 ; ... ; xn ]
%
% @param x0:
%       变量初始值，列向量
%       x0 = [a1 ; a2 ; ... ; an]

    min_sk = double( input('s（k）范数的精度:') );
    min_fkk = double( input('F（x(k)）范数的精度:') );
    num = input('迭代次数:');
    
    % 直接用自带函数求雅克比矩阵:
    jacobi = jacobian(F,x);
    
    % 开始求解前的初值:
    Bk = inv( double( subs(jacobi,x,x0) ) );
    Fk = double( subs(F,x,x0) );
    % 循环完全按照流程来
    for k = 1:num
        sk = -Bk*Fk;
        x0 = x0 + sk;
        Fkk = double( subs(F,x,x0) );
        if norm(sk) < min_sk || norm(Fkk) < min_fkk
            break;
        end
        yk = Fkk - Fk;
        Bkk = Bk + (sk-Bk*yk)*(sk')*Bk/( sk'*Bk*yk );  % 不同方法改这里表达式即可
        Fk = Fkk;
        Bk = Bkk;
    end
    
    if k < num
        x_result = x0;
        fprintf('精度已达要求，迭代提前结束!\n');
        fprintf('%d次迭代后, 近似解为:\n',k);
    else
        x_result = x0;
        fprintf('迭代次数已达上限!\n');
        fprintf("近似解为:\n",k);
    end

fprintf('F(x)结果为:%f\n',double( subs(F,x,x0) ));
fprintf('sk范数:%f\n',norm(sk));
fprintf('fkk范数:%f\n',norm(Fkk));

end

function x = pre_seidel(dF,F,n,varargin)
% @illustrate:
%       线性方程组预处理后的"万能收敛"高斯-赛德尔迭代法
%       无视各种问题，只要有解，一定能解出。
%
% TDOT:原理未知
%

    % 预处理 
    F = dF'*F;
    dF = dF'*dF;
    % 下面是正常的高斯-赛德尔操作:
    D = diag(diag(dF));
    L = tril(dF,-1);     % 向左下移一格的下三角阵;
    U = triu(dF,1);      % 向右上移一格的上三角阵;
    B2 = -(D+L)\U;   % 赛德尔迭代矩阵, 用来计算和判断收敛与否;
    g2 = (D+L)\F;
    
    radius = max(abs(eig(B2)));  % 特征值有可能为复数, abs取绝对值 + 取模
    fprintf('第%d次求解线性方程组, 赛德尔迭代矩阵谱半径为(越小越好): %.4f\n',n,radius);
    
    % 迭代计算部分:
    x = zeros(length(F),1);  % 初始迭代4x1的矩阵(列矩阵): 初始值全为0
    
    if lenght(varargin)==0
        error = 0.0001; % 统一用一个精度即可,可修改
    else
        error = varargin(1);
    end
    
    count = 0;      % 迭代计数器
    while 1
        tmp = B2*x + g2;
        if max(abs(tmp - x)) < error
            break;
        end
        x = tmp;
        count = count + 1;
    end

end

