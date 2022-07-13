function GreyPredictive(x0,mode)
%%%%%%%%%%%%%%%%%%%   灰色预测模型  %%%%%%%%%%%%%%%%%%%%%%%
% 适用范围：    
%       该模型使用的不是原始数据的序列，而是
%       生成的数据序列。核心体系是Grey Model.即对原始
%       数据作累加生成（或其他处理生成）得到近似的指数
%       规律再进行建模的方法。
% 优点：       
%       在处理较少的特征值数据，不需要数据的样本
%       空间足够大，就能解决历史数据少、序列的完整性以
%       及可靠性低的问题，能将无规律的原始数据进行生成
%       得到规律较强的生成序列。
% 缺点：       
%       只适用于中短期的预测，只适合近似于指数增
%       长的预测。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @illustrate:
%       灰度预测的模型均为微分方程模型
%       无返回值，结果在命令窗打印
%
% @param x0:
% @param mode:预测方法
%       1 -> GM1_1
%       2 -> GM2_1
%       3 -> DGM2_1
%       4 -> Verhulst
%

    switch(mode)
        case 1
            disp("GM1_1");
            GM1_1(x0);
        case 2
            disp("GM2_1");
            GM2_1(x0);
        case 3
            disp("DGM2_1");
            DGM2_1(x0);
        case 4
            disp("Verhulst");
            Verhulst(x0);
        otherwise
            disp("输入参数错误");
    end

end

function GM1_1(x0)
% @illustrate:
%       GM(1,1)表示模型是 1 阶微分方程，且只含 1 个变量的灰色模型
%       GM(1,1)模型适用于具有较强指数规律的序列，只能描述单调的变化过程，
%       对于非单调的摆动发展序列或有饱和的 S 形序列，可以考虑建立 GM(2,1)，DGM 和 Verhulst 模型
%
% @function: 
%       该函数会进行数据的检验，未落入可容覆盖内时无法使用
%       命令窗打印结果：
%           1.微分方程的解 x(t)
%           2.残差
%           3.相对误差
%           4.级比偏差值
%       其中x(t)是对1次累加序列进行预测，k=t+1
%       实践预测值 prediction(k+1) = x(k+1) - x(k) = x(t)-x(t-1)
%
% @param x0: 
%       观察序列，列向量，nx1
%
    n = length(x0);
    % 数据检验
    lamda = x0(1:n-1)./x0(2:n); 
    if min(lamda)<exp(-2/(n+1))||max(lamda)>exp(2/(n+1))
        disp('未落入可容覆盖内')
        disp('请取适当常数使y=y+c落入可容覆盖内')
        return
    end
    % 模型求解
    y1 = cumsum(x0); % 累加
    
    Y = x0(2:n);
    B = [0.5*(y1(1:n-1)+y1(2:n)),ones(n-1,1)];
    u=B\Y ; % 拟合参数u=[a,b]
    
    %%%% dsolve（）解微分方程
    syms x(t);
    x = dsolve(diff(x,t)+u(1)*x==u(2),x(0)==x0(1)); %求微分方程的符号解
    %%%%
    xt = vpa(x,6); %以小数格式显示微分方程的解
    disp("微分方程的解 x(t)= "+ xt);
    disp("prediction(k+1) = x(k+1) - x(k) = x(t)-x(t-1) ");
    disp("");
    
    y_p=subs(x,t,0:n-1); %求已知数据的预测值
    y_p=double(y_p); % 符号数转换成数值类型，否则无法作差分运算
    
    % 差分运算，还原数据
    x0_p=[x0(1),diff(y_p)]; % 预测值
    
    epsilon=x0'-x0_p; % 计算残差
    disp("残差 = "+ epsilon);
    delta=abs(epsilon./x0'); % 计算相对误差
    disp("相对误差 = "+ delta);
    rho=1-(1-0.5*u(1))/(1+0.5*u(1))*lamda'; % 计算级比偏差值，u(1)=a
    disp("级比偏差值 = "+ rho);
end

function GM2_1(x0)
% @illustrate:
%       模型功能？？？原理？？？
%
% @function: 
%       命令窗打印结果：
%           1.微分方程的解 x(t)
%           2.残差
%           3.相对误差
%       其中x(t)是对1次累加序列进行预测，k=t+1
%       实践预测值 prediction(k+1) = x(k+1) - x(k) = x(t)-x(t-1)
%
% @param x0: 
%       观察序列，列向量，nx1
%   
    x0 = x0';
    n=length(x0); 
    
    x1=cumsum(x0); %计算1次累加序列
    z=0.5*(x1(2:end)+x1(1:end-1))'; %计算均值生成序列
    
    Y=diff(x0)'; %计算1次累减序列
    B=[-x0(2:end)', -z ,ones(n-1,1)];
    u=B\Y ;%最小二乘法拟合参数
    
    syms x(t)
    x=dsolve(diff(x,2)+u(1)*diff(x)+u(2)*x==u(3),x(0)==x1(1),x(5)==x1(6));
    %求符号解
    xt=vpa(x,6) ;%显示小数形式的符号解
    disp("微分方程的解 x(t)= "+ xt);
    disp("prediction(k+1) = x(k+1) - x(k) = x(t)-x(t-1) ");
    disp("");
    
    y_p=subs(x,t,0:n-1); %求已知数据点1次累加序列的预测值
    y_p=double(y_p) ;%符号数转换成数值类型，否则无法作差分运算
    
    x0_P = [y_p(1),diff(y_p)]; % 求已知数据点的预测值
    epsilon = x0-x0_P; % 求残差
    disp("残差 = "+ epsilon);
    delta = abs(epsilon./x0); % 求相对误
    disp("相对误差 = "+ delta);
end

function DGM2_1(x0)
% @illustrate:
%       模型功能？？？原理？？？
%
% @function: 
%       命令窗打印结果：
%           1.微分方程的解 x(t)
%           2.残差
%           3.相对误差
%       其中x(t)是对1次累加序列进行预测，k=t+1
%       实践预测值 prediction(k+1) = x(k+1) - x(k) = x(t)-x(t-1)
%
% @param x0: 
%       观察序列，列向量，nx1
%
    x0 = x0';
    n=length(x0);
    
    Y=diff(x0)'; %求1次累减序列，即1阶向前差分
    B=[-x0(2:end)',ones(n-1,1)];
    u=B\Y ;%最小二乘法拟合参数
    
    syms x(t); % t=0 -> k=1
    dx = diff(x); %定义一阶和二阶导数
    d2x = diff(x,2); 
    
    %求二阶微分方程符号解
    x = dsolve(d2x+u(1)*dx==u(2),x(0)==x0(1),dx(0)==x0(1));
    xt = vpa(x,6); %显示小数形式的符号解
    disp("微分方程的解 x(t)= "+ xt);
    disp("prediction(k+1) = x(k+1) - x(k) = x(t)-x(t-1) ");
    disp("");
    
    y_p = subs(x,t,0:n-1); % 求已知数据点1次累加序列的预测值
    y_p = double(y_p); % 符号数转换成数值类型，否则无法作差分运算
    
    x0_P = [y_p(1),diff(y_p)]; % 求已知数据点的预测值
    epsilon = x0-x0_P; % 求残差
    disp("残差 = "+ epsilon);
    delta = abs(epsilon./x0); % 求相对误
    disp("相对误差 = "+ delta);
    
end

function Verhulst(x0)
% @illustrate:
%       Verhulst 模型主要用来描述具有饱和状态的过程，
%       即 S 形过程，常用于人口预测、生物生长、繁殖预测及产品经济寿命预测等。
%
% @function: 
%       命令窗打印结果：
%           1.微分方程的解 x(t)
%           2.残差
%           3.相对误差
%       其中x(t)是对1次累加序列进行预测，k=t+1
%       实践预测值 prediction(k+1) = x(k+1) - x(k) = x(t)-x(t-1)
%
% @param x0: 
%       观察序列，列向量，nx1
%
    x0 = x0';
    n = length(x0);
    
    x1 = cumsum(x0); %求1次累加序列
    z = 0.5*(x1(2:n)+x1(1:n-1)); %求x1的均值生成序列
    
    Y = x0(2:end)';
    B = [-z',z'.^2];
    u = B\Y; %估计参数a,b的值
    
    %%%% dsolve（）解微分方程
    syms x(t);
    x = dsolve(diff(x)+u(1)*x==u(2)*x^2,x(0)==x0(1)); %求符号解
    %%%%
    % vpa(x,d) uses at least d significant digits, instead of the value of digits.
    xt = vpa(x,6); % 调整为小数形式
    disp("微分方程的解 x(t)= "+ xt);
    disp("prediction(k+1) = x(k+1) - x(k) = x(t)-x(t-1) ");
    disp("");
    
    % subs（x,t,a）将 x中的 t换成 a
    y_p = subs(x,t,0:n-1); % 求已知数据点1次累加序列的预测值
    y_p = double(y_p);% 符号数转换成数值类型，否则无法作差分运算
    
    % diff差分后会减少一项，用y_p(1)补上
    x0_p = [y_p(1),diff(y_p)] ;% 求已知数据点的预测值
    
    epsilon = x0-x0_p; % 求残差
    disp("残差 = "+ epsilon);
    delta = abs(epsilon./x0); %求相对误差
    disp("相对误差 = "+ delta);
    
    % writematrix([x0',x0_p',epsilon',delta'], 'data15_5.xlsx')
end



