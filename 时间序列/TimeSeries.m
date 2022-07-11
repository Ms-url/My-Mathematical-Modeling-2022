function [prediction,params] = TimeSeries(y,mode)
    % function
    % @param y:观测序列，列向量
    % @param mode：预测方法
    % @return
    % 
    
    prediction = MovingAverage(y);
    [prediction,params] = ExponentialSmoothing(y);
    [prediction,params] = DifferenceExponentialSmoothing(y);

end

function prediction = MovingAverage(y)
    
    %%%%%%%%   简单移动平均法  %%%%%%%%
    % 适用范围：
    %       当预测目标的基本趋势是在某一水平上下波动时，可用一次简单移动平均方法建立预测模型
    %       简单移动平均法只适合做近期预测，而且是预测目标的发展趋势变化不大的情况。
    % 优缺点：
    %       如果目标的发展趋势存在其它的变化，采用简单移动平均法就会产生较大的预测偏差和滞后。
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % function：简单移动平均法
    %       输入观察序列y，计算平移项 N 由2到y长度2/3 的预测结果，依据均方差最小原则，输出最优解
    %
    % @param y:观测序列，列向量nx1
    %
    % @return prediction：
    %       行向量1x3
    %           N(平移项数), 均方差， y(t+1)
    %
    len=length(y);
    S = []; % 均方差
    output = [];
    for N = 2:floor(len*2/3)
        count=0;
        for i=N:len-1
            Mt = sum(y(i-N+1:i))/N;
            count = count + (y(i+1)-Mt)^2 ;
        end
        S = [S;N,(count/(len-N))^0.5];
        Mtout = sum(y(len-N+1:len)/N);
        output = [output;N,Mtout];
    end
    
    % 绘图比较
    figure(99)
    title('简单移动平均法标准差')
    plot(S(:,1),S(:,2),'-o');
    xlabel('N');ylabel('S');
    
    [x,y]=find(S==min(S(:,2))); 
    prediction = [N,min(S(:,2)),output(x,y)];  % 注意prediction 的结构
    fprintf("N = %d\nprediction = %f",output(x,1),prediction);
end

function [prediction,params] = ExponentialSmoothing(y)
    %%%%%%%    指数平滑法   %%%%%%%%%
    % 适用范围：
    %       指数平滑法根据平滑次数的不同，又分为一次指数平滑法、二次指数平滑法和三次指数平滑法
    % 优缺点：
    %       当发展趋势存在其它的变化，采用简单移动平均法就会产生较大的预测偏差和滞后。
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % @function 指数平滑法
    %       输入观察序列y，计算alpha 取0.1、0.2、...0.9、1 的结果，依据均方差最小原则，输出最优解
    %
    % @param y：
    %       观测序列，列向量nx1
    %
    % @return prediction：
    %       数据矩阵 prediction = 3x3 
    %       输出一、二、三次指数模型最优解
    %           alpha1，均方差1，y1（t+1）
    %           alpha2，均方差2，y2（t+1）
    %           alpha3，均方差3，y3（t+1）
    %
    % @return params:
    %       元组 params
    %       输出二、三次指数模型的系数
    %           a2,b2
    %           a3,b3,c3
    %       对应模型：
    %       y(t+T) = a2 + b2*T;
    %       y(t+T) = a3 + b3*T + c3*T^2;
    
    output = zeros([3,length(y)]); % 第t+1项
    S = zeros([3,length(y)]); % 均方差
    St_end = zeros(3,10); % 第t个平滑项
    len=length(y);
    for alpha = 0.1:0.1:1
        St1 = zeros([1,len(y)]); % 一次平滑项
        St2 = zeros([1,len(y)]); % 二次平滑项
        St3 = zeros([1,len(y)]); % 三次平滑项
        St1_0 = sum(y(1:4))/4; 
        St2_0 = St1_0;
        St3_0 = St2_0;  % 初始值，取前几项均值
        for i = 1:len(y)
            if i==1
                St1(i) = alpha*y(i)+(1-alpha)*St1_0;
                St2(i) = alpha*St1(i)+(1-alpha)*St2_0;
                St3(i) = alpha*St2(i)+(1-alpha)*St3_0;
            else
                St1(i) = alpha*y(i)+(1-alpha)*St1(i-1);
                St2(i) = alpha*St1(i)+(1-alpha)*St2(i-1);
                St3(i) = alpha*St2(i)+(1-alpha)*St3(i-1);
            end
        end
        yt1 = St1(end);
        yt2 = (1+1/(1-alpha))*St1(end)-St2(end)/(1-alpha);
        yt3 = ((3-3*alpha+alpha^2)/(1-alpha)^2)*St1(end)-((3-alpha)/(1-alpha)^2)*St2(end)+St3(end)/(1-alpha)^2;
        output(:,alpha*10)=[yt1;yt2;yt3];
        
        S(1,alpha*10) = (sum((St1-y').^2)/len)^0.5;
        S(2,alpha*10) = (sum((St2-y').^2)/len)^0.5;
        S(3,alpha*10) = (sum((St3-y').^2)/len)^0.5;
  
        St_end(:,alpha*10)= [St1;St2;St3];
    end
    
    a2 = 2*St_end(:,1) - St_end(:,2);
    b2 = alpha/(1-alpha)*(St_end(:,1) - St_end(:,2));
    param2 = [a2;b2];
    
    a3 = 3*St_end(:,1) - 3*St_end(:,2) + St_end(:,3);
    b3 = alpha/(2*(1-alpha)^2) * ( (5-6*alpha)*St_end(:,1)- 2*(5-4*alpha)*St_end(:,1) +(4-3*alpha)*St_end(:,1));
    c3 = alpha^2/(2*(1-alpha)^2) * (St_end(:,1) - 2*St_end(:,2) + St_end(:,3));
    param3 = [a3;b3;c3];
    
    % 绘图比较
    figure(98)
    title('指数平滑法标准差');hold on;
    for i=1:3
        plot(0.1:0.1:1,S(i,:),'-o');
    end
    legend('一次','二次','三次')
    xlabel('alpha');ylabel('S');
    
    prediction = zeros([3,3]); 
    % 数据矩阵
    % prediction = 3x3 
    % 输出一、二、三次指数模型最优解
    %   alpha1，均方差1，y1（t+1）
    %   alpha2，均方差2，y2（t+1）
    %   alpha3，均方差3，y3（t+1）
    % 
    for i=1:3
        [x,y]=find(S==min(S(i,:)));
        prediction(i,1) = x/10;
        prediction(i,2) = min(S(i,:));
        prediction(i,3) = output(x,y);
    end
    params = {param2(:,prediction(2,1)*10)';param3(:,prediction(3,1)*10)'};
    % 元组
    % params
    %   a2，b2
    %   a3, b3, c3
    %
    
end

function DifferenceExponentialSmoothing(y)
%
%
%


end




