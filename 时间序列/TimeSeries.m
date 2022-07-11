function [prediction, params, W] = TimeSeries(y,mode,varargin)
    % function
<<<<<<< HEAD
    %    ����
    %
    % @param y:�۲����У�������
    % @param mode��Ԥ�ⷽ��
    %       1 -> ���ƶ�ƽ����
    %       2 -> 
    %
=======
    % @param y:观测序列，列向量
    % @param mode：预测方法
>>>>>>> f92c044ecdee66bf80d6b2c4905686b79504e842
    % @return
    % 
    prediction = 0;
    params = 0;
    W = 0;
    N = varargin(1);
    switch(mode)
        case 1
            disp("���ƶ�ƽ����");
            prediction = MovingAverage(y);
        case 2
            disp("ָ��ƽ����");
            [prediction,params] = ExponentialSmoothing(y);
        case 3
            disp("���ָ��ƽ����");
            prediction = DifferenceExponentialSmoothing(y);
        case 4
            disp("����Ӧ�˲���");
            [prediction, W] = AdaptiveFiltering(y,N);
        otherwise
             disp("�����������");
    end
    
end

function prediction = MovingAverage(y)
    
<<<<<<< HEAD
%%%%%%%%%%%%   ���ƶ�ƽ����  %%%%%%%%%%%%
%     ���÷�Χ��
%           ��Ԥ��Ŀ��Ļ�����������ĳһˮƽ���²���ʱ������һ�μ��ƶ�ƽ����������Ԥ��ģ��
%           ���ƶ�ƽ����ֻ�ʺ�������Ԥ�⣬������Ԥ��Ŀ��ķ�չ���Ʊ仯����������
%     ��ȱ�㣺
%           ���Ŀ��ķ�չ���ƴ��������ı仯�����ü��ƶ�ƽ�����ͻ�����ϴ��Ԥ��ƫ����ͺ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % function�����ƶ�ƽ����
    %       ����۲�����y������ƽ���� N ��2��y����2/3 ��Ԥ���������ݾ�������Сԭ��������Ž�
=======
    %%%%%%%%   简单移动平均法  %%%%%%%%
    % 适用范围：
    %       当预测目标的基本趋势是在某一水平上下波动时，可用一次简单移动平均方法建立预测模型
    %       简单移动平均法只适合做近期预测，而且是预测目标的发展趋势变化不大的情况。
    % 优缺点：
    %       如果目标的发展趋势存在其它的变化，采用简单移动平均法就会产生较大的预测偏差和滞后。
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % function：简单移动平均法
    %       输入观察序列y，计算平移项 N 由2到y长度2/3 的预测结果，依据均方差最小原则，输出最优解
>>>>>>> f92c044ecdee66bf80d6b2c4905686b79504e842
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
<<<<<<< HEAD
%%%%%%%%%%%%    ָ��ƽ����   %%%%%%%%%%%%
%     ���÷�Χ��
%           ָ��ƽ��������ƽ�������Ĳ�ͬ���ַ�Ϊһ��ָ��ƽ����������ָ��ƽ����������ָ��ƽ����
%     ��ȱ�㣺
%           ����չ���ƴ��������ı仯�����ü��ƶ�ƽ�����ͻ�����ϴ��Ԥ��ƫ����ͺ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % @function ָ��ƽ����
    %       ����۲�����y������alpha ȡ0.1��0.2��...0.9��1 �Ľ�������ݾ�������Сԭ��������Ž�
=======
    %%%%%%%    指数平滑法   %%%%%%%%%
    % 适用范围：
    %       指数平滑法根据平滑次数的不同，又分为一次指数平滑法、二次指数平滑法和三次指数平滑法
    % 优缺点：
    %       当发展趋势存在其它的变化，采用简单移动平均法就会产生较大的预测偏差和滞后。
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % @function 指数平滑法
    %       输入观察序列y，计算alpha 取0.1、0.2、...0.9、1 的结果，依据均方差最小原则，输出最优解
>>>>>>> f92c044ecdee66bf80d6b2c4905686b79504e842
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
<<<<<<< HEAD
    for j = 1:10
        alpha = j/10;
        St1 = zeros([1,len(y)]); % һ��ƽ����
        St2 = zeros([1,len(y)]); % ����ƽ����
        St3 = zeros([1,len(y)]); % ����ƽ����
=======
    for alpha = 0.1:0.1:1
        St1 = zeros([1,len(y)]); % 一次平滑项
        St2 = zeros([1,len(y)]); % 二次平滑项
        St3 = zeros([1,len(y)]); % 三次平滑项
>>>>>>> f92c044ecdee66bf80d6b2c4905686b79504e842
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
        output(:,j)=[yt1;yt2;yt3]; % ÿһ�ж�Ӧһ��alpha
        
        S(1,j) = (sum((St1-y').^2)/len)^0.5;
        S(2,j) = (sum((St2-y').^2)/len)^0.5;
        S(3,j) = (sum((St3-y').^2)/len)^0.5;
  
        St_end(:,j)= [St1;St2;St3];% ÿһ�ж�Ӧһ��alpha
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
        prediction(i,1) = x/10; %%%%%%%%%%%%%%%% ��չ����
        prediction(i,2) = min(S(i,:));
        prediction(i,3) = output(x,y);
    end
<<<<<<< HEAD
    params = {param2(:,prediction(2,1)*10)';param3(:,prediction(3,1)*10)'}; %%%%%%%%%%%%% ��չ����
    % Ԫ��
=======
    params = {param2(:,prediction(2,1)*10)';param3(:,prediction(3,1)*10)'};
    % 元组
>>>>>>> f92c044ecdee66bf80d6b2c4905686b79504e842
    % params
    %   a2，b2
    %   a3, b3, c3
    %
    
end

<<<<<<< HEAD
function prediction = DifferenceExponentialSmoothing(y)
%%%%%%%%%%%%%%   ���ָ��ƽ����  %%%%%%%%%%%%%%%
%     ��ʱ�����еı䶯����ֱ������ʱ����һ��ָ��ƽ�������
%     ���ͺ�ƫ���ԭ���������ݲ�����ģ��Ҫ����ˣ�����Ҳ���Դ����ݱ任�ĽǶ���
%     ���ǸĽ���ʩ����������ָ��ƽ������ǰ�ȶ�������һЩ�����ϵĴ���ʹ֮���ʺ���
%     һ��ָ��ƽ��ģ�ͣ��Ժ��ٶ��������������ϵķ��ش���ʹ֮�ָ�Ϊԭ��������̬��
%     ��ַ����Ǹı����ݱ䶯���Ƶļ��׷�����
%
%     ��ʱ�����г�ֱ������ʱ��������һ�ײ��ָ��ƽ��ģ����Ԥ�⡣
%     ��ʱ�����г��ֶ�����������ʱ�����ö��ײ��ָ��ƽ��ģ����Ԥ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % function �����ָ��ƽ����
    %       ��������Ƚ�
    %       
    % @param y: �۲����У���������nx1
    %
    % @return prediction��
    %
    % @return params��
    %
    len = length(y);
    % һ�ײ��
    d_y = zeros([10,len]);
    d_y1 = zeros([10,len+1]);
    % ���ײ��
    d2_y = zeros([10,len]);
    d2_y2 = zeros([10,len+1]);
    % Ԥ��ֵ
    yt1 = zeros([10,len+1]);
    yt2 = zeros([10,len+1]);
    % Ԥ��ֵ
    % ÿһ��Ϊ��Ӧ��һ��alpha�µ�һ��Ԥ��ֵ
    %
    for i = 1:10
        alpha = i/10;
        for t = 2:len
            d_y(i,t) = y(t)-y(t-1);
            d2_y(i,t) = d_y(i,t) - d_y(i,t-1);
            
            d_y1(i,t+1) = alpha*d_y(i,t)+(1-alpha)*d_y1(i,t);
            d2_y2(i,t+1) = alpha*d2_y(i,t)+(1-alpha)*d2_y2(i,t);
            
            yt1(i,t+1) = d_y1(i,t+1) + y(t);
            yt2(i,t+1) = d2_y2(i,t+1) + d_y(i,t) + y(t);
        end
    end
    S1 = (yt1(:,3:end-1)-y(3:end)').^2./(len-2);
    S2 = (yt2(:,3:end-1)-y(3:end)').^2./(len-2);
    S1 = sum(S1'); % �����������
    S2 = sum(S2'); % �����������
    
    % ��ͼ�Ƚ�
    figure(97);
    hold on;
    for i = 1:2
        plot(1:10,S1,'-o');
    end
    legend('һ�β��','���β��')
    
    [~,y1]=find(S1==min(S1));
    [~,y2]=find(S2==min(S2));
    prediction = [y1/10,min(S1),yt1(y1,end);y2/10,min(S2),yt2(y2,end)];
    % ���ݾ���2x3
    %   alpha1, ������1��y(t+1)
    %   alpha2, ������2��y(t+1)
    
end
=======
function DifferenceExponentialSmoothing(y)
%
%
>>>>>>> f92c044ecdee66bf80d6b2c4905686b79504e842

function [predication, W] = AdaptiveFiltering(y,N)
%%%%%%%%%%  ����Ӧ�˲��� %%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   N=2ʱ
%   y(t+1) = w(1)*y(t) + w(2)*y(t-1)

    if N==0
        N=2; % Ȩ�ظ���
    end
    len=length(y); 
    k=0.9; % ѧϰ����
    
    diff=10000;
    W = ones(1,N)/N;
    while abs(diff)>0.0001
        diff=[];
        for j=N+1:len-1
            predication(j)= W * y(j-1:-1:j-N)';
            d = abs(y(j)-predication(j));
            diff=[diff,d];
            W = W + 2 * k * d *y(j-1:-1:j-N);
        end
        diff=max(diff);
    end
    fprintf("����ֵ = %f",diff);
end



