function [prediction, params, W] = TimeSeries(y,mode,varargin)
    % function
<<<<<<< HEAD
    %    包含
    %
    % @param y:观测序列，列向量
    % @param mode：预测方法
    %       1 -> 简单移动平均法
    %       2 -> 
    %
=======
    % @param y:瑙傛祴搴忓垪锛屽垪鍚戦噺
    % @param mode锛氶娴嬫柟娉�
>>>>>>> f92c044ecdee66bf80d6b2c4905686b79504e842
    % @return
    % 
    prediction = 0;
    params = 0;
    W = 0;
    N = varargin(1);
    switch(mode)
        case 1
            disp("简单移动平均法");
            prediction = MovingAverage(y);
        case 2
            disp("指数平滑法");
            [prediction,params] = ExponentialSmoothing(y);
        case 3
            disp("差分指数平滑法");
            prediction = DifferenceExponentialSmoothing(y);
        case 4
            disp("自适应滤波法");
            [prediction, W] = AdaptiveFiltering(y,N);
        otherwise
             disp("输入参数错误");
    end
    
end

function prediction = MovingAverage(y)
    
<<<<<<< HEAD
%%%%%%%%%%%%   简单移动平均法  %%%%%%%%%%%%
%     适用范围：
%           当预测目标的基本趋势是在某一水平上下波动时，可用一次简单移动平均方法建立预测模型
%           简单移动平均法只适合做近期预测，而且是预测目标的发展趋势变化不大的情况。
%     优缺点：
%           如果目标的发展趋势存在其它的变化，采用简单移动平均法就会产生较大的预测偏差和滞后。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % function：简单移动平均法
    %       输入观察序列y，计算平移项 N 由2到y长度2/3 的预测结果，依据均方差最小原则，输出最优解
=======
    %%%%%%%%   绠�鍗曠Щ鍔ㄥ钩鍧囨硶  %%%%%%%%
    % 閫傜敤鑼冨洿锛�
    %       褰撻娴嬬洰鏍囩殑鍩烘湰瓒嬪娍鏄湪鏌愪竴姘村钩涓婁笅娉㈠姩鏃讹紝鍙敤涓�娆＄畝鍗曠Щ鍔ㄥ钩鍧囨柟娉曞缓绔嬮娴嬫ā鍨�
    %       绠�鍗曠Щ鍔ㄥ钩鍧囨硶鍙�傚悎鍋氳繎鏈熼娴嬶紝鑰屼笖鏄娴嬬洰鏍囩殑鍙戝睍瓒嬪娍鍙樺寲涓嶅ぇ鐨勬儏鍐点��
    % 浼樼己鐐癸細
    %       濡傛灉鐩爣鐨勫彂灞曡秼鍔垮瓨鍦ㄥ叾瀹冪殑鍙樺寲锛岄噰鐢ㄧ畝鍗曠Щ鍔ㄥ钩鍧囨硶灏变細浜х敓杈冨ぇ鐨勯娴嬪亸宸拰婊炲悗銆�
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % function锛氱畝鍗曠Щ鍔ㄥ钩鍧囨硶
    %       杈撳叆瑙傚療搴忓垪y锛岃绠楀钩绉婚」 N 鐢�2鍒皔闀垮害2/3 鐨勯娴嬬粨鏋滐紝渚濇嵁鍧囨柟宸渶灏忓師鍒欙紝杈撳嚭鏈�浼樿В
>>>>>>> f92c044ecdee66bf80d6b2c4905686b79504e842
    %
    % @param y:瑙傛祴搴忓垪锛屽垪鍚戦噺nx1
    %
    % @return prediction锛�
    %       琛屽悜閲�1x3
    %           N(骞崇Щ椤规暟), 鍧囨柟宸紝 y(t+1)
    %
    len=length(y);
    S = []; % 鍧囨柟宸�
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
    
    % 缁樺浘姣旇緝
    figure(99)
    title('绠�鍗曠Щ鍔ㄥ钩鍧囨硶鏍囧噯宸�')
    plot(S(:,1),S(:,2),'-o');
    xlabel('N');ylabel('S');
    
    [x,y]=find(S==min(S(:,2))); 
    prediction = [N,min(S(:,2)),output(x,y)];  % 娉ㄦ剰prediction 鐨勭粨鏋�
    fprintf("N = %d\nprediction = %f",output(x,1),prediction);
end

function [prediction,params] = ExponentialSmoothing(y)
<<<<<<< HEAD
%%%%%%%%%%%%    指数平滑法   %%%%%%%%%%%%
%     适用范围：
%           指数平滑法根据平滑次数的不同，又分为一次指数平滑法、二次指数平滑法和三次指数平滑法
%     优缺点：
%           当发展趋势存在其它的变化，采用简单移动平均法就会产生较大的预测偏差和滞后。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % @function 指数平滑法
    %       输入观察序列y，计算alpha 取0.1、0.2、...0.9、1 的结果，依据均方差最小原则，输出最优解
=======
    %%%%%%%    鎸囨暟骞虫粦娉�   %%%%%%%%%
    % 閫傜敤鑼冨洿锛�
    %       鎸囨暟骞虫粦娉曟牴鎹钩婊戞鏁扮殑涓嶅悓锛屽張鍒嗕负涓�娆℃寚鏁板钩婊戞硶銆佷簩娆℃寚鏁板钩婊戞硶鍜屼笁娆℃寚鏁板钩婊戞硶
    % 浼樼己鐐癸細
    %       褰撳彂灞曡秼鍔垮瓨鍦ㄥ叾瀹冪殑鍙樺寲锛岄噰鐢ㄧ畝鍗曠Щ鍔ㄥ钩鍧囨硶灏变細浜х敓杈冨ぇ鐨勯娴嬪亸宸拰婊炲悗銆�
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % @function 鎸囨暟骞虫粦娉�
    %       杈撳叆瑙傚療搴忓垪y锛岃绠梐lpha 鍙�0.1銆�0.2銆�...0.9銆�1 鐨勭粨鏋滐紝渚濇嵁鍧囨柟宸渶灏忓師鍒欙紝杈撳嚭鏈�浼樿В
>>>>>>> f92c044ecdee66bf80d6b2c4905686b79504e842
    %
    % @param y锛�
    %       瑙傛祴搴忓垪锛屽垪鍚戦噺nx1
    %
    % @return prediction锛�
    %       鏁版嵁鐭╅樀 prediction = 3x3 
    %       杈撳嚭涓�銆佷簩銆佷笁娆℃寚鏁版ā鍨嬫渶浼樿В
    %           alpha1锛屽潎鏂瑰樊1锛寉1锛坱+1锛�
    %           alpha2锛屽潎鏂瑰樊2锛寉2锛坱+1锛�
    %           alpha3锛屽潎鏂瑰樊3锛寉3锛坱+1锛�
    %
    % @return params:
    %       鍏冪粍 params
    %       杈撳嚭浜屻�佷笁娆℃寚鏁版ā鍨嬬殑绯绘暟
    %           a2,b2
    %           a3,b3,c3
    %       瀵瑰簲妯″瀷锛�
    %       y(t+T) = a2 + b2*T;
    %       y(t+T) = a3 + b3*T + c3*T^2;
    
    output = zeros([3,length(y)]); % 绗瑃+1椤�
    S = zeros([3,length(y)]); % 鍧囨柟宸�
    St_end = zeros(3,10); % 绗瑃涓钩婊戦」
    len=length(y);
<<<<<<< HEAD
    for j = 1:10
        alpha = j/10;
        St1 = zeros([1,len(y)]); % 一次平滑项
        St2 = zeros([1,len(y)]); % 二次平滑项
        St3 = zeros([1,len(y)]); % 三次平滑项
=======
    for alpha = 0.1:0.1:1
        St1 = zeros([1,len(y)]); % 涓�娆″钩婊戦」
        St2 = zeros([1,len(y)]); % 浜屾骞虫粦椤�
        St3 = zeros([1,len(y)]); % 涓夋骞虫粦椤�
>>>>>>> f92c044ecdee66bf80d6b2c4905686b79504e842
        St1_0 = sum(y(1:4))/4; 
        St2_0 = St1_0;
        St3_0 = St2_0;  % 鍒濆鍊硷紝鍙栧墠鍑犻」鍧囧��
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
        output(:,j)=[yt1;yt2;yt3]; % 每一列对应一个alpha
        
        S(1,j) = (sum((St1-y').^2)/len)^0.5;
        S(2,j) = (sum((St2-y').^2)/len)^0.5;
        S(3,j) = (sum((St3-y').^2)/len)^0.5;
  
        St_end(:,j)= [St1;St2;St3];% 每一列对应一个alpha
    end
    
    a2 = 2*St_end(:,1) - St_end(:,2);
    b2 = alpha/(1-alpha)*(St_end(:,1) - St_end(:,2));
    param2 = [a2;b2];
    
    a3 = 3*St_end(:,1) - 3*St_end(:,2) + St_end(:,3);
    b3 = alpha/(2*(1-alpha)^2) * ( (5-6*alpha)*St_end(:,1)- 2*(5-4*alpha)*St_end(:,1) +(4-3*alpha)*St_end(:,1));
    c3 = alpha^2/(2*(1-alpha)^2) * (St_end(:,1) - 2*St_end(:,2) + St_end(:,3));
    param3 = [a3;b3;c3];
    
    % 缁樺浘姣旇緝
    figure(98)
    title('鎸囨暟骞虫粦娉曟爣鍑嗗樊');hold on;
    for i=1:3
        plot(0.1:0.1:1,S(i,:),'-o');
    end
    legend('涓�娆�','浜屾','涓夋')
    xlabel('alpha');ylabel('S');
    
    prediction = zeros([3,3]); 
    % 鏁版嵁鐭╅樀
    % prediction = 3x3 
    % 杈撳嚭涓�銆佷簩銆佷笁娆℃寚鏁版ā鍨嬫渶浼樿В
    %   alpha1锛屽潎鏂瑰樊1锛寉1锛坱+1锛�
    %   alpha2锛屽潎鏂瑰樊2锛寉2锛坱+1锛�
    %   alpha3锛屽潎鏂瑰樊3锛寉3锛坱+1锛�
    % 
    for i=1:3
        [x,y]=find(S==min(S(i,:)));
        prediction(i,1) = x/10; %%%%%%%%%%%%%%%% 扩展限制
        prediction(i,2) = min(S(i,:));
        prediction(i,3) = output(x,y);
    end
<<<<<<< HEAD
    params = {param2(:,prediction(2,1)*10)';param3(:,prediction(3,1)*10)'}; %%%%%%%%%%%%% 扩展限制
    % 元组
=======
    params = {param2(:,prediction(2,1)*10)';param3(:,prediction(3,1)*10)'};
    % 鍏冪粍
>>>>>>> f92c044ecdee66bf80d6b2c4905686b79504e842
    % params
    %   a2锛宐2
    %   a3, b3, c3
    %
    
end

<<<<<<< HEAD
function prediction = DifferenceExponentialSmoothing(y)
%%%%%%%%%%%%%%   差分指数平滑法  %%%%%%%%%%%%%%%
%     当时间序列的变动具有直线趋势时，用一次指数平滑法会出
%     现滞后偏差，其原因在于数据不满足模型要求。因此，我们也可以从数据变换的角度来
%     考虑改进措施，即在运用指数平滑法以前先对数据作一些技术上的处理，使之能适合于
%     一次指数平滑模型，以后再对输出结果作技术上的返回处理，使之恢复为原变量的形态。
%     差分方法是改变数据变动趋势的简易方法。
%
%     当时间序列呈直线增加时，可运用一阶差分指数平滑模型来预测。
%     当时间序列呈现二次曲线增长时，可用二阶差分指数平滑模型来预测
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % function ：差分指数平滑法
    %       输出最优先解
    %       
    % @param y: 观察序列，列项量，nx1
    %
    % @return prediction：
    %
    % @return params：
    %
    len = length(y);
    % 一阶差分
    d_y = zeros([10,len]);
    d_y1 = zeros([10,len+1]);
    % 二阶差分
    d2_y = zeros([10,len]);
    d2_y2 = zeros([10,len+1]);
    % 预测值
    yt1 = zeros([10,len+1]);
    yt2 = zeros([10,len+1]);
    % 预测值
    % 每一行为对应的一个alpha下的一组预测值
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
    S1 = sum(S1'); % 均方差，行向量
    S2 = sum(S2'); % 均方差，行向量
    
    % 绘图比较
    figure(97);
    hold on;
    for i = 1:2
        plot(1:10,S1,'-o');
    end
    legend('一次差分','二次差分')
    
    [~,y1]=find(S1==min(S1));
    [~,y2]=find(S2==min(S2));
    prediction = [y1/10,min(S1),yt1(y1,end);y2/10,min(S2),yt2(y2,end)];
    % 数据矩阵2x3
    %   alpha1, 均方差1，y(t+1)
    %   alpha2, 均方差2，y(t+1)
    
end
=======
function DifferenceExponentialSmoothing(y)
%
%
>>>>>>> f92c044ecdee66bf80d6b2c4905686b79504e842

function [predication, W] = AdaptiveFiltering(y,N)
%%%%%%%%%%  自适应滤波法 %%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   N=2时
%   y(t+1) = w(1)*y(t) + w(2)*y(t-1)

    if N==0
        N=2; % 权重个数
    end
    len=length(y); 
    k=0.9; % 学习常数
    
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
    fprintf("最大差值 = %f",diff);
end



