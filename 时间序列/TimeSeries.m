function [prediction,params] = TimeSeries(y,mode)
    % function
    % @param y:�۲����У�������
    % @param mode��Ԥ�ⷽ��
    % @return
    % 
    
    prediction = MovingAverage(y);
    [prediction,params] = ExponentialSmoothing(y);
    [prediction,params] = DifferenceExponentialSmoothing(y);

end

function prediction = MovingAverage(y)
    
    %%%%%%%%   ���ƶ�ƽ����  %%%%%%%%
    % ���÷�Χ��
    %       ��Ԥ��Ŀ��Ļ�����������ĳһˮƽ���²���ʱ������һ�μ��ƶ�ƽ����������Ԥ��ģ��
    %       ���ƶ�ƽ����ֻ�ʺ�������Ԥ�⣬������Ԥ��Ŀ��ķ�չ���Ʊ仯����������
    % ��ȱ�㣺
    %       ���Ŀ��ķ�չ���ƴ��������ı仯�����ü��ƶ�ƽ�����ͻ�����ϴ��Ԥ��ƫ����ͺ�
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % function�����ƶ�ƽ����
    %       ����۲�����y������ƽ���� N ��2��y����2/3 ��Ԥ���������ݾ�������Сԭ��������Ž�
    %
    % @param y:�۲����У�������nx1
    %
    % @return prediction��
    %       ������1x3
    %           N(ƽ������), ����� y(t+1)
    %
    len=length(y);
    S = []; % ������
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
    
    % ��ͼ�Ƚ�
    figure(99)
    title('���ƶ�ƽ������׼��')
    plot(S(:,1),S(:,2),'-o');
    xlabel('N');ylabel('S');
    
    [x,y]=find(S==min(S(:,2))); 
    prediction = [N,min(S(:,2)),output(x,y)];  % ע��prediction �Ľṹ
    fprintf("N = %d\nprediction = %f",output(x,1),prediction);
end

function [prediction,params] = ExponentialSmoothing(y)
    %%%%%%%    ָ��ƽ����   %%%%%%%%%
    % ���÷�Χ��
    %       ָ��ƽ��������ƽ�������Ĳ�ͬ���ַ�Ϊһ��ָ��ƽ����������ָ��ƽ����������ָ��ƽ����
    % ��ȱ�㣺
    %       ����չ���ƴ��������ı仯�����ü��ƶ�ƽ�����ͻ�����ϴ��Ԥ��ƫ����ͺ�
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % @function ָ��ƽ����
    %       ����۲�����y������alpha ȡ0.1��0.2��...0.9��1 �Ľ�������ݾ�������Сԭ��������Ž�
    %
    % @param y��
    %       �۲����У�������nx1
    %
    % @return prediction��
    %       ���ݾ��� prediction = 3x3 
    %       ���һ����������ָ��ģ�����Ž�
    %           alpha1��������1��y1��t+1��
    %           alpha2��������2��y2��t+1��
    %           alpha3��������3��y3��t+1��
    %
    % @return params:
    %       Ԫ�� params
    %       �����������ָ��ģ�͵�ϵ��
    %           a2,b2
    %           a3,b3,c3
    %       ��Ӧģ�ͣ�
    %       y(t+T) = a2 + b2*T;
    %       y(t+T) = a3 + b3*T + c3*T^2;
    
    output = zeros([3,length(y)]); % ��t+1��
    S = zeros([3,length(y)]); % ������
    St_end = zeros(3,10); % ��t��ƽ����
    len=length(y);
    for alpha = 0.1:0.1:1
        St1 = zeros([1,len(y)]); % һ��ƽ����
        St2 = zeros([1,len(y)]); % ����ƽ����
        St3 = zeros([1,len(y)]); % ����ƽ����
        St1_0 = sum(y(1:4))/4; 
        St2_0 = St1_0;
        St3_0 = St2_0;  % ��ʼֵ��ȡǰ�����ֵ
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
    
    % ��ͼ�Ƚ�
    figure(98)
    title('ָ��ƽ������׼��');hold on;
    for i=1:3
        plot(0.1:0.1:1,S(i,:),'-o');
    end
    legend('һ��','����','����')
    xlabel('alpha');ylabel('S');
    
    prediction = zeros([3,3]); 
    % ���ݾ���
    % prediction = 3x3 
    % ���һ����������ָ��ģ�����Ž�
    %   alpha1��������1��y1��t+1��
    %   alpha2��������2��y2��t+1��
    %   alpha3��������3��y3��t+1��
    % 
    for i=1:3
        [x,y]=find(S==min(S(i,:)));
        prediction(i,1) = x/10;
        prediction(i,2) = min(S(i,:));
        prediction(i,3) = output(x,y);
    end
    params = {param2(:,prediction(2,1)*10)';param3(:,prediction(3,1)*10)'};
    % Ԫ��
    % params
    %   a2��b2
    %   a3, b3, c3
    %
    
end

function DifferenceExponentialSmoothing(y)




end




