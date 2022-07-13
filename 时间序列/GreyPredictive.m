function GreyPredictive(x0,mode)
%%%%%%%%%%%%%%%%%%%   ��ɫԤ��ģ��  %%%%%%%%%%%%%%%%%%%%%%%
% ���÷�Χ��    
%       ��ģ��ʹ�õĲ���ԭʼ���ݵ����У�����
%       ���ɵ��������С�������ϵ��Grey Model.����ԭʼ
%       �������ۼ����ɣ��������������ɣ��õ����Ƶ�ָ��
%       �����ٽ��н�ģ�ķ�����
% �ŵ㣺       
%       �ڴ�����ٵ�����ֵ���ݣ�����Ҫ���ݵ�����
%       �ռ��㹻�󣬾��ܽ����ʷ�����١����е���������
%       ���ɿ��Ե͵����⣬�ܽ��޹��ɵ�ԭʼ���ݽ�������
%       �õ����ɽ�ǿ���������С�
% ȱ�㣺       
%       ֻ�������ж��ڵ�Ԥ�⣬ֻ�ʺϽ�����ָ����
%       ����Ԥ�⡣
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @illustrate:
%       �Ҷ�Ԥ���ģ�;�Ϊ΢�ַ���ģ��
%       �޷���ֵ������������ӡ
%
% @param x0:
% @param mode:Ԥ�ⷽ��
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
            disp("�����������");
    end

end

function GM1_1(x0)
% @illustrate:
%       GM(1,1)��ʾģ���� 1 ��΢�ַ��̣���ֻ�� 1 �������Ļ�ɫģ��
%       GM(1,1)ģ�������ھ��н�ǿָ�����ɵ����У�ֻ�����������ı仯���̣�
%       ���ڷǵ����İڶ���չ���л��б��͵� S �����У����Կ��ǽ��� GM(2,1)��DGM �� Verhulst ģ��
%
% @function: 
%       �ú�����������ݵļ��飬δ������ݸ�����ʱ�޷�ʹ��
%       �����ӡ�����
%           1.΢�ַ��̵Ľ� x(t)
%           2.�в�
%           3.������
%           4.����ƫ��ֵ
%       ����x(t)�Ƕ�1���ۼ����н���Ԥ�⣬k=t+1
%       ʵ��Ԥ��ֵ prediction(k+1) = x(k+1) - x(k) = x(t)-x(t-1)
%
% @param x0: 
%       �۲����У���������nx1
%
    n = length(x0);
    % ���ݼ���
    lamda = x0(1:n-1)./x0(2:n); 
    if min(lamda)<exp(-2/(n+1))||max(lamda)>exp(2/(n+1))
        disp('δ������ݸ�����')
        disp('��ȡ�ʵ�����ʹy=y+c������ݸ�����')
        return
    end
    % ģ�����
    y1 = cumsum(x0); % �ۼ�
    
    Y = x0(2:n);
    B = [0.5*(y1(1:n-1)+y1(2:n)),ones(n-1,1)];
    u=B\Y ; % ��ϲ���u=[a,b]
    
    %%%% dsolve������΢�ַ���
    syms x(t);
    x = dsolve(diff(x,t)+u(1)*x==u(2),x(0)==x0(1)); %��΢�ַ��̵ķ��Ž�
    %%%%
    xt = vpa(x,6); %��С����ʽ��ʾ΢�ַ��̵Ľ�
    disp("΢�ַ��̵Ľ� x(t)= "+ xt);
    disp("prediction(k+1) = x(k+1) - x(k) = x(t)-x(t-1) ");
    disp("");
    
    y_p=subs(x,t,0:n-1); %����֪���ݵ�Ԥ��ֵ
    y_p=double(y_p); % ������ת������ֵ���ͣ������޷����������
    
    % ������㣬��ԭ����
    x0_p=[x0(1),diff(y_p)]; % Ԥ��ֵ
    
    epsilon=x0'-x0_p; % ����в�
    disp("�в� = "+ epsilon);
    delta=abs(epsilon./x0'); % ����������
    disp("������ = "+ delta);
    rho=1-(1-0.5*u(1))/(1+0.5*u(1))*lamda'; % ���㼶��ƫ��ֵ��u(1)=a
    disp("����ƫ��ֵ = "+ rho);
end

function GM2_1(x0)
% @illustrate:
%       ģ�͹��ܣ�����ԭ������
%
% @function: 
%       �����ӡ�����
%           1.΢�ַ��̵Ľ� x(t)
%           2.�в�
%           3.������
%       ����x(t)�Ƕ�1���ۼ����н���Ԥ�⣬k=t+1
%       ʵ��Ԥ��ֵ prediction(k+1) = x(k+1) - x(k) = x(t)-x(t-1)
%
% @param x0: 
%       �۲����У���������nx1
%   
    x0 = x0';
    n=length(x0); 
    
    x1=cumsum(x0); %����1���ۼ�����
    z=0.5*(x1(2:end)+x1(1:end-1))'; %�����ֵ��������
    
    Y=diff(x0)'; %����1���ۼ�����
    B=[-x0(2:end)', -z ,ones(n-1,1)];
    u=B\Y ;%��С���˷���ϲ���
    
    syms x(t)
    x=dsolve(diff(x,2)+u(1)*diff(x)+u(2)*x==u(3),x(0)==x1(1),x(5)==x1(6));
    %����Ž�
    xt=vpa(x,6) ;%��ʾС����ʽ�ķ��Ž�
    disp("΢�ַ��̵Ľ� x(t)= "+ xt);
    disp("prediction(k+1) = x(k+1) - x(k) = x(t)-x(t-1) ");
    disp("");
    
    y_p=subs(x,t,0:n-1); %����֪���ݵ�1���ۼ����е�Ԥ��ֵ
    y_p=double(y_p) ;%������ת������ֵ���ͣ������޷����������
    
    x0_P = [y_p(1),diff(y_p)]; % ����֪���ݵ��Ԥ��ֵ
    epsilon = x0-x0_P; % ��в�
    disp("�в� = "+ epsilon);
    delta = abs(epsilon./x0); % �������
    disp("������ = "+ delta);
end

function DGM2_1(x0)
% @illustrate:
%       ģ�͹��ܣ�����ԭ������
%
% @function: 
%       �����ӡ�����
%           1.΢�ַ��̵Ľ� x(t)
%           2.�в�
%           3.������
%       ����x(t)�Ƕ�1���ۼ����н���Ԥ�⣬k=t+1
%       ʵ��Ԥ��ֵ prediction(k+1) = x(k+1) - x(k) = x(t)-x(t-1)
%
% @param x0: 
%       �۲����У���������nx1
%
    x0 = x0';
    n=length(x0);
    
    Y=diff(x0)'; %��1���ۼ����У���1����ǰ���
    B=[-x0(2:end)',ones(n-1,1)];
    u=B\Y ;%��С���˷���ϲ���
    
    syms x(t); % t=0 -> k=1
    dx = diff(x); %����һ�׺Ͷ��׵���
    d2x = diff(x,2); 
    
    %�����΢�ַ��̷��Ž�
    x = dsolve(d2x+u(1)*dx==u(2),x(0)==x0(1),dx(0)==x0(1));
    xt = vpa(x,6); %��ʾС����ʽ�ķ��Ž�
    disp("΢�ַ��̵Ľ� x(t)= "+ xt);
    disp("prediction(k+1) = x(k+1) - x(k) = x(t)-x(t-1) ");
    disp("");
    
    y_p = subs(x,t,0:n-1); % ����֪���ݵ�1���ۼ����е�Ԥ��ֵ
    y_p = double(y_p); % ������ת������ֵ���ͣ������޷����������
    
    x0_P = [y_p(1),diff(y_p)]; % ����֪���ݵ��Ԥ��ֵ
    epsilon = x0-x0_P; % ��в�
    disp("�в� = "+ epsilon);
    delta = abs(epsilon./x0); % �������
    disp("������ = "+ delta);
    
end

function Verhulst(x0)
% @illustrate:
%       Verhulst ģ����Ҫ�����������б���״̬�Ĺ��̣�
%       �� S �ι��̣��������˿�Ԥ�⡢������������ֳԤ�⼰��Ʒ��������Ԥ��ȡ�
%
% @function: 
%       �����ӡ�����
%           1.΢�ַ��̵Ľ� x(t)
%           2.�в�
%           3.������
%       ����x(t)�Ƕ�1���ۼ����н���Ԥ�⣬k=t+1
%       ʵ��Ԥ��ֵ prediction(k+1) = x(k+1) - x(k) = x(t)-x(t-1)
%
% @param x0: 
%       �۲����У���������nx1
%
    x0 = x0';
    n = length(x0);
    
    x1 = cumsum(x0); %��1���ۼ�����
    z = 0.5*(x1(2:n)+x1(1:n-1)); %��x1�ľ�ֵ��������
    
    Y = x0(2:end)';
    B = [-z',z'.^2];
    u = B\Y; %���Ʋ���a,b��ֵ
    
    %%%% dsolve������΢�ַ���
    syms x(t);
    x = dsolve(diff(x)+u(1)*x==u(2)*x^2,x(0)==x0(1)); %����Ž�
    %%%%
    % vpa(x,d) uses at least d significant digits, instead of the value of digits.
    xt = vpa(x,6); % ����ΪС����ʽ
    disp("΢�ַ��̵Ľ� x(t)= "+ xt);
    disp("prediction(k+1) = x(k+1) - x(k) = x(t)-x(t-1) ");
    disp("");
    
    % subs��x,t,a���� x�е� t���� a
    y_p = subs(x,t,0:n-1); % ����֪���ݵ�1���ۼ����е�Ԥ��ֵ
    y_p = double(y_p);% ������ת������ֵ���ͣ������޷����������
    
    % diff��ֺ�����һ���y_p(1)����
    x0_p = [y_p(1),diff(y_p)] ;% ����֪���ݵ��Ԥ��ֵ
    
    epsilon = x0-x0_p; % ��в�
    disp("�в� = "+ epsilon);
    delta = abs(epsilon./x0); %��������
    disp("������ = "+ delta);
    
    % writematrix([x0',x0_p',epsilon',delta'], 'data15_5.xlsx')
end



