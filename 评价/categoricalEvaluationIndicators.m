function [cMatrix,Precision,Recall,F1_score]=categoricalEvaluationIndicators(x,y,label)
% 
% function�������������ָ��
% param x��Ԥ�����ݣ�������
% param y����ʵ���ݣ�������
% param label����ǩ����������������
% 
% return Precision����׼�ʣ����ȣ�
% return Recall���ٻ��� (��ȫ��)
% return F1_score��
%

   
    n = length(label);
    
    cMatrix = confusionMatrix(x,y,label);
    % ��������ĵ�i�еĺͣ���������Ԥ��Ϊlabel��i������������j�еĺͣ���������ʵ��Ϊlabel��j����������
    forcast_sample = sum(transpose(cMatrix));
    ture_sample = sum(cMatrix);
    TP = diag(cMatrix)';
    
    % ��׼�ʣ����ȣ��� Precision=TP/(TP+FP)
    % �ٻ��� (��ȫ��)��Recall=TP/( TP+ FN)
    % F1 score= 2P*R / (P+R)��
    Precision = TP./forcast_sample;
    Recall = TP./ture_sample;
    F1_score = 2* Precision* Recall/(Precision + Recall);
    
%     figure(99)
%     % ROC���ߣ���������ʣ�FPR�ʣ�������Ϊ�����ʣ�TPR�ʣ���
%     for i=1:n
%         subplot(ceil(n/3),min(3,n),i);
%         FP = forcast_sample(i)-TP(i);
%         TPR = TP(i)/ ture_sample(i);
%         FPR = FP / sum(ture_sample)-ture_sample(i);
%         plot(FPR,TPR);
%     end

%      https://blog.csdn.net/Manketon/article/details/44587829
    
    
end


function cMatrix = confusionMatrix(x,y,label)
% function ����ȡ��������
% param x: Ԥ������, ������
% param y����ʵ���ݣ�������
% param label����ǩ����������������
% illustrate��
%       ��������ĵ�i�еĺͣ���������Ԥ��Ϊlabel��i����������
%                ��j�еĺͣ���������ʵ��Ϊlabel��j����������
% 
    n = length(label);
    cMatrix = zeros(n,n);
  
    for i = 1:n
       forecast_t = y(x==label(i)); % Ԥ��Ϊlabel(i)��������Ӧ��ʵ������
       for j=1:n
          temp = sum(forecast_t==label(j)); % ͳ��Ԥ��Ϊlabel(i)�£�ʵ��������label��j��������
          cMatrix(i,j) = temp;
       end       
    end    

end

% ���۽����ָ��
%
% ���۷�������
%       ��׼�ȡ��������󡢾�׼�ʡ��ٻ��ʡ�F1 Score��ROC���ߣ�AUCֵ��
% ���ۻع�����
%       MSE��RMSE��MAE��R Squared������ R Squared
% 

% ��������:NxN
%   ��Ԫ        
%       TP    FP
%       FN    TN
% 
%       T: Ԥ����ȷ�� 
%       P: Ԥ��Ϊ����������

% ��׼�ʣ����ȣ���
%       Precision=TP/(TP+FP)

% �ٻ��� (��ȫ��)��
%       Recall=TP/( TP+ FN)

% F1 score= 2P*R / (P+R)��
%       PΪ��׼�ʣ�RΪ��ȫ�ʡ�F1score�����ڲ�׼�����ȫ��֮��Ѱ��һ��ƽ��㡣

% ROC���ߣ�
%       ��������ʣ�FPR��������Ϊ�����ʣ�TPR����

% AUCֵ��
%       ROC���������֮��Χ�������

% ����ϵ����
%       �ж�����������ϵ��Ӧ����60%�������ģ�͡�
%       ԭ������ϵ���������ڷ������⣬�����ֱ�Ӵ�AUC�еõ����乫ʽΪ��
%       Gini �� 2 * AUC �� 1

