function [cMatrix,Precision,Recall,F1_score]=categoricalEvaluationIndicators(x,y,label)
% 
% function：计算分类评价指标
% param x：预测数据，列向量
% param y：真实数据，列向量
% param label：标签，列向量，纯数字
% 
% return Precision：查准率（精度）
% return Recall：召回率 (查全率)
% return F1_score：
%

   
    n = length(label);
    
    cMatrix = confusionMatrix(x,y,label);
    % 混淆矩阵的第i行的和，等于所有预测为label（i）的样本，第j列的和，等于所有实际为label（j）的样本。
    forcast_sample = sum(transpose(cMatrix));
    ture_sample = sum(cMatrix);
    TP = diag(cMatrix)';
    
    % 查准率（精度）： Precision=TP/(TP+FP)
    % 召回率 (查全率)：Recall=TP/( TP+ FN)
    % F1 score= 2P*R / (P+R)。
    Precision = TP./forcast_sample;
    Recall = TP./ture_sample;
    F1_score = 2* Precision* Recall/(Precision + Recall);
    
%     figure(99)
%     % ROC曲线：横轴假正率（FPR率），纵轴为真正率（TPR率）。
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
% function ：获取混淆矩阵
% param x: 预测数据, 列向量
% param y：真实数据，列向量
% param label：标签，列向量，纯数字
% illustrate：
%       混淆矩阵的第i行的和，等于所有预测为label（i）的样本，
%                第j列的和，等于所有实际为label（j）的样本。
% 
    n = length(label);
    cMatrix = zeros(n,n);
  
    for i = 1:n
       forecast_t = y(x==label(i)); % 预测为label(i)的样本对应的实际样本
       for j=1:n
          temp = sum(forecast_t==label(j)); % 统计预测为label(i)下，实际样本中label（j）的数量
          cMatrix(i,j) = temp;
       end       
    end    

end

% 评价结果的指标
%
% 评价分类结果：
%       精准度、混淆矩阵、精准率、召回率、F1 Score、ROC曲线，AUC值等
% 评价回归结果：
%       MSE、RMSE、MAE、R Squared，调整 R Squared
% 

% 混淆矩阵:NxN
%   二元        
%       TP    FP
%       FN    TN
% 
%       T: 预测正确， 
%       P: 预测为正向样本，

% 查准率（精度）：
%       Precision=TP/(TP+FP)

% 召回率 (查全率)：
%       Recall=TP/( TP+ FN)

% F1 score= 2P*R / (P+R)。
%       P为查准率，R为查全率。F1score就是在查准率与查全率之间寻找一个平衡点。

% ROC曲线：
%       横轴假正率（FPR），纵轴为真正率（TPR）。

% AUC值：
%       ROC曲线与横轴之间围的面积。

% 基尼系数：
%       判定方法：基尼系数应大于60%，就算好模型。
%       原理：基尼系数经常用于分类问题，其可以直接从AUC中得到。其公式为：
%       Gini ＝ 2 * AUC － 1

