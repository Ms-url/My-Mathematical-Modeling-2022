
% 评价结果的指标
%
% 评价分类结果：
%       精准度、混淆矩阵、精准率、召回率、F1 Score、ROC曲线，AUC值等
% 评价回归结果：
%       MSE、RMSE、MAE、R Squared，调整 R Squared
% 

% 混淆矩阵:NxN
%           
%       TP    FP
%       FN    TN
% 
%       T: 预测正确， 
%       P: 预测为正向样本，

% 查准率（精度）：
%       Precision=TP/(TP+FP)

% 召回率 (查全率)：
%       Recall=TP/(TP+FN)

% F1 score=2P*R/(P+R)。
%       P为查准率，R为查全率。F1score就是在查准率与查全率之间寻找一个平衡点。

% ROC曲线：
%       横轴假正率（FP率），纵轴为真正率（TP率）。

% AUC值：
%       ROC曲线与横轴之间围的面积。

% 基尼系数：
%       判定方法：基尼系数应大于60%，就算好模型。
%       原理：基尼系数经常用于分类问题，其可以直接从AUC中得到。其公式为：
%       Gini ＝ 2 * AUC － 1














