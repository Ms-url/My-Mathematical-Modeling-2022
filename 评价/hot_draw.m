
% 绘制相关性热力图
data_1 = data_mean./sum(data_mean.^2).^0.5;
 
relevance = corr(data_1(:,3:8));

figure(77)
string_name = txt(3:8);
xvalues = string_name;
yvalues = string_name;
h1 = heatmap(xvalues , yvalues, relevance, 'FontSize',10,'Colormap',hot);
h1.Title = '相关性热力图';
colorbar




