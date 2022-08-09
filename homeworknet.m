% Solve an Autoregression Problem with External Input with a NARX Neural Network
% Script generated by Neural Time Series app
% Created 04-Aug-2022 14:05:01
%
% This script assumes these variables are defined:
%
%   year - input time series.
%   data - feedback time series.

% 使用保存的网络模型（net）运行该脚本，即可进行一步预测

X = tonndata(year,true,false);
T = tonndata(data,true,false);

nets = removedelay(net2);
nets.name = [net2.name ' - Predict One Step Ahead'];
view(nets)
[xs,xis,ais,ts] = preparets(nets,X,{},T);
ys = nets(xs,xis,ais);
stepAheadPerformance = perform(nets,ts,ys);

