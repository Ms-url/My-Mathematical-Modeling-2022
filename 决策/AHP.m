function AHP()
% @illustrate:
%       层次分析法
%       AHP—一种将定性问题定量化、系统化、层次化的综合评判方法
%
% 步骤：
%   ① 分析系统中各元素之间的关系，建立系统的递阶层析结构；
%   ② 对同一层次的各元素关于上一层次中某一准则的重要性进行两两比较，构造两两比较判断矩阵；
%   ③ 由判断矩阵计算被比较元素对于该准则的相对权重；
%   ④ 计算各层元素对系统目标的合成权重，并进行排序
% 
% 分值选择（1-9）
% 采用1～9的比例表度的依据是：
%   ①心理学的实验表明，大多数人对不同事物在相同属性上的差别的分辨能力在5～9，采用1～9的标度反映了大多数人的判断能力；
%   ②大量的社会调查表明，1～9的比例标度早已被人们所熟悉和采用；
%   ③科学考察和实践表明，1～9的比例标度已完全能区分引起人们感觉差别的事物的各种属性
% 
% 一致性检查
%       CI = (lamda_max - n)/(n - 1) % lamda 为特征向量
% 只要 CI < 0.1就可认为判断矩阵具有满意的一致性。
disp('该函数的子函数包含AHP辅助函数：')
disp('    get_W()')
disp('    ConsistencyDetection()')

end

function [W,lamda_max] = get_W(a)
% @function：
%       输入分值向量a（行向量），构建比值矩阵A，
%       返回A的最大特征值lamda_max即其对应的特征向量W（规范化后）
%       当一致性检测未通过时，返回 空
% @illustrate：
%       AHP需要计算各层最大特征值对应的特征向量，
%       该函数用于辅助AHP完整过程的编写
%
% @param a: 分值向量，行向量，1xn
% @return W：规范化后的，比值矩阵最大特征值对应的特征向量
% @return lamda_max：比值矩阵的最大特征值
%
    len = length(a);
    A = zeros([len,len])+a;
    A = a'./A; % 比值矩阵
    [W,lamda] = eig(A); % 求矩阵A的全部特征值，构成对角阵lamda，并求A的特征向量构成W的列向量。
    [lamda_max ,col] = max(max(lamda));
    W = W(:,col); % 最大特征值对应的特征向量
    W = W./sum(W); % 规范化
    
    % 一致性检查
    CI = (lamda_max - len)/(len - 1);
    if CI < 0.1
        disp('一致性检查通过');
    else
        disp('error：一致性检查未通过')
        disp("CI = "+ CI);
        W = [];
        lamda_max = [];
    end
    
end

function flag = ConsistencyDetection(lamda , n)
%
% @function：一致性检查
% @illustrate：该函数可单独对某一层进行一致性检查，用于辅助AHP完整过程的编写
% @param lamda: 该层最大特征值
% @param n: 该层分值向量长度
% @return flag: 是否通过标志，通过为1，未通过为0
%
    CI = (lamda - n)/(n - 1);
    if CI < 0.1
        disp('一致性检测通过');
        flag = 1;
    else
        disp('error：一致性检测未通过')
        disp("CI = "+ CI);
        flag = 0;
    end
        
end

