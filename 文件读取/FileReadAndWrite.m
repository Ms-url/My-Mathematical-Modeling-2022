%%%%%%%%%%%%%%%%%% excel文件读写
%%% 文件读取
%   xlsread()
    [num,txt,raw] = xlsread(fliename,sheet);

%%% 文件写入
%   xlswrite()
    xlswrite(fliename);
    
%%% csv
    csvread();
    
%%%%%%%%%%%%%%%%%% 记事本文件读写
%   load() 包含字符串时一般不能使用此函数
%   textread() 该函数也可读取 .dat  .m  .csv 的文件
% 
%%% 指针读写
%   fopen()
%%% 扫描函数 
%   fscanf()

    fid = fopen(filename); % 获取文件指针
    a = scanf(fid,'%s/n'); % 扫描一行字符赋值给a
    fclose(fid); % 关闭文件

%%% 文件流打印函数 
%   fprintf()

    fid = fopen(filename,'wt'); % 'wt'表示数据写入命令
    fprintf(fid,"Strings"); % 写入
    fclose(fid); % 关闭文件
