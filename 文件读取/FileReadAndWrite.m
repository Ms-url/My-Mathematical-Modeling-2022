%%%%%%%%%%%%%%%%%% excel�ļ���д
%%% �ļ���ȡ
%   xlsread()
    [num,txt,raw] = xlsread(fliename,sheet);

%%% �ļ�д��
%   xlswrite()
    xlswrite(fliename);
    
%%% csv
    csvread();
    
%%%%%%%%%%%%%%%%%% ���±��ļ���д
%   load() �����ַ���ʱһ�㲻��ʹ�ô˺���
%   textread() �ú���Ҳ�ɶ�ȡ .dat  .m  .csv ���ļ�
% 
%%% ָ���д
%   fopen()
%%% ɨ�躯�� 
%   fscanf()

    fid = fopen(filename); % ��ȡ�ļ�ָ��
    a = scanf(fid,'%s/n'); % ɨ��һ���ַ���ֵ��a
    fclose(fid); % �ر��ļ�

%%% �ļ�����ӡ���� 
%   fprintf()

    fid = fopen(filename,'wt'); % 'wt'��ʾ����д������
    fprintf(fid,"Strings"); % д��
    fclose(fid); % �ر��ļ�
