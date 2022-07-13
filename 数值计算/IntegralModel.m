% int
%       �÷�1: ��ʽ�� int��fun,x,a,b��
%       ���ܣ� ���㶨����
%       �÷�2: ��ʽ�� int(f,x)
%       ���ܣ� ���㲻������
%       ע�� ʹ��int����֮ǰ�� ����syms����x�Ƿ��ű���
% 
% trapz ���������η���
%       ��ʽ�� I=trapz(x,y)
%       ���ܣ� ��ȡ�����֣� �����ڱ�����������ɢ����
%           eg�� 
%               ac=@(x)sin(x)./x %��@���뺯�����
%               x1=pi/4:pi/50:pi/2;
%               y1=ac(x1);
%               s1=trapz(x1,y1)
% 
% quad �����ڱ䲽������ɭ����
%       ��ʽ�� [I,n]=quad(��fname��,a,b,Tol,trace)
%       ���У� ��fname�Ǳ���������
%       a,b�ǻ���������
%       Tol�Ǿ��ȿ���ֵ�� ʡȴʱȡ0.001
%       Trace:�����Ƿ���ʾչ�ֻ��ֹ��̣� ȡ0��չ��
%       I�� ����ֵ
%       n�� �����������ô���
%           eg�� 
%               ac=@(x)sin(x)./x
%               s=quad(ac,pi/4,pi/2)
% 
% integral1
%       ��ʽ�� q = integral(fun,xmin,xmax)
%       �÷�1�� �������
%           �������� f(x)=exp(-x2)(lnx)2��
%           **fun = @(x) exp(-x.2).*log(x).2;
%           ���� x=0 �� x=Inf �Ļ��֡�
%           q = integral(fun,0,Inf)
%           �� q = 1.9475
%       �÷�2�� ����������
%           ��������һ������ c �ĺ��� f(x)=1/(x3-2x-c)��
%           fun = @(x,c) 1./(x.^3-2x-c);
%           �� c=5 ʱ�� ����� x=0 �� x=2 �Ļ��֡�
%           q = integral(@(x)fun(x,5),0,2)
%           �� q = -0.4605
%       �÷�3�� ����������
%           ��������ֵ���� f(x)=[sinx,sin2x,sin3x,sin4x,sin5x] �� ���� x=0 �� x=1 �Ļ��֡� ָ��
%           ��ArrayValued��,true �Ա��������ֵ������ֵ�����Ļ��֡�
%           fun = @(x)sin((1:5)*x);
%           q = integral(fun,0,1,��ArrayValued��,true) %true��ʾ��������������ֵ����
%           �� q = 1��5
%           0.4597 0.7081 0.6633 0.4134 0.1433
% 
% �߾���Lobatto���ַ�
%       ��ʽ�� z = quadl(Fun,a,b)
% 
% ����ӦGauss-Kronrod��ֵ����
%       ��ʽ�� z = quadgk(Fun,a,b)
%       ���ܣ� �����ڸ߾��Ⱥ�����ֵ���֣� �Լ�������ֵ����
%       ���ַ�ʸ����
% 
% ����Ӧsimpson��ֵ����
%       ��ʽ�� z = quadv(Fun,a,b)
%       ���ܣ� ����������
%           eg��
%               F=@(x,n)1./((1:n)+x.^2);
%               quadv(@(x)F(x,6),0,1)
% 
% dblquad ��ֵ���ػ���
%       ��ʽ�� I=dblquad(f,a,b,c,d,tol,trace)��
%       ���ܣ� ��f(x,y)��[a,b]��[c,d]�����ϵĶ��ػ���
%           eg.
%               f=@(x,y)exp(-x.2/2).*sin(x.2+y)
%               I=dblquad(f,-2,2,-1,1)
% 
% integral2 ��ֵ���ػ���
%       ��ʽ�� q = integral2(fun,xmin,xmax,ymin,ymax)
%       ���ܣ� ��ƽ������ xmin �� x �� xmax �� ymin(x) �� y �� ymax(x) �ϱƽ����� z = fun(x,y) �Ļ��֡�
%       �÷�1�� ���������������������ڱ߽紦����
%           eg��
%               �˺����� x �� y Ϊ��ʱδ���塣
%               ��������λ�ڻ��ֱ߽���ʱ�� integral2 ��������ѡ�
%               ��������������
%               fun = @(x,y) 1./( sqrt(x + y) .* (1 + x* + y).^2 )
%               �� 0��x��1 �� 0��y��1-x �޶������������������֡�
%               ymax = @(x) 1 - x;
%               q = integral2(fun,0,1,0,ymax)
%               ��q = 0.2854
%       �÷�2: ��������ػ���
%           eg��
%               fun = @(x,y) 1./( sqrt(x + y) . (1 + x + y).^2 ); %��x y ��ʾrcos rsin
%               polarfun = @(theta,r) fun(r.*cos(theta),r.*sin(theta)).*r;
%               Ϊ r �����޶���һ��������
%               rmax = @(theta) 1./(sin(theta) + cos(theta));
%               �� 0�ܦȡܦ�/2 �� 0��r��rmax �޶������������֡�
%               q = integral2(polarfun,0,pi/2,0,rmax)
% 
% ntegral3 ��ֵ���ػ���
%       ��ʽ�� q = integral3(fun,xmin,xmax,ymin,ymax,zmin,zmax)
%       �÷�1�� �������޷�Χ�����ػ���
%           eg.
%               fun = @(x,y,z) y.*sin(x)+z.*cos(x)
%               q = integral3(fun,0,pi,0,1,-1,1)
%               �� q = 2.0000
%       �÷�2�� �ڵѿ��������жԵ�λ����������
%           eg.
%               fun = @(x,y,z) x.*cos(y) + x.^2.*cos(z)
%               %���ַ�Χ
%               xmin = -1;
%               xmax = 1;
%               ymin = @(x)-sqrt(1 - x.^2);
%               ymax = @(x) sqrt(1 - x.^2);
%               zmin = @(x,y)-sqrt(1 - x.^2 - y.^2);
%               zmax = @(x,y) sqrt(1 - x.^2 - y.^2);
%               ʹ�� ��tiled�� �������㶨���֡�
%               q = integral3(fun,xmin,xmax,ymin,ymax,zmin,zmax,��Method��,��tiled��)
