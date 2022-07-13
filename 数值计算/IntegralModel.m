% int
%       用法1: 格式： int（fun,x,a,b）
%       功能： 计算定积分
%       用法2: 格式： int(f,x)
%       功能： 计算不定积分
%       注： 使用int函数之前， 先用syms声明x是符号变量
% 
% trapz （利用梯形法）
%       格式： I=trapz(x,y)
%       功能： 求取定积分， 适用于被积函数是离散数据
%           eg： 
%               ac=@(x)sin(x)./x %用@引入函数句柄
%               x1=pi/4:pi/50:pi/2;
%               y1=ac(x1);
%               s1=trapz(x1,y1)
% 
% quad （基于变步长辛普森法）
%       格式： [I,n]=quad(‘fname’,a,b,Tol,trace)
%       其中： ‘fname是被积函数名
%       a,b是积分上下限
%       Tol是精度控制值， 省却时取0.001
%       Trace:控制是否显示展现积分过程， 取0不展现
%       I： 积分值
%       n： 被积函数调用次数
%           eg： 
%               ac=@(x)sin(x)./x
%               s=quad(ac,pi/4,pi/2)
% 
% integral1
%       格式： q = integral(fun,xmin,xmax)
%       用法1： 广义积分
%           创建函数 f(x)=exp(-x2)(lnx)2。
%           **fun = @(x) exp(-x.2).*log(x).2;
%           计算 x=0 至 x=Inf 的积分。
%           q = integral(fun,0,Inf)
%           》 q = 1.9475
%       用法2： 参数化函数
%           创建带有一个参数 c 的函数 f(x)=1/(x3-2x-c)。
%           fun = @(x,c) 1./(x.^3-2x-c);
%           在 c=5 时， 计算从 x=0 至 x=2 的积分。
%           q = integral(@(x)fun(x,5),0,2)
%           》 q = -0.4605
%       用法3： 向量化积分
%           创建向量值函数 f(x)=[sinx,sin2x,sin3x,sin4x,sin5x] ， 并求 x=0 到 x=1 的积分。 指定
%           ‘ArrayValued’,true 以便计算数组值或向量值函数的积分。
%           fun = @(x)sin((1:5)*x);
%           q = integral(fun,0,1,‘ArrayValued’,true) %true表示被积函数是数组值函数
%           》 q = 1×5
%           0.4597 0.7081 0.6633 0.4134 0.1433
% 
% 高精度Lobatto积分法
%       格式： z = quadl(Fun,a,b)
% 
% 自适应Gauss-Kronrod数值积分
%       格式： z = quadgk(Fun,a,b)
%       功能： 适用于高精度和震荡数值积分， 以及广义数值积分
%       积分法矢量化
% 
% 自适应simpson数值积分
%       格式： z = quadv(Fun,a,b)
%       功能： 向量化积分
%           eg：
%               F=@(x,n)1./((1:n)+x.^2);
%               quadv(@(x)F(x,6),0,1)
% 
% dblquad 数值二重积分
%       格式： I=dblquad(f,a,b,c,d,tol,trace)，
%       功能： 求f(x,y)在[a,b]×[c,d]区域上的二重积分
%           eg.
%               f=@(x,y)exp(-x.2/2).*sin(x.2+y)
%               I=dblquad(f,-2,2,-1,1)
% 
% integral2 数值二重积分
%       格式： q = integral2(fun,xmin,xmax,ymin,ymax)
%       功能： 在平面区域 xmin ≤ x ≤ xmax 和 ymin(x) ≤ y ≤ ymax(x) 上逼近函数 z = fun(x,y) 的积分。
%       用法1： 将三角形区域与奇异性在边界处集成
%           eg：
%               此函数在 x 和 y 为零时未定义。
%               当奇异性位于积分边界上时， integral2 的性能最佳。
%               创建匿名函数。
%               fun = @(x,y) 1./( sqrt(x + y) .* (1 + x* + y).^2 )
%               对 0≤x≤1 和 0≤y≤1-x 限定的三角形区域计算积分。
%               ymax = @(x) 1 - x;
%               q = integral2(fun,0,1,0,ymax)
%               》q = 0.2854
%       用法2: 极坐标二重积分
%           eg：
%               fun = @(x,y) 1./( sqrt(x + y) . (1 + x + y).^2 ); %用x y 表示rcos rsin
%               polarfun = @(theta,r) fun(r.*cos(theta),r.*sin(theta)).*r;
%               为 r 的上限定义一个函数。
%               rmax = @(theta) 1./(sin(theta) + cos(theta));
%               对 0≤θ≤π/2 和 0≤r≤rmax 限定的区域计算积分。
%               q = integral2(polarfun,0,pi/2,0,rmax)
% 
% ntegral3 数值三重积分
%       格式： q = integral3(fun,xmin,xmax,ymin,ymax,zmin,zmax)
%       用法1： 具有有限范围的三重积分
%           eg.
%               fun = @(x,y,z) y.*sin(x)+z.*cos(x)
%               q = integral3(fun,0,pi,0,1,-1,1)
%               》 q = 2.0000
%       用法2： 在笛卡尔坐标中对单位球面计算积分
%           eg.
%               fun = @(x,y,z) x.*cos(y) + x.^2.*cos(z)
%               %积分范围
%               xmin = -1;
%               xmax = 1;
%               ymin = @(x)-sqrt(1 - x.^2);
%               ymax = @(x) sqrt(1 - x.^2);
%               zmin = @(x,y)-sqrt(1 - x.^2 - y.^2);
%               zmax = @(x,y) sqrt(1 - x.^2 - y.^2);
%               使用 ‘tiled’ 方法计算定积分。
%               q = integral3(fun,xmin,xmax,ymin,ymax,zmin,zmax,‘Method’,‘tiled’)
