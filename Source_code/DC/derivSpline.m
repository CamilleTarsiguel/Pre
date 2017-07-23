%% symbolic derivatives for B-splines
syms t ax bx cx dx ay by cy dy f_

%f=t*(t*(t*ax + bx)+cx)+dx;
f=t^3*ax + t^2*bx + t*cx + dx;
%  f_=diff(f,t)


%% Elin
syms t ax bx cx dx ay by cy dy f fx_ fy_ speed frR nf hory k

syms c1x c2x c3x c4x c1y c2y c3y c4y t1 t2 t3 t4 t5 t6 curt

% fx_=c1x; fy_=c1y;
fx_=-((((c2x*(curt - t1) - c1x*(curt - t4))*(curt - t4))/(t1 - t4) - ((c3x*(curt - t2) - c2x*(curt - t5))*(curt - t2))/(t2 - t5))/(t2 - t4) - (((c3x*(curt - t2) - c2x*(curt - t5))*(curt - t5))/(t2 - t5) - ((c4x*(curt - t3) - c3x*(curt - t6))*(curt - t3))/(t3 - t6))/(t3 - t5) + ((curt - t4)*((c2x*(curt - t1) - c1x*(curt - t4))/(t1 - t4) - (c3x*(curt - t2) - c2x*(curt - t5))/(t2 - t5) - ((c1x - c2x)*(curt - t4))/(t1 - t4) + ((c2x - c3x)*(curt - t2))/(t2 - t5)))/(t2 - t4) - ((curt - t3)*((c3x*(curt - t2) - c2x*(curt - t5))/(t2 - t5) - (c4x*(curt - t3) - c3x*(curt - t6))/(t3 - t6) - ((c2x - c3x)*(curt - t5))/(t2 - t5) + ((c3x - c4x)*(curt - t3))/(t3 - t6)))/(t3 - t5))/(t3 - t4);
fy_=-((((c2y*(curt - t1) - c1y*(curt - t4))*(curt - t4))/(t1 - t4) - ((c3y*(curt - t2) - c2y*(curt - t5))*(curt - t2))/(t2 - t5))/(t2 - t4) - (((c3y*(curt - t2) - c2y*(curt - t5))*(curt - t5))/(t2 - t5) - ((c4y*(curt - t3) - c3y*(curt - t6))*(curt - t3))/(t3 - t6))/(t3 - t5) + ((curt - t4)*((c2y*(curt - t1) - c1y*(curt - t4))/(t1 - t4) - (c3y*(curt - t2) - c2y*(curt - t5))/(t2 - t5) - ((c1y - c2y)*(curt - t4))/(t1 - t4) + ((c2y - c3y)*(curt - t2))/(t2 - t5)))/(t2 - t4) - ((curt - t3)*((c3y*(curt - t2) - c2y*(curt - t5))/(t2 - t5) - (c4y*(curt - t3) - c3y*(curt - t6))/(t3 - t6) - ((c2y - c3y)*(curt - t5))/(t2 - t5) + ((c3y - c4y)*(curt - t3))/(t3 - t6)))/(t3 - t5))/(t3 - t4);
fy =-(((((c2y*(curt - t1) - c1y*(curt - t4))*(curt - t4))/(t1 - t4) - ((c3y*(curt - t2) - c2y*(curt - t5))*(curt - t2))/(t2 - t5))*(curt - t4))/(t2 - t4) - ((((c3y*(curt - t2) - c2y*(curt - t5))*(curt - t5))/(t2 - t5) - ((c4y*(curt - t3) - c3y*(curt - t6))*(curt - t3))/(t3 - t6))*(curt - t3))/(t3 - t5))/(t3 - t4);

% obj=(speed - (frR*(px_^2 + py_^2)^(1/2))/nf)^2;
obj=(sqrt((fx_)^2 + (fy_)^2)/nf*frR-speed)^2;
obj=( (1/sqrt((fy-hory)^2+k) * sqrt((fx_)^2 + (fy_)^2) )/nf*frR)^2; % ETH

varnames={'c1x','c2x','c3x','c4x','c1y','c2y','c3y','c4y'};

%% Eang bspline
syms t ax bx cx dx ay by cy dy f fx_ fy_ fx__ fy__ frR k epsil hory
syms c1x c2x c3x c4x c1y c2y c3y c4y t1 t2 t3 t4 t5 t6 curt

% fx_ = 3*t*t*ax + 2*t*bx + cx;
% fy_ = 3*t*t*ay + 2*t*by + cy;
% fx__= 6*t*ax + 2*bx;
% fy__= 6*t*ay + 2*by;
fy =-(((((c2y*(curt - t1) - c1y*(curt - t4))*(curt - t4))/(t1 - t4) - ((c3y*(curt - t2) - c2y*(curt - t5))*(curt - t2))/(t2 - t5))*(curt - t4))/(t2 - t4) - ((((c3y*(curt - t2) - c2y*(curt - t5))*(curt - t5))/(t2 - t5) - ((c4y*(curt - t3) - c3y*(curt - t6))*(curt - t3))/(t3 - t6))*(curt - t3))/(t3 - t5))/(t3 - t4);
fx_=-((((c2x*(curt - t1) - c1x*(curt - t4))*(curt - t4))/(t1 - t4) - ((c3x*(curt - t2) - c2x*(curt - t5))*(curt - t2))/(t2 - t5))/(t2 - t4) - (((c3x*(curt - t2) - c2x*(curt - t5))*(curt - t5))/(t2 - t5) - ((c4x*(curt - t3) - c3x*(curt - t6))*(curt - t3))/(t3 - t6))/(t3 - t5) + ((curt - t4)*((c2x*(curt - t1) - c1x*(curt - t4))/(t1 - t4) - (c3x*(curt - t2) - c2x*(curt - t5))/(t2 - t5) - ((c1x - c2x)*(curt - t4))/(t1 - t4) + ((c2x - c3x)*(curt - t2))/(t2 - t5)))/(t2 - t4) - ((curt - t3)*((c3x*(curt - t2) - c2x*(curt - t5))/(t2 - t5) - (c4x*(curt - t3) - c3x*(curt - t6))/(t3 - t6) - ((c2x - c3x)*(curt - t5))/(t2 - t5) + ((c3x - c4x)*(curt - t3))/(t3 - t6)))/(t3 - t5))/(t3 - t4);
fy_=-((((c2y*(curt - t1) - c1y*(curt - t4))*(curt - t4))/(t1 - t4) - ((c3y*(curt - t2) - c2y*(curt - t5))*(curt - t2))/(t2 - t5))/(t2 - t4) - (((c3y*(curt - t2) - c2y*(curt - t5))*(curt - t5))/(t2 - t5) - ((c4y*(curt - t3) - c3y*(curt - t6))*(curt - t3))/(t3 - t6))/(t3 - t5) + ((curt - t4)*((c2y*(curt - t1) - c1y*(curt - t4))/(t1 - t4) - (c3y*(curt - t2) - c2y*(curt - t5))/(t2 - t5) - ((c1y - c2y)*(curt - t4))/(t1 - t4) + ((c2y - c3y)*(curt - t2))/(t2 - t5)))/(t2 - t4) - ((curt - t3)*((c3y*(curt - t2) - c2y*(curt - t5))/(t2 - t5) - (c4y*(curt - t3) - c3y*(curt - t6))/(t3 - t6) - ((c2y - c3y)*(curt - t5))/(t2 - t5) + ((c3y - c4y)*(curt - t3))/(t3 - t6)))/(t3 - t5))/(t3 - t4);
fx__ = -((2*((c2x*(curt - t1) - c1x*(curt - t4))/(t1 - t4) - (c3x*(curt - t2) - c2x*(curt - t5))/(t2 - t5) - ((c1x - c2x)*(curt - t4))/(t1 - t4) + ((c2x - c3x)*(curt - t2))/(t2 - t5)))/(t2 - t4) - (2*((c3x*(curt - t2) - c2x*(curt - t5))/(t2 - t5) - (c4x*(curt - t3) - c3x*(curt - t6))/(t3 - t6) - ((c2x - c3x)*(curt - t5))/(t2 - t5) + ((c3x - c4x)*(curt - t3))/(t3 - t6)))/(t3 - t5) - (((2*(c1x - c2x))/(t1 - t4) - (2*(c2x - c3x))/(t2 - t5))*(curt - t4))/(t2 - t4) + (((2*(c2x - c3x))/(t2 - t5) - (2*(c3x - c4x))/(t3 - t6))*(curt - t3))/(t3 - t5))/(t3 - t4);
fy__ = -((2*((c2y*(curt - t1) - c1y*(curt - t4))/(t1 - t4) - (c3y*(curt - t2) - c2y*(curt - t5))/(t2 - t5) - ((c1y - c2y)*(curt - t4))/(t1 - t4) + ((c2y - c3y)*(curt - t2))/(t2 - t5)))/(t2 - t4) - (2*((c3y*(curt - t2) - c2y*(curt - t5))/(t2 - t5) - (c4y*(curt - t3) - c3y*(curt - t6))/(t3 - t6) - ((c2y - c3y)*(curt - t5))/(t2 - t5) + ((c3y - c4y)*(curt - t3))/(t3 - t6)))/(t3 - t5) - (((2*(c1y - c2y))/(t1 - t4) - (2*(c2y - c3y))/(t2 - t5))*(curt - t4))/(t2 - t4) + (((2*(c2y - c3y))/(t2 - t5) - (2*(c3y - c4y))/(t3 - t6))*(curt - t3))/(t3 - t5))/(t3 - t4);


nominator=(fx_*fy__ - fy_*fx__);

% simple curvature (pseudo huber)
obj=frR * k^2 * (sqrt(1+(nominator/k)^2)-1)/(fx_^2 + fy_^2);

% simple Charbonnier
obj = frR * sqrt(epsil+nominator^2) / (fx_^2 + fy_^2 + epsil);

yscale=(1/sqrt((fy-hory)^2+k));

% lorentzian penalty
% obj = -log(1/(obj^2+1));
obj = -log(1/((yscale*obj)^2+1)); % ETH yscaling
varnames={'c1x','c2x','c3x','c4x','c1y','c2y','c3y','c4y'};

%% Eseg
syms t ax bx cx dx ay by cy dy a2x b2x c2x d2x d2y f_ T
obj = (T^3*ax + T^2*bx + T*cx + dx - d2x)^2 + (T^3*ay + T^2*by + T*cy + dy - d2y)^2;
varnames={'ax','bx','cx','dx','d2x','d2y'};

%% Eseg der
syms t ax bx cx dx ay by cy dy c2x c2y T
obj = (3*T^2*ax + 2*T*bx + cx - c2x)^2 + (3*T^2*ay + 2*T*by + cy - c2y)^2;
varnames={'ax','bx','cx','ay','by','cy','c2x','c2y'};

%% Eseg2
syms c1x c2x c3x c4x c1y c2y c3y c4y t1 t2 t3 t4 t5 t6 curt
syms t1_2 t2_2 t3_2 t4_2 t5_2 t6_2 curt_2

fx_=-((((c2x*(curt - t1) - c1x*(curt - t4))*(curt - t4))/(t1 - t4) - ((c3x*(curt - t2) - c2x*(curt - t5))*(curt - t2))/(t2 - t5))/(t2 - t4) - (((c3x*(curt - t2) - c2x*(curt - t5))*(curt - t5))/(t2 - t5) - ((c4x*(curt - t3) - c3x*(curt - t6))*(curt - t3))/(t3 - t6))/(t3 - t5) + ((curt - t4)*((c2x*(curt - t1) - c1x*(curt - t4))/(t1 - t4) - (c3x*(curt - t2) - c2x*(curt - t5))/(t2 - t5) - ((c1x - c2x)*(curt - t4))/(t1 - t4) + ((c2x - c3x)*(curt - t2))/(t2 - t5)))/(t2 - t4) - ((curt - t3)*((c3x*(curt - t2) - c2x*(curt - t5))/(t2 - t5) - (c4x*(curt - t3) - c3x*(curt - t6))/(t3 - t6) - ((c2x - c3x)*(curt - t5))/(t2 - t5) + ((c3x - c4x)*(curt - t3))/(t3 - t6)))/(t3 - t5))/(t3 - t4);
fy_=-((((c2y*(curt - t1) - c1y*(curt - t4))*(curt - t4))/(t1 - t4) - ((c3y*(curt - t2) - c2y*(curt - t5))*(curt - t2))/(t2 - t5))/(t2 - t4) - (((c3y*(curt - t2) - c2y*(curt - t5))*(curt - t5))/(t2 - t5) - ((c4y*(curt - t3) - c3y*(curt - t6))*(curt - t3))/(t3 - t6))/(t3 - t5) + ((curt - t4)*((c2y*(curt - t1) - c1y*(curt - t4))/(t1 - t4) - (c3y*(curt - t2) - c2y*(curt - t5))/(t2 - t5) - ((c1y - c2y)*(curt - t4))/(t1 - t4) + ((c2y - c3y)*(curt - t2))/(t2 - t5)))/(t2 - t4) - ((curt - t3)*((c3y*(curt - t2) - c2y*(curt - t5))/(t2 - t5) - (c4y*(curt - t3) - c3y*(curt - t6))/(t3 - t6) - ((c2y - c3y)*(curt - t5))/(t2 - t5) + ((c3y - c4y)*(curt - t3))/(t3 - t6)))/(t3 - t5))/(t3 - t4);
fx_2=-((((c2x*(curt_2 - t1) - c1x*(curt_2 - t4))*(curt_2 - t4))/(t1 - t4) - ((c3x*(curt_2 - t2) - c2x*(curt_2 - t5))*(curt_2 - t2))/(t2 - t5))/(t2 - t4) - (((c3x*(curt_2 - t2) - c2x*(curt_2 - t5))*(curt_2 - t5))/(t2 - t5) - ((c4x*(curt_2 - t3) - c3x*(curt_2 - t6))*(curt_2 - t3))/(t3 - t6))/(t3 - t5) + ((curt_2 - t4)*((c2x*(curt_2 - t1) - c1x*(curt_2 - t4))/(t1 - t4) - (c3x*(curt_2 - t2) - c2x*(curt_2 - t5))/(t2 - t5) - ((c1x - c2x)*(curt_2 - t4))/(t1 - t4) + ((c2x - c3x)*(curt_2 - t2))/(t2 - t5)))/(t2 - t4) - ((curt_2 - t3)*((c3x*(curt_2 - t2) - c2x*(curt_2 - t5))/(t2 - t5) - (c4x*(curt_2 - t3) - c3x*(curt_2 - t6))/(t3 - t6) - ((c2x - c3x)*(curt_2 - t5))/(t2 - t5) + ((c3x - c4x)*(curt_2 - t3))/(t3 - t6)))/(t3 - t5))/(t3 - t4);
fy_2=-((((c2y*(curt_2 - t1) - c1y*(curt_2 - t4))*(curt_2 - t4))/(t1 - t4) - ((c3y*(curt_2 - t2) - c2y*(curt_2 - t5))*(curt_2 - t2))/(t2 - t5))/(t2 - t4) - (((c3y*(curt_2 - t2) - c2y*(curt_2 - t5))*(curt_2 - t5))/(t2 - t5) - ((c4y*(curt_2 - t3) - c3y*(curt_2 - t6))*(curt_2 - t3))/(t3 - t6))/(t3 - t5) + ((curt_2 - t4)*((c2y*(curt_2 - t1) - c1y*(curt_2 - t4))/(t1 - t4) - (c3y*(curt_2 - t2) - c2y*(curt_2 - t5))/(t2 - t5) - ((c1y - c2y)*(curt_2 - t4))/(t1 - t4) + ((c2y - c3y)*(curt_2 - t2))/(t2 - t5)))/(t2 - t4) - ((curt_2 - t3)*((c3y*(curt_2 - t2) - c2y*(curt_2 - t5))/(t2 - t5) - (c4y*(curt_2 - t3) - c3y*(curt_2 - t6))/(t3 - t6) - ((c2y - c3y)*(curt_2 - t5))/(t2 - t5) + ((c3y - c4y)*(curt_2 - t3))/(t3 - t6)))/(t3 - t5))/(t3 - t4);

obj=((fx_-fx_2)^2 + (fy_-fy_2)^2);

varnames={'c1x','c2x','c3x','c4x','c1y','c2y','c3y','c4y'};


%% Eexc
syms t1 t2 ax bx cx dx ay by cy dy a2x b2x c2x d2x a2y b2y c2y d2y sigscale siga sigb
syms c1x c2x c3x c4x c1y c2y c3y c4y t1 t2 t3 t4 t5 t6 curt
syms c1x_2 c2x_2 c3x_2 c4x_2 c1y_2 c2y_2 c3y_2 c4y_2 t1_2 t2_2 t3_2 t4_2 t5_2 t6_2


px = t1^3*ax + t1^2*bx + t1*cx + dx; py = t1^3*ay + t1^2*by + t1*cy + dy;
p2x = t2^3*a2x + t2^2*b2x + t2*c2x + d2x; p2y = t2^3*a2y + t2^2*b2y + t2*c2y + d2y;

px=-(((((c2x*(curt - t1) - c1x*(curt - t4))*(curt - t4))/(t1 - t4) - ((c3x*(curt - t2) - c2x*(curt - t5))*(curt - t2))/(t2 - t5))*(curt - t4))/(t2 - t4) - ((((c3x*(curt - t2) - c2x*(curt - t5))*(curt - t5))/(t2 - t5) - ((c4x*(curt - t3) - c3x*(curt - t6))*(curt - t3))/(t3 - t6))*(curt - t3))/(t3 - t5))/(t3 - t4);
py=-(((((c2y*(curt - t1) - c1y*(curt - t4))*(curt - t4))/(t1 - t4) - ((c3y*(curt - t2) - c2y*(curt - t5))*(curt - t2))/(t2 - t5))*(curt - t4))/(t2 - t4) - ((((c3y*(curt - t2) - c2y*(curt - t5))*(curt - t5))/(t2 - t5) - ((c4y*(curt - t3) - c3y*(curt - t6))*(curt - t3))/(t3 - t6))*(curt - t3))/(t3 - t5))/(t3 - t4);
p2x=-(((((c2x_2*(curt - t1_2) - c1x_2*(curt - t4_2))*(curt - t4_2))/(t1_2 - t4_2) - ((c3x_2*(curt - t2_2) - c2x_2*(curt - t5_2))*(curt - t2_2))/(t2_2 - t5_2))*(curt - t4_2))/(t2_2 - t4_2) - ((((c3x_2*(curt - t2_2) - c2x_2*(curt - t5_2))*(curt - t5_2))/(t2_2 - t5_2) - ((c4x_2*(curt - t3_2) - c3x_2*(curt - t6_2))*(curt - t3_2))/(t3_2 - t6_2))*(curt - t3_2))/(t3_2 - t5_2))/(t3_2 - t4_2);
p2y=-(((((c2y_2*(curt - t1_2) - c1y_2*(curt - t4_2))*(curt - t4_2))/(t1_2 - t4_2) - ((c3y_2*(curt - t2_2) - c2y_2*(curt - t5_2))*(curt - t2_2))/(t2_2 - t5_2))*(curt - t4_2))/(t2_2 - t4_2) - ((((c3y_2*(curt - t2_2) - c2y_2*(curt - t5_2))*(curt - t5_2))/(t2_2 - t5_2) - ((c4y_2*(curt - t3_2) - c3y_2*(curt - t6_2))*(curt - t3_2))/(t3_2 - t6_2))*(curt - t3_2))/(t3_2 - t5_2))/(t3_2 - t4_2);


obj = sqrt((px-p2x)^2 + (py-p2y)^2);
obj = 1-1./(1+exp(-siga*obj+sigb));
varnames={'c1x','c2x','c3x','c4x','c1y','c2y','c3y','c4y', ...
    'c1x_2','c2x_2','c3x_2','c4x_2','c1y_2','c2y_2','c3y_2','c4y_2'};

%% Edat
syms t ax bx cx dx ay by cy dy Dx Dy dist deltasq sc nf k ksq

dist=sqrt((t^3*ax + t^2*bx + t*cx + dx - Dx)^2 + (t^3*ay + t^2*by + t*cy + dy - Dy)^2);

% L2^2
obj = sc * (dist/nf)^2;

% Pseudo-huber d^2 * (sqrt(1+(x/d)^2)-1)
obj = sc * k * (sqrt(1 + (dist/nf)^2/ksq) - 1);

% simple Charbonnier
obj = sc * sqrt(k + (dist/nf)^2);
varnames={'ax','bx','cx','dx','ay','by','cy','dy'};

%% Edatsp
syms t curt ax bx cx dx ay by cy dy Dx Dy dist deltasq sc nf k ksq
syms c1x c2x c3x c4x c1y c2y c3y c4y t1 t2 t3 t4 t5 t6

t1=t1-curt;t2=t2-curt;t3=t3-curt;
t4=t4-curt;t5=t5-curt;t6=t6-curt;

c1x=(t4*c1x-t1*c2x)/(t4-t1);
c2x=(t5*c2x-t2*c3x)/(t5-t2);
c3x=(t6*c3x-t3*c4x)/(t6-t3);
c1x=(t4*c1x-t2*c2x)/(t4-t2);
c2x=(t5*c2x-t3*c3x)/(t5-t3);
c1x=(t4*c1x-t3*c2x)/(t4-t3);

c1y=(t4*c1y-t1*c2y)/(t4-t1);
c2y=(t5*c2y-t2*c3y)/(t5-t2);
c3y=(t6*c3y-t3*c4y)/(t6-t3);
c1y=(t4*c1y-t2*c2y)/(t4-t2);
c2y=(t5*c2y-t3*c3y)/(t5-t3);
c1y=(t4*c1y-t3*c2y)/(t4-t3);


% dist=sqrt((t^3*ax + t^2*bx + t*cx + dx - Dx)^2 + (t^3*ay + t^2*by + t*cy + dy - Dy)^2);
dist = sqrt((c1x-Dx)^2 + (c1y-Dy)^2);

% L2^2
obj = sc * (dist/nf)^2;

% Pseudo-huber d^2 * (sqrt(1+(x/d)^2)-1)
obj = sc * k * (sqrt(1 + (dist/nf)^2/ksq) - 1);

% simple Charbonnier
obj = sc * sqrt(k + (dist/nf)^2);
varnames={'c1x','c2x','c3x','c4x','c1y','c2y','c3y','c4y'};


%% Efid
syms t ax bx cx dx ay by cy dy Dx Dy dist deltasq sc nf k ksq siga sigb

dist=sqrt((t^3*ax + t^2*bx + t*cx + dx - Dx)^2 + (t^3*ay + t^2*by + t*cy + dy - Dy)^2);
% Pseudo-huber d^2 * (sqrt(1+(x/d)^2)-1)
dist = k * (sqrt(1 + dist^2/ksq) - 1);

% Simple Charbonnier
dist = sqrt(k + dist^2);


% L2^2
obj = 1./(1+exp(-siga*dist+sigb));

varnames={'ax','bx','cx','dx','ay','by','cy','dy'};

%% Eper start
syms t ax bx cx dx ay by cy dy x0 x1 x2 x3 x4 y0 y1 y2 y3 y4 px py nf k
px = dx; py=dy;

% signed distance
dist = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt((x1-x0)^2 + (y1-y0)^2) );

dl = ((y1-y2)*px + (x2-x1)*py + (x1*y2-x2*y1) ) / (sqrt((x2-x1)^2 + (y2-y1)^2) );
dt = ((y2-y3)*px + (x3-x2)*py + (x2*y3-x3*y2) ) / (sqrt((x3-x2)^2 + (y3-y2)^2) );
dr = ((y3-y4)*px + (x4-x3)*py + (x3*y4-x4*y3) ) / (sqrt((x4-x3)^2 + (y4-y3)^2) );
db = ((y4-y1)*px + (x1-x4)*py + (x4*y1-x1*y4) ) / (sqrt((x1-x4)^2 + (y1-y4)^2) );

% pseudo huber
dl = sqrt(1 + dl^2)-1; dl = dl / nf;
dt = sqrt(1 + dt^2)-1; dt = dt / nf;
dr = sqrt(1 + dr^2)-1; dr = dr / nf;
db = sqrt(1 + db^2)-1; db = db / nf;

dist = dist / nf;

% Tukey
dl=(k^2/6)*(1-(1-(dl/k)^2)^3); dl=k*6/k^2*dl;
dt=(k^2/6)*(1-(1-(dt/k)^2)^3); dt=k*6/k^2*dt;
dr=(k^2/6)*(1-(1-(dr/k)^2)^3); dr=k*6/k^2*dr;
db=(k^2/6)*(1-(1-(db/k)^2)^3); db=k*6/k^2*db;

obj=(k^2/6)*(1-(1-(dist/k)^2)^3); obj=obj*6/k;
varnames={'dx','dy'};

%% Eper end
syms t ax bx cx dx ay by cy dy x0 x1 x2 x3 x4 y0 y1 y2 y3 y4 px py nf k
syms c1x c2x c3x c4x c1y c2y c3y c4y t1 t2 t3 t4 t5 t6 curt

% px = t^3*ax + t^2*bx + t*cx + dx; py = t^3*ay + t^2*by + t*cy + dy;

px=-(((((c2x*(curt - t1) - c1x*(curt - t4))*(curt - t4))/(t1 - t4) - ((c3x*(curt - t2) - c2x*(curt - t5))*(curt - t2))/(t2 - t5))*(curt - t4))/(t2 - t4) - ((((c3x*(curt - t2) - c2x*(curt - t5))*(curt - t5))/(t2 - t5) - ((c4x*(curt - t3) - c3x*(curt - t6))*(curt - t3))/(t3 - t6))*(curt - t3))/(t3 - t5))/(t3 - t4);
py=-(((((c2y*(curt - t1) - c1y*(curt - t4))*(curt - t4))/(t1 - t4) - ((c3y*(curt - t2) - c2y*(curt - t5))*(curt - t2))/(t2 - t5))*(curt - t4))/(t2 - t4) - ((((c3y*(curt - t2) - c2y*(curt - t5))*(curt - t5))/(t2 - t5) - ((c4y*(curt - t3) - c3y*(curt - t6))*(curt - t3))/(t3 - t6))*(curt - t3))/(t3 - t5))/(t3 - t4);
% signed distance
dist = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt((x1-x0)^2 + (y1-y0)^2) );
dist = dist / nf;
obj=(k^2/6)*(1-(1-(dist/k)^2)^3); obj=obj*6/k;

varnames={'c1x','c2x','c3x','c4x','c1y','c2y','c3y','c4y'};



%% Matlab code

objstr=char(obj);
objstr=optimizeString(objstr);
fprintf('obj=%s;\n',objstr);
fprintf('fx=fx+obj;\n');

eper=0;
eexc=0;
%
for allvars=1:length(varnames)
    derivestr=sprintf('dobj = diff(obj,''%s'');',varnames{allvars});
    eval(derivestr);
    if eper
        codestring = sprintf('derl(%i) = %s;',allvars,char(dobj));
    elseif eexc
        if allvars<=8
            codestring = sprintf('dfx(svpos1(%i)) = dfx(svpos1(%i)) + %s;',allvars,allvars,char(dobj));
        else
            codestring = sprintf('dfx(svpos2(%i)) = dfx(svpos2(%i)) + %s;',allvars-8,allvars-8,char(dobj));
        end
    else
%         if allvars<=4
            codestring = sprintf('dfx(svpos(%i)) = dfx(svpos(%i)) + %s;',allvars,allvars,char(dobj));
%         else
%             codestring = sprintf('dfx(%i+sn + os) = dfx(%i+sn + os) + %s;',allvars-4,allvars-4,char(dobj));
%         end
    end
    %     codestring = sprintf('tmpdfx(%i,det) = %s;',allvars,char(dobj));
    codestring = optimizeString(codestring);
    disp(codestring);
end

%% C code

objstr=char(obj);
% fprintf('obj=%s;\n',char(objstr));
fprintf('obj=%s\n',ccode(obj));
objstr=optimizeString(objstr);
fprintf('*fx+=obj;\n');

eper=0;

for allvars=1:length(varnames)
    derivestr=sprintf('dobj = diff(obj,''%s'');',varnames{allvars});
    eval(derivestr);
%     codestring = sprintf('dfx[%i + os] += %s;',allvars-1,char(dobj));
        codestring = sprintf('dfx[svpos[%d]] += %s',allvars-1,ccode(dobj));
    if eper
        codestring = sprintf('dderl[%d] = %s',allvars-1,ccode(dobj));
    end
    codestring = optimizeString(codestring);
    disp(codestring);
end
%%
% dobjax=diff(obj,ax);
% dobjbx=diff(obj,bx);
% dobjcx=diff(obj,cx);
% dobjdx=diff(obj,dx);

% diff(obj,ay)
% diff(obj,by)
% diff(obj,cy)
% diff(obj,dy)
%%
% f=(2*t^2 + t)^(1/2)
% diff(f,t)