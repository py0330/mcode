function s = bezier3_arc2s(p0,p1,p2,arc)
%   give bezier3_length under condition:
%   
%   p0 p1 p1 p2 是控制点
%
%   p0-p1  和  p1-p2 的距离相同

p0=p0(:);
p1=p1(:);
p2=p2(:);

dis = norm(p1-p0);
theta = acos((p2-p1)'*(p1-p0)/dis/dis);

arc = arc/dis;

a = 12-6*cos(theta/2);
b = -a;
c = 3;

d = (a-6)/pi;
e = (a-6)/pi*0.5;
f = -(a-6)/pi*0.5;

% 全弧长 
max_arc = (1/3*a + 1/2*b - d/pi*cos(pi) + e/2/pi*sin(2*pi) + (f+c) + d/pi);

% initial guess
s = arc / max_arc;

% newton-raphson method
v  =(1/3*a*s.*s.*s + 1/2*b*s.*s - d/pi*cos(s*pi) + e/2/pi*sin(s*2*pi) + (f+c)*s + d/pi);
while(abs(v-arc) > 1e-10)
    dv = a*s.*s + b*s + c + d*sin(s*pi) + e*cos(s*2*pi) + f;
    s  = s - (v-arc)/dv;
    v  =(1/3*a*s.*s.*s + 1/2*b*s.*s - d/pi*cos(s*pi) + e/2/pi*sin(s*2*pi) + (f+c)*s + d/pi);
end


end








