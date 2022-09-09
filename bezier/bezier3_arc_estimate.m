function len = bezier3_arc_estimate(p0,p1,p2,s)
%   give bezier3_length under condition:
%   
%   p0 p1 p1 p2 是控制点
%
%   p0-p1  和  p1-p2 的距离相同

p0=p0(:);
p1=p1(:);
p2=p2(:);

dis = norm(p1-p0);
theta = asin((p2-p1)'*(p1-p0)/dis/dis);


x2 = sin(theta) + 1;
y2 = cos(theta);

p0 = [0,0,0]';
p1 = [1,0,0]';
p2 = [x2,y2,0]';

% 数值求解
ds=0.0001;
s=0:ds:s;
[p,dp] = bezier3(p0,p1,p1,p2,s);
len = sum(sqrt(sum(diff(p')'.*diff(p')')))*dis
% len = sum(sqrt(sum(dp.*dp)))*ds*dis



end








