function [p, dp] = bezier3(p0,p1,p2,p3,s)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明
% B1 = p0 * (1-s) + p1 * s;
% B2 = p1 * (1-s) + p2 * s;
% 
% % p = B1 .* (1-s) +  B2 .* s;

%p = (p0 + p2 - 2*p1)*(s.*s) + 2*(p1-p0)*s + p0;
%dp = 2*(p0 + p2 - 2*p1)*s + 2*(p1-p0);

%p = (1-s).*bezier2(p0,p1,p2,s)+ s.*bezier2(p1,p2,p3,s);
%dp = s;

p0=p0(:);
p1=p1(:);
p2=p2(:);
p3=p3(:);

p  = (3*p1 - p0 - 3*p2 + p3)*s.^3 + (3*p0 - 6*p1 + 3*p2)*s.^2 + (3*p1 - 3*p0)*s + p0;
dp = 3*p1 - 3*p0 - 3*(p0 - 3*p1 + 3*p2 - p3)*s.^2 + 2*(3*p0 - 6*p1 + 3*p2)*s;

end