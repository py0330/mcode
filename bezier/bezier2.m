function [p, dp] = bezier2(p0,p1,p2,s)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明
% B1 = p0 * (1-s) + p1 * s;
% B2 = p1 * (1-s) + p2 * s;
% 
% % p = B1 .* (1-s) +  B2 .* s;
p0=p0(:);
p1=p1(:);
p2=p2(:);

p = (p0 + p2 - 2*p1)*(s.*s) + 2*(p1-p0)*s + p0;
dp = 2*(p0 + p2 - 2*p1)*s + 2*(p1-p0);



end