function [arc, darc, ddarc] = s_estimate_bezier3_arc(darc0, d2arc0, darc1, d2arc1, darc50, s)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明

% 对于 theta 较大时，以下式子估算较为准确
% 使用1元3次方程 a s^3 + b s^2 + c s + d  模拟在s=0和s=1处的dp
% [ 0 0 0 1 ] [ a ] = [ darc0  ] 
% | 1 1 1 1 | | b |   | darc1  |
% | 0 0 1 0 | | c |   | d2arc0 |
% [ 3 2 1 0 ] [ d ]   [ d2arc1 ]
% 
% 系数矩阵的逆为：
% [  2 -2  1  1 ]
% | -3  3 -2 -1 |
% |  0  0  1  0 | 
% [  1  0  0  0 ]
%
% 因而：
% a =  2*darc0 - 2*darc1 +   d2arc0 + d2arc1;
% b = -3*darc0 + 3*darc1 - 2*d2arc0 - d2arc1;
% c =                        d2arc0         ;
% d =    darc0                              ;
%
% 后续需要修正d，因此先计算前三项。

a =  2*darc0 - 2*darc1 +   d2arc0 + d2arc1;
b = -3*darc0 + 3*darc1 - 2*d2arc0 - d2arc1;
c =                        d2arc0         ;
d =    darc0                              ;

ddarcA = a.*s.^2.*3 + b.*s   .*2 + c                ;
darcA  = a.*s.^3    + b.*s.^2    + c.*s       + d   ;
arcA   = a.*s.^4./4 + b.*s.^3./3 + c.*s.^2./2 + d.*s;

% 对于theta较小时，以下式子估算较为准确
% 用经验公式 atan2(x, h) *i 来模拟 d2arc, 其中 x = s - 0.5
% 对其积分，可得：
% darc = (x*atan2(x, h) - h/2*log(x*x + h*h))*i + j
%      = x*d2darc - (h/2*log(x*x + h*h))*i + j;
% 【参考】https://zh.m.wikipedia.org/zh-hant/%E5%8F%8D%E4%B8%89%E8%A7%92%E5%87%BD%E6%95%B0%E7%A7%AF%E5%88%86%E8%A1%A8
%
% 进一步积分，可得：
% arc = i*( (x*x+h*h)/2*atan2(x, h)-h*x/2 - (x*log(x*x + h*h)+2*h*atan2(x,h)-2*x)*h/2 ) + k
%     = k*( (x*x-h*h)*i*atan2(x, h))/2 + (h*i*x)/2 - (h*i*x*log(h^2 + x^2))/2
% 【参考】https://www.symbolab.com/solver/step-by-step/%5Cint%20ln%5Cleft(x%5E%7B2%7D%2Bc%5E%7B2%7D%5Cright)%20dx?or=input

h = ((darc50/(max(darc0 + darc1, eps)))*4*pi)^1.5 * 0.05;   % 经验公式
i = (abs(d2arc1) + abs(d2arc0)) /2/ atan2(1.0, 2*h); % 经验公式
j = -(0.5*atan2(0.5, h) - h/2*log(0.25 + h*h))*i + darc0; % darc的常数项
k = -i*(h/4 - (h*(1 + 2*h*atan2(-1/2, h) - log(h^2 + 1/4)/2))/2 + atan2(-1/2, h)*(h^2/2 + 1/8)) + 0.5*j; % arc 的常数项

x = s - 0.5;
lx = log(x.*x + h*h);%可能为nan
at = atan2(x, h);

ddarcB = i*at;
darcB  = i*x.*at - h*i/2*lx + j;
arcB   = k - i*at*h^2/2 + i*x.*x.*at/2 - h*i/2*x.*lx + x*i*h/2+ x*j;

% 上述两方程，无法保证s = 0.5时，darc的正确性，因此需修正
% 使用 cos 函数，在不改变0，1处的darc的前提下，将s=0.5处的darc修正到正确值
% 修正函数为 (1-cos(2 pi s))/2 * darc_error_at_0.5
% 而 darc_error_at_0.5 = dp50 - 0.125.*a - 0.25.*b - 0.5.*c - darc0;
% 
% 将上述等式化简，可得修正之后的一次项d，以及cos函数的系数 e：

ratio = 4*darc50/(darc0+darc1);
ratio = ratio*ratio*ratio*ratio;
ratio = max(ratio, 0);
ratio = min(ratio, 1);

lh = max(log(h), -100);

darcA50 = 0.125.*a + 0.25.*b + 0.5.*c + d;
darcB50 = -lh*h*i + j;
darcE50 = darc50 - (ratio * darcA50   + (1-ratio)*darcB50);

arc   = ratio*arcA    + (1-ratio)*arcB   + s*darcE50/2 - sin(2*pi*s)/2*darcE50/2/pi;
darc  = ratio*darcA   + (1-ratio)*darcB  + (1-cos(2*pi*s))/2 * darcE50;
ddarc = ratio*ddarcA  + (1-ratio)*ddarcB + sin(2*pi*s)*pi*darcE50;



A = ((a*ratio)/4);
B = ((b*ratio)/3);
C = (c*ratio)/2;
D = darc50/2 - (ratio - 1)*(j + (h*i)/2) + d*ratio - (ratio*(a/8 + b/4 + c/2 + d))/2 + ((j - h*i*lh)*(ratio - 1))/2;
E = h*i/2*(ratio - 1);
F = h*i/2*(1-ratio)/2;
G = i*(1-ratio)/2;
H = i*(h^2/2 - 1/8)*(ratio - 1);
I = -(darc50 - ratio*(a/8 + b/4 + c/2 + d) + (j - h*i*lh)*(ratio - 1))/(4*pi);

X = -(F*log(0.25+h*h) + H * atan2(-0.5,h)); % 将 s = 0 带入，arc应该为0
Y = D + (ratio - 1)*((h*i)/2);
Z = c*ratio;

lx = max(log(x.*x + h*h), -100); % log(eps/1e28) ~= -100

% lx = log(x.*x + h*h);%可能为nan
at = atan2(x, h);

A
B
C
D
E
F
G
H
I
X
Y
Z


arc  = A*s.^4 + B*s.^3 + C*s.^2 + D*s + E*lx.*s + F*lx + G*at.*s.^2 - G*at.*s + H*at + I.*sin(2*pi*s) + X;
darc = 4*A*s.^3 + 3*B*s.^2 + 2*C*s + 2*G*at.*s - G*at + E*lx + I*2*pi*cos(2*pi*s) + Y;
ddarc = 12*A*s.^2 + 6*B*s + 2*G*at + 4*I*sin(2*pi*s) + Z;
end