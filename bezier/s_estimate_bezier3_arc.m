function [arc, darc, ddarc] = s_estimate_bezier3_arc(darc0, d2arc0, darc1, d2arc1, darc50, s)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明

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

% 上述一元三次方程，无法保证s = 0.5时，darc的正确性，因此需修正
% 使用 cos 函数，在不改变0，1处的darc的前提下，将s=0.5处的darc修正到正确值
% 修正函数为 (1-cos(2 pi s))/2 * darc_error_at_0.5
% 而 darc_error_at_0.5 = dp50 - 0.125.*a - 0.25.*b - 0.5.*c - darc0;
% 
% 将上述等式化简，可得修正之后的一次项d，以及cos函数的系数 e：
d =  0.5*darc50 - 0.0625*a - 0.125*b - 0.25*c + 0.5*darc0;
e = -0.5*darc50 + 0.0625*a + 0.125*b + 0.25*c + 0.5*darc0;

darc  = a.*s.^3    + b.*s.^2    + c.*s       + d    + e*cos(2*pi*s);
arc   = a.*s.^4./4 + b.*s.^3./3 + c.*s.^2./2 + d.*s + e*sin(2*pi*s)/2/pi;
ddarc = a.*s.^2.*3 + b.*s   .*2 + c                 - e*sin(2*pi*s)*2*pi;
end