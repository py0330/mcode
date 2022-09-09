function q = s_q_dot_q(q1,q2)
%UNTITLED6 此处提供此函数的摘要
%   此处提供详细说明
q1 = q1(:);
q2 = q2(:);

v1 = q1(1:3);
v2 = q2(1:3);
s1 = q1(4);
s2 = q2(4);

q = [cross(v1,v2)+s1*v2+s2*v1;s1*s2-v1'*v2];


end