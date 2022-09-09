function [q,dq] = slerp(q0,q1,s)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明

q0=q0(:);
q1=q1(:);

theta = acos(q0'*q1);

q  = zeros(4,length(s));
dq = zeros(4,length(s));
for i=1:length(s)
    if(theta > 1e-10)
        q(:,i) =(sin((1-s(i))*theta)*q0 + sin(s(i)*theta)*q1) / sin(theta);
        dq(:,i)=(-theta * cos((1-s(i))*theta)*q0 + theta * cos(s(i)*theta)*q1) / sin(theta);
    else
        q(:,i) =((1-s(i))*q0 + s(i)*q1);
        dq(:,i)=(-cos((1-s(i))*theta)*q0 + cos(s(i)*theta)*q1);
    end
end

end