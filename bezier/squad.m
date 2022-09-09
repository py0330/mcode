function [q,dq,darc] = squad(q0,q1,q2,s)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明
% q1 need to be [0,0,0,1]'

q0=q0(:);
q2=q2(:);

% using slerp to compute
q=zeros(4,length(s));
for i=1:length(s)
%     q(:,i)=slerp(slerp(q0,q2,s(i)),q1,2*s(i)*(1-s(i)));
    qa=slerp(slerp(q0,q1,s(i)),q1,s(i));
    qb=slerp(q1,slerp(q1,q2,s(i)),s(i));
    q(:,i) = slerp(qa,qb,s(i));
end

dq=[]
darc=[]
return;
% using exp to compute

[qs,dqs] = slerp(q0,q2,s);
q = zeros(4,length(s));
for i=1:length(s)
    h=s(i);
    q(:,i) = qs(:,i);
    t = acos(q(4,i));
    v = q(1:3,i)/sin(t);
    q(:,i)=[sin((1-2*h*(1-h))*t)*v;cos((1-2*h*(1-h))*t)];
end

thetas = acos(qs(4,:));
dthetas = -dqs(4,:)./sin(thetas);

vs = qs(1:3,:)./sin(thetas);
dvs = (dqs(1:3,:)-cos(thetas).*dthetas.*vs)./sin(thetas);


dvs_norm_square = (sum(dqs.*dqs) - dthetas.*dthetas)./sin(thetas).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = 1-2*s.*(1-s);
dr = -2+4*s;

drtheta = dr.*thetas + r.*dthetas;



dq = [cos(r.*thetas).*drtheta.*vs+sin(r.*thetas).*dvs;
    -sin(r.*thetas).*drtheta];

darc = 2*sqrt(drtheta.^2 + dvs_norm_square .* sin(r.*thetas).^2);

% theta = 
% dq = darc; 


end

