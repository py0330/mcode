function [q,dq,ddq] = s_blend_quaternion_bezier3(q0, q1, q2, s)
% 使用3阶bezier曲线，控制点为：[p0,p1,p1,p2]
q0 = q0(:);
q1 = q1(:);
q2 = q2(:);

% 将 q1 用单位四元数代替
q0 = s_q_dot_q(s_q_inv(q1), q0);
q2 = s_q_dot_q(s_q_inv(q1), q2);

% 求取角度与转轴
theta0 = acos(q0(4));
v0 = q0(1:3)/sin(theta0);

theta2 = acos(q2(4));
v2 = q2(1:3)/sin(theta2);

% 求 qa
qa  = [sin(theta0*((1-s).^2)).*v0;cos(theta0*((1-s).^2))];
dqa = theta0*(2*s-2).*[cos(theta0*((1-s).^2)).*v0;-sin(theta0*((1-s).^2))];
ddqa = theta0*2.*[cos(theta0*((1-s).^2)).*v0;-sin(theta0*((1-s).^2))]...
    +theta0*(2*s-2).*theta0.*(2*s-2).*[-sin(theta0*((1-s).^2)).*v0;-cos(theta0*((1-s).^2))];

% qb1
qb1  = [sin(theta0*(-(1-s).^2)).*v0;cos(theta0*(-(1-s).^2))];
dqb1 = -theta0*(2*s-2).*[cos(theta0*((1-s).^2)).*v0;sin(theta0*((1-s).^2))];
ddqb1 = -theta0*2.*[cos(theta0*((1-s).^2)).*v0;sin(theta0*((1-s).^2))]...
    -theta0*(2*s-2).*theta0.*(2*s-2).*[-sin(theta0*((1-s).^2)).*v0;cos(theta0*((1-s).^2))];

% qb2
qb2  = [sin(theta2*(s.^2)).*v2;cos(theta2*(s.^2))];
dqb2 = 2*theta2*s.*[cos(theta2*(s.^2)).*v2;-sin(theta2*(s.^2))];
ddqb2 = 2*theta2.*[cos(theta2*(s.^2)).*v2;-sin(theta2*(s.^2))]...
    + 2*theta2*s.*2*theta2.*s.*[-sin(theta2*(s.^2)).*v2;-cos(theta2*(s.^2))];

% qb3
qb3 = [ qb1(2,:).*qb2(3,:) - qb2(2,:).*qb1(3,:) + qb1(4,:).*qb2(1,:) + qb2(4,:).*qb1(1,:)
        qb2(1,:).*qb1(3,:) - qb1(1,:).*qb2(3,:) + qb1(4,:).*qb2(2,:) + qb2(4,:).*qb1(2,:)
        qb1(1,:).*qb2(2,:) - qb2(1,:).*qb1(2,:) + qb1(4,:).*qb2(3,:) + qb2(4,:).*qb1(3,:)
        qb1(4,:).*qb2(4,:) - qb1(1,:).*qb2(1,:) - qb1(2,:).*qb2(2,:) - qb1(3,:).*qb2(3,:)
        ];
dqb3 = [
    qb1(2,:).*dqb2(3,:)+dqb1(2,:).*qb2(3,:) - qb2(2,:).*dqb1(3,:)-dqb2(2,:).*qb1(3,:) + qb1(4,:).*dqb2(1,:)+dqb1(4,:).*qb2(1,:) + qb2(4,:).*dqb1(1,:)+dqb2(4,:).*qb1(1,:)
    qb2(1,:).*dqb1(3,:)+dqb2(1,:).*qb1(3,:) - qb1(1,:).*dqb2(3,:)-dqb1(1,:).*qb2(3,:) + qb1(4,:).*dqb2(2,:)+dqb1(4,:).*qb2(2,:) + qb2(4,:).*dqb1(2,:)+dqb2(4,:).*qb1(2,:)
    qb1(1,:).*dqb2(2,:)+dqb1(1,:).*qb2(2,:) - qb2(1,:).*dqb1(2,:)-dqb2(1,:).*qb1(2,:) + qb1(4,:).*dqb2(3,:)+dqb1(4,:).*qb2(3,:) + qb2(4,:).*dqb1(3,:)+dqb2(4,:).*qb1(3,:)
    qb1(4,:).*dqb2(4,:)+dqb1(4,:).*qb2(4,:) - qb1(1,:).*dqb2(1,:)-dqb1(1,:).*qb2(1,:) - qb1(2,:).*dqb2(2,:)-dqb1(2,:).*qb2(2,:) - qb1(3,:).*dqb2(3,:)-dqb1(3,:).*qb2(3,:)
];

ddqb3 = [
    qb1(2,:).*ddqb2(3,:)+2*dqb1(2,:).*dqb2(3,:)+ddqb1(2,:).*qb2(3,:) - qb2(2,:).*ddqb1(3,:)-2*dqb2(2,:).*dqb1(3,:)-ddqb2(2,:).*qb1(3,:) + qb1(4,:).*ddqb2(1,:)+2*dqb1(4,:).*dqb2(1,:)+ddqb1(4,:).*qb2(1,:) + qb2(4,:).*ddqb1(1,:)+2*dqb2(4,:).*dqb1(1,:)+ddqb2(4,:).*qb1(1,:)
    qb2(1,:).*ddqb1(3,:)+2*dqb2(1,:).*dqb1(3,:)+ddqb2(1,:).*qb1(3,:) - qb1(1,:).*ddqb2(3,:)-2*dqb1(1,:).*dqb2(3,:)-ddqb1(1,:).*qb2(3,:) + qb1(4,:).*ddqb2(2,:)+2*dqb1(4,:).*dqb2(2,:)+ddqb1(4,:).*qb2(2,:) + qb2(4,:).*ddqb1(2,:)+2*dqb2(4,:).*dqb1(2,:)+ddqb2(4,:).*qb1(2,:)
    qb1(1,:).*ddqb2(2,:)+2*dqb1(1,:).*dqb2(2,:)+ddqb1(1,:).*qb2(2,:) - qb2(1,:).*ddqb1(2,:)-2*dqb2(1,:).*dqb1(2,:)-ddqb2(1,:).*qb1(2,:) + qb1(4,:).*ddqb2(3,:)+2*dqb1(4,:).*dqb2(3,:)+ddqb1(4,:).*qb2(3,:) + qb2(4,:).*ddqb1(3,:)+2*dqb2(4,:).*dqb1(3,:)+ddqb2(4,:).*qb1(3,:)
    qb1(4,:).*ddqb2(4,:)+2*dqb1(4,:).*dqb2(4,:)+ddqb1(4,:).*qb2(4,:) - qb1(1,:).*ddqb2(1,:)-2*dqb1(1,:).*dqb2(1,:)-ddqb1(1,:).*qb2(1,:) - qb1(2,:).*ddqb2(2,:)-2*dqb1(2,:).*dqb2(2,:)-ddqb1(2,:).*qb2(2,:) - qb1(3,:).*ddqb2(3,:)-2*dqb1(3,:).*dqb2(3,:)-ddqb1(3,:).*qb2(3,:)
];

% qb
theta_b3 = zeros(1,size(qb3,2));
vb3 = zeros(3,size(qb3,2));
dtheta_b3 = zeros(1,size(qb3,2));
dvb3 = zeros(3,size(qb3,2));
ddtheta_b3 = zeros(1,size(qb3,2));
ddvb3 = zeros(3,size(qb3,2));

for i=1:size(qb3,2)
[theta_b3(:,i),vb3(:,i),dtheta_b3(:,i),dvb3(:,i),ddtheta_b3(:,i),ddvb3(:,i)] = ...
    s_q_to_theta_v(qb3(:,i),dqb3(:,i),ddqb3(:,i));
end
qb  = [sin(theta_b3.*s).*vb3;cos(theta_b3.*s)];
dqb = [
    cos(s.*theta_b3).*(theta_b3+s.*dtheta_b3).*vb3+sin(s.*theta_b3).*dvb3
    -sin(s.*theta_b3).*(theta_b3+s.*dtheta_b3)
];

ddqb_1 = -sin(s.*theta_b3).*(theta_b3 + s.*dtheta_b3).*(theta_b3+s.*dtheta_b3).*vb3 ...
        + cos(s.*theta_b3).*(dtheta_b3+ dtheta_b3 + s.*ddtheta_b3).*vb3 ...
        + cos(s.*theta_b3).*(theta_b3+s.*dtheta_b3).*dvb3 ...
        ...
        + cos(s.*theta_b3).*(theta_b3 + s.*dtheta_b3).*dvb3 ...
        + sin(s.*theta_b3).*ddvb3;

ddqb_2 = -cos(s.*theta_b3).*(theta_b3 + s.*dtheta_b3).*(theta_b3+s.*dtheta_b3) ...
        -sin(s.*theta_b3).*(dtheta_b3 + dtheta_b3 + s.*ddtheta_b3);

ddqb = [
    ddqb_1
    ddqb_2
    ];

% 最后计算 q
q = [ 
    qa(2,:).*qb(3,:) - qb(2,:).*qa(3,:) + qa(4,:).*qb(1,:) + qb(4,:).*qa(1,:)
    qb(1,:).*qa(3,:) - qa(1,:).*qb(3,:) + qa(4,:).*qb(2,:) + qb(4,:).*qa(2,:)
    qa(1,:).*qb(2,:) - qb(1,:).*qa(2,:) + qa(4,:).*qb(3,:) + qb(4,:).*qa(3,:)
    qa(4,:).*qb(4,:) - qa(1,:).*qb(1,:) - qa(2,:).*qb(2,:) - qa(3,:).*qb(3,:)
];
dq = [
    qa(2,:).*dqb(3,:)+dqa(2,:).*qb(3,:) - qb(2,:).*dqa(3,:)-dqb(2,:).*qa(3,:) + qa(4,:).*dqb(1,:)+dqa(4,:).*qb(1,:) + qb(4,:).*dqa(1,:)+dqb(4,:).*qa(1,:)
    qb(1,:).*dqa(3,:)+dqb(1,:).*qa(3,:) - qa(1,:).*dqb(3,:)-dqa(1,:).*qb(3,:) + qa(4,:).*dqb(2,:)+dqa(4,:).*qb(2,:) + qb(4,:).*dqa(2,:)+dqb(4,:).*qa(2,:)
    qa(1,:).*dqb(2,:)+dqa(1,:).*qb(2,:) - qb(1,:).*dqa(2,:)-dqb(1,:).*qa(2,:) + qa(4,:).*dqb(3,:)+dqa(4,:).*qb(3,:) + qb(4,:).*dqa(3,:)+dqb(4,:).*qa(3,:)
    qa(4,:).*dqb(4,:)+dqa(4,:).*qb(4,:) - qa(1,:).*dqb(1,:)-dqa(1,:).*qb(1,:) - qa(2,:).*dqb(2,:)-dqa(2,:).*qb(2,:) - qa(3,:).*dqb(3,:)-dqa(3,:).*qb(3,:)
];
ddq = [
    qa(2,:).*ddqb(3,:)+2*dqa(2,:).*dqb(3,:)+ddqa(2,:).*qb(3,:) - qb(2,:).*ddqa(3,:)-2*dqb(2,:).*dqa(3,:)-ddqb(2,:).*qa(3,:) + qa(4,:).*ddqb(1,:)+2*dqa(4,:).*dqb(1,:)+ddqa(4,:).*qb(1,:) + qb(4,:).*ddqa(1,:)+2*dqb(4,:).*dqa(1,:)+ddqb(4,:).*qa(1,:)
    qb(1,:).*ddqa(3,:)+2*dqb(1,:).*dqa(3,:)+ddqb(1,:).*qa(3,:) - qa(1,:).*ddqb(3,:)-2*dqa(1,:).*dqb(3,:)-ddqa(1,:).*qb(3,:) + qa(4,:).*ddqb(2,:)+2*dqa(4,:).*dqb(2,:)+ddqa(4,:).*qb(2,:) + qb(4,:).*ddqa(2,:)+2*dqb(4,:).*dqa(2,:)+ddqb(4,:).*qa(2,:)
    qa(1,:).*ddqb(2,:)+2*dqa(1,:).*dqb(2,:)+ddqa(1,:).*qb(2,:) - qb(1,:).*ddqa(2,:)-2*dqb(1,:).*dqa(2,:)-ddqb(1,:).*qa(2,:) + qa(4,:).*ddqb(3,:)+2*dqa(4,:).*dqb(3,:)+ddqa(4,:).*qb(3,:) + qb(4,:).*ddqa(3,:)+2*dqb(4,:).*dqa(3,:)+ddqb(4,:).*qa(3,:)
    qa(4,:).*ddqb(4,:)+2*dqa(4,:).*dqb(4,:)+ddqa(4,:).*qb(4,:) - qa(1,:).*ddqb(1,:)-2*dqa(1,:).*dqb(1,:)-ddqa(1,:).*qb(1,:) - qa(2,:).*ddqb(2,:)-2*dqa(2,:).*dqb(2,:)-ddqa(2,:).*qb(2,:) - qa(3,:).*ddqb(3,:)-2*dqa(3,:).*dqb(3,:)-ddqa(3,:).*qb(3,:)
];

for i = 1:size(q,2)
    q(:,i) = s_q_dot_q(q1, q(:,i));
    dq(:,i) = s_q_dot_q(q1, dq(:,i));
    ddq(:,i) = s_q_dot_q(q1, ddq(:,i));
end


end