%% test squad
theta = pi-0.8;
dis = pi/6;

q0 = [[1,0,0]*sin(dis/2),cos(dis/2)];
q1 = [0,0,0,1];
q2 = [[cos(theta),sin(theta),0]*sin(dis/2),cos(dis/2)];

ds = 0.1;
s = 0:ds:1;
[q,dq,darc] = squad(q0,q1,q2,s);
q;
% [q,dq] = slerp(q0,q1,s);
% q(:,end)-q(:,end-1)

hold on
axis equal

plot3(q(1,:)',q(2,:)',q(4,:)')

% plot(s,q)

% plot(s',darc)
% darc = sqrt(sum(dq.*dq))*2;
% plot(s',darc)
% darc = sqrt(sum(diff(q')'.*diff(q')'))*2/ds;
% plot(s(2:end),darc)

% plot(s',dq(4,:)')
% plot(s(2:end)',diff(q(4,:)')/ds)



qa = s_q_dot_q(s_q_inv(q0),q2);
qb = s_q_dot_q(s_q_inv(q0),q1);

v1 = qa(1:3)*acos(qa(4))/sin(acos(qa(4)));
v2 = qb(1:3)*acos(qb(4))/sin(acos(qb(4)));

norm(v1+2*v2)*2




% 
% % 数值积分1：基于贝塞尔曲线计算
% len1 = sum(sqrt(sum(diff(p')'.*diff(p')')))
% 
% % 数值积分2：基于贝塞尔曲线导数计算
% len2 = sum(sqrt(sum(dp.*dp)))*ds
% 
% % 数值积分3：
% len3 = bezier2_length(p0,p1,p2,1)
%%
clear;
s=0:0.001:1;

hold on
axis equal
plot(s,2*s.*(1-s))
plot(s,s)
%% generate
dis = pi/6;
theta = dis/2;
v = [1,0,0]';
q = [v*sin(theta);cos(theta)];

dq_norm = 1.5;
dq = rand(4,1);
% dq = dq/norm(dq);
dq = dq - (q'*dq)*q;
dq = dq/norm(dq) * dq_norm;
%% test

clc
dq=[        -0.432946633204834
         0.192657817822397
                         0
        0.0848368181230961
];

q = [         0.202970072417778
        0.0250782661768659
                         0
         0.978863744485556];

theta = acos(q(4));
dtheta = -dq(4)/sin(theta);
v = q(1:3)/sin(theta);
dv = (dq(1:3)-cos(theta)*dtheta*v)/sin(theta);

s=sin(theta);
c=cos(theta);

% dv = (sin(theta) - cos(theta)*dtheta)/sin(theta) * v;
% dq2 = [cos(theta)*dtheta*v+sin(theta)*dv;-sin(theta)*dtheta];

% dq-dq2





dq - [c*dtheta*v+s*dv;-s*dtheta]

norm(s*dv) - norm(dq-[c*dtheta*v; -s*dtheta])

norm(s*dv)^2 - (norm(dq)^2-   2*dq'*[c*dtheta*v; -s*dtheta]  +dtheta^2)

% || dv ||^2 = (|| dq ||^2 - dtheta^2) / s^2
norm(dv)^2 - (norm(dq)^2  - dtheta^2) / s^2

dv'*v


%%
theta = pi-0.8;
dis = pi/6;

q0 = [[1,0,0]*sin(dis/2),cos(dis/2)]';
q1 = [0,0,0,1]';
q2 = [[cos(theta),sin(theta),0]*sin(dis/2),cos(dis/2)]';

ds = 0.001;
s = 0:0.001:1;

theta0 = acos(q0(4));
v0 = q0(1:3)/sin(theta0);

qa  = [sin(theta0*((1-s).^2)).*v0;cos(theta0*((1-s).^2))];
dqa = theta0*(2*s-2).*[cos(theta0*((1-s).^2)).*v0;-sin(theta0*((1-s).^2))];
ddqa = theta0*2.*[cos(theta0*((1-s).^2)).*v0;-sin(theta0*((1-s).^2))]...
    +theta0*(2*s-2).*(2*s-2).*[-sin(theta0*((1-s).^2)).*v0;-cos(theta0*((1-s).^2))];

size(diff(qa')')
size(dqa(:,2:end))
subplot(2,3,1)
hold on
plot(s,dqa)
plot(s(2:end),diff(qa')'/ds)
hold off
subplot(2,3,2)
plot(s(2:end),dqa(:,2:end)-diff(qa')'/ds)
subplot(2,3,3)
plot(s(2:end),ddqa(:,2:end)-diff(dqa')'/ds)
%%
theta2 = acos(q2(4));
v2 = q2(1:3)/sin(theta2);

qb1  = [sin(theta0*(-(1-s).^2)).*v0;cos(theta0*(-(1-s).^2))];
dqb1 = -theta0*(2*s-2).*[cos(theta0*((1-s).^2)).*v0;sin(theta0*((1-s).^2))];

subplot(2,2,3);
hold on
plot(s,dqb1);
plot(s(2:end),diff(qb1')'/ds);
hold off
subplot(2,2,4);
plot(s(2:end),dqb1(:,2:end)-diff(qb1')'/ds);
% qba = s_q_dot_q()

qb2  = [sin(theta2*(s.^2)).*v2;cos(theta2*(s.^2))];
dqb2 = 2*theta2*s.*[cos(theta2*(s.^2)).*v2;-sin(theta2*(s.^2))];

subplot(2,2,3);
plot(s,dqb2);
hold on
plot(s(2:end),diff(qb2')'/ds);
hold off
subplot(2,2,4);
plot(s(2:end),dqb2(:,2:end)-diff(qb2')'/ds);

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

% 验证dq3
subplot(2,2,3);
plot(s,dqb3);
hold on
plot(s(2:end),diff(qb3')'/ds);
hold off
subplot(2,2,4);
plot(s(2:end),dqb3(:,2:end)-diff(qb3')'/ds);

% qb32=zeros(size(qb1));
% for i=1:size(qb1,2)
% qb32(:,i) = s_q_dot_q(qb1(:,i),qb2(:,i));
% end
% qb3 - qb32(1:4,:)

% q&dq to v theta dv dtheta
theta_b3 = zeros(1,size(qb3,2));
vb3 = zeros(3,size(qb3,2));
dtheta_b3 = zeros(1,size(qb3,2));
dvb3 = zeros(3,size(qb3,2));
for i=1:size(qb3,2)
[theta_b3(:,i),vb3(:,i),dtheta_b3(:,i),dvb3(:,i)] = s_q_to_theta_v(qb3(:,i),dqb3(:,i));
end

% 计算qb 和 dqb
qb  = [sin(theta_b3.*s).*vb3;cos(theta_b3.*s)];
dqb = [
    cos(s.*theta_b3).*(theta_b3+s.*dtheta_b3).*vb3+sin(s.*theta_b3).*dvb3
    -sin(s.*theta_b3).*(theta_b3+s.*dtheta_b3)
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



% q = qa.*qb;
% dq = dqa.*qb+qa.*dqb;



subplot(2,2,3);
plot(s,dq);
hold on
plot(s(2:end),diff(q')'/ds);
hold off
subplot(2,2,4);
plot(s(2:end),dq(:,2:end)-diff(q')'/ds);




