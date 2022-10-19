%% init data
theta = 0;
dis = pi/3;

q0 = [[1,0,0]*sin(dis/2),cos(dis/2)]';
q1 = [0,0,0,1]';
q2 = [[cos(theta),sin(theta),0]*sin(dis/2),cos(dis/2)]';

ds = 0.0001;
s = 0:ds:1;

theta0 = acos(q0(4));
v0 = q0(1:3)/sin(theta0);

theta2 = acos(q2(4));
v2 = q2(1:3)/sin(theta2);

debug = 1;

%% qa
qa  = [sin(theta0*((1-s).^2)).*v0;cos(theta0*((1-s).^2))];
dqa = theta0*(2*s-2).*[cos(theta0*((1-s).^2)).*v0;-sin(theta0*((1-s).^2))];
ddqa = theta0*2.*[cos(theta0*((1-s).^2)).*v0;-sin(theta0*((1-s).^2))]...
    +theta0*(2*s-2).*theta0.*(2*s-2).*[-sin(theta0*((1-s).^2)).*v0;-cos(theta0*((1-s).^2))];

if(debug)
    subplot(2,3,1);
    hold on
    plot(s,dqa)
    plot(s(2:end),diff(qa')'/ds)
    hold off
    subplot(2,3,2)
    plot(s(2:end),dqa(:,2:end)-diff(qa')'/ds)
    subplot(2,3,3)
    plot(s(2:end),ddqa(:,2:end)-diff(dqa')'/ds)
end
%% qb1
qb1  = [sin(theta0*(-(1-s).^2)).*v0;cos(theta0*(-(1-s).^2))];
dqb1 = -theta0*(2*s-2).*[cos(theta0*((1-s).^2)).*v0;sin(theta0*((1-s).^2))];
ddqb1 = -theta0*2.*[cos(theta0*((1-s).^2)).*v0;sin(theta0*((1-s).^2))]...
    -theta0*(2*s-2).*theta0.*(2*s-2).*[-sin(theta0*((1-s).^2)).*v0;cos(theta0*((1-s).^2))];

if(debug)
    subplot(2,3,4);
    hold on
    plot(s,dqb1);
    plot(s(2:end),diff(qb1')'/ds);
    hold off
    subplot(2,3,5);
    plot(s(2:end),dqb1(:,2:end)-diff(qb1')'/ds);
    subplot(2,3,6)
    plot(s(2:end),ddqb1(:,2:end)-diff(dqb1')'/ds)
end

%% qb2
qb2  = [sin(theta2*(s.^2)).*v2;cos(theta2*(s.^2))];
dqb2 = 2*theta2*s.*[cos(theta2*(s.^2)).*v2;-sin(theta2*(s.^2))];
ddqb2 = 2*theta2.*[cos(theta2*(s.^2)).*v2;-sin(theta2*(s.^2))]...
    + 2*theta2*s.*2*theta2.*s.*[-sin(theta2*(s.^2)).*v2;-cos(theta2*(s.^2))];

if(debug)
    subplot(2,3,4);
    plot(s,dqb2);
    hold on
    plot(s(2:end),diff(qb2')'/ds);
    hold off
    subplot(2,3,5);
    plot(s(2:end),dqb2(:,2:end)-diff(qb2')'/ds);
    subplot(2,3,6);
    plot(s(2:end),ddqb2(:,2:end)-diff(dqb2')'/ds);
end
%% qb3
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
if(debug)
% 验证dq3
subplot(2,3,4);
plot(s,dqb3);
hold on
plot(s(2:end),diff(qb3')'/ds);
hold off
subplot(2,3,5);
plot(s(2:end),dqb3(:,2:end)-diff(qb3')'/ds);
subplot(2,3,6)
plot(s(2:end),ddqb3(:,2:end)-diff(dqb3')'/ds);
% qb32=zeros(size(qb1));
% for i=1:size(qb1,2)
% qb32(:,i) = s_q_dot_q(qb1(:,i),qb2(:,i));
% end
% qb3 - qb32(1:4,:)
end
%% qb
% 计算qb 和 dqb

% q&dq to v theta dv dtheta
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

if(debug)
    subplot(2,3,2);
    plot(s,dvb3);
    hold on
    plot(s(2:end),diff(vb3')'/ds);
    hold off
    subplot(2,3,2);
    plot(s(2:end),dvb3(:,2:end)-diff(vb3')'/ds);
    subplot(2,3,3)
    plot(s(2:end),ddvb3(:,2:end)-diff(dvb3')'/ds);
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

if(debug)
    subplot(2,3,4);
    plot(s,dqb);
    hold on
    plot(s(2:end),diff(qb')'/ds);
    hold off
    subplot(2,3,5);
    plot(s(2:end),dqb(:,2:end)-diff(qb')'/ds);
    subplot(2,3,6)
    plot(s(2:end),ddqb(:,2:end)-diff(dqb')'/ds);
end
%% q
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
if(debug)
    subplot(2,3,4);
    plot(s,dq);
    hold on
    plot(s(2:end),diff(q')'/ds);
    hold off
    subplot(2,3,5);
    plot(s(2:end),dq(:,2:end)-diff(q')'/ds);
    subplot(2,3,6)
    plot(s(2:end),ddq(:,2:end)-diff(dq')'/ds);
end

%% compute arc
darc = sqrt(sum(dq.*dq));
ddarc = sum(dq.*ddq)./darc;

[arc_est, darc_est, ddarc_est] = s_estimate_bezier3_arc(darc(1), ddarc(1), darc(end), ddarc(end), darc((length(darc) + 1)/2), s);

subplot(1,2,1)
hold on
axis equal
plot(s,darc)
plot(s,ddarc)
plot(s,darc_est)
plot(s,ddarc_est)
% plot(s,2*a.*s+b)

subplot(1,2,2)
hold on
axis equal
plot(s(2:end)-ds/2,sqrt(sum(diff(q')'.*diff(q')'))/ds)
% plot(s,ddarc)

plot(s,darc)
plot(s,ddarc)
% plot(s,a.*s.*s+b.*s+c+ d*sin(s*pi) + e*cos(s*2*pi) + f)
% % plot(s,a.*s.*s+b.*s+c+ d*sin(s*pi) + e*cos(s*2*pi) + f + (c/7)*(1-cos(4*pi*s))/2)
% plot(s,2*a.*s+b+d*cos(s*pi)*pi - e*sin(s*2*pi)*2*pi)









