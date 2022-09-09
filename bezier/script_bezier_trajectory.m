%% test multi ee, make data
m  = 100;% point num
n  = 5; % ee num

v0 = zeros(1,n);
p0 = zeros(1,n);

pos      = rand(m,n) * 5;
pos(7,1) = 0.0;

vel  = ones(m,n) * 5;
acc  = ones(m,n) * 10;
jerk = ones(m,n) * 100;


max_vb = zeros(m,n);
% max_vb(1:4,:)   = rand(4,n)*4.0;
max_vb(1:4,:)   = 5.0;
max_vb(5,:)     = 0;
max_vb(6:end,:) = 5.0;
max_vb(end,:)   = 0;
for i = 2:size(pos,1)
    pos(i,:) = pos(i-1,:) + pos(i,:);
end
%% compute
p0     = zeros(1, n);
v0     = zeros(1, n);

[T,Ta,Tb,vb, real_v] = s_make_s_curve_multiple(p0,v0,pos,max_vb,vel,acc,jerk);
%%
total_t = sum(T);
dt = 0.01;
t = 0:dt:total_t;
p = s_s_curve_multiple(t,p0,v0,pos,vb,real_v,acc,jerk,T,Ta,Tb);
v = diff(p)/dt;
a = diff(v)/dt;

% for j=1:n
    subplot(1,3,1);
    plot(t,p);
%     axis equal
    for i=1:m
        line([sum(T(1:i)) sum(T(1:i))],[0 1.1*max(pos(end,:))],'linestyle','--', 'Color','k');
    end
    subplot(1,3,2);
    plot(t(2:end)-dt/2,v);
%     axis equal
    for i=1:m
        line([sum(T(1:i)) sum(T(1:i))],[-1.1*max(vel(end,:)) 1.1*max(vel(end,:))],'linestyle','--', 'Color','k');
    end

    subplot(1,3,3);
    plot(t(3:end)-dt/2*3,a);
%     axis equal
    for i=1:m
        line([sum(T(1:i)) sum(T(1:i))],[-1.1*max(acc(end,:)) 1.1*max(acc(end,:))],'linestyle',':', 'Color','k');
    end
% end

% s_s_curve_multiple(2.09,p0,v0,pos,vb,real_v,acc,jerk,T,Ta,Tb)
% s_s_curve_multiple(2.1,p0,v0,pos,vb,real_v,acc,jerk,T,Ta,Tb)