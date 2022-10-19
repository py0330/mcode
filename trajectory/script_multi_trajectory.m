%% test multi ee, make data
clear;
m  = 10;% point num
n  = 5; % ee num

v0 = zeros(1,n);
p0 = zeros(1,n);

pos      = rand(m,n) * 10;
pos(7,1) = 0.0;

max_v  = (rand(1,1) + 0.01)*10;
max_a  = (rand(1,1) + 0.01)*10;
max_j  = (rand(1,1) + 0.01)*10;

vel  = ones(m,n) * max_v;
acc  = ones(m,n) * max_a;
jerk = ones(m,n) * max_j;

max_vb = zeros(m,n);
max_vb(1:4,:)   = min(rand(4,n)*0.5, max_v);
max_vb(5,:)     = 0;
max_vb(6:end,:) = min(rand(m-5,n) + 0.5,ones(m-5,n)) * max_v;
max_vb(end,:)   = 0;
for i = 2:size(pos,1)
    pos(i,:) = pos(i-1,:) + pos(i,:);
end
% jerk
%% compute
rows = 1:m;
if(rows(1) == 1)
p0     = zeros(1, n);
v0     = zeros(1, n);
else
p0     = pos(rows(1)-1,:);
v0     = zeros(1, n);
end
tic
[T,Ta,Tb,vb, real_v, mode] = s_make_s_curve_multiple(p0,v0,pos(rows,:),max_vb(rows,:),vel(rows,:),acc(rows,:),jerk(rows,:));
toc
%%
cols = 1:n;

total_t = sum(T);
dt = 0.01;
t = 0:dt:total_t;
p = s_s_curve_multiple(t,p0(:,cols),v0(:,cols),pos(rows,cols),vb(:,cols),real_v(:,cols),acc(rows,cols),jerk(rows,cols),T,Ta(:,cols),Tb(:,cols),mode(:,cols));
v = diff(p)/dt;
a = diff(v)/dt;
j = diff(a)/dt;

if(max(max(v)) > max(max(vel)) + 1e-10)
    max(max(v)) - max(max(vel))
    [value, pos_at]=max(v)
    error('v error')
end
if(max(max(a)) > max(max(acc)) + 1e-7)
    max(max(a))
    error('a error')
end

%%
% for j=1:n
    subplot(1,4,1);
    plot(t,p);
%     axis equal
    for i=1:m
        line([sum(T(1:i)) sum(T(1:i))],[0 1.1*max(pos(end,:))],'linestyle','--', 'Color','k');
    end
    subplot(1,4,2);
    plot(t(2:end)-dt/2,v);
%     axis equal
    for i=1:m
        line([sum(T(1:i)) sum(T(1:i))],[-1.1*max(vel(end,:)) 1.1*max(vel(end,:))],'linestyle','--', 'Color','k');
    end

    subplot(1,4,3);
    plot(t(3:end)-dt/2*3,a);
%     axis equal
    for i=1:m
%         line([sum(T(1:i)) sum(T(1:i))],[-1.1*max(acc(end,:)) 1.1*max(acc(end,:))],'linestyle',':', 'Color','k');
    end

    subplot(1,4,4);
    plot(t(4:end)-dt/2*5,j);
%     axis equal
    for i=1:m
        line([sum(T(1:i)) sum(T(1:i))],[-1.1*max(acc(end,:)) 1.1*max(acc(end,:))],'linestyle',':', 'Color','k');
    end
% end

