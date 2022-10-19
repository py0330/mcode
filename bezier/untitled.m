%%
syms pa va Ta a j si t_ vb Tacc vc pb T Tb
p_ = pa + (va + vc) / 2.0 * Ta - (vc * (Ta - t_) - si / 6.0 * j * (Ta - t_) * (Ta - t_) * (Ta - t_));
p_ = pa + va * t_ + si / 6.0 * j * t_ * t_ * t_;
p_ = pa + (va + vc) / 2.0 * Ta - (vc * (Ta - t_) - si / 6.0 * j * (Ta - t_) * (Ta - t_) * (Ta - t_));
p_ = pa + (va + vc) / 2.0 * Ta + vc * (t_ - Ta);
p_ = pb - vb * (T - t_) + si / 6.0 * j * (T - t_) * (T - t_) * (T - t_);
p_ = pb - vb * a / j + si / 6.0 * a * a * a / j / j - (vb - si / 2.0 * a * a / j) * (T - t_ - a / j) + si * a / 2.0 * (T - t_ - a / j) * (T - t_ - a / j);
p_ = pb - (vc + vb) / 2.0 * Tb + (vc * (Tb - T + t_) + si / 6.0 * j * (Tb - T + t_) * (Tb - T + t_) * (Tb - T + t_));
p_ = pb - vb * (T - t_) + si / 6.0 * j * (T - t_) * (T - t_) * (T - t_);
p_ = pb - (vc + vb) / 2.0 * Tb + vc * (Tb - T + t_) + si / 6.0 * j * (Tb - T + t_) * (Tb - T + t_) * (Tb - T + t_);
p_ = pa + va * t_ + si / 6.0 * j * t_ * t_ * t_;
v_ = diff(p_, t_)
v_ = expand(diff(p_, t_))
a_ = diff(v_, t_)
a_ = expand(diff(v_, t_))
%%
theta0 = 0.1;
q0 = [sin(theta0)*[1,0,0],cos(theta0)]';
q1 = [0,0,0,1]';
q2 = [sin(theta0)*[1,0,0],cos(theta0)]';


[q,dq,ddq] = s_blend_quaternion_bezier3(q0,q1,q2,0.5)

%%
syms theta vx vy vz dtheta t dvx dvy dvz d2vx d2vy d2vz

% theta = dtheta*t
% 
% vx = dvx*t
% vy = dvy*t
% vz = dvz*t

q = [sin(theta)*vx; sin(theta)*vy; sin(theta)*vz; cos(theta)]
dq = diff(q,theta) * dtheta + diff(q,vx) * dvx;
%%
clear
syms theta v dtheta dv d2theta d2v t

% theta = dtheta*t
% 
% vx = dvx*t
% vy = dvy*t
% vz = dvz*t

q = [sin(theta)*v; cos(theta)]
dq = diff(q,theta) * dtheta + diff(q,v) * dv
d2q = diff(dq,theta)*dtheta + diff(dq,v)*dv + diff(dq,dtheta)*d2theta + diff(dq,dv)*d2v

