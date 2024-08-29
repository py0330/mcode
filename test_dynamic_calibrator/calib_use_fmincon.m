%% load data
A = dlmread("data_after\A.txt");
x = dlmread("data_after\x.txt");
b = dlmread("data_after\b.txt");

%% compute using moore-penrose inverse for verification
x1 = pinv(A)*b

%% compute solution without constraint
G = 0.5*A'*A;
f = -b'*A;
fun = @(x)x'*G*x + f*x + 0.5*b'*b;

x2 = fmincon(fun, x)

% test
fun(x1)-fun(x2)

%% add constraints for mass and frictions
CI = zeros(2,26)
ci = zeros(2,1);
CI(1,1) = -1;
CI(2,11) = -1;

x3 = fmincon(fun, x2, CI, ci);
[fun(x1),fun(x2),fun(x3)]

%% add constraints for mass and frictions
CI = zeros(8,26);
ci = zeros(8,1);
CI(1,1) = -1;
CI(2,11) = -1;
CI(3,21) = -1;
CI(4,22) = -1;
CI(5,23) = -1;
CI(6,24) = -1;
CI(7,25) = -1;
CI(8,26) = -1;
% CI(9,21) = 1;
% CI(10,22) = 1;
% CI(11,23) = 1;
% CI(12,24) = 1;
% CI(13,25) = 1;
% CI(14,26) = 1;
% ci(9) = 0.1;
% ci(10) = 0.1;
% ci(11) = 0.1;
% ci(12) = 0.1;
% ci(13) = 0.1;
% ci(14) = 0.1;


x4 = fmincon(fun, x2, CI, ci);
[fun(x1),fun(x2),fun(x3),fun(x4)]


%% add constraints for rotate inerter
cons_id = [4,5,6,8,9,10,14,15,16,18,19,20];
cons_v  = [0,0,0,0,0, 0, 0, 0, 0, 0, 0, 0];

CE = zeros(length(cons_id),26);
ce = zeros(length(cons_id),1);

for i=1:length(cons_id)
    CE(i,cons_id(i)) = 1;
end


% CE(2,26) = 1;
% CE(3,23) = 1;
% CE(4,24) = 1;
% CE(3,25) = 1;
% CE(4,26) = 1;

%% add constraints for inertias

x_init = [1.469, 0, 0, 0,0,0,7667.511 * 1e-6,0,0,0,...
    2.285,0,0,0,0,0,9103.438 * 1e-6,0,0,0,...
    0,0,1002 * 1e-3 * 1e-4,0,0,1002 * 1e-3 * 1e-4]'

x_init2 = [0, 0, 0, 0,0,0,0,0,0,0,...
    2,0,0,0,0,0,0,0,0,0,...
    0,0,0,0,0,0]'
x5 = fmincon(fun, x_init, CI, ci,CE,ce,[],[],@calib_constraint);
[fun(x1),fun(x2),fun(x3),fun(x4),fun(x5)]



% calib_constraint(x5)


% l1.m = 1.469;
% l1.cx = 0.009 * 1e-3;
% l1.cy = 165.468 * 1e-3;
% l1.cz = -19.239 * 1e-3;
% l1.Ixx = 6995.572 * 1e-6;
% l1.Iyy = 1357.599 * 1e-6;
% l1.Izz = 7667.511 * 1e-6;
% l1.Ixy = 0.387 * 1e-6;
% l1.Ixz = -0.332 * 1e-6;
% l1.Iyz = 467.024 * 1e-6;
% l1.ki = 1002 * 1e-3 * 1e-4;
% 
% l2.m = 2.285;
% l2.cx = 0.001 * 1e-3;
% l2.cy = 166.209 * 1e-3;
% l2.cz = 1.981 * 1e-3;
% l2.Ixx = 9145.552 * 1e-6;
% l2.Iyy = 1637.069 * 1e-6;
% l2.Izz = 9103.438 * 1e-6;
% l2.Ixy = -0.316 * 1e-6;
% l2.Ixz = -0.042 * 1e-6;
% l2.Iyz = 534.887 * 1e-6;
% l2.ki = 1002 * 1e-3 * 1e-4;
x5(1:10)'
x5(11:20)'
x5(21:end)'



%%

x_init = [1.469, 0, 0, 0,0,0,7667.511 * 1e-6,0,0,0,...
    2.285,0,0,0,0,0,9103.438 * 1e-6,0,0,0,...
    0,0,1002 * 1e-3 * 1e-4,0,0,1002 * 1e-3 * 1e-4]'

x_init2 = [0, 0, 0, 0,0,0,0,0,0,0,...
    0,0,0,0,0,0,0,0,0,0,...
    0,0,0,0,0,0]'

x_init3 = [0, 0, 0, 0,0,0,0,0,0,0,...
    0,0.5,0,0,0,0,0.1,0,0,0,...
    0,0,0,0,0,0]'
[fun(x_init),fun(x_init2),fun(x_init3)]

