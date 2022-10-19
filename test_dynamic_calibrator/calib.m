%% clb
m = dlmread('C:/Users/py033/Desktop/calib/calib_dyn_par1.txt');

%% clb
for i=1:6
    subplot(2,3,i)
    hold on
%     plot(m(:,1 + (i-1) * 3))
    plot(diff(m(:,1 + (i-1) * 3))*1000)
%     
    plot(m(:,2 + (i-1) * 3))
%      plot(m(:,3 + (i-1) * 3)/1000)
end

%%
subplot(4,1,1)
plot(m(:,7))
hold on
subplot(4,1,2)
hold on
plot(diff(m(:,7))*1000)
plot(m(:,8)*1.414*1e4)
subplot(4,1,3)
plot(-m(:,9)*0.1454)
subplot(4,1,4)
plot(m(:,4) + m(:,7))
% plot(m(:,)*1.414*1e4)
%% clb result
pos1 = dlmread('C:/Users/py033/Desktop/calib/data_after/pos0.txt');
pos2 = dlmread('C:/Users/py033/Desktop/calib/data_after/pos1.txt');
pos3 = dlmread('C:/Users/py033/Desktop/calib/data_after/pos2.txt');
pos4 = dlmread('C:/Users/py033/Desktop/calib/data_after/pos3.txt');
pos5 = dlmread('C:/Users/py033/Desktop/calib/data_after/pos4.txt');
pos6 = dlmread('C:/Users/py033/Desktop/calib/data_after/pos5.txt');

acc1 = dlmread('C:/Users/py033/Desktop/calib/data_after/acc0.txt');
acc2 = dlmread('C:/Users/py033/Desktop/calib/data_after/acc1.txt');
acc3 = dlmread('C:/Users/py033/Desktop/calib/data_after/acc2.txt');
acc4 = dlmread('C:/Users/py033/Desktop/calib/data_after/acc3.txt');
acc5 = dlmread('C:/Users/py033/Desktop/calib/data_after/acc4.txt');
acc6 = dlmread('C:/Users/py033/Desktop/calib/data_after/acc5.txt');

vel1 = dlmread('C:/Users/py033/Desktop/calib/data_after/vel0.txt');
vel2 = dlmread('C:/Users/py033/Desktop/calib/data_after/vel1.txt');
vel3 = dlmread('C:/Users/py033/Desktop/calib/data_after/vel2.txt');
vel4 = dlmread('C:/Users/py033/Desktop/calib/data_after/vel3.txt');
vel5 = dlmread('C:/Users/py033/Desktop/calib/data_after/vel4.txt');
vel6 = dlmread('C:/Users/py033/Desktop/calib/data_after/vel5.txt');

fce1 = dlmread('C:/Users/py033/Desktop/calib/data_after/fce0.txt');
fce2 = dlmread('C:/Users/py033/Desktop/calib/data_after/fce1.txt');
fce3 = dlmread('C:/Users/py033/Desktop/calib/data_after/fce2.txt');
fce4 = dlmread('C:/Users/py033/Desktop/calib/data_after/fce3.txt');
fce5 = dlmread('C:/Users/py033/Desktop/calib/data_after/fce4.txt');
fce6 = dlmread('C:/Users/py033/Desktop/calib/data_after/fce5.txt');
f1 = dlmread('C:/Users/py033/Desktop/calib/data_after/f0.txt');
f2 = dlmread('C:/Users/py033/Desktop/calib/data_after/f1.txt');
f3 = dlmread('C:/Users/py033/Desktop/calib/data_after/f2.txt');
f4 = dlmread('C:/Users/py033/Desktop/calib/data_after/f3.txt');
f5 = dlmread('C:/Users/py033/Desktop/calib/data_after/f4.txt');
f6 = dlmread('C:/Users/py033/Desktop/calib/data_after/f5.txt');
%% clb result
subplot(2,3,1)
hold on
% plot(f1)
% plot(fce1)
plot(vel1)

subplot(2,3,2)
hold on
% plot(f2)
% plot(fce2)
plot(vel2)

subplot(2,3,3)
hold on
% plot(f3)
% plot(fce3)
plot(vel3)

subplot(2,3,4)
hold on
% plot(f4)
% plot(fce4)
plot(vel4)

subplot(2,3,5)
hold on
% plot(f5)
% plot(fce5)
plot(vel5)

subplot(2,3,6)
hold on
% plot(f6)
% plot(fce6)
plot(vel6)


%% 
pnts=[-0.1444   0.0024   -0.1046
0.3309   0.1494   -0.1546
-0.2527   0.3793   -0.0256
0.2226   0.5263   -0.0755
-0.0226   -0.1263   0.6755
0.4527   0.0207   0.6256
-0.1309   0.2506   0.7546
0.3444   0.3976   0.7046


0.0400   0.1760   0.1319
0.1390   0.1859   0.1219
0.0340   0.3077   0.2034
0.1330   0.3176   0.1935
0.0670   0.0824   0.3065
0.1660   0.0923   0.2966
0.0610   0.2141   0.3781
0.1600   0.2240   0.3681]

pnts1 = pnts(1:8,:);

pnts2 = pnts(9:end,:);

hold on
axis equal
axis on

scatter3(pnts1(6,1), pnts1(6,2), pnts1(6,3))
vertex=pnts1;
facet=[1 2 4 3;1 2 6 5;1 3 7 5;2 4 8 6;3 4 8 7;5 6 8 7];
% color=[0;0;0;0;1;1;1;1];
patch('Vertices',vertex,'Faces',facet,'FaceVertexCData',color,'FaceColor','interp','FaceAlpha',0.5);

scatter3(pnts2(:,1), pnts2(:,2), pnts2(:,3))
vertex=pnts2;
facet=[1 2 4 3;1 2 6 5;1 3 7 5;2 4 8 6;3 4 8 7;5 6 8 7];
% color=[0;0;0;0;1;1;1;1];
patch('Vertices',vertex,'Faces',facet,'FaceVertexCData',color,'FaceColor','interp','FaceAlpha',0.5);

%%
pnts2 = [-0.0332   0.3747   0.7546
0.1747   0.3956   0.5394
0.1202   0.6976   0.9342
0.3281   0.7184   0.7189
0.2719   0.0816   1.0211
0.4798   0.1024   0.8058
0.4253   0.4044   1.2006
0.6332   0.4253   0.9854];

sphere_center = [0.4,0.4,0.8];
radius = 0.1;

hold on
axis equal
axis on

[X,Y,Z] = sphere;

X2 = X * radius;
Y2 = Y * radius;
Z2 = Z * radius;

surf(X2+sphere_center(1),Y2+sphere_center(2),Z2+sphere_center(3))

scatter3(pnts2(:,1), pnts2(:,2), pnts2(:,3))
vertex=pnts2;
facet=[1 2 4 3;1 2 6 5;1 3 7 5;2 4 8 6;3 4 8 7;5 6 8 7];
% color=[0;0;0;0;1;1;1;1];
patch('Vertices',vertex,'Faces',facet,'FaceVertexCData',color,'FaceColor','interp','FaceAlpha',0.5);

%%

sphere1_center = [0.3,0.3,0.3];
radius1 = 0.1;
sphere2_center = [0.4,0.44,0.2];
radius2 = 0.3;

hold on
axis equal
axis on

[X,Y,Z] = sphere;

X2 = X * radius1;
Y2 = Y * radius1;
Z2 = Z * radius1;

surf(X2+sphere1_center(1),Y2+sphere1_center(2),Z2+sphere1_center(3))

[X,Y,Z] = sphere;

X3 = X * radius2;
Y3 = Y * radius2;
Z3 = Z * radius2;

surf(X3+sphere2_center(1),Y3+sphere2_center(2),Z3+sphere2_center(3))