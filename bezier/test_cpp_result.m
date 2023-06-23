%%
m = dlmread('C:\Users\py033\Desktop\test_data\pes.txt');
m_v = dlmread('C:\Users\py033\Desktop\test_data\vpes.txt');
m_a = dlmread('C:\Users\py033\Desktop\test_data\apes.txt');
plot(m(:,:))

%%
ee_b = 1;
ee_e = 2;
% ee_b = 10;
% ee_e = 12;

subplot(3,1,1)
hold on
plot(diff(m(1:end,ee_b:ee_e))*1e3)
plot(m_v(1:end,ee_b:ee_e))
subplot(3,1,2)
hold on
% plot(diff(diff(m(1:end,ee_b:ee_e)))*1e6)
plot(m_a(1:end,ee_b:ee_e))
subplot(3,1,3)
hold on
plot(diff(diff(diff(m(1:end,ee_b:ee_e))))*1e9)

%%
b = 27000;
% e = size(m,1);
e = 32000;
plot3(m(b:e,10),m(b:e,11),m(b:e,12))
axis equal

%%
subplot(1,2,1)
plot(diff(m(:,1:3)))
subplot(1,2,2)
plot(diff(diff(m(:,1:3))))
% plot(diff(diff(m(:,1:3))))