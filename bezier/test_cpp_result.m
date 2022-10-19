%%
m = dlmread('C:\Users\py033\Desktop\test_data\pes.txt');
m_v = dlmread('C:\Users\py033\Desktop\test_data\vpes.txt');
m_a = dlmread('C:\Users\py033\Desktop\test_data\apes.txt');
plot(m(:,:))

%%
ee_b = 3;
ee_e = 5;
% ee_b = 1;
% ee_e = size(m,2);

subplot(3,1,1)
hold on
plot(diff(m(1:end,ee_b:ee_e))*1000)
plot(m_v(1:end,ee_b:ee_e))
subplot(3,1,2)
hold on
plot(diff(diff(m(1:end,ee_b:ee_e)))*1e6)
plot(m_a(1:end,ee_b:ee_e))
subplot(3,1,3)
hold on
plot(diff(diff(diff(m(1:end,ee_b:ee_e)))))

%%
b = 20000;
% e = size(m,1);
e = 22000;
plot3(m(b:e,3),m(b:e,4),m(b:e,5))
axis equal