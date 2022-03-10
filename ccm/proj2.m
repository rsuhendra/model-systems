clc
clear all

%% numerical solver for lorenz system
T = 100;
dt = 0.01;
sigma = 10;
beta = 8/3;
rho = 28;

L = T/dt + 1;
x = zeros(L,1); y = zeros(L,1); z = zeros(L,1); t = zeros(L,1);
x(1) = 5;
for i = 1:L
    dx = sigma * (y(i) - x(i)) + 0 * randn;
    dy = x(i) * (rho - z(i)) - y(i) + 0 * randn;
    dz = x(i) * y(i) - beta * z(i) + 0 * randn;
    
    x(i+1) = x(i) + dx*dt;
    y(i+1) = y(i) + dy*dt;
    z(i+1) = z(i) + dz*dt;
    t(i+1) = t(i) + dt;
end

plot3(x,y,z)
xlabel('x')
ylabel('y')
zlabel('z')

%% assymetric coupled system from paper
T = 100000;
x = zeros(T,1); y = zeros(T,1);
x(1) = .51; y(1) = .51;

rx = 3.8;
ry = 3.5;
b_yx = 0.1;
b_xy = 0.0002;

for i = 1:T
    x(i+1) = x(i) * (rx - rx * x(i) - b_xy * y(i)) + 0.0 * randn;
    y(i+1) = y(i) * (ry - ry * y(i) - b_yx * x(i)) + 0.0 * randn;
end
figure(1)
plot(x)
hold on
plot(y)
hold off


%% for eq1 in paper
[y_truth,y_hat] = CCM(x,y,3,4,5,500,5705);
[x_truth,x_hat] = CCM(y,x,3,4,5,500,5705);

close all
subplot(1,2,1)
plot(y_truth,y_hat,'.')
title(['x predict y, corr=' num2str(corr(y_truth, y_hat))])
subplot(1,2,2)
plot(x_truth, x_hat,'.')
title(['y predict x, corr=' num2str(corr(x_truth, x_hat))])

%% for lorenz
[y_truth,y_hat] = CCM(x,y,100,3,15,3500,10000);
[x_truth,x_hat] = CCM(y,x,100,3,15,3500,10000);
[z_truth,z_hat] = CCM(x,z,100,3,15,3500,10000);

close all
subplot(1,3,1)
plot(y_truth,y_hat,'.')
title(['x predict y, corr=' num2str(corr(y_truth, y_hat))])
subplot(1,3,2)
plot(x_truth, x_hat,'.')
title(['y predict x, corr=' num2str(corr(x_truth, x_hat))])
subplot(1,3,3)
plot(z_truth, z_hat,'.')
title(['x predict z, corr=' num2str(corr(z_truth, z_hat))])
%% sensitivity for eq1 dependence on tau
TAU = 1:2:30;
res_xy = zeros(1,length(TAU));
res_yx = zeros(1,length(TAU));
for i = 1:length(TAU)
    i
    [y_truth,y_hat] = CCM(x,y,TAU(i),4,5,500,5705);
    [x_truth,x_hat] = CCM(y,x,TAU(i),4,5,500,5705);
    
    res_xy(i) = corr(y_truth, y_hat);
    res_yx(i) = corr(x_truth, x_hat);
end

figure(1)
subplot(1,2,1)
plot(TAU, res_xy, 'LineWidth',3)
xlabel('tau')
ylabel('x predict y')
subplot(1,2,2)
plot(TAU, res_yx, 'LineWidth',3)
xlabel('tau')
ylabel('y predict x')
%% sensitivity eq1 dependence on N-dim
NDIM = 2:2:30;
res_xy = zeros(1,length(NDIM));
res_yx = zeros(1,length(NDIM));
for i = 1:length(NDIM)
    i
    [y_truth,y_hat] = CCM(x,y,3,NDIM(i),5,500,5705);
    [x_truth,x_hat] = CCM(y,x,3,NDIM(i),5,500,5705);
    
    res_xy(i) = corr(y_truth, y_hat);
    res_yx(i) = corr(x_truth, x_hat);
end

figure(1)
subplot(1,2,1)
plot(NDIM, res_xy, 'LineWidth',3)
xlabel('n dim')
ylabel('x predict y')
subplot(1,2,2)
plot(NDIM, res_yx, 'LineWidth',3)
xlabel('n dim')
ylabel('y predict x')
%% functions for CCM
function [Y_truth, Y_hat] = CCM(X, Y, tau, N_dim, K_neighbour, st, ed)
time_coor = st:ed;

M_x = zeros(N_dim,length(st:ed)); %rows refers to time dimension
Y_truth = Y(time_coor);

for i = 1:N_dim
    st1 = st - (i-1) * tau;
    ed1 = ed - (i-1) * tau;
    M_x(i,:) = X(st1:ed1);
end


distMx = squareform(pdist(M_x'));
[~,W] = size(M_x);
Y_hat = zeros(W,1);

for i = 1:W
dist = distMx(i,:); dist(i) = max(dist) + 1;
[~,order] = sort(dist,'ascend');
u = dist(order(1:K_neighbour))/min(dist);
u = exp(-u);
w = u/sum(u);

Y_predict = Y_truth(order(1:K_neighbour))' * w';
Y_hat(i) = Y_predict;

end
end