clc
clear all
%% implement tent map
M = 100; % size of lattice
N = 300; % length of time series
dat = zeros(N,M); % entire data matrix
dat(1,:) = rand(1,M); % random initialization
epsilon = 0.2; % coupling strength
% make it periodic just in case
for i = 1:(N-1)
    for j = 1:M

        if j == 1
            inside = (1-epsilon) * dat(i,j) + epsilon * dat(i,M);
        
        else
        
        inside = epsilon * dat(i,j-1) + (1-epsilon) * dat(i,j);
        end
        if inside < 0.5
            fx = 2 * inside;
        elseif inside >= 0.5
            fx = 2-2*inside;
        end
        dat(i+1,j) = fx;
    end
end
plot(dat)
%% implement ulam map
close all
M = 100; % size of lattice
N = 50000; % length of time series
dat = zeros(N,M); % entire data matrix
dat(1,:) = rand(1,M); % random initialization
epsilon = 1; % coupling strength
for i = 1:(N-1)
    for j = 1:M
        if j == 1
            inside = (1-epsilon) * dat(i,j) + epsilon * dat(i,M);
        
        else
        
        inside = epsilon * dat(i,j-1) + (1-epsilon) * dat(i,j);
        end
        fx = 2 - inside^2;
        dat(i+1,j) = fx;
    end
end
plot(dat)

%% 
r = 0.01;
var1 = dat(:,4); % var2 causes var1
var2 = dat(:,9);
N_states = 50;

var1_states = linspace(min(var1), max(var1), N_states);
var2_states = linspace(min(var2), max(var2), N_states);

[xnp1, xn, yn] = sort3(var1,var2);
T = 0;

for i = 1:N_states
    xnp1_i = var1_states(i);
    for j = 1:N_states
        xn_i = var1_states(j);
        for k = 1:N_states
            yn_i = var2_states(k);
            pijk = joint3(xnp1, xn, yn, xnp1_i, xn_i, yn_i,r);
            pj = joint1(xn, xn_i, r);
            pjk = joint2(xn, yn, xn_i, yn_i, r);
            pij = joint2(xnp1, xn, xnp1_i, xn_i,r);
            T_here = pijk*log2(pijk * pj / (pjk * pij));
            %[pijk, pij, pjk, pj]
            if isnan(T_here) 
                T_here = 0;
            elseif isinf(T_here)
                T_here = 0;
                
            end
            T = T_here + T;
        end
    end
end
T

%%
% set r to 0.05 seem to work, 
% longer time series seem to reduce non-related edges
test = 1:10;
output = zeros(1,length(test));
test_case = 5;
for i = 1:length(test)
    if test(i) ~= test_case
        output(i) = transferEnt2(dat(:,test(i)), dat(:,test_case), 0.05, 50);
    end
    i
end
%%
figure
plot(test(test~=test_case),output(test~=test_case), 'linewidth', 3)
hold on
plot(test(test~=test_case),output(test~=test_case), 'x', 'MarkerSize', 5,'LineWidth',3)
hold off
xlabel(['transfer entropy between said box and box', num2str(test_case)])
ylabel('transfer entropy')
%% recreate figure 1
eps = 0.001:0.005:0.05;
M = 10; N = 50000;
x1 = 3; x2 = 2;
res = zeros(1,length(eps));
for i = 1:length(eps)
    i
    newdat = tent_map(M,N,eps(i));
    newdat = newdat(45000:end,:);
    T = transferEnt2(newdat(:,x1), newdat(:,x2), 0.05, 50);
    res(i) = T;
end

%%

a = 6.5;
x = linspace(0,0.05,100);
theo = (a*x).^2 ./ log(2) + res(1) ;
close all
plot(eps,res, 'x', 'LineWidth', 3)
hold on
plot(x,theo, 'LineWidth', 1)
legend('observation','theory')
hold off
xlabel(' $\epsilon$', 'Interpreter',"latex");
ylabel('transfer entropy')

%% recreate figure2
eps = [0.01,0.1,0.18,0.19,0.2,0.3,0.4,0.5,0.6,0.7,0.81,0.82,0.83,0.9,0.95];
%eps = [0.01 0.5 0.9]
M = 10; N = 100000;
x1 = 3; x2 = 2;
res12 = zeros(1,length(eps));
res21 = zeros(1,length(eps));
for i = 1:length(eps)
    i
    newdat = ulam_map(M,N,eps(i));
    newdat = newdat(75000:end,:);
    T = transferEnt2(newdat(:,x1), newdat(:,x2), 0.1, 50);
    res12(i) = T;
    T = transferEnt2(newdat(:,x2), newdat(:,x1), 0.1, 50);
    res21(i) = T;
end

%%
close all
plot(eps,res12, 'linewidth',3)
hold on
plot(eps,res21, 'linewidth',3)
hold off
legend('$T_{21}$','$T_{12}$', 'Interpreter',"latex")
hold off
xlabel(' $\epsilon$', 'Interpreter',"latex");
ylabel('transfer entropy')
%% functions
function T = transferEnt2(var1, var2, r, N_states)

var1_states = linspace(min(var1), max(var1), N_states);
var2_states = linspace(min(var2), max(var2), N_states);

[xnp1, xn, yn] = sort3(var1,var2); %var2 causes var1
T = 0;
pijk = zeros(N_states,N_states,N_states);
pj = zeros(1,N_states);
pjk = zeros(N_states,N_states);
pij = zeros(N_states,N_states);

for i = 1:N_states
    xnp1_i = var1_states(i);
    for j = 1:N_states
        xn_i = var1_states(j);
        for k = 1:N_states
            yn_i = var2_states(k);
            pijk(i,j,k) = joint3(xnp1, xn, yn, xnp1_i, xn_i, yn_i,r);
            pj(j) = joint1(xn, xn_i, r);
            pjk(j,k) = joint2(xn, yn, xn_i, yn_i, r);
            pij(i,j) = joint2(xnp1, xn, xnp1_i, xn_i,r);
            
          
        end
    end
end
pijk = pijk/sum(pijk,'all');
pj = pj/sum(pj,'all');
pjk = pjk/sum(pjk,'all');
pij = pij/sum(pij,'all');

for i = 1:N_states
    for j = 1:N_states
        for k = 1:N_states
            T_here = pijk(i,j,k)*log2(pijk(i,j,k) * pj(j) / (pjk(j,k) * pij(i,j)));
            
            if isnan(T_here) 
                T_here = 0;
            elseif isinf(T_here)
                T_here = 0;
                
            end
            T = T_here + T;
          
        end
    end
end
end
function p = joint3(xnp1, xn, yn, xnp1_i, xn_i, yn_i, r)
% [xnp1, xn, yn] are data, [xnp1_i, xn_i, yn_i] are location at which
% probability is to be computed
N = length(yn);
res = zeros(N,1);
for i = 1:N
    %res(i) = r - NORM([xnp1_i - xnp1(i); xn_i - xn(i); yn_i - yn(i)]);
    res(i) = ((r-abs(xnp1_i-xnp1(i)))>0) * ((r-abs(xn_i-xn(i)))>0) * ((r-abs(yn_i-yn(i)))>0);
end
p = 1/N * sum(res > 0);
end

function p = joint2(x1, x2, x1_i, x2_i, r)
N = length(x1);
res = zeros(N,1);
for i = 1:N
    %res(i) = r - NORM([xnp1_i - xnp1(i); xn_i - xn(i); yn_i - yn(i)]);
    res(i) = ((r-abs(x1_i-x1(i)))>0) * ((r-abs(x2_i-x2(i)))>0); 
end
p = 1/N * sum(res > 0);
end

function p = joint1(x1, x1_i, r)
N = length(x1);
res = zeros(N,1);
for i = 1:N
    res(i) = r - abs(x1(i)-x1_i);
end
p = 1/N * sum(res > 0);
end

function [xnp1, xn, yn] = sort3(x,y)
L = length(y);
xnp1 = x(2:L);
xn = x(1:(L-1));
yn = y(1:(L-1));
end

function T = transferEnt(var1, var2, r, N_states)

var1_states = linspace(min(var1), max(var1), N_states);
var2_states = linspace(min(var2), max(var2), N_states);

[xnp1, xn, yn] = sort3(var1,var2); %var2 causes var1
T = 0;

for i = 1:N_states
    xnp1_i = var1_states(i);
    for j = 1:N_states
        xn_i = var1_states(j);
        for k = 1:N_states
            yn_i = var2_states(k);
            pijk = joint3(xnp1, xn, yn, xnp1_i, xn_i, yn_i,r);
            pj = joint1(xn, xn_i, r);
            pjk = joint2(xn, yn, xn_i, yn_i, r);
            pij = joint2(xnp1, xn, xnp1_i, xn_i,r);
            T_here = pijk*log2(pijk * pj / (pjk * pij));
            %[pijk, pij, pjk, pj]
            if isnan(T_here) 
                T_here = 0;
            elseif isinf(T_here)
                T_here = 0;
                
            end
            T = T_here + T;
        end
    end
end
end

function dat = tent_map(M,N,epsilon)
dat = zeros(N,M); % entire data matrix
dat(1,:) = rand(1,M); % random initialization
% make it periodic just in case
for i = 1:(N-1)
    for j = 1:M

        if j == 1
            inside = (1-epsilon) * dat(i,j) + epsilon * dat(i,M);
        
        else
        
        inside = epsilon * dat(i,j-1) + (1-epsilon) * dat(i,j);
        end
        if inside < 0.5
            fx = 2 * inside;
        elseif inside >= 0.5
            fx = 2-2*inside;
        end
        dat(i+1,j) = fx;
    end
end
end

function dat = ulam_map(M,N,epsilon)
dat = zeros(N,M); % entire data matrix
dat(1,:) = rand(1,M); % random initialization
for i = 1:(N-1)
    for j = 1:M
        if j == 1
            inside = (1-epsilon) * dat(i,j) + epsilon * dat(i,M);
        
        else
        
        inside = epsilon * dat(i,j-1) + (1-epsilon) * dat(i,j);
        end
        fx = 2 - inside^2;
        dat(i+1,j) = fx;
    end
end
end