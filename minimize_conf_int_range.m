% This script demonstrate how to minimize confidence interval range as
% found using profile likelihoods.
% This accompanies the paper "Optimal experiment design for practical 
% parameter identifiability and model discrimination" by 
% Liu, Baker, Maini, 2025, and illustrates the steps to create Fig.2 and 4.

%% the default controls
u_K0=200;
tau0=10;
tau=10;
warning('off','MATLAB:ode45:IntegrationTolNotMet');

%% one-dimensional parameter scan for the control parameters
u_K0s=linspace(50,800,101);
ff=@(u_K0) logistic_bd_uk_rrange(u_K0,tau0,tau);
widths_uK0=arrayfun(ff,u_K0s);

tau0s=linspace(0,15,101);
ff=@(tau0) logistic_bd_uk_rrange(u_K0,tau0,tau);
widths_tau0=arrayfun(ff,tau0s);

taus=linspace(0,15,101);
ff=@(tau) logistic_bd_uk_rrange(u_K0,tau0,tau);
widths_tau=arrayfun(ff,taus);

%% plot one-dimensional parameter scan

fig1=figure;
plot(u_K0s,widths_uK0);
xlim([0,800]);
xticks(0:400:800);
yticks(0:2);
xlabel('$u_{\textrm{max}}$','Interpreter','latex');
ylabel('$\Delta r$','Interpreter','latex');

fig2=figure;
plot(tau0s,widths_tau0);
xlim([0,15]);
ylim([0,2]);
yticks(0:2);
xlabel('$\tau_0$','Interpreter','latex');
ylabel('$\Delta r$','Interpreter','latex');

fig3=figure;
plot(taus,widths_tau);
xlim([0,15]);
ylim([0,2]);
yticks(0:2);
xlabel('$\tau$','Interpreter','latex');
ylabel('$\Delta r$','Interpreter','latex');

%% 2D parameter scan for the control parameters

tau0ss=linspace(0,25,41);
tauss=linspace(1,25,41);
[X,Y]=meshgrid(tau0ss,tauss);
Z=arrayfun(@(tau0,tau)logistic_bd_uk_rrange(u_K0,tau0,tau), X,Y);
fig4=figure;
surf(X,Y,Z,'LineStyle','none');
xlabel('tau0');
ylabel('tau');
zlabel('r range');
clim([0,2]);
view(0,90);
colorbar;
title('rrange vs tau,tau0');
