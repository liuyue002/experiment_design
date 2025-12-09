% This script visualize local parameter sensitivity for the logistic model.
% This accompanies the paper "Optimal experiment design for practical 
% parameter identifiability and model discrimination" by 
% Liu, Baker, Maini, 2025, and illustrates the steps to create Fig.3.

%% Finite difference with respect to K with exact model solution

C0=100;
nt=100;
T=25;
t=linspace(0,T,nt);
params=[0.45,0.15,1,3900];
dp=0.001;
sol1 = sol_richards(t,params,C0);
solK = sol_richards(t,[params(1:3),params(4)+dp],C0);
diffK=(solK-sol1)./dp;
solr = sol_richards(t,[params(1)+dp,params(2:4)],C0);
diffr=(solr-sol1)./dp;
sold = sol_richards(t,[params(1),params(2)+dp,params(3:4)],C0);
diffd=(sold-sol1)./dp;


%%
figK=figure;
plot(t,diffK);
xlabel('$t$','Interpreter','latex');
ylabel('$\phi_{K}$','Interpreter','latex');
xlim([0,25]);
xticks(0:5:25);
xtickangle(0);
ytickformat('%.1f');

figr=figure;
plot(t,diffr);
xlabel('$t$','Interpreter','latex');
ylabel('$\phi_{r}$','Interpreter','latex');
%title('sensitivity');
xlim([0,25]);
xticks(0:5:25);
xtickangle(0);
ylim([0,10000]);
yticks(0:5000:10000);

figd=figure;
plot(t,diffd);
xlabel('$t$','Interpreter','latex');
ylabel('$\phi_{\delta}$','Interpreter','latex');
xlim([0,25]);
xticks(0:5:25);
xtickangle(0);
ylim([-12000,0]);
yticks(-12000:4000:0);
