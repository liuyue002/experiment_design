% This script demonstrate how to compute profile likelihoods for an ODE
% model, given some control parameters.
% This accompanies the paper "Optimal experiment design for practical 
% parameter identifiability and model discrimination" by 
% Liu, Baker, Maini, 2025, and illustrates the steps to create Fig.1.

% Consider the Richards model: dC/dt=r*C*(1-(C/K)^gamma) - delta*C
% parameters: [r,d,gamma,K]
% The logistic model is with gamma=1

%% Illustrate what typical solutions to Logistic model looks like
fig=figure;
tt=linspace(0,25,50);
plot(tt,sol_richards(tt,[0.3,0,1,2600],100));
xlabel('t');
ylabel('C');
title('Typical solution to logistic model');

%% Generate synthetic data
% The synthetic data is generated with ground truth parameter values:
% r=0.45, delta=0.15, gamma=1, K=3900
rng(0);
C0=100;
nt=101;
T=25;
t=linspace(0,T,nt);
params=[0.45,0.15,1,3900];

% specify the control functions u_K(t), u_d(t), u_r(t)
% the values here corresponds to the case of Fig.1(d) in the paper
% vary u_K0, tau0, and tau to recreate other subfigures
u_K0=400;
tau0=10;
tau=10;
u_d=@(t)0;
u_r=@(t)0;
u_K=@(t) ((t>tau0)&(t<(tau0+tau)))*u_K0;
clean_data=sol_richards_control(t,params,C0,u_d,u_K,u_r);
sigma=20;
noisy_data=clean_data+randn(size(clean_data))*sigma;

figure;
hold on
plot(t,clean_data);
plot(t,noisy_data);
plot(t,u_K(t));
xlabel('t');
legend('Clean data','Noisy data','u_K');

%% MLE for logistic model
fixed=[0,0,1,0]; % gamma=1 is fixed for logistic model
fixed_param_val=params;
numeric_params={T,nt,u_d,u_K,u_r};
lb_opt=[0,0,0,0];
ub_opt=[5,5,9,10000];
opt.lb=lb_opt;
opt.ub=ub_opt;
opt.logging=true;
opt.alg=2;
[mle,mle_sigma,max_l] = optimize_likelihood_general(@richards_ctrl_sq_err,fixed_param_val,fixed,noisy_data,numeric_params,C0,opt);
optimal_param_vals=fixed_param_val';
optimal_param_vals(fixed==0)=mle;

fprintf("The MLE is r=%.3f, d=%.3f, K=%.0f\n",mle);
fprintf("Whereas the ground truth is r=%.3f, d=%.3f, K=%.0f\n",params([1,2,4]));

%% MLE solution

mle_soln=sol_richards_control(t,optimal_param_vals,C0,u_d,u_K,u_r);
figure;
hold on;
plot(t,clean_data);
plot(t,noisy_data);
plot(t,mle_soln);
hold off
xlabel('t');
legend('Clean data','Noisy data','MLE solution');

%% univariate profile likelihood for logistic model with control

fixed_param_val=params;
fixed=[0,0,1,0];
num_params=size(fixed_param_val,2);
num_free_params=sum(1-fixed);
lb=[ 0.10, 0.00, 1,   600];
ub=[ 1.00, 1.00, 1, 10000];
opt.ub=[10,10,10,99000];
opt.alg=1;
param_names={'r','\delta','\gamma','K'};

numpts=41;
param_vals=zeros(num_params,numpts);
max_ls=zeros(num_params,numpts);
minimizers=cell(num_params,numpts);
% add the global optimum to the list of param vals
param_vals=[param_vals,optimal_param_vals];

for param=1:num_params
    if fixed(param)
        continue;
    end
    param_vals(param,1:numpts)=linspace(lb(param),ub(param),numpts);
    param_vals(param,:)=sort(param_vals(param,:));
    [~,mle_idx]=min(abs(param_vals(param,:)-optimal_param_vals(param)));
    fixed_params=fixed;
    fixed_params(param)=1;
    initial=fixed_param_val;
    minimizers{param,mle_idx}=optimal_param_vals(fixed_params==0);
    max_ls(param,mle_idx)=max_l;
    for i=mle_idx+1:numpts+1
        fprintf('Optimizing for %s=%.3f\n',param_names{param},param_vals(param,i));
        initial(fixed_params==0)=minimizers{param,i-1};
        initial(param)=param_vals(param,i);
        [minimizer,~,max_ls(param,i)] = optimize_likelihood_general(@richards_ctrl_sq_err,initial,fixed_params,noisy_data,numeric_params,C0,opt);
        minimizers{param,i}=minimizer;
    end
    for i=mle_idx-1:-1:1
        fprintf('Optimizing for %s=%.3f\n',param_names{param},param_vals(param,i));
        initial(fixed_params==0)=minimizers{param,i+1};
        initial(param)=param_vals(param,i);
        [minimizer,~,max_ls(param,i)] = optimize_likelihood_general(@richards_ctrl_sq_err,initial,fixed_params,noisy_data,numeric_params,C0,opt);
        minimizers{param,i}=minimizer;
    end
end

%% plot univariate profile likelihood (this recreates first 3 columns of Fig.1)
true_params=params;
fig=figure('Position',[100 100 1200 400],'color','w');
tl=tiledlayout(1,num_free_params);
sgtitle(sprintf('$u_{max}=%.0f,\\tau_0=%.0f,\\tau=%.0f$',u_K0,tau0,tau),'interpreter','latex','fontSize',30);
free_param_count=0;
zs = cell(num_params,1);
conf_interval=nan(num_params,1);
for param=1:num_params
    if fixed(param)
        continue;
    end
    free_param_count = free_param_count+1;
    %subplot(1,num_free_params,free_param_count);
    nexttile;
    hold on;
    xx=param_vals(param,:);
    yy=max_ls(param,:)-max(max_ls(param,:));
    plot(xx,yy);
    plot([0,1e10],[-1.92,-1.92]);
    plot([true_params(param),true_params(param)],[-10,10],'--r'); % true max
    plot([optimal_param_vals(param),optimal_param_vals(param)],[-10,10],'--g'); % MLE
    xlabel(param_names{param},'Interpreter','latex');
    ylabel('$l$','Interpreter','latex');
    axis('square');
    xlim([min(param_vals(param,:)),max(param_vals(param,:))]);
    ylim([-2.5,0]);
    if param == 1 || param == 2
        xtickformat('%.1f');
    end
    ytickformat('%.1f');
    hold off;

    zs{param}=interp_zero(xx,yy+1.92);
    if size(zs{param},2) == 2
        conf_interval(param)=zs{param}(2)-zs{param}(1);
        fprintf('95%% Confidence interval for param %s is: (intercept at -1.92)\n',param_names{param});
        fprintf('width=%.4f: [%.4f,%.4f]\n',conf_interval(param),zs{param}(1),zs{param}(2));
    elseif size(zs{param},2) == 1 && strcmp(param_names{param},'n')
        conf_interval(param)=zs{param}(1);
        fprintf('95%% Confidence interval for param %s is: (intercept at -1.92)\n',param_names{param});
        fprintf('width=%.4f: [0,%.4f]\n',conf_interval(param),zs{param}(1));
    elseif size(zs{param},2) == 1 && strcmp(param_names{param},'gamma')
        conf_interval(param)=Inf;
        fprintf('95%% Confidence interval for param %s is: (intercept at -1.92)\n',param_names{param});
        fprintf('width=Inf: [%.4f,Inf]\n',zs{param}(1));
    else
        fprintf('Do not have 2 intercepts for param %s, they are:\n',param_names{param});
        disp(zs{param});
    end
end

%% bivariate profile likelihoods 

% choose which 2 params to loop over, valid choices are 1/2, 1/4, 2/4
param1=2; param2=4;

fixed=[0,0,1,0];
fixed_param_val=params;
numeric_params={T,nt,u_d,u_K,u_r};
lb_opt=[0,0,0,0];
ub_opt=[5,5,9,50000];
opt.lb=lb_opt;
opt.ub=ub_opt;
opt.logging=true;
opt.alg=1;

numpts=201;
lb=[0, 0, 1, 0]; 
ub=[1, 1, 1, 10000];
param_names={'$r$','$\delta$','$\gamma$','$K$'};
fixed_params=fixed;
fixed_params(param1)=1;
fixed_params(param2)=1;
num_free_params=sum(1-fixed);
p1s=linspace(lb(param1),ub(param1),numpts);
p2s=linspace(lb(param2),ub(param2),numpts);
ls=zeros(numpts,numpts);
minimizers=cell(numpts,numpts);

for i=1:numpts
    for j=1:numpts
        initial=fixed_param_val;
        initial(param1)=p1s(i);
        initial(param2)=p2s(j);
        shortcut=0;
        % if I know the likelihod is very small, skip the calculation
        if param1==1 && param2==2
            if abs(initial(param2)-(initial(param1)-0.3))>0.03
                shortcut=1;
            end
        end
        if param1==1 && param2==4
            if abs(initial(param2)-8666.66*initial(param1))>500
                shortcut=1;
            end
        end
        if param1==2 && param2==4
            if abs(initial(param2)-(8666.66*initial(param1)+2600))>500
                shortcut=1;
            end
        end
        if shortcut
            fprintf('shortcut for %s=%.3f,%s=%.3f\n',param_names{param1},initial(param1),param_names{param2},initial(param2));
            ls(i,j)=-Inf;
        else
            fprintf('Optimizing for %s=%.3f,%s=%.3f\n',param_names{param1},initial(param1),param_names{param2},initial(param2));
            if num_free_params==0
                % not doing this case
            else
                % optimize over parameters other than r and D
                [minimizer,~,ls(i,j)] = optimize_likelihood_general(@richards_ctrl_sq_err,initial,fixed_params,noisy_data,numeric_params,C0,opt);
                initial(fixed==0)=minimizer;
                minimizers{i,j}=initial;
            end
        end
    end
end

%% plot bivariate profile likelihoods (this recreates columns 4-6 in Fig.1)

fig=figure;
hold on;
imagesc(ls'-max(ls,[],'all'),[-20,0]); % need transpose + reverse axis to make it right
colorbar;
set(gca,'YDir','normal');
xlabel(param_names{param1},Interpreter="latex");
ylabel(param_names{param2},Interpreter="latex");

[~,I] = max(ls',[],'all','linear');
[ix, iy] = ind2sub(size(ls'),I);
plot(iy,ix,'g*','MarkerSize',20);
plot(params(param1),params(param2),'r*','MarkerSize',20);