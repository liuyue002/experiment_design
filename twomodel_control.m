% This script demonstrate how to use optimal control to maximize the
% difference between two models.
% This accompanies the paper "Optimal experiment design for practical 
% parameter identifiability and model discrimination" by 
% Liu, Baker, Maini, 2025, and illustrates the steps to create Fig.5. The
% steps for Fig.6-11 are very similar.

%% specify the two models
f = @(C,r,d,gamma,K,uk) r*C.*(1-(C./(K-uk)).^gamma)-d*C;
dfdC= @(C,r,d,gamma,K,uk) r*(1-(1+gamma)*(C./(K-uk)).^gamma)-d;

r1=0.45;
d1=0.15;
gamma1=1;
K1=3900;

r2=0.30;
d2=0;
gamma2=1;
K2=2600;

C0=100;
T=25;
upts=100;
odeopts = odeset('RelTol',1e-4,'AbsTol',1e-4,'MaxStep',T/upts);
alpha=0.03; % weight of control cost
omega=0.1; % control update rate
uklim = 1200; %upper bound for uk
bangbang=false; % whether using continuous or bang-bang control

if bangbang
    H = @(C1,C2,lambda1,lambda2,u) (C1-C2).^2-alpha*u + lambda1*f(C1,r1,d1,gamma1,K1,u) + lambda2*f(C2,r2,d2,gamma2,K2,u);
else
    H = @(C1,C2,lambda1,lambda2,u) (C1-C2).^2-alpha*u.^2 + lambda1*f(C1,r1,d1,gamma1,K1,u) + lambda2*f(C2,r2,d2,gamma2,K2,u);
end
fminconopts=optimoptions(@fmincon,'Algorithm','interior-point','Display','none');
optiminits=[0,uklim/3,uklim*2/3,uklim];
ukvar=optimvar("uk",LowerBound=0,UpperBound=uklim);
ms=MultiStart("Display","off");

%% figure and initializations

uk=@(t) 100; % initial guess
tfine=linspace(0,T,upts);
uknum=arrayfun(uk,tfine);
Lambinit=[0;0];
Jold=nan;
maxiter=1000;
Js=zeros(maxiter,1);
Jmodels=zeros(maxiter,1);
Jcontrols=zeros(maxiter,1);

fig=figure('visible','on','Position',[50,200,1700,600]);
tiles=tiledlayout(1,3);
ax1=nexttile;
hold on
c1plot=plot(ax1,tfine,zeros(size(tfine)));
c2plot=plot(ax1,tfine,zeros(size(tfine)));
cdiffplot=plot(ax1,tfine,zeros(size(tfine)));
uplot=plot(ax1,tfine,zeros(size(tfine)));
hold off
xlabel('t');
legend('C_1','C_2','|C_1-C_2|','u');
xlim([0,T]);
ylim([0,2600]);
xtitle=title('Forward');
ax2=nexttile;
hold on
l1plot=plot(ax2,tfine,zeros(size(tfine)));
l2plot=plot(ax2,tfine,zeros(size(tfine)));
hold off
xlabel('t');
legend('\lambda_1','\lambda_2');
xlim([0,T]);
ltitle=title('Backward');
ax3=nexttile;
uknewplot=plot(ax3,tfine,zeros(size(tfine)));
xlabel('t');
ylabel('uknew');
xlim([0,T]);
ylim([0,uklim+100]);
title('u update');
figtitle=sgtitle('title');
tiles.Padding="tight";
tiles.TileSpacing="tight";

%% forward-backward sweep


for iter=1:maxiter
    [Js(iter),Jcontrols(iter),Jmodels(iter),C,Lambda]=J_uk(C0,T,tfine,uknum,r1,d1,gamma1,K1,r2,d2,gamma2,K2,alpha,bangbang);

    if iter==1
        Jgain = nan;
    else
        Jgain = -(Js(iter)-Js(iter-1));
    end
    fprintf('Iteration %d: Jcontrol=%.5f, Jmodel=%.5f, J=%.5f, gain = %.5f\n',iter,Jcontrols(iter),Jmodels(iter),Js(iter),Jgain);

    if Jgain<0
        omega=omega/2;
        fprintf('reduced omega to %f\n',omega);
    end

    % update control
    uknumnew=zeros(size(uknum));
    for i=1:upts
        objective=@(u) H(C(i,1),C(i,2),Lambda(i,1),Lambda(i,2),u);
        if bangbang
            if objective(uklim)>objective(0)
                uknumnew(i)=uklim;
            else
                uknumnew(i)=0;
            end
        else
            prob = optimproblem(Objective=objective(ukvar), ObjectiveSense='maximize');
            sol = solve(prob, optimvalues(prob,uk=optiminits), ms,Options=fminconopts);
            uknumnew(i)=sol.uk;
        end
    end

    c1plot.YData=C(:,1);
    c2plot.YData=C(:,2);
    cdiffplot.YData=abs(C(:,1)-C(:,2));
    uplot.YData=uknum;
    l1plot.YData=Lambda(:,1);
    l2plot.YData=Lambda(:,2);
    uknewplot.YData=uknumnew;
    figtitle.String=['Iteration=',num2str(iter),', Jcontrol=',num2str(Jcontrols(iter),'%.02f'),', Jmodel=',num2str(Jmodels(iter),'%.02f'),', Jtotal=',num2str(Js(iter),'%.02f'),];
    drawnow;


    if abs(Jgain)<0.01
        fprintf('Converged.\n');
        break;
    end

    uknum= uknum*(1-omega) +uknumnew*omega;
    uk=@(tt) interp1(tfine,uknum,tt)';
end

%% history of J
figJ = figure;
hold on;
plot(Js(1:iter));
xlabel('iter');
ylabel('J');

