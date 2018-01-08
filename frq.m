function frq

%%% 3-variable FRQ model for the circadian clock of Neurospora
%%% Didier Gonze
%%% Created: 25/1/2008
%%% Created: 29/12/2010

clc;

%%% Verbosity (display comments):

verbo=0; 

%%% Number of variable and initial conditions:

nbvar=3; 
xini=ones(1,nbvar)/10;

%%% Time parameters:

trans=100;
tend=96;
tstep=0.1;

%%% Task:

integration(verbo,xini,trans,tend,tstep);



%====================================================================
% Integration
%====================================================================

function output=integration(v,x0,trans,tend,tstep);

if v==1
    fprintf('Let s go...\n');
    tic
end

[t,x] = run(x0,trans,tend,tstep);

if v==1
    fprintf('Integration finished!...\n');
    toc
end

set(figure(1),'Position', [400 400 500 300]);  
clf;

plot(t,x(:,1:3));           
xlabel('Time (h)','fontsize',18);
ylabel('Concentration (nM)','fontsize',18);
xlim([0 tend]);
%ylim([0 0.0005]);
set(gca,'xtick',[0:12:tend],'fontsize',14);
legend('M','P_C','P_N');
set(findobj(gca,'Type','line'),'LineWidth',2);


%%% Computation of the period

xt=x(:,1); % variable test pour calculer la periode
xp=[xt(2:end); 0]; % variable test shiftee
xmean=mean(xt);  % seuil

k=find((xt(:)<xmean) & (xp(:)>xmean));

tp=t(k);

periods=[];
for i=1:length(tp)-1
    periods=[periods; tp(i+1)-tp(i)];
end
    
period=mean(periods);      % mean of the period (peak-to-peak intervals)
stdperiod=std(periods);    % std dev of the period

fprintf('Mean(Periods) = %g \n', period);
fprintf('Std(Periods) = %g \n', stdperiod);



% ============================================================================================
% Run
% ============================================================================================

function [t,x]=run(x0,trans,tend,tstep)

ttrans = [0:tstep:trans];
tspan = [0:tstep:tend];

option = [];
%option = odeset('RelTol', 1e-5);
option=odeset('OutputS',[1:3],'OutputF','odeplot');

if trans > 0 
    [t x] = ode45(@dxdt,ttrans,x0,option);
    x0=x(end,:);
end

[t x] = ode45(@dxdt,tspan,x0,option);


% ============================================================================================
% dxdt
% ============================================================================================

function y = dxdt(t,x)


%%% parameter set 1

vs=1.6; ki=1; n=4;
vm=0.505; km=0.5;
ks=0.5;
vd=1.4; kd=0.13;
pk1=0.5; pk2=0.6;

%%% variables

m=x(1);
pc=x(2);
pn=x(3);

%%% equations

y = [
    vs*ki^n/(ki^n+pn^n)-vm*m/(km+m);    % dm/dt
    ks*m-vd*pc/(kd+pc)-pk1*pc+pk2*pn;   % dpc/dt
    pk1*pc-pk2*pn;                      % dpn/dt
] ; 
