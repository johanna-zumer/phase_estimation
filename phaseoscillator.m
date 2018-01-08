function phaseoscillator
clc;

fprintf('The Phase Oscillator\n');

par=param(1);
ini=[0.1;0.1];

trans=0; tend=100; tstep=0.1;

integration(par,ini,trans,tend,tstep);

% ==========  Integration =================================

function integration(par,x0,trans,tend,tstep)

[t,x] = run(x0,par,trans,tend,tstep,[],[]);

figure(1)
clf;

subplot(1,2,1)
plot(t,x);

subplot(1,2,2)
plot(x(:,1),x(:,2));


% ==========  dxdt =================================    

function output = dxdt(t,x,par,forcing,pulse)

output = myequa(x,par);

% ==========  Equations =================================    

function output = myequa(x,par)

g = par(1);
k = par(2);
a = par(3);
c = par(4);

r2=x(1)^2+x(2)^2;

output = [-k*x(2)+g*(x(1)-c)*(a^2-r2);
          g*x(2)*(a^2-r2)+k*(x(1)-c)];

% ==========  Parameters ================================= 

function par = param(set)
g=1; k=1; a=1; c=0;
par=[g,k,a,c];

% ==========  Run ===================================== 

function [t,x]=run(x0,par,trans,tend,tstep,forcing,pulse)

ttrans = [0:tstep:trans];
tspan = [0:tstep:tend];

option = odeset('RelTol', 1e-5);

if ttrans(end)>0 
    sol = ode45(@dxdt,ttrans,x0,option,par,forcing,pulse);
    x0 = sol.y(:,end);
end

[t x] = ode45(@dxdt,tspan,x0,option,par,forcing,pulse);