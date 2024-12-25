%This code solves the Singleton model using the z-tran toolbox

% Date: 7/3/2023
% Contact: adamsjonathan@ufl.edu

% Dependencies: z-tran toolbox Han, Tan & Wu (2022): Analytic Policy Function Iteration.

clear
close all
clc

% Parameters - baseline calibration in Nimark (2017)
beta = 0.95;
rho_theta =.9;
sigma_u = .05;
sigma_eps = 1;
sigma_eta = .1;

% Variables
xp = 1;

% Index signal dimensions (Han et al refer to these as "shocks")
sd = 1; %"d" refers to theta + epsilon in Nimark's notation
sz = 2; %the exogenous private signal z = theta + eta_i

% Index shock dimensions (Han et al refer to these as "innovations")
eu = 1; %u shocks to the AR(1) theta
eeps = 2; %aggregate transient shocks
eeta = 3; %idiosyncratic noise shocks

% Dimensions
nx = 1; %one endogneous variable (p)
ns = 2; %two signals
ne = 3; %three fundamental shocks

% Model
Ax0 = zeros(nx);
As0 = zeros(nx,ns);
Bx0 = zeros(nx);
Bx1 = zeros(nx);
C1 = zeros(ns);
D0 = zeros(ns,ne);
D1 = zeros(ns,ne);
V = zeros(ne);

Ax0(1,xp) = -1;
As0(1,sd) = -1;
Bx1(1,xp) = beta;

%exogenous process for theta + epsilon:
C1(sd,sd) = rho_theta;
D0(sd,eu) = 1;
D0(sd,eeps) = 1;
D1(sd,eeps) = -rho_theta;
%exogenous process for theta + eta_i:
C1(sz,sz) = rho_theta;
D0(sz,eu) = 1;
D0(sz,eeta) = 1;
D1(sz,eeta) = -rho_theta;


V(eeps,eeps) = sigma_eps^2;
V(eu,eu) = sigma_u^2;
V(eeta,eeta) = sigma_eta^2;

agg = {xp,[eu eeps]};     % {var_id, inn_id}
sig = {
    1, xp, sz, true     % {eqn_id, end_id, exo_id, ave_ex}
    };

% Solve model
tic
m = ztran(nx,ns);
m.Ax = {Ax0};
m.As = {As0};
m.Bx = {Bx0,Bx1};
m.C = {C1};
m.D = {D0,D1};
m.V = V;
m.agg = agg;
m.sig = sig;
m = solve(m,'nit',{10,500,10});%Solve with default settings
ztran_singleton_time = toc

% IRF
PlotT = 40;
FullT = 300;
var_lab = {'$p$'};
inn_lab = {'$u$','$\epsilon$','$\eta_{i}$'};

%Save each IRF:
zirfs=NaN(3,FullT);

figure('Name','IRF')
for i = 1:ne
    imp = zeros(ne,FullT);
    imp(i,1) = sqrt(V(i,i)); %This plots the impulse to a one s.d. shock
    res = irf(m.sol,imp);
    zirfs(i,:)=res;
    if ismember(i,setdiff(1:ne,agg{2}))
        res(agg{1},:) = 0;
    end
    for j = 1:nx
        subplot(ne,nx,(i-1)*nx+j)
        plot(1:PlotT,res(j,1:PlotT),'linewidth',1.2)
        xlim([1 PlotT])
        xticks([1 round(PlotT/2) PlotT])
        if i==1
            title(var_lab{j},'Interpreter','latex','FontSize',12)
        end
        if j==1
            ylabel(inn_lab{i},'Interpreter','latex','FontSize',12)
        end
    end
end

save('ztran_singleton.mat','zirfs');
save('ztran_singleton_time','ztran_singleton_time');
