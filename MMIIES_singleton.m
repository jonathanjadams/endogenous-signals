%This code solves the Singleton model as a component of MMIIES

% Date: 12/19/2024
% Contact: adamsjonathanj@gmail.com

% Dependencies: SIGNAL_OP Package of MATLAB functions, Version 1.2+ (add to path), downloadable from jonathanjadams.com
% Matrix of results from the Nimark algorithm: 'nimark_pirfs.mat'
% Matrices of results from Han et al (2023) algorithm: 'ztran_singleton.mat', 'ztran_singleton_time.mat'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Set Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tau_space = [(10:10:100)'; 300]; %vector of approximation lengths
Tau_times = NaN*Tau_space;
Tau_IRFs = zeros(3,max(Tau_space),length(Tau_space));


calculate_full_info = 0; %set to 1 to calculate full information solution too

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Define the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Set dimensions (most of the code is written for general dimensions)
nc=1; %number of control dimensions
ns=0; %number of state dimensions
control_vars = ([1]==1)'; %identifies which entries in choice vector X are controls (versus states)
n=nc+ns; %number of equilibrium equations
m_z=2; %observed signal length (p,z)
m_eps=3; %fundamental shock length [theta,epsilon,eta] = (aggr fundamental, aggr noise, idio noise)

%labels:
vartitles={'f'}; %endogenous variables in X
shocktitles={'theta','epsilon','eta'}; %Shock names


%define model parameters
beta=.95;

%shock persistences
rho_theta=.9;

%Shock vector [theta,eta,epsilon]
%Shock AR coefficient:
RHO=diag([rho_theta 0 0]);
%Set aggregation matrix
aggr_shocks = [1 1 0]; %indicates aggregate shocks
P_G=diag(aggr_shocks); 
%Set shock variance matrix:
sigma_u = .05;
sigma_eps = 1;
sigma_eta = .1;
Sigma =  diag([sigma_u sigma_eps sigma_eta]).^2; 



%Set matrices of economic model 
%Endogenous Vector Coefficients: (X = [forecast]')
BX0=-1;  
BX1=0;

%Innovation Coefficients: (Z=[p, z]')
BA0=[0 0]; 
BA1=[1 0];

for ttt = 1:length(Tau_space)
T = Tau_space(ttt); %truncation order


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Define exogenous polynomials S_X (n x m) and A (n x n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set matrices of information feedback.  (X = f -> Z=[p, z]')
G_0=[beta; 0];  %the coefficient matrix for current aggregate variables

%Construct the polynomial G
G_sig = zeros(size(G_0,1),size(G_0,2),T); %If A were noncausal, this polynomial would have to be 2T-1
G_sig(:,:,1)=G_0;

%Set matrices of determining the exogenous signal component
%3 shocks:
S_X0 = [-1 -1 0; 1 0 1];
S_X1 = [0 0 0; 0 0 0];


    
%Define the process for exogenous random (but not i.i.d) variables
Shocks_sig = NaN(m_eps,m_eps,T);
for j = 1:T
Shocks_sig(:,:,j)=RHO^(j-1);
end

%Calculate S_X as the component that depends on exogenous variables, times the process for those variables
S_X_noshocks_sig = zeros(m_z,m_eps,T);
S_X_noshocks_sig(:,:,1)=S_X0;
S_X_noshocks_sig(:,:,2)=S_X1;
S_X_sig = smulti(S_X_noshocks_sig,Shocks_sig,T);


%%
tic
[X_shock_sig, X_sig_full, IFRnorm, S_sig] = soi_solve(BX0,BX1,BA0,BA1,nc,P_G,G_sig,S_X_sig,Sigma,calculate_full_info);
toc;
Tau_times(ttt)=toc;
Tau_IRFs(:,1:T,ttt) = squeeze(S_sig(1,:,:));

end


%% Analyze Truncation Sensitivity 

%calculate errors for price IRFs:
for ttt=1:length(Tau_space)-1
    for jj = 1:3
        Errors(ttt,jj) = norm(Tau_IRFs(jj,:,ttt)-Tau_IRFs(jj,:,end));
    end
end
Errors_Aggr = sum(Errors(:,1:2),2);
Times_plot = Tau_times(1:end-1); %/Tau_times(end-1); %if you want to normalize



%Tau_plot_indices = [1 3 5 7 9 11];
Tau_plot_indices = [1 2 3 4 5 11];
plotN = 40;

load('nimark_pirfs.mat')
Nimark_IRF1 = pirfs(1,1:plotN+1);
Nimark_IRF2 = pirfs(2,1:plotN+1);
load('ztran_singleton.mat');
Ztran_IRF1 = zirfs(1,1:plotN+1);
Ztran_IRF2 = zirfs(2,1:plotN+1);

%Calculate ztran accuracy vs tau=300 IRFs:
sigmavec = [sigma_u sigma_eps sigma_eta]; %Ztran irfs are saved wrt one s.d. impulses
for jj = 1:2
   Ztran_errors(jj) =  norm( zirfs(jj,:,end)/sigmavec(jj)- Tau_IRFs(jj,:,end) );
end
load('ztran_singleton_time.mat');

%how long does it take SOI to achieve this error?
soi_time = interp1(Errors_Aggr,Times_plot,sum(Ztran_errors),'linear','extrap');
soi_ztran_time_ratio = round(soi_time/ztran_singleton_time,1);
display(strcat('Signal Operator Iteration solves the model as accurately as z-tran in ',sprintf(' %g ',soi_ztran_time_ratio),'x time.'))


%close all
fig7 = figure(7);
hold on
plot(Tau_space(1:end-1),Errors_Aggr/max(Errors_Aggr),'LineWidth',2)
plot(Tau_space(1:end-1),Times_plot/max(Times_plot),'--','LineWidth',2)
hold off
legend('Solution Error','Computation Time','Location','E')
xlabel('Truncation Length \tau','FontSize',12,'FontName', 'AvantGarde');
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XMinorTick'  , 'off'      , ...
  'LineWidth'   , 1         )
if plot_saving==1
saveas(gcf,'graphs/nimark_soi_tau_comparisons.png')
end



fig8 = figure(8);
hold on
for ii=1:length(Tau_plot_indices)
    Tau_index = Tau_plot_indices(ii);
    txt = ['\tau = ',num2str(Tau_space(Tau_index))];
    plot(0:plotN,sigma_u*Tau_IRFs(1,1:plotN+1,Tau_index),'--','Color', (1-ii/length(Tau_plot_indices))*[1 1 1],'LineWidth',2,'DisplayName',txt)
end
%plot(0:plotN,Tau_IRFs(1,1:plotN+1,end),'Color', 'r','LineWidth',2,'DisplayName','Max \tau')
plot(0:plotN,Nimark_IRF1,'Color', 'r','LineWidth',2,'DisplayName','Nimark Solution')
plot(0:plotN,Ztran_IRF1,'b:','LineWidth',2,'DisplayName','Han et al Solution')
legend show
legend('Location','SE')
hold off
xlabel('Periods','FontSize',12,'FontName', 'AvantGarde');
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XMinorTick'  , 'on'      , ...
  'LineWidth'   , 1         )
if plot_saving==1
saveas(gcf,'graphs/nimark_soi_theta.png')
end


%% Analyze IFR sensitivity

T=50;
maxiter = 50; %we are going to measure errors after 50 iterations

beta_space = [.95:.01:1.25];
beta_times = NaN*beta_space;
beta_errors = NaN(2,length(beta_space));
beta_IRFs = zeros(3,T,length(beta_space));

%calculate for a couple persistences
rho_theta_large = .99;
rho_theta_huge = .999;
rho_vec = [rho_theta rho_theta_large];
%rho_vec = [rho_theta rho_theta_large rho_theta_huge];

for rrr = 1:length(rho_vec)
for ttt = 1:length(beta_space)
beta = beta_space(ttt); %discount factor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Define exogenous polynomials S_X (n x m) and A (n x n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set matrices of information feedback.  (X = f -> Z=[p, z]')
G_0=[beta; 0];  %the coefficient matrix for current aggregate variables

%Construct the polynomial G
G_sig = zeros(size(G_0,1),size(G_0,2),T); %If A were noncausal, this polynomial would have to be 2T-1
G_sig(:,:,1)=G_0;

%Set matrices of determining the exogenous signal component
%3 shocks:
S_X0 = [-1 -1 0; 1 0 1];
S_X1 = [0 0 0; 0 0 0];


RHO=diag([rho_vec(rrr) 0 0]);

%Define the process for exogenous random (but not i.i.d) variables
Shocks_sig = NaN(m_eps,m_eps,T);
for j = 1:T
Shocks_sig(:,:,j)=RHO^(j-1);
end

%Calculate S_X as the component that depends on exogenous variables, times the process for those variables
S_X_noshocks_sig = zeros(m_z,m_eps,T);
S_X_noshocks_sig(:,:,1)=S_X0;
S_X_noshocks_sig(:,:,2)=S_X1;
S_X_sig = smulti(S_X_noshocks_sig,Shocks_sig,T);


%
tic
[X_shock_sig, X_sig_full, IFRnorm, S_sig, A_sig, X_sig, relativeerror] = soi_solve(BX0,BX1,BA0,BA1,nc,P_G,G_sig,S_X_sig,Sigma,calculate_full_info,maxiter);
toc;
%beta_times(ttt)=toc;
%beta_IRFs(:,1:T,ttt) = squeeze(S_sig(1,:,:));
beta_errors(rrr,ttt)=relativeerror;

end
end


%plot sensitivity results

yaxislimits = [0 2];

%FIRE solution: 1/(1-beta*rho) requires beta*rho<1
firesupvec = 1./rho_vec;

fig9 = figure(9);
hold on
plot(beta_space,beta_errors(1,:),'LineWidth',2)
plot(beta_space,beta_errors(2,:),'--','LineWidth',2)
legend('\rho = 0.90','\rho = 0.99','Location','NW')
ylim(yaxislimits)
hold off
xlabel('Information Feedback (\beta)','FontSize',12,'FontName', 'AvantGarde');
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XMinorTick'  , 'on'      , ...
  'LineWidth'   , 1         )
if plot_saving==1
saveas(gcf,'graphs/nimark_ifr_errors.png')
end
