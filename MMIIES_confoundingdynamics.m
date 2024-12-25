%This code solves the Confounding Dynamics Asset Pricing model as a component of MMIIES

% Date: 12/19/2024
% Contact: adamsjonathanj@gmail.com

% Dependencies: SIGNAL_OP Package of MATLAB functions, Version 1.2+ (add to path), downloadable from jonathanjadams.com
% cdi_discriminant


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Set Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

calculate_full_info = 0; %set to 1 to calculate full information solution too

Delta = .0001; %perturbation size

N = 200; %truncation length
vecT = 200; %vector length for convolution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Calculate Feedback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha_vec = linspace(1,2,N);
beta_vec = linspace(.001,2,N);

[alpha_grid, beta_grid] = meshgrid(alpha_vec,beta_vec);

normgrid = NaN(N);

Frechetgrid = NaN(N);


for aa=1:N
    for bb = 1:N
        alpha = alpha_grid(aa,bb);
        beta = beta_grid(aa,bb);
        theta = 1/alpha; %B_CD parameter
        xi = (Delta/beta + theta); %B_Delta parameter
        B_CD_denom = (-theta).^[0:(vecT-1)];
        B_CD_numer = [theta 1];
        B_CD = conv(B_CD_numer,B_CD_denom);  %check: B_CD(1:end-1)*B_CD(2:end)' should be zero, norm(B_CD) should be one
        B_Delta_denom = (-xi).^[0:(vecT-1)];
        B_Delta_numer = [xi 1];
        B_Delta = conv(B_Delta_numer,B_Delta_denom);  %check: B_Delta(1:end-1)*B_Delta(2:end)' should be zero, norm(B_Delta) should be one
        
        Dp = (Delta*alpha + beta)*B_Delta - beta*B_CD;
        Frechetgrid(aa,bb) = norm(Dp)/Delta;
    end
end

%choose contour lines:
lines = [2.415 2.45 2.5 2.6:.2:3 3:.5:4 4 5 10]; 

close all
fig1 = figure(1);
contour(alpha_grid,beta_grid,Frechetgrid,lines,'ShowText','on')
xlabel('\alpha','FontSize',12,'FontName', 'AvantGarde');
ylabel('\beta','FontSize',12,'FontName', 'AvantGarde');
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XMinorTick'  , 'off'      , ...
  'LineWidth'   , 1         )
if plot_saving==1
saveas(gcf,'graphs/cd_model_instability.png')
end

%numerical minimum:
min(min(Frechetgrid))



%% Solve the confounding dynamics equilibrium

T = 100;


%Set dimensions (most of the code is written for general dimensions)
nc=1; %number of control dimensions
ns=0; %number of state dimensions
control_vars = ([1]==1)'; %identifies which entries in choice vector X are controls (versus states)
n=nc+ns; %number of equilibrium equations
m_z=2; %observed signal length (p,z)
m_eps=2; %fundamental shock length [theta,epsilon,eta] = (aggr fundamental, aggr noise, idio noise)


%define model parameters
beta = .99;
alpha = 2;

%Set aggregation matrix
aggr_shocks = [1 0]; %indicates aggregate shocks
P_G=diag(aggr_shocks); 
%Set shock variance matrix:
Sigma =  diag([1 1]).^2; 

%Set matrices of economic model 
%Endogenous Vector Coefficients: (X = [forecast]')
BX0=-1;  
BX1=0;

%Innovation Coefficients: (Z=[z, p]')
BA0=[0 0]; 
BA1=[beta 0];

%Set matrices of information feedback.  (X = p -> Z=[z, p]')
G_0=[0; 1];  %the coefficient matrix for current aggregate variables

%Construct the polynomial G
G_sig = zeros(size(G_0,1),size(G_0,2),T); %If A were noncausal, this polynomial would have to be 2T-1
G_sig(:,:,1)=G_0;

%Set matrices of determining the exogenous signal component
%2 shocks:
S_X0 = [ 1 1; 0 0];
S_X1 = [alpha 0; 0 0];
   
%Calculate S_X as the component that depends on exogenous variables
S_X_sig = zeros(m_z,m_eps,T);
S_X_sig(:,:,1)=S_X0;
S_X_sig(:,:,2)=S_X1;

[X_shock_sig, X_sig_full, IFRnorm, S_sig] = soi_solve(BX0,BX1,BA0,BA1,nc,P_G,G_sig,S_X_sig,Sigma,calculate_full_info);


% Next, calculate divergence from a perturbation of the CD equilibrium
% Find a perturbed signal:

Delta = .01; %controls perturbation size

perturb = Delta * (alpha^-1).^[0:1:(T-1)];
perturb=permute(perturb,[1 3 2]);

S_CD = S_sig;

S_Delta = S_sig;
S_Delta(2,1,:)=S_sig(2,1,:)+perturb;

%Also find the full information equilibrium by starting with FIRE as the
%initial guess:
calculate_full_info=1;
[X_shock_sig, X_sig_full, IFRnorm, S_sig] = soi_solve(BX0,BX1,BA0,BA1,nc,P_G,G_sig,S_X_sig,Sigma,calculate_full_info);

S_FIRE = S_sig;

% In the following, I have unpacked the soi_solve code to start from the perturbed guess,
% and to track IRFs from each iteration:
S_array = [];

%Initial guess is Confounding Dynamics Perturbed!
S_sig = S_Delta;

%dimension variables
ns = size(BX0,2)-nc;
T = size(S_X_sig,3); %truncation length


%set computational parameters:
tolerance=1e-6;  
error_upperbound = 1e5; %set an upper bound to exit loop if algorithm diverges
maxiter=500;


%Construct the Theta, Xi operators from the model ingredients:
[Theta, Xi] = theta_xi(BX0,BX1,nc,T);

%Is Information Feedback Regularity satisfied?
IFRnorm = ifrnorm(Theta,Xi,BA0,BA1,nc,G_sig,T);

%Solve
%reset error, iteration counter:
error=100; iter=1;

while error>tolerance && iter<=maxiter && error<error_upperbound
%Steps 1 and 2: Find the Wold Representation
[A_sig, W_sig] = wold(S_sig,T,Sigma);

%Step 3: define A tilde= [(BA1 L^-1 +BA0)A]_+
    QA_tilde_sig =  smulti(BA0,A_sig,T) + smulti(BA1,lagmulti(A_sig,-1),T);

    %Step 4: Find the policy function
    XiA_noncausal = smulti(Xi,QA_tilde_sig,T);
    XiA_causal = XiA_noncausal(:,:,T:end);
    X_sig = smulti(Theta,XiA_causal,T);

%Step 5: Calculate the next signal polynomial

if size(G_sig,3)==T %FOR CAUSAL G (standard)
    S_next_sig = S_X_sig + smulti(G_sig,smulti(X_sig,smulti(W_sig,P_G,T),T),T);
end
if size(G_sig,3)==2*T-1 %FOR NONCAUSAL G 
    GX_sig = smulti(G_sig,X_sig,T);
    S_next_sig = S_X_sig + smulti(GX_sig(:,:,T:2*T-1),smulti(W_sig,P_G,T),T); 
end

%Calculate difference:
error=norm(s2t(S_next_sig,T)-s2t(S_sig,T));
display(strcat('DI iteration ',sprintf(' %g ',iter),': error ',sprintf(' %g ',error)))

    if iter==1
        initialerror = error;
    end

    
%log iteration into array
S_array(:,:,:,iter)=S_sig;   

%update:
S_sig=S_next_sig;
iter=iter+1;
end
iter=iter-1; %role back one iteration if the loop ended


%Now plot the results:

plotT=10;
price_array = S_array(2,1,1:plotT,:);
iterations = [1 20 30 33 35:1:40]; %only plot some iterations (first 20 diverge slowly)

fig2 = figure(2);
hold on
plot(1:plotT,squeeze(price_array(:,:,:,1)),':','Color',	[0 1 0],'LineWidth',3)
plot(1:plotT,squeeze(price_array(:,:,:,end)),'--','Color',	[0 0 0],'LineWidth',2)
for ii = 1:length(iterations)
    jj = iterations(ii);
    reldistance = norm(squeeze(price_array(:,:,:,jj)-price_array(:,:,:,end)))/norm(squeeze(price_array(:,:,:,1)-price_array(:,:,:,end)));
    plot(1:plotT,squeeze(price_array(:,:,:,jj)),'Color',	[1-reldistance 0 reldistance],'LineWidth',1)
end
plot(1:plotT,squeeze(price_array(:,:,:,end)),'--','Color',	[0 0 0],'LineWidth',2)
plot(1:plotT,squeeze(price_array(:,:,:,1)),':','Color',	[0 1 0],'LineWidth',3)
hold off
l=legend('Confounding Dynamics (perturbed)','Full Information (final iteration)','Location','NE');
set(l,'FontSize',20);
xlabel('Shock Lag','FontSize',12,'FontName', 'AvantGarde');
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XMinorTick'  , 'on'      , ...
  'LineWidth'   , 1         )

if plot_saving==1
saveas(gcf,'graphs/cd_model_divergence.png')
end

%% Modify the model to feature idiosyncratic noise

%Baseline parameterization:
alpha = 2;
theta = 1/alpha;
beta = .1;
tau_idio = (.25)^2;
tau_y = tau_idio;
tau_v = tau_idio;

%Define function for calculating SIC
SIC_rhs_cdi = @(alpha,theta,beta,tau_y) (1+beta)*(1-beta)^3/(4*beta*(sqrt(1+alpha^2) + theta/sqrt(tau_y))^2);


% vary the tau_v parameter:
tausqrtvec = linspace(.1, 5, 1000);
tau_v_vec = tausqrtvec.^2;
ROOTS_vec = zeros(3,length(tau_v_vec));
multi_vec = ROOTS_vec==0;
SIC_vec = tau_v_vec==0;
for tt = 1:length(tau_v_vec)
    tau_v_tt = tau_v_vec(tt);
    %calculate solutions to the cubic equation:
    ROOTS_raw = roots([tau_v_tt, 0, (tau_y + 1 - beta*tau_v_tt), -beta*tau_y]);
    [realsorted realindex] = sort(real(ROOTS_raw),'descend'); %order by real component
    ROOTS_vec(:,tt) = ROOTS_raw(realindex); %store roots
    for rrrr = 1:3
        multi_vec(rrrr,tt) = isreal(ROOTS_vec(rrrr,tt));
    end
    %is SIC satisfied?:
    SIC_vec(tt) = SIC_rhs_cdi(alpha,theta,beta,tau_y) > max(tau_y,tau_v_tt); 
end

%find the min and max tau_v such that SIC holds:
tausqrtSIC_min = 0; %is necessarily zero
tausqrtSIC_max = max(tausqrtvec(SIC_vec==1));



%make the bifurcation diagram
%close all
fig3 = figure(3);
hold on
plot(tausqrtvec(multi_vec(1,:)), ROOTS_vec(1,multi_vec(1,:)))
plot(tausqrtvec(multi_vec(2,:)), ROOTS_vec(2,multi_vec(2,:)))
plot(tausqrtvec(multi_vec(3,:)), ROOTS_vec(3,multi_vec(3,:)))
taulim = xlim;
blim = ylim;
hold off

close(fig3)
fig3 = figure(3);
hold on
patch([tausqrtSIC_min tausqrtSIC_min tausqrtSIC_max tausqrtSIC_max], [min(blim) max(blim) max(blim) min(blim)], [0.8 0.8 0.8],'EdgeColor','none')
plot(tausqrtvec(multi_vec(1,:)), ROOTS_vec(1,multi_vec(1,:)),'LineWidth',2)
plot(tausqrtvec(multi_vec(2,:)), ROOTS_vec(2,multi_vec(2,:)),'LineWidth',2)
plot(tausqrtvec(multi_vec(3,:)), ROOTS_vec(3,multi_vec(3,:)),'LineWidth',2)
hold off
legend('SIC Satisfied','Solution 1','Solution 2','Solution 3','Location','SW')
ylabel('Solution Coefficient','FontSize',12,'FontName', 'AvantGarde');
xlabel({'$\sqrt{\tau_v}$'},'Interpreter','latex','FontSize',12);
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XMinorTick'  , 'off'      , ...
  'LineWidth'   , 1         )
if plot_saving==1
saveas(gcf,'graphs/cdi_bifurcation.png')
end

%Next, vary both beta and tau_v


% Define the range of x- and y-coordinates.
beta_min = 0;
beta_max = 1;
tau_v_sqrt_min = 0;
tau_v_sqrt_max = 3;

% Create a grid of x- and y-coordinates.
Ngrid = 200;
betavec = linspace(beta_min, beta_max, Ngrid);
tausqrtvec = linspace(tau_v_sqrt_min, tau_v_sqrt_max, Ngrid);
discrim = NaN(Ngrid);
sicdiff = NaN(Ngrid);


% Evaluate the function at each point in the grid.
for jbeta = 1:Ngrid
    for jtau = 1:Ngrid
        discrim(jbeta,jtau) = cdi_discriminant(tau_y,tausqrtvec(jtau)^2,betavec(jbeta));
    end
end

%Calculate the SIC difference for beta/tau_v pairs:
for jbeta = 1:Ngrid
    for jtau = 1:Ngrid
        sicdiff(jbeta,jtau) = SIC_rhs_cdi(alpha,theta,betavec(jbeta),tau_y) - max(tau_y,tausqrtvec(jtau)^2);
    end
end

%Create the indeterminacy region plot
fig4 = figure(4);
hold on
contourf(tausqrtvec,betavec,discrim,[.001 .001],'LineWidth',2) %note: tau is x-axis, and COLUMNS of discrim
contourf(tausqrtvec,betavec,sicdiff,[0 0],'LineWidth',2) %note: tau is x-axis, and COLUMNS of sicdiff
text(.4,.11,'\leftarrow SIC Satisfied','FontSize',20,'FontName', 'AvantGarde','Color','blue');
text(2,.7,'Multiplicity','FontSize',20,'FontName', 'AvantGarde','Color','blue');
hold off
ylabel('\beta','FontSize',12,'FontName', 'AvantGarde');
xlabel({'$\sqrt{\tau_v}$'},'Interpreter','latex','FontSize',12);
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XMinorTick'  , 'off'      , ...
  'LineWidth'   , 1         )
if plot_saving==1
saveas(gcf,'graphs/cdi_multiplicity.png')
end

