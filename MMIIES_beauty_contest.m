%This code solves for and graphs the beauty contest indeterminacy regions as a component of MMIIES

% Date: 12/19/2024
% Contact: adamsjonathanj@gmail.com

% Dependencies: 
% beauty_contest_discriminant.m

%Baseline parameterization:
tau_u = .3;
tau_v = .15;
varphi = .8;
alpha = .98;

% Define the range of x- and y-coordinates.
alpha_min = 0;
alpha_max = 1;
varphi_min = 0;
varphi_max = 1;

% Create a grid of x- and y-coordinates.
Ngrid = 400;
alphvec = linspace(alpha_min, alpha_max, Ngrid);
varphivec = linspace(varphi_min, varphi_max, Ngrid);
discrim = NaN(Ngrid);
irfnorm = NaN(Ngrid);


% Evaluate the function at each point in the grid.
for ja = 1:Ngrid
    for jv = 1:Ngrid
        discrim(ja,jv) = beauty_contest_discriminant(tau_u,tau_v, alphvec(ja), varphivec(jv));
    end
end


%Toeplitz components for calculating IFR:
NN=10;
vC = 4 *ones(NN, 1); vC(1) = 2; vC(2) = 3;
vR = 4 *ones(NN, 1); vR(1) = 2; vR(2) = 1;

%Calculate the IFR norms for alpha/varphi pairs:
for ja = 1:Ngrid
    for jv = 1:Ngrid
        irf_oper_blocks = {[0 alphvec(ja); 0 0] , diag([varphivec(jv) alphvec(ja)]), [0 0; varphivec(jv) 0], zeros(2)}; %blocks
        irf_oper = cell2mat(irf_oper_blocks(toeplitz(vC, vR)));
        irfnorm(ja,jv) = norm(irf_oper);
    end
end

%Create a tau grid
Mgrid = 40;
tauuvec = linspace(0, 2, Mgrid);
tauvvec = linspace(0, 2, Mgrid);
minnorm = NaN(Mgrid);
maxdiscnorm1 = NaN(Mgrid);
for iu = 1:Mgrid
    for iv = 1:Mgrid
        for ja = 1:Ngrid
            for jv = 1:Ngrid
                discrim_tau(ja,jv) = beauty_contest_discriminant(tauuvec(iu),tauvvec(iv), alphvec(ja), varphivec(jv));
            end
        end
        minnorm(iu,iv)=min(min((discrim_tau>0).*irfnorm+(discrim_tau<=0)*999));
        maxdiscnorm1(iu,iv)=max(max((irfnorm<=1).*discrim_tau-(irfnorm>1)*999));
    end
end



%Create the discriminant contours plot
fig5 = figure(5);
contour(tauuvec,tauvvec,maxdiscnorm1,'ShowText','on','LineWidth',2)
xlabel('\tau_U','FontSize',12,'FontName', 'AvantGarde');
ylabel('\tau_V','FontSize',12,'FontName', 'AvantGarde');
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XMinorTick'  , 'off'      , ...
  'LineWidth'   , 1         )
if plot_saving==1
saveas(gcf,'graphs/beauty_contest_discriminant.png')
end


%Create the indeterminacy region plot
fig6 = figure(6);
hold on
contourf(alphvec,varphivec,discrim',[0 0],'LineWidth',2)
contourf(alphvec,varphivec,-irfnorm',[-1 -1],'LineWidth',2) 
text(.1,.25,'IFR Satisfied','FontSize',20,'FontName', 'AvantGarde','Color','yellow');
text(.57,.83,'Multiplicity \rightarrow','FontSize',20,'FontName', 'AvantGarde','Color','blue');
hold off
xlabel('\alpha','FontSize',12,'FontName', 'AvantGarde');
ylabel({'$\varphi$'},'Interpreter','latex','FontSize',12);
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XMinorTick'  , 'off'      , ...
  'LineWidth'   , 1         )
if plot_saving==1
saveas(gcf,'graphs/beauty_contest_indeterminacy.png')
end


