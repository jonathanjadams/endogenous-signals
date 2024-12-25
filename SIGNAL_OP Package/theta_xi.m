function [Theta, Xi] = theta_xi(BX0,BX1,nc,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Preliminary variable assignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dimension variables
ns = size(BX0,2)-nc;
n = ns+nc;
control_vars = [zeros(ns,1);ones(nc,1)]==1; %logical vector identifying control dimensions, assumed to be ordered after states
%T = size(S_X_sig,3); %truncation length


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Generate the policy operators (Xi and Theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%USE GENERALIZED SCHUR TO REDUCE SYSTEM
[AA,BB,Q,Z] = qz(BX1,BX0); 
eigs=ordeig(AA,BB).^-1; %output from ordeig(AA,BB) is INVERTED from eig(BX1^-1*BX0), i.e. explosive are <1 before inversion

[T1schur,T0schur,Qprime,Zprime]=ordqz(AA,BB,Q,Z,'udo'); %Matlab spits out T1schur = Qschur*BX1*Zschur and T0schur = Qschur*BX0*Zschur 
%to be consistent with the paper's notation:
Qschur=Qprime'; 
Zschur=Zprime';
%if you want to check: BX0 == Qschur*T0schur*Zschur; BX1 == Qschur*T1schur*Zschur

%check to ensure number of state variables corresponds to number of stable
%generalized eigenvalues:
if sum(abs(eigs)>1)~= nc
    display('Incorrect number of explosive eigenvalues: model not exactly determined')
end

T1_SS = T1schur(1:ns,1:ns); 
T1_SC = T1schur(1:ns,ns+1:end);
T1_CC = T1schur(ns+1:end,ns+1:end);
T0_SS = T0schur(1:ns,1:ns);
T0_SC = T0schur(1:ns,ns+1:end);
T0_CC = T0schur(ns+1:end,ns+1:end);
%inverting once saves some computation time for large models
T1_SS_inv = T1_SS^-1; 
T0_CC_inv = T0_CC^-1; 


%Define relevant submatrices (Equation 7)
Z_CC = Zschur(ns+1:end, control_vars);
Z_CC_inv = Z_CC^-1;
Z_SC = Zschur(1:ns, control_vars);

%Xi and Theta operators as defined in MMIIES
Xi_C = zeros(nc,nc,2*T-1);
for j=0:(T-1)
    Xi_C(:,:,T-j) =  -(- T0_CC_inv*T1_CC)^j*T0_CC_inv ; %NEGATIVE INSIDE OR OUTSIDE?? both!
end
Xi_op = zeros(nc+n,n,2*T-1);
%final nc rows determine the unstable dimensions (ZX)_{C}:
Xi_op(n+1:end,ns+1:end,:)=Xi_C;
%first nc rows correspond to the forecast error (ZX)_{C,0}:
Xi_op(1:nc,ns+1:end,:)=lagmulti(Xi_C,-1);
%remaining ns rows correspond to the stable dimensions (ZX)_{S}:
Xi_op(nc+1:nc+ns,1:ns,T) = eye(ns,ns);
%Take a look: Xi_op(:,:,T-2:T+1)
Xi=smulti(Xi_op,Qschur',T);

B_S_inv =  zeros(ns,ns,T);
for jj = 1:T
B_S_inv(:,:,jj)=(-T1_SS_inv*T0_SS)^(jj-1); %needs that negative!
end
Theta_1 = lagmulti(smulti(-B_S_inv,Z_SC*Z_CC_inv+T1_SS_inv*T1_SC,T),1);
Theta_2 = lagmulti(smulti(-B_S_inv,T1_SS_inv,T),1);
Theta_3 = smulti(B_S_inv,Z_SC*Z_CC_inv,T) - lagmulti(smulti(B_S_inv,T1_SS_inv*T0_SC,T),1);
Theta_op = zeros(n,nc+n,T);
Theta_op(1:ns,1:nc,:) = Theta_1;
Theta_op(1:ns,nc+1:n,:) = Theta_2;
Theta_op(1:ns,n+1:end,:) = Theta_3;
Theta_op(ns+1:end,n+1:end,1) = eye(nc,nc);
%Take a look: Theta_op(:,:,1:3)
Theta = smulti(Zschur',Theta_op,T);


end