function frechetnorm = soifrechet(S_sig,Theta,Xi,BA0,BA1,nc,G_sig,P_G,T) %this function solves a dispersed information model with Signal Operator Iteration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Preliminary variable assignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dimension variables
n_zeta = size(Xi,1); %n_zeta =  n + ns;
m_z = size(S_sig,1);
m_eps = size(S_sig,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Component operator construction 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Construct Signal-derived operators
PS = s2t(S_sig,T)'*(s2t(S_sig,T)*s2t(S_sig,T)')^-1*s2t(S_sig,T); %why are the primes where they are? s2t() converts S_sig to an upper triangular block operator wihtout block transposing
MS = eye(m_eps*T)-PS;%m_eps x m_eps blocks
S_Linv = s2t(S_sig,T)'*(s2t(S_sig,T)*s2t(S_sig,T)')^-1; 

%Construct zeta operators
zeta = smulti(Xi,BA0,T) + lagmulti(smulti(Xi,BA1,T),-1); %zeta = Xi(L)*(BA0 + BA1 L^-1)
zetaS = smulti(zeta,S_sig,T); %zeta * S
zetaS=zetaS(:,:,T:end); %Only keep the causal component
zetaS_MS = t2s(s2t(zetaS,T)*MS,T,1); %orthogonalize with respect to signals

%Construct Hankel operator
zetaS_MS_kron = NaN(m_z*size(zetaS_MS,1),m_z*size(zetaS_MS,2),T); 
for jj = 1:T
zetaS_MS_kron(:,:,jj)=kron(eye(m_z),zetaS_MS(:,:,jj)); 
end
L_HzetaS_MS=s2t(perm(n_zeta,m_z),T)*s2h(zetaS_MS_kron,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Calculate Subspace Coefficient Q_MS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R_TzetaS_SLinv = block_r(s2t(zetaS,T)*S_Linv,m_eps,T); 
Q_MS = R_TzetaS_SLinv + block_l(S_Linv,n_zeta,T)*L_HzetaS_MS;
term2 = block_l(P_G,m_z,T)*block_r(smulti(G_sig,Theta,T),m_eps,T)*Q_MS*block_l(MS,m_z,T); 

[outputtext,IFRnorm] = evalc('ifrnorm(Theta,Xi,BA0,BA1,nc,G_sig,T);');

frechetnorm = sqrt(IFRnorm^2,norm(term2)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Display Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if frechetnorm <1
    display(strcat('Solution signal operator is stable with Frechet derivative norm',sprintf(' %g',frechetnorm)))
else
    display(strcat('Solution signal operator is unstable with Frechet derivative norm',sprintf(' %g',frechetnorm)))
end



end