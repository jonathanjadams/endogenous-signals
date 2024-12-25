function [X_shock_sig, X_sig_full, IFRnorm, S_sig, A_sig, X_sig, relativeerror] = soi_solve(BX0,BX1,BA0,BA1,nc,P_G,G_sig,S_X_sig,Sigma,calculate_full_info,maxiter) %this function solves a dispersed information model with Signal Operator Iteration


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Preliminary variable assignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dimension variables
ns = size(BX0,2)-nc;
%n = ns+nc;
%m_z = size(BA0,2);
T = size(S_X_sig,3); %truncation length


%set computational parameters:
tolerance=1e-6;  
error_upperbound = 1e5; %set an upper bound to exit loop if algorithm diverges
if nargin<11
    maxiter=500;
end

%Construct the Theta, Xi operators from the model ingredients:
[Theta, Xi] = theta_xi(BX0,BX1,nc,T);

%Is Information Feedback Regularity satisfied?
IFRnorm = ifrnorm(Theta,Xi,BA0,BA1,nc,G_sig,T);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Apply S.O.I. to solve the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%Calculate the Full Information solution
if calculate_full_info==1
%now, your Wold representation is S!  

%initial guess:
S_sig_full = S_X_sig;

%set error, iteration counter:
error=100; iter=1;
    
while error>tolerance && iter<=maxiter && error<error_upperbound
    %Steps 1 and 2: Find the Wold Representation
    %[A_sig, W_sig] = wold(S_sig,T,Sigma);
    A_sig_full = S_sig_full; 

    %Step 3: define A tilde= [(BA1 L^-1 +BA0)A]_+
    %A_tilde_sig =  smulti(Qschur*BZ0,A_sig,T) + smulti(Qschur*BZ1,lagmulti(A_sig,-1),T);
    QA_tilde_sig_full =  smulti(BA0,A_sig_full,T) + smulti(BA1,lagmulti(A_sig_full,-1),T);
    
    %Step 4: Find the policy function
    XiA_noncausal_full = smulti(Xi,QA_tilde_sig_full,T);
    XiA_causal_full = XiA_noncausal_full(:,:,T:end);
    X_sig_full = smulti(Theta,XiA_causal_full,T);
    
    %Step 5: Calculate the next signal polynomial
    if size(G_sig,3)==T %FOR CAUSAL G (standard)
        S_next_sig_full = S_X_sig + smulti(G_sig,smulti(X_sig_full,P_G,T),T);
    end
    if size(G_sig,3)==2*T-1 %FOR NONCAUSAL G 
        GX_sig_full = smulti(G_sig,X_sig_full,T);
        S_next_sig_full = S_X_sig + smulti(GX_sig_full(:,:,T:2*T-1),P_G,T);
    end
    %Calculate difference:
    error=norm(s2t(S_next_sig_full,T)-s2t(S_sig_full,T));
    display(strcat('FI iteration ',sprintf(' %g ',iter),': error ',sprintf(' %g ',error)))

    %update:
    S_sig_full= S_next_sig_full;
    iter=iter+1;
    end
  iter=iter-1; %role back one iteration if the loop ended
  
if error<=tolerance
    display(' ')
    display(strcat('Full Information algorithm converged in',sprintf(' %g ',iter),' iterations, with final error',sprintf(' %g. ',error)))
else
    if error>=error_upperbound
        display(' ')
        display(strcat('Full Information algorithm diverged in',sprintf(' %g ',iter),' iterations, with final error',sprintf(' %g. ',error)))
    else
        display(' ')
        display(strcat('Full Information algorithm reached iteration limit of',sprintf(' %g ',maxiter),', with final error',sprintf(' %g. ',error)))   

    end
end


%Initial guess: full info solution
S_sig = S_sig_full;

else
    X_sig_full = NaN*S_X_sig; %to output if you did not calculate it
    S_sig = S_X_sig; %shock sig
end


%Calculate the DISPERSED INFORMATION solution

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
    

%update:
S_sig=S_next_sig;
iter=iter+1;
end
iter=iter-1; %role back one iteration if the loop ended

if error<=tolerance
    display(' ')
    display(strcat('Dispersed Information algorithm converged in',sprintf(' %g ',iter),' iterations, with final error',sprintf(' %g. ',error)))
else
    if error>=error_upperbound
        display(' ')
        display(strcat('Dispersed Information algorithm diverged in',sprintf(' %g ',iter),' iterations, with final error',sprintf(' %g. ',error)))
    else
        display(' ')
        display(strcat('Dispersed Information algorithm reached iteration limit of',sprintf(' %g ',maxiter),', with final error',sprintf(' %g. ',error)))   

    end
end


relativeerror = error/initialerror; %record the relative error compared to error given initial guess

% Construct IRFS 
% X_sig gives the responses of the X vector (consumption and capital) to innovations in each signal.
% We also want to know the response to the fundamental shocks (aggregate)
X_shock_sig = smulti(X_sig,W_sig,T);


end


%   SIGNAL_OP Package: Version 1.2
%   Jonathan J. Adams
%   Contact: adamsjonathan@ufl.edu