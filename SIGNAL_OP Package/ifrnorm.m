function IFRnorm = ifrnorm(Theta,Xi,BA0,BA1,nc,G_sig,T) %this function solves a dispersed information model with Signal Operator Iteration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Preliminary variable assignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dimension variables
ns = size(Theta,1)-nc;
n = ns+nc;
m_z = size(BA0,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Construct IFR Operator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BA_terms = zeros(n,m_z,2*T-1);
BA_terms(:,:,T)=BA0; BA_terms(:,:,T-1)=BA1;

IFR_operator_causal = smulti(G_sig,Theta,T);
IFR_operator_noncausal = smulti(Xi,BA_terms,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Calculate IFR Norm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IFRnorm = norm(s2t(IFR_operator_causal,T)*s2t(IFR_operator_noncausal,T));

if IFRnorm<1
display(' ')
display(strcat('Information Feedback Regularity satisfied: operator norm = ',sprintf(' %g ',IFRnorm)))
else
display(' ')
display(strcat('Information Feedback Regularity not satisfied: operator norm = ',sprintf(' %g ',IFRnorm)))
end

end