function [A, W] = wold(S_sig,N,Sigma,method)
if nargin<3
    Sigma=eye(size(S_sig,2)); %If shock variance is unspecified, set to identity matrix
end
if nargin<4
    method = 'yule walker';  %Cholesky version is slightly less accurate and slightly slower
end

%Drop zero rows from S_sig:
rowsums = sum(sum(abs(S_sig),3),2)~=0;
S=S_sig(rowsums~=0,:,:);

%Rewrite S to depend on unit variance white noise shocks
if issymmetric(Sigma)
    if all(eig(Sigma)>0) % if positive definite
       C_eps = chol(Sigma)'; %this Cholesky is transposed so that C_eps*C_eps' = Sigma
    else      
        Sigma_nonzero=Sigma((all(Sigma==0)==0),(all(Sigma==0)==0));
        C_eps = 0*Sigma;
        C_eps((all(Sigma==0)==0),(all(Sigma==0)==0))=chol(Sigma_nonzero)'; %this Cholesky is transposed so that C_eps*C_eps' = Sigma
    end
else
    disp('Sigma is not symmetric')
    stop
end
Swn = smulti(S,C_eps,N); 

if size(S,3)==N %if S is causal... otherwise it will have 2N-1 array dimension
    Gamma=t2s(s2t(Swn,N)*s2t(Swn,N)',N,1); %selects the correct entries for Gamma
    Gamma_T = s2t(Gamma,N); %transforms back to Toeplitz
else %%%  
    Snc=Swn;
    Snc_Gamma=s2t(Snc,N)*s2t(Snc,N)';
    Snc_Gamma_v = v2s(Snc_Gamma(:,(end-size(S,1)+1):end),N);
    Snc_Gamma = cat(3,Snc_Gamma_v(:,:,1:N-1),flip(Snc_Gamma_v,3));
    Gamma_T = s2t(Snc_Gamma,N);
    display('This non-causal case has not been updated since I switched to upper triangular form. If you need to Wold a noncausal signal, check this code.')
end

if strcmp(method , 'cholesky')
%You cannot directly Cholesky decompose Gamma, because it will give you LT * UT
ACw_t=chol(Gamma_T^-1)^-1;
Cw_inv = ACw_t(1:size(S,1), 1:size(S,1))^(-1); %Cholesky of the W variance
ACw = t2s(ACw_t,N);
Z = smulti(ACw,Cw_inv,N);
Zinv = inv_signal(Z,N);
end

if strcmp(method , 'yule walker')
%Here's the Yule Walker solution
R = v2s(s2v(Gamma(:,:,N+1:end),N-1,0)*s2t(Gamma(:,:,2:end-1),N-1)^-1,N-1,0);
Zinv = cat(3,eye(size(R,1)),-R);
Z = t2s(s2t(Zinv,N)^-1,N);
end

 
W = smulti (Zinv, S , N);  %W process that depends on shocks with variance matrix Sigma: W * Sigma * W' = Sigma_W (block diagonal)

%reintroduce any zero rows:
A = zeros(size(S_sig,1),size(W,1),N);
A(rowsums~=0,:,:) = Z;

end


%   SIGNAL_OP Package: Version 1.1
%   Jonathan J. Adams
%   Contact: adamsjonathan@ufl.edu