
function S = t2s(T,N,noncausal)
if nargin<3
    noncausal=0;
end
%causal signal process = upper triangular toeplitz matrix
    V = T(1:size(T,1)/N,:); %Select the top block row of the Toeplitz matrix
    S = v2s(V,N);  
if noncausal==1
    %This takes the first block column, after the first row:
    S = cat(3,v2s( T(1+size(T,1)/N:end,1:size(T,2)/N),N-1,1),S);  
end

%for noncausal == 2 option:
%this version uses the final row of the matrix rather than the first column to identify the non-causal terms
%and final column to identify causal terms
if noncausal==2 
    V = T(end+1-size(T,1)/N:end,:); %The bottom block row of the Toeplitz matrix
    Snc = v2s(V,N);
    V=T(:,end+1-size(T,2)/N:end,:);
    S=v2s(V,N,1);
    S = cat(3,Snc(:,:,1:N-1),S);  
end

%both methods correctly invert s2t, i.e. t2s(s2t(S,N),N,k) = S for k=1 and k=2

end


%   SIGNAL_OP Package: Version 1.1
%   Jonathan J. Adams
%   Contact: adamsjonathan@ufl.edu