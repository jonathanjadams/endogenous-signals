function block_r_S = block_r(S,k,N)

if size(S,3)==1 && size(S,1) >=N %Toeplitz matrix as input 
        
    T=S;
    block_r_S = kron(T,eye(k)); %notice that if S is m x n, output is km x kn
    
else  %signal array as input

    block_r_S_sig = NaN(size(S,1)*k,size(S,2)*k,size(S,3));

    for jj=1:size(S,3)
        block_r_S_sig(:,:,jj) = kron(S(:,:,jj),eye(k));
    end
        block_r_S=s2t(block_r_S_sig,N);

end

end

%   SIGNAL_OP Package: Version 1.1
%   Jonathan J. Adams
%   Contact: adamsjonathan@ufl.edu