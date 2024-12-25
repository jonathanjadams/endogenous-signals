function block_l_S = block_l(S,k,N)

if size(S,3)==1 && size(S,1) >=N %Toeplitz matrix as input
     T=S;
    block_l_S = kron(eye(k),T);  %notice that if S is m x n, output is km x kn
    
else %signal array as input

    block_l_S_sig = NaN(k*size(S,2),k*size(S,1),size(S,3));

    for jj=1:size(S,3)
        block_l_S_sig(:,:,jj) = kron(eye(k),S(:,:,jj)');
    end
        block_l_S=s2t(block_l_S_sig,N);

end
    
    
end

%   SIGNAL_OP Package: Version 1.1
%   Jonathan J. Adams
%   Contact: adamsjonathan@ufl.edu