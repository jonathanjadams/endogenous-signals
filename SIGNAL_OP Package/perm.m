function P = perm(m,n)

%create m x n permuation matrix: for m x n matrix A, such that Perm*vec(A')=vec(A)

P = zeros(m*n); 
for ii=1:m
    for jj = 1:n
        %A_{ii,jj} is entry m*(jj-1)+ii in vec(A)
        %Need to map entry n*(ii-1)+jj to m*(jj-1)+ii
        P(m*(jj-1)+ii,n*(ii-1)+jj)=1;
    end
end


%   Reference:
%   H. V. Henderson and S.R. Searle The vec-permutation matrix,
%   the vec operator and Kronecker products: A review Linear and
%   Multilinear Algebra, 9 (1981), pp. 271-288.

end

