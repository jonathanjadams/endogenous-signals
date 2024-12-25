function H = s2h(S,N)

if size(S,3)==N %if S is causal; otherwise it will have 2N-1 array dimension
Snc = cat(3,0*S(:,:,1:end-1),S); %adds non-causal zero terms
else
    if size(S,3)== 2*N-1 %if S is non-causal, use the whole thing
    Snc=S;
    else
        if size(S,3)== 1 %if S is just a matrix, assume it's the only non-zero matrix in the array; place in 0 position
            Snc=cat(3,zeros([size(S) N-1]),S,zeros([size(S) N-1]));
        else
            display('error: S does nto have an acceptable number of array dimensions (N, 2N-1, or 1)')
            Snc=0;
        end
    end
end


CellS=mat2cell(Snc,size(Snc,1),size(Snc,2),ones(2*N-1,1));  %Note that cell indexing is the reverse of the array indexing
CellS(:,:,2*N)=mat2cell(zeros(size(Snc,1),size(Snc,2)),size(Snc,1),size(Snc,2)); %makes 0 matrix the final entry
H=cell2mat(CellS(hankel(N:1:(2*N-1),[2*N-1, 2*N*ones(1,N-1)]))); %Creates the nN x mN matrix; causal nonzero entries appear in the upper left triangle.


end


%   SIGNAL_OP Package: Version 1.0
%   Jonathan J. Adams
%   Contact: adamsjonathan@ufl.edu
