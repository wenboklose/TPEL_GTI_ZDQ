% Function "sortloci" returns sorted set of eigenvalues.
%
% Call: [sorted_eigenvalues]= sortloci(eigenvalues)


function [sorted]= sortloci(eigen)

n=length(eigen);
mirror=eigen;

% sorting algorithm for continuity
dif1=abs(eigen(1,1:n-1)-eigen(1,2:n));
dif2=abs(eigen(1,1:n-1)-eigen(2,2:n));
for k=1:n-1
    if dif1(k)>dif2(k)
        eigen(1,k+1:n)=mirror(2,k+1:n);
        eigen(2,k+1:n)=mirror(1,k+1:n);
        mirror=eigen;
    end
end
sorted=eigen;



