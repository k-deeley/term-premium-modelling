function Q = vech(A)
[n1,n2] = size(A);
if n1 ~= n2
    error('A is not square')
end
Q       = zeros(n1*(n1+1)/2,1);
idx     = 0;
for i1=1:n1
    for i2=i1:n1
        idx = idx + 1;
        Q(idx,1) = A(i1,i2);
    end
end

end

