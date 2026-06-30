function y=nmod(A,N)

%NMOD N-modulus of an integer.
%   Y=NMOD(A,N) reduces the integer A into 1->N, mod N.
%
%   E.g. N=8. 
%
%   A   ... -2 -1 0 1 2 3 4 5 6 7 8 9 10 ...
%   Y   ...  6  7 8 1 2 3 4 5 6 7 8 1 2  ...
%



y = mod(A-1-N,N)+1;
