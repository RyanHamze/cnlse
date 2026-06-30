function y=fastexp(x)

%FASTEXP Calculates exp( i * x ) quickly
%   Y=FASTEXP(X) 
%   It's a fast computation of the expression:
%       Y = EXP( i * X )


y=complex(cos(x),sin(x));
