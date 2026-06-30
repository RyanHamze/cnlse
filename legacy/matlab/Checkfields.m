function checkfields(x,allfields)

%CHECKFIELDS check for valid input fields
%   CHECKFIELDS(X,ALLFIELDS) checks if the struct variable X has valid
%   fields. ALLFIELDS is a cell variable containing all possible valid
%   fields of X. Each element of the cell is a string. The check is 
%   case insensitive.
%
%   Note: Use CHECKFIELDS only for low level functions. If a function calls
%   a subfunction, CHECKFIELDS should be run only in the nested subfunction.
%   Do not abuse the function use. Beware.
%

xf = fieldnames(x);
Lx = length(xf);

for k=1:Lx
    chk = strcmp(xf{k},allfields);
    if ~any(chk)
        error('Non existing field %s',xf{k});
    end
end
