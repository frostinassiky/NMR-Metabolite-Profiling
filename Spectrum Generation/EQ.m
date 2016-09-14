function C = EQ( A,B )
%EQ Summary of this function goes here
%   Detailed explanation goes here
    
    tol = eps(20);
    C = ( abs(A-B) < tol); 
   
end

