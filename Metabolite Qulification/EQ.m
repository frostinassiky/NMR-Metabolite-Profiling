function C = EQ( A,B )
%EQ Summary of this function goes here
%   Detailed explanation goes here
    
    tol = 1e-5;
    C = ( abs(A-B) < tol); 
   
end

