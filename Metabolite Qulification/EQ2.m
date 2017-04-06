function C = EQ2( A,B )
%EQ Summary of this function goes here
%   Detailed explanation goes here
    
    tol = abs((A(1)-A(2)))/2;
    C = ( abs(A-B) < tol); 
   
end

