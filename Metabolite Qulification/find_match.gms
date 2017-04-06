$include GAMS_data.txt
*display S
binary variable x(i,j);
positive variable concentration(m);    concentration.up(m)=maxConcentration;
positive variable height(i);          height.up(i)= maxConcentration*pks_theory(i);
positive variable c1(i,j);
variable obj;
*-------------------------------------------------------------------------------
equation eq1; eq1(i) .. height(i) =e= sum(m$K(m,i), concentration(m)*pks_theory(i));
equation eq2; eq2(i,j)$S(i,j) .. c1(i,j) =e=  exp( -sqr(height(i)-pks_sample(j))/sqr(sigma)/2    );
equation eq3; eq3(i) .. sum(j$S(i,j),x(i,j)) =l= 1;
equation eq4; eq4(j)$(sum(i$S(i,j),1) ge 1) .. sum(i$S(i,j),x(i,j)) =l= 1;
equation eq5; eq5 .. obj =e= sum((i,j)$S(i,j), x(i,j)*c1(i,j)*c2(i,j));


model test /all/;
option minlp=baron;
option optcr=0;
solve test max obj use minlp;

*-------------------------------------------------------------------------------
*display x.l
parameter c(i,j), error(i,j);
loop((i,j)$S(i,j),
         c(i,j) = x.l(i,j)*c1.l(i,j)*c2(i,j);
         error(i,j) = x.l(i,j)* ( height.l(i) - pks_sample(j))/pks_sample(j);
);
*display c;
*display error;

*-------------------------------------------------------------------------------
file solFile /"GAMS_output.txt"/;
put solFile;
put '#'/;
loop((i,j)$(S(i,j) and x.l(i,j) ge 0.99 and c(i,j) ge 0.6),
   put ord(i):0:0', 'ord(j):0:0/;
);
put '#'/;
loop(m$(concentration.l(m) ge 0.0001),
   put ord(m):0:0', 'concentration.l(m):0:8/;
);
putclose;
