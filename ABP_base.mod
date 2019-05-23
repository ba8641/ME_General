/*
	ABPTSS, TCS, and TPF MODEL ELEMENTS (model file) and BASE FORMULATION
	-  based on Proano et al, Omega 40(2012):54-63
	 - modular implementation
	by Ruben A. Proano
	rpmeie@rit.edu
	ROCHESTER INSTITUTE OF TECHNOLOGY
	
	NOTE:
	Changes with respect to ABP_base.mod version Feb 2017:
	- Eliminates set A2 
	- Redefines set B2:
	    - previously 'set B2{B} within P;'
	    - current set B2{P} within B;'
	
    - Adds additional elements needed for referring
      upper and lower bounds on X and Y, as well as 
      the declaration of variable W needed for McCormick 
      envelop linearization
      
    - Change ElasticityNL constraint to modify the 
      scalar in the RHS to ensure that price is small
      if bundle demand equals its maximum upper bound (Xup).
	
	See additional comments at the end of the file
	
	07.05.2017
	version I.02

*/


/* General elements of the ABP formulation */ 

param n_antigens >=0;
#param n_markets  >=0;

param n_bundles >=0;
param n_producers >=0;

set B:= 1..n_bundles; 
set A;
set M_ALL:=1..194;
set M_DROP within M_ALL default {13..194};
set M:= M_ALL diff M_DROP;
#set M:= 1..n_markets;
set P:= 1..n_producers;

set A1{B} within A;
#set A2{M} within A; 
set B1{A} within B;
set B2{P} within B;

param R{B,M}; 
param l{M};
param C{B};
param d{A,M};
param D{B,M};
param S{B};
param k{B} within P;		#Isn't this basically B2?
param phi; 
param name{B} symbolic; 
param gni_p{M};
param theta;
param producing{B} default 0;	#Displays 1 if the bundle is being purchased by any entity; 0 otherwise.;

var X{B,M} >=0 default 0;
var Y{B,M}>=0 default 0; 
var g{B} binary default 0;


/* Elements used in SupperAdditivity constraint */
set Q{B} default {}; 
set QQ:= union {b in B} Q[b];
set N{QQ} within B;  

/* Elements needed for implementation of McCormic Linearization  6/12/2017*/
var W{B,M} >=0 default 0;
param Xlow{b in B, m in M}>=0 default 0;
param Xup{b in B,m in M} default min(S[b], D[b,m]*l[m]);
param Ylow{b in B, m in M} >=0 default 0;
param Yup{b in B, m in M} default theta*R[b,m];

var loss{P} binary default 0;	#Becomes 1 when the manufacturer has a loss.
#var loss{B} binary default 0;	#Becomes 1 when the bundle has a loss.





maximize TSS: 
sum{b in B, m in M} R[b,m]*X[b,m]- sum{b in B}C[b]*g[b];
#sum{b in B, m in M} R[b,m]*X[b,m]- sum{b in B}C[b]*g[b] - sum{b in B}pen[b]*1e9; 

#maximize TCS: sum{b in B, m in M} (R[b,m]- Y[b,m])*X[b,m];
#maximize TCS: sum{b in B, m in M} (R[b,m]- Y[b,m])*X[b,m] - sum{p in P, b in B: k[b]==p and producing[b]==1}(C[b]-sum{m in M}X[b,m]*Y[b,m]) ;	#Needs some way to keep prices from inflating. Should only discount when <0.
maximize TCS: sum{b in B, m in M} (R[b,m]- Y[b,m])*X[b,m] - 10e10*sum{p in P}loss[p];
#maximize TCS: sum{b in B, m in M} (R[b,m]- Y[b,m])*X[b,m] - 10e10*sum{b in B}loss[b];

#maximize TPF: sum{b in B, m in M}Y[b,m]*X[b,m] - sum{b in B} C[b]*g[b];
#maximize TPF: sum{b in B, m in M}Y[b,m]*X[b,m] - sum{b in B} C[b]*g[b] - sum{p in P, b in B: k[b]==p and producing[b]==1}(C[b]-sum{m in M}X[b,m]*Y[b,m]);
maximize TPF: sum{b in B, m in M}Y[b,m]*X[b,m] - sum{b in B} C[b]*g[b] - 10e10*sum{p in P}loss[p];
#maximize TPF: sum{b in B, m in M}Y[b,m]*X[b,m] - sum{b in B} C[b]*g[b] - 10e10*sum{b in B}loss[b];

subject to SuperAdditivity{b in B, m in M, q in Q[b]: card(Q[b])>=1}:  
	R[b,m]- Y[b,m] >= sum{t in N[q]}(R[t,m]*g[t] - Y[t,m]);

/* If Q[b] = 3 6;
subject to SuperAdditivity{b in B, m in M: card(Q[b])>=1}:  
	R[b,m]- Y[b,m] >= sum{t in Q[b]}(R[t,m]*g[t] - Y[t,m]);

/*subject to SuperAdditivity{b in B, m in M, q in Q[b]: card(Q[b])>=1}:
    Y[b,m]+(1-g[b])*(sum{t in N[q]}R[t,m])>= sum{t in N[q]}Y[t,m];
*/

subject to AntigenDemand{a in A, m in M}: sum{b in B1[a]}X[b,m] <= d[a,m]*l[m];

subject to Capacity{b in B}: sum{m in M} X[b,m] <= S[b]*g[b];

subject to RecoverAnnuity{ b in B}: sum{m in M} Y[b,m]*X[b,m] >= C[b]*g[b];

# subject to RecoverAnnuity{ b in B}: sum{m in M} Y[b,m]*X[b,m] >= C[b]*(1 - pen[b]);
# subject to PEN{ b in B}:pen[B] <= interest_rate;

subject to RecoverAnnuityTotal{ p in P}: sum{ m in M, b in M: k[b]==p}Y[b,m]*X[b,m] >=sum{b in B:k[b]==p} C[b]*g[b];	

subject to Elasticity{b in B, m in M}: 
	X[b,m]<= Xup[b,m]*(1-(Y[b,m]/(Yup[b,m]))^phi); 


/* Linear Elasticity constraint to be used in the LP problems ABPTCS and ABPTPF, that consider Y as the only variable of interest*/
subject to ElasticityNLY{b in B, m in M }:
 	Y[b,m] <= Yup[b,m]*(1.00000000001- X[b,m]/Xup[b,m])^(1/phi);   

subject to SegmentedPrices{b in B, m in M, n in M: m<>n and gni_p[m] >=gni_p[n]}: 	#Should only be active when buying the thing, or will conflict with minPrice
	Y[b,m]>=Y[b,n];  

subject to ReservationPriceLim{b in B, m in M}: 
	Y[b,m] <= theta*R[b,m]*g[b];

subject to DefineLoss{p in P}:
	sum{b in B: producing[b]=1 and k[b]=p}(C[b] - sum{m in M}X[b,m]*Y[b,m]) <= 10e15*loss[p];

#subject to DefineLoss{b in B}:
#	C[b] - sum{m in M}X[b,m]*Y[b,m] <= 10e15*loss[b];



