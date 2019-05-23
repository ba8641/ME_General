/*
	ABPTSS_Lin and ABPTSS_H MODEL ELEMENTS -  based on Proano et al, Omega 40(2012):54-63
	LINEARIZED VERSION - modular implementation
	by Ruben A. Proano
	rpmeie@rit.edu
	ROCHESTER INSTITUTE OF TECHNOLOGY
	
	NOTE:
	# 1- Linearization of constraint (2) - subadditivity on value for a combination vaccine
	# 2- Linearization of constraint (6) - demand elasticity restriction
	# 3- It assumes that for each bundle id, there is only one manufacturer (e.i.,there could be several 
	#   bundles with the same antigen composition but with different and unique bundle ids.
	# 4- Instead of using a set of producers to own the bundles, each bundle owns a producer via a parameter
		
	07.05.2017
	version I.01
    
    ############################################################
    # INCLUDE ABP_base.mod PRIOR THIS FILE SO THAT ALL GENERAL #
    # ABP model elements are properly loaded                   #
    ############################################################
    
	See additional comments at the end of the file
*/




/*
    General ABP model elements loaded in ABP_base.mod
*/


/* Elements needed for linearization of the NL elasticity constraints
 - number of line segments in which the allocation function is partition
 - auxiliary price variable
 - price breakpoints at which linear segments intersect
 - rates (slopes) of each linear segment */

param n_linsegments >=0;
var y{B, M} >=0;	
param bkpt{b in B, m in M, bp in 1..n_linsegments-1}:=bp*R[b,m]/n_linsegments;
param rate{ b in B, m in M, s in 1..n_linsegments}:= n_linsegments*D[b,m]*l[m]*(((s-1)/n_linsegments)^phi-(s/n_linsegments)^phi)/R[b,m];



param lo{a in A,m in M} default l[m];

/* Elements needed to keep track of allocations during the iterative 
stage of the ABP_heuristic 
   - parameter to account for total antigen demand in a market left to be met
   - parameter used for storing bundle allocations in each market segment
   - parameter used for storing information of bundles that must be produced
   - tabu set to store bundles that have been used in the iterative process
   - parameters used for defining the range of bundles in the Sorted bundle set
     that would go in the ABP iteration
*/
param numiterations default card(B);
param dH{a in A, m in M};
param Xo{B,M} default 0;
param go{B} binary default 0;

set LL default {}; 

# parameters used to define the indices of bundles included in the 
# sorted list of bundles arranged by decreasing fitness factor
 
param loo default 1;  
param hii default card(B);

# Fitness feature to sort all bundles, and lowest and highest fitness 
# among all bundles.

param bfit{b in B} default(sum{m in M}R[b,m]*D[b,m])/log(C[b]);   # bundle fitness
param LO:= min {b in B} bfit[b]; 
param HI:= max {b in B} bfit[b];

# creating an ordered set of bundles based on their fitness

set temps ordered by reversed [LO,HI] := setof{b in B} bfit[b];
set Bsort ordered:= setof{kk in temps, b in B: bfit[b]=kk} b;

# deficit variable
var deficit{A,M} >=0;

maximize TSSH_com: 
sum{b in setof{z in loo..hii} member(z, Bsort), m in M:b not in LL} R[b,m]*X[b,m] 
- sum{b in setof{z in loo..hii} member(z, Bsort):b not in LL}C[b]*g[b]; 

maximize TSSH_mono: 
sum{b in setof{z in loo..hii} member(z, Bsort), m in M:b not in LL} R[b,m]*X[b,m] 
- sum{b in setof{z in loo..hii} member(z, Bsort):b not in LL}C[b]*g[b]
-1000*sum{a in A, m in M} deficit[a,m]; 

subject to AntigenDemandH{a in A, m in M}: 
sum{b in B1[a]: b in setof{z in loo..hii} member(z, Bsort) and b not in LL}
X[b,m] <= dH[a,m];

subject to AntigenDemandH_mono{a in A, m in M}: 
	sum{b in B1[a]: b in setof{z in loo..hii} member(z, Bsort)and b not in LL}
	X[b,m] + deficit[a,m]<= dH[a,m];

subject to CapacityH{b in setof{z in loo..hii} member(z, Bsort): b not in LL}: 
	sum{m in M} X[b,m] <= S[b]*g[b];

subject to RecoverAnnuityH{p in P, b in setof{z in loo..hii} member(z, Bsort): 
k[b]==p and b not in LL}: 	sum{m in M} Y[b,m]*X[b,m] >= C[b]*g[b];

# Linearized Elasticity constraint (ABPTSS_Lin, ABPTSS_LH_com, ABPTSS_LH_mono)

subject to ElasticityL{b in B, m in M}: 
	X[b,m] <= <<{bp in 1..n_linsegments-1} bkpt[b,m,bp] ;{s in 1..n_linsegments} rate[b,m,s] >>(Y[b,m], theta*R[b,m]);


subject to ElasticityLH{b in setof{z in loo..hii} member(z, Bsort), m in M: b not in LL}: 
	X[b,m] <= <<{bp in 1..n_linsegments-1} bkpt[b,m,bp] ;{s in 1..n_linsegments} rate[b,m,s] >>(y[b,m], R[b,m]);
	
subject to ElasticityHH{b in setof{z in loo..hii} member(z, Bsort), m in M:b not in LL}: 
	X[b,m]<= D[b,m]*l[m]*(1-(Y[b,m]/R[b,m])^phi); 

subject to SegmentedPricesH{b in setof{z in loo..hii} member(z, Bsort), m in M, n in M: gni_p[m] >=gni_p[n] and b not in LL}: 
	Y[b,m]>=Y[b,n];  

subject to ReservationPriceLimH{b in setof{z in loo..hii} member(z, Bsort), m in M: b not in LL}: 
	Y[b,m] <= R[b,m]*g[b];

# Constraints needed for ABP_heuristic to fixed partially obtained allocations

subject to PartialProcurementH{b in B, m in M:Xo[b,m]>=1}: X[b,m]=Xo[b,m];

subject to PartialSetUpH{b in B:go[b]>=1}: g[b]>=go[b];


# Constraints needed for implementing McCormick envelop linearization on the RecoverAnnuity

subject to RecoverAnnuityLMcMain{p in P, b in B: k[b]==p}: sum{m in M} W[b,m]>= C[b]*g[b];
subject to RecoverAnnuityLMc1{b in B, m in M}: W[b,m] >= Xup[b,m]*Y[b,m]+ X[b,m]*Yup[b,m]-Xup[b,m]*Yup[b,m];
subject to RecoverAnnuityLMc2{b in B, m in M}: W[b,m] >= Xlow[b,m]*Y[b,m]+ X[b,m]*Ylow[b,m]-Xlow[b,m]*Ylow[b,m];
subject to RecoverAnnuityLMc3{b in B, m in M}: W[b,m] <= Xup[b,m]*Y[b,m]+ X[b,m]*Ylow[b,m]-Xup[b,m]*Ylow[b,m];
subject to RecoverAnnuityLMc4{b in B, m in M}: W[b,m] <= Xlow[b,m]*Y[b,m]+ X[b,m]*Yup[b,m]-Xlow[b,m]*Yup[b,m];
subject to RecoverAnnuityLMc5{b in B, m in M}: Xlow[b,m]<=X[b,m]<=Xup[b,m];
subject to RecoverAnnuityLMc6{b in B, m in M}: Ylow[b,m]<=Y[b,m]<=Yup[b,m];





