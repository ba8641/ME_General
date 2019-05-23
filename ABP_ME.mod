/*
	ABPTSSME, TCSME, and TPFME MODEL ELEMENTS (model file) and BASE FORMULATION
	-  based on Proano et al, Omega 40(2012):54-63
	- modular implementation
	- adapted from ABP_base.mod, by Dr.Proano.
	by Bruno Alves Maciel
	ba8641@rit.edu
	ROCHESTER INSTITUTE OF TECHNOLOGY
	
	NOTE:
	
	07.11.2017
	version I.01

	See additional comments at the end of the file
*/


/* General elements of the ABPME formulation */ 

set T{B} within M;	#Defines the target market

param alpha;	#How much we penalize not buying the entire supply of a bundle

param minPrice{B};

param X2{b in B, m in M} default 0;

param Cl{B}; 	#Keeps track of how much of the annuity costs remain to be paid.
param Countries{1..194,M_ALL} default 0;

#Parameters with the weigths of each part of the objective function
param WOri;	#Weight of the original portion of the TSS, which tries to maximize vaccine allocations taking care not to raise costs too much.
param WSat;	# Weight of the new portion of the model, which sees how much of the demand is satisfied.

param targetGNI default 1036;
param lmGNI default 4085;
param hmGNI default 12615;


maximize TSSME: 
#WOri*(sum{b in B, m in T[b]} R[b,m]*X[b,m]- sum{b in B}C[b]*sum{m in T[b]}X[b,m]/(S[b]+1e-10)) - WSat*(sum{a in A, m in M}(d[a,m]*l[m] - sum{b in B1[a]}X[b,m]))- sum{b in B}g[b]; #Added demand
WOri*(sum{b in B, m in T[b]} R[b,m]*X[b,m]- sum{b in B}C[b]*sum{m in T[b]}X[b,m]/(S[b]+1e-10)) - WSat*(sum{a in A, m in M: gni_p[m]<targetGNI}(d[a,m]*l[m] - sum{b in B1[a]}X[b,m]))- sum{b in B}g[b]; #Tweight-TS
#WOri*(sum{b in B, m in T[b]} R[b,m]*X[b,m]- sum{b in B}C[b]*sum{m in T[b]}X[b,m]/(S[b]+1e-10)) - WSat*100*(sum{a in A, m in M: gni_p[m]<targetGNI}(d[a,m]*l[m] - sum{b in B1[a]}X[b,m])) - WSat*(sum{a in A, m in M: gni_p[m]>targetGNI}(d[a,m]*l[m] - sum{b in B1[a]}X[b,m]))- sum{b in B}g[b]; #Tweight-TSG

maximize TCSME: sum{b in B, m in T[b]} (R[b,m]- Y[b,m])*X[b,m];
#maximize TCSME: sum{b in B, m in T[b]: gni_p[m] < targetGNI} (R[b,m]- Y[b,m])*X[b,m] + .01*sum{b in B, m in T[b]: gni_p[m] >= targetGNI} (R[b,m]- Y[b,m])*X[b,m];	#TWeight-TC

maximize TPFME: sum{b in B, m in T[b]}Y[b,m]*X[b,m] - sum{b in B} C[b]*g[b]*sum{m in T[b]}X[b,m]/(S[b]+1e-10);

subject to SuperAdditivity2{b in B, m in T[b], q in Q[b]: card(Q[b])>=1}:  
	Y[b,m] + (1-g[b])*sum{t in N[q]}R[t,m] >= sum{t in N[q]}(Y[t,m]);

subject to CapacityME{b in B}: sum{m in T[b]} X[b,m] <= S[b]*g[b];

subject to GenericPrice{b in B, m in T[b]}: Y[b,m] >= g[b]*minPrice[b];	#Should only count when X[b,m]>0

subject to GenericPrice2{b in B, m in T[b]: X2[b,m]>0}: Y[b,m] >= g[b]*minPrice[b];	#Should only count when X[b,m]>0

subject to PriceElasticity{b in B, m in T[b]:S[b]>0}:
	Y[b,m] >= ((S[b]*g[b]-X[b,m])*alpha+S[b]*g[b])/S[b] * (C[b]/S[b]);		#Look for better elasticity relationships?	min when X = supply; 

subject to ReservationPriceLimME{b in B, m in T[b]}: 
	Y[b,m] <= theta*R[b,m]*g[b];

subject to TargetedSale{b in B}:
	sum{m in T[b]}X[b,m] >= sum{m in M}X[b,m];

subject to MaxBDose{b in B, m in T[b]}:
	X[b,m] <= D[b,m]*l[m];
