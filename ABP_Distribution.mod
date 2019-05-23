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

set E = 1..num_ent;	#Entities

param s{M,M};	#Similarity between markets (by geography, for example)

param chi;	#Minimum accepable similarity for 2 markets in an entity.

param zeta;	#Maximum acceptable gni difference wihin an entiy.


var Z{E,M} binary;
var SZ;

maximize Distribution: 
	SZ;

subject to Belonging{m in M}:  
	sum{e in E}Z[e,m] = 1;

subject to MinSZ{e in E}:
	SZ <= sum{m in M, b in B}R[b,m]*Z[e,m];

/*subject to Similarity{e in E, m1 in M, m2 in M}:
	W[e,m1]*s[m1,m2]+W[e,m2]*s[m1,m2] >= chi;

subject to GNIRes: 
	abs(W[e,m1]*gni_p[m1]-E[e,m2]*gni_p[m2])>=zeta;	#Should only count when X[b,m]>0*/

#kEEP A GAP BETWEEN THE ip AND THE INITIAL GAP. <- There might be a different.


