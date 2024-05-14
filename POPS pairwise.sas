libname mvp 'C:\Users\u0118563\OneDrive - KU Leuven\Projecten\MV probit\Data\POPS';

/*independence*/
data binary;
set mvp.popslda;
bil100=BIL/100;
run;
proc logistic data=binary ;
  model ABIL1 = COV CGM bil100 /link=probit;
run;
proc logistic data=binary ;
  model ABIL2 = COV CGM bil100 /link=probit;
run;
proc logistic data=binary ;
  model ABIL3 = COV CGM bil100 /link=probit;
run;


/*pseudo-likelihood: complete cases*/
data POPS;
set mvp.popslda ;
response123=100*ABIL1+10*ABIL2+ABIL3;
WHERE NOT MISSING (BIL) & NOT MISSING (COV) & NOT MISSING (CGM) ;
run; 

data pops;
set pops;
where not missing(response123);
run;
**ABIL1 and ABIL2**;
proc nlmixed data=POPS qpoints=5 maxiter=10000 HESSIAN SUBGRADIENT=grad_12
maxfunc=10000 technique=quanew;
parms
int1	2.03
beta_ns1	-1.14
beta_cm1	-0.59
beta_b1	0
int2	2.19
beta_ns2	-1.25
beta_cm2	-0.55
beta_b2	0
r12	0.72
;
eta1 = int1+beta_ns1*COV+beta_cm1*CGM+beta_b1*BIL/100;
eta2 = int2+beta_ns2*COV+beta_cm2*CGM+beta_b2*BIL/100;
if ABIL1=0 and ABIL2=0 then do;
lik = probbnrm(eta1,eta2,r12);
ll_12=log(lik);
end;
if ABIL1=0 and ABIL2=1 then do;
lik = cdf('NORMAL',(eta1))-probbnrm(eta1,eta2,r12);
ll_12=log(lik);
end;
if ABIL1=1 and ABIL2=0 then do;
lik = cdf('NORMAL',(eta2))-probbnrm(eta1,eta2,r12);
ll_12=log(lik);
end;
if ABIL1=1 and ABIL2=1 then do;
lik = 1-
(cdf('NORMAL',(eta1))-probbnrm(eta1,eta2,r12))-
(cdf('NORMAL',(eta2))-probbnrm(eta1,eta2,r12))-
probbnrm(eta1,eta2,r12);
ll_12=log(lik);
end;
ll=(ll_12);
model response123 ~ general(ll);
ods output hessian=hess_12 parameterestimates=parms_12;
run;
**ABIL1 and ABIL3**;
proc nlmixed data=POPS qpoints=5 maxiter=10000 HESSIAN SUBGRADIENT=grad_13
maxfunc=10000 technique=quanew;
parms
int1	2.03
beta_ns1	-1.14
beta_cm1	-0.59
beta_b1	0
int3	1.84
beta_ns3	-0.92
beta_cm3	-0.5
beta_b3	0
r13	0.82
;
eta1 = int1+beta_ns1*COV+beta_cm1*CGM+beta_b1*BIL/100;
eta3 = int3+beta_ns3*COV+beta_cm3*CGM+beta_b3*BIL/100;
if ABIL1=0 and ABIL3=0 then do;
lik = probbnrm(eta1,eta3,r13);
ll_13=log(lik);
end;
if ABIL1=0 and ABIL3=1 then do;
lik= cdf('NORMAL',(eta1))-probbnrm(eta1,eta3,r13);
ll_13=log(lik);
end;
if ABIL1=1 and ABIL3=0 then do;
lik = cdf('NORMAL',(eta3))-probbnrm(eta1,eta3,r13);
ll_13=log(lik);
end;
if ABIL1=1 and ABIL3=1 then do;
lik = 1-
(cdf('NORMAL',(eta1))-probbnrm(eta1,eta3,r13))-
(cdf('NORMAL',(eta3))-probbnrm(eta1,eta3,r13))-
probbnrm(eta1,eta3,r13);
ll_13=log(lik);
end;
ll=ll_13;
model response123 ~ general(ll);
ods output hessian=hess_13 parameterestimates=parms_13;
run;
**ABIL2 and ABIL3**;
proc nlmixed data=POPS qpoints=5 maxiter=10000 HESSIAN SUBGRADIENT=grad_23
maxfunc=10000 technique=quanew;
parms
int2	2.19
beta_ns2	-1.25
beta_cm2	-0.55
beta_b2	0
int3	1.84
beta_ns3	-0.92
beta_cm3	-0.5
beta_b3	0
r23	0.7304
;
eta2 = int2+beta_ns2*COV+beta_cm2*CGM+beta_b2*BIL/100;
eta3 = int3+beta_ns3*COV+beta_cm3*CGM+beta_b3*BIL/100;
if ABIL2=0 and ABIL3=0 then do;
lik = probbnrm(eta2,eta3,r23);
ll_23 = log(lik);
end;
if ABIL2=0 and ABIL3=1 then do;
lik = cdf('NORMAL',(eta2))-probbnrm(eta2,eta3,r23);
ll_23 = log(lik);
end;
if ABIL2=1 and ABIL3=0 then do;
lik = cdf('NORMAL',(eta3))-probbnrm(eta2,eta3,r23);
ll_23 = log(lik);
end;
if ABIL2=1 and ABIL3=1 then do;
lik = 1-
(cdf('NORMAL',(eta2))-probbnrm(eta2,eta3,r23))-
(cdf('NORMAL',(eta3))-probbnrm(eta2,eta3,r23))-
probbnrm(eta2,eta3,r23);
ll_23 = log(lik);
end;
ll=(ll_23);
model response123 ~ general(ll);
ods output hessian=hess_23 parameterestimates=parms_23;
run;

*macro for the combining the pairwise estimates estimator;
%macro combine;
/*drop unnessecary columns in the hessians, gradients, and parameter estimates*/
data H_1;
set hess_12;
drop row Parameter;
run;
data G_1;
set grad_12;
drop observation;
run;
data H_2;
set hess_13;
drop row Parameter;
run;
data G_2;
set grad_13;
drop observation;
run;
data H_3;
set hess_23;
drop row Parameter;
run;
data G_3;
set grad_23;
drop observation;
run;
data P_1;
set parms_12;
estimate_12=estimate; 
keep Parameter estimate_12;
run;
data P_2;
set parms_13;
estimate_13=estimate; 
keep Parameter estimate_13;
run;
data P_3;
set parms_23;
estimate_23=estimate; 
keep Parameter estimate_23;
run;

/*combine the hessians and gradients*/
proc iml;
use H_1; 
read all into H_1; 
close H_1;
use H_2; 
read all into H_2; 
close H_2;
use H_3;
read all into H_3; 
close H_3;
H_comb=block(H_1, H_2, H_3);
create H from H_comb;
append from H_comb;
use G_1; 
read all into G_1; 
close G_1;use G_2; 
read all into G_2; 
close G_2;
use G_3; 
read all into G_3;
close G_3;
G_comb=t(G_1)//t(G_2)//t(G_3);
create G from G_comb;
append from G_comb;
quit; 

**create J and K**;
proc iml;
use G; 
read all into G; 
close G;
nsubjects=ncol(G);
K_1=(G)*t(G);
K=K_1#1/nsubjects;
create K from K;
append from K;
use H; 
read all into H; 
close H;
J=H#1/nsubjects;
create J from J; 
append from J;
quit;

*estimate sigma**;
proc iml;
use J; 
read all into J; 
close J;
use K; 
read all into K;
close K;
use G; 
read all into G; 
close G;
nsubjects=ncol(G);
Sigma=inv(J)*K*inv(J);
Sigma0=Sigma#1/nsubjects;
create Sigma0 from Sigma0; 
append from Sigma0;
var_pm = (vecdiag(Sigma0)); 
create var_pm from var_pm;append from var_pm;
quit; 

**create a list with the variables names in the same ordering as the sigma matrix**;
data parms;
set parms_12(keep=parameter) parms_13(keep=parameter) parms_23(keep=parameter);
run;

**calculate variance of the parameters**;
data result_se;
merge parms var_pm;
run; 

proc sort data=result_se;
by parameter;
run;

proc means data=result_se mean noprint;
var COL1;
by parameter;
output out=variances_pm;
run;

data variances_pm;
set variances_pm;
SE=sqrt(col1);
where _stat_='MEAN';
keep parameter SE;
run;


**calculate parameter estimates**;
proc sort data = P_1;
by Parameter;
run;
proc sort data = P_2;
by Parameter;
run;
proc sort data = P_3;
by Parameter;
run;

data parameters;
merge P_1 P_2 P_3;
by Parameter;
run;

data result;
set parameters;
estimate= mean(estimate_12,estimate_13,estimate_23);
keep parameter estimate;
run;


/*combine the results*/
proc sort data=variances_pm;
by parameter;
run;

proc sort data=result;
by parameter;
run;

data final;
merge result variances_pm;
by parameter;
run;
proc print data=final;
run;
%mend combine;
%combine;

PROC SORT DATA= + FINAL


/*pseudolikelihood:complete pairs*/;
data POPS;
set mvp.popslda ;
response123=-1000000;
WHERE NOT MISSING (BIL) & NOT MISSING (COV) & NOT MISSING (CGM) &
NOT (missing(ABIL1) & missing(abil2) & missing(abil3));
run; 
**ABIL1 and ABIL2**;
proc nlmixed data=POPS qpoints=5 maxiter=10000 HESSIAN SUBGRADIENT=grad_12
maxfunc=10000 technique=quanew;
parms
int1	2.03
beta_ns1	-1.14
beta_cm1	-0.59
beta_b1	0
int2	2.19
beta_ns2	-1.25
beta_cm2	-0.55
beta_b2	0
r12	0.72
;
eta1 = int1+beta_ns1*COV+beta_cm1*CGM+beta_b1*BIL/100;
eta2 = int2+beta_ns2*COV+beta_cm2*CGM+beta_b2*BIL/100;
if ABIL1=0 and ABIL2=0 then do;
lik = probbnrm(eta1,eta2,r12);
ll_12=log(lik);
end;
if ABIL1=0 and ABIL2=1 then do;
lik = cdf('NORMAL',(eta1))-probbnrm(eta1,eta2,r12);
ll_12=log(lik);
end;
if ABIL1=1 and ABIL2=0 then do;
lik = cdf('NORMAL',(eta2))-probbnrm(eta1,eta2,r12);
ll_12=log(lik);
end;
if ABIL1=1 and ABIL2=1 then do;
lik = 1-
(cdf('NORMAL',(eta1))-probbnrm(eta1,eta2,r12))-
(cdf('NORMAL',(eta2))-probbnrm(eta1,eta2,r12))-
probbnrm(eta1,eta2,r12);
ll_12=log(lik);
end;
if ABIL1=. and ABIL2=0 then do;
lik = CDF('NORMAL',eta2);
ll_12=0;
end;
if ABIL1=. and ABIL2=1 then do;
lik = 1-CDF('NORMAL',eta2);
ll_12=0;
end;
if ABIL1=0 and ABIL2=. then do;
lik = CDF('NORMAL',eta1);
ll_12=0;
end;
if ABIL1=1 and ABIL2=. then do;
lik = 1-CDF('NORMAL',eta1);
ll_12=0;
end;
model response123 ~ general(ll_12);
ods output hessian=hess_12 parameterestimates=parms_12;
run;
**ABIL1 and ABIL3**;
proc nlmixed data=POPS qpoints=5 maxiter=10000 HESSIAN SUBGRADIENT=grad_13
maxfunc=10000 technique=quanew;
parms
int1	2.03
beta_ns1	-1.14
beta_cm1	-0.59
beta_b1	0
int3	1.84
beta_ns3	-0.92
beta_cm3	-0.5
beta_b3	0
r13	0.82
;
eta1 = int1+beta_ns1*COV+beta_cm1*CGM+beta_b1*BIL/100;
eta3 = int3+beta_ns3*COV+beta_cm3*CGM+beta_b3*BIL/100;
if ABIL1=0 and ABIL3=0 then do;
lik = probbnrm(eta1,eta3,r13);
ll_13=log(lik);
end;
if ABIL1=0 and ABIL3=1 then do;
lik= cdf('NORMAL',(eta1))-probbnrm(eta1,eta3,r13);
ll_13=log(lik);
end;
if ABIL1=1 and ABIL3=0 then do;
lik = cdf('NORMAL',(eta3))-probbnrm(eta1,eta3,r13);
ll_13=log(lik);
end;
if ABIL1=1 and ABIL3=1 then do;
lik = 1-
(cdf('NORMAL',(eta1))-probbnrm(eta1,eta3,r13))-
(cdf('NORMAL',(eta3))-probbnrm(eta1,eta3,r13))-
probbnrm(eta1,eta3,r13);
ll_13=log(lik);
end;
/* observations with missing values*/
if ABIL1=. and ABIL3=0 then do;
lik = CDF('NORMAL',eta3);
ll_13=0;
end;
if ABIL1=. and ABIL3=1 then do;
lik= 1-CDF('NORMAL',eta3);
ll_13 = 0;
end;
if ABIL1=0 and ABIL3=. then do;
lik = CDF('NORMAL',eta1);
ll_13=0;
end;
if ABIL1=1 and ABIL3=. then do;
lik = 1-CDF('NORMAL',eta1);
ll_13=0;
end;
model response123 ~ general(ll_13);
ods output hessian=hess_13 parameterestimates=parms_13;
run;
**ABIL2 and ABIL3**;
proc nlmixed data=POPS qpoints=5 maxiter=10000 HESSIAN SUBGRADIENT=grad_23
maxfunc=10000 technique=quanew;
parms
int2	2.19
beta_ns2	-1.25
beta_cm2	-0.55
beta_b2	0
int3	1.84
beta_ns3	-0.92
beta_cm3	-0.5
beta_b3	0
r23	0.7304
;
eta2 = int2+beta_ns2*COV+beta_cm2*CGM+beta_b2*BIL/100;
eta3 = int3+beta_ns3*COV+beta_cm3*CGM+beta_b3*BIL/100;
if ABIL2=0 and ABIL3=0 then do;
lik = probbnrm(eta2,eta3,r23);
ll_23 = log(lik);
end;
if ABIL2=0 and ABIL3=1 then do;
lik = cdf('NORMAL',(eta2))-probbnrm(eta2,eta3,r23);
ll_23 = log(lik);
end;
if ABIL2=1 and ABIL3=0 then do;
lik = cdf('NORMAL',(eta3))-probbnrm(eta2,eta3,r23);
ll_23 = log(lik);
end;
if ABIL2=1 and ABIL3=1 then do;
lik = 1-
(cdf('NORMAL',(eta2))-probbnrm(eta2,eta3,r23))-
(cdf('NORMAL',(eta3))-probbnrm(eta2,eta3,r23))-
probbnrm(eta2,eta3,r23);
ll_23 = log(lik);
end;
/* observations with missing values*/
if ABIL2=. and ABIL3=0 then do;
lik = CDF('NORMAL',eta3);
ll_23 = 0;
end;
if ABIL2=. and ABIL3=1 then do;
lik = 1-CDF('NORMAL',eta3);
ll_23 = 0;
end;
if ABIL2=0 and ABIL3=. then do;
lik = CDF('NORMAL',eta2);
ll_23 = 0;
end;
if ABIL2=1 and ABIL3=. then do;
lik = 1-CDF('NORMAL',eta2);
ll_23 = 0;
end;
model response123 ~ general(ll_23);
ods output hessian=hess_23 parameterestimates=parms_23;
run;
%combine;



/*pseudolikelihood:available cases*/;
data POPS;
set mvp.popslda ;
response123=-1000000;
WHERE NOT MISSING (BIL) & NOT MISSING (COV) & NOT MISSING (CGM) &
NOT (missing(ABIL1) & missing(abil2) & missing(abil3));
run; 

**ABIL1 and ABIL2**;
proc nlmixed data=POPS qpoints=5 maxiter=10000 HESSIAN SUBGRADIENT=grad_av
maxfunc=10000 technique=quanew;
parms
int1	2.03
beta_ns1	-1.14
beta_cm1	-0.59
beta_b1	0
int2	2.19
beta_ns2	-1.25
beta_cm2	-0.55
beta_b2	0
r12	0.72
;
eta1 = int1+beta_ns1*COV+beta_cm1*CGM+beta_b1*BIL/100;
eta2 = int2+beta_ns2*COV+beta_cm2*CGM+beta_b2*BIL/100;
if ABIL1=0 and ABIL2=0 then do;
lik = probbnrm(eta1,eta2,r12);
ll_12=log(lik);
end;
if ABIL1=0 and ABIL2=1 then do;
lik = cdf('NORMAL',(eta1))-probbnrm(eta1,eta2,r12);
ll_12=log(lik);
end;
if ABIL1=1 and ABIL2=0 then do;
lik = cdf('NORMAL',(eta2))-probbnrm(eta1,eta2,r12);
ll_12=log(lik);
end;
if ABIL1=1 and ABIL2=1 then do;
lik = 1-
(cdf('NORMAL',(eta1))-probbnrm(eta1,eta2,r12))-
(cdf('NORMAL',(eta2))-probbnrm(eta1,eta2,r12))-
probbnrm(eta1,eta2,r12);
ll_12=log(lik);
end;
/* observations with missing values*/
if ABIL1=. and ABIL2=0 then do;
lik = CDF('NORMAL',eta2);
ll_12=log(lik);
end;
if ABIL1=. and ABIL2=1 then do;
lik = 1-CDF('NORMAL',eta2);
ll_12=log(lik);
end;
if ABIL1=0 and ABIL2=. then do;
lik = CDF('NORMAL',eta1);
ll_12=log(lik);
end;
if ABIL1=1 and ABIL2=. then do;
lik = 1-CDF('NORMAL',eta1);
ll_12=log(lik);
end;
model response123 ~ general(ll_12);
ods output hessian=hess_12 parameterestimates=parms_12;
run;
**ABIL1 and ABIL3**;
proc nlmixed data=POPS qpoints=5 maxiter=10000 HESSIAN SUBGRADIENT=grad_av
maxfunc=10000 technique=quanew;
parms
int1	2.03
beta_ns1	-1.14
beta_cm1	-0.59
beta_b1	0
int3	1.84
beta_ns3	-0.92
beta_cm3	-0.5
beta_b3	0
r13	0.82
;
eta1 = int1+beta_ns1*COV+beta_cm1*CGM+beta_b1*BIL/100;
eta3 = int3+beta_ns3*COV+beta_cm3*CGM+beta_b3*BIL/100;
if ABIL1=0 and ABIL3=0 then do;
lik = probbnrm(eta1,eta3,r13);
ll_13=log(lik);
end;
if ABIL1=0 and ABIL3=1 then do;
lik= cdf('NORMAL',(eta1))-probbnrm(eta1,eta3,r13);
ll_13=log(lik);
end;
if ABIL1=1 and ABIL3=0 then do;
lik = cdf('NORMAL',(eta3))-probbnrm(eta1,eta3,r13);
ll_13=log(lik);
end;
if ABIL1=1 and ABIL3=1 then do;
lik = 1-
(cdf('NORMAL',(eta1))-probbnrm(eta1,eta3,r13))-
(cdf('NORMAL',(eta3))-probbnrm(eta1,eta3,r13))-
probbnrm(eta1,eta3,r13);
ll_13=log(lik);
end;
/* observations with missing values*/
if ABIL1=. and ABIL3=0 then do;
lik = CDF('NORMAL',eta3);
ll_13=log(lik);
end;
if ABIL1=. and ABIL3=1 then do;
lik= 1-CDF('NORMAL',eta3);
ll_13 = log(lik);
end;
if ABIL1=0 and ABIL3=. then do;
lik = CDF('NORMAL',eta1);
ll_13=log(lik);
end;
if ABIL1=1 and ABIL3=. then do;
lik = 1-CDF('NORMAL',eta1);
ll_13=log(lik);
end;
model response123 ~ general(ll_13);
ods output hessian=hess_13 parameterestimates=parms_13;
run;
**ABIL2 and ABIL3**;
proc nlmixed data=POPS qpoints=5 maxiter=10000 HESSIAN SUBGRADIENT=grad_23
maxfunc=10000 technique=quanew;
parms
int2	2.19
beta_ns2	-1.25
beta_cm2	-0.55
beta_b2	0
int3	1.84
beta_ns3	-0.92
beta_cm3	-0.5
beta_b3	0
r23	0.7304
;
eta2 = int2+beta_ns2*COV+beta_cm2*CGM+beta_b2*BIL/100;
eta3 = int3+beta_ns3*COV+beta_cm3*CGM+beta_b3*BIL/100;
if ABIL2=0 and ABIL3=0 then do;
lik = probbnrm(eta2,eta3,r23);
ll_23 = log(lik);
end;
if ABIL2=0 and ABIL3=1 then do;
lik = cdf('NORMAL',(eta2))-probbnrm(eta2,eta3,r23);
ll_23 = log(lik);
end;
if ABIL2=1 and ABIL3=0 then do;
lik = cdf('NORMAL',(eta3))-probbnrm(eta2,eta3,r23);
ll_23 = log(lik);
end;
if ABIL2=1 and ABIL3=1 then do;
lik = 1-
(cdf('NORMAL',(eta2))-probbnrm(eta2,eta3,r23))-
(cdf('NORMAL',(eta3))-probbnrm(eta2,eta3,r23))-
probbnrm(eta2,eta3,r23);
ll_23 = log(lik);
end;
/* observations with missing values*/
if ABIL2=. and ABIL3=0 then do;
lik = CDF('NORMAL',eta3);
ll_23 = log(lik);
end;
if ABIL2=. and ABIL3=1 then do;
lik = 1-CDF('NORMAL',eta3);
ll_23 = log(lik);
end;
if ABIL2=0 and ABIL3=. then do;
lik = CDF('NORMAL',eta2);
ll_23 = log(lik);
end;
if ABIL2=1 and ABIL3=. then do;
lik = 1-CDF('NORMAL',eta2);
ll_23 = log(lik);
end;
model response123 ~ general(ll_23);
ods output hessian=hess_23 parameterestimates=parms_23;
run;

%combine;




/*Pseudolikelihood after multiple imputation*/
/*imputation step*/
proc mi data=mvp.popslda out=mvp.popslda_imp seed=123 nimpute=10;
class ABIL1 ABIL2 ABIL3 COV CGM;
fcs discrim(ABIL1 ABIL2 ABIL3 COV CGM);
fcs reg(BIL);
var ABIL1 ABIL2 ABIL3 COV CGM BIL;
run;

/*sequence of parameters*/
data parms_seq;
input Parameter $;
datalines;
int1
beta_ns1
beta_cm1
beta_b1
int2
beta_ns2
beta_cm2
beta_b2
r12
int3
beta_ns3
beta_cm3
beta_b3
r13
r23
;run;

data nlcov;
run;
data parms_imput;
run;
/*analysis step*/
data mvp.popslda_imp;
set mvp.popslda_imp;
response123=-50000;
run;
proc sort data=mvp.popslda_imp; by _imputation_;run;

%macro pairwise_imput(i=);
proc nlmixed data=mvp.popslda_imp qpoints=5 maxiter=10000 HESSIAN SUBGRADIENT=grad_12
maxfunc=10000 technique=quanew;
parms
int1	2.03
beta_ns1	-1.14
beta_cm1	-0.59
beta_b1	0
int2	2.19
beta_ns2	-1.25
beta_cm2	-0.55
beta_b2	0
r12	0.72
;
eta1 = int1+beta_ns1*COV+beta_cm1*CGM+beta_b1*BIL/100;
eta2 = int2+beta_ns2*COV+beta_cm2*CGM+beta_b2*BIL/100;
if ABIL1=0 and ABIL2=0 then do;
lik = probbnrm(eta1,eta2,r12);
ll_12=log(lik);
end;
if ABIL1=0 and ABIL2=1 then do;
lik = cdf('NORMAL',(eta1))-probbnrm(eta1,eta2,r12);
ll_12=log(lik);
end;
if ABIL1=1 and ABIL2=0 then do;
lik = cdf('NORMAL',(eta2))-probbnrm(eta1,eta2,r12);
ll_12=log(lik);
end;
if ABIL1=1 and ABIL2=1 then do;
lik = 1-
(cdf('NORMAL',(eta1))-probbnrm(eta1,eta2,r12))-
(cdf('NORMAL',(eta2))-probbnrm(eta1,eta2,r12))-
probbnrm(eta1,eta2,r12);
ll_12=log(lik);
end;
ll=(ll_12);
model response123 ~ general(ll);
ods output hessian=hess_12 parameterestimates=parms_12;
where _imputation_=&i;
run;
proc nlmixed data=mvp.popslda_imp qpoints=5 maxiter=10000 HESSIAN SUBGRADIENT=grad_13
maxfunc=10000 technique=quanew;
parms
int1	2.03
beta_ns1	-1.14
beta_cm1	-0.59
beta_b1	0
int3	1.84
beta_ns3	-0.92
beta_cm3	-0.5
beta_b3	0
r13	0.82
;
eta1 = int1+beta_ns1*COV+beta_cm1*CGM+beta_b1*BIL/100;
eta3 = int3+beta_ns3*COV+beta_cm3*CGM+beta_b3*BIL/100;
if ABIL1=0 and ABIL3=0 then do;
lik = probbnrm(eta1,eta3,r13);
ll_13=log(lik);
end;
if ABIL1=0 and ABIL3=1 then do;
lik= cdf('NORMAL',(eta1))-probbnrm(eta1,eta3,r13);
ll_13=log(lik);
end;
if ABIL1=1 and ABIL3=0 then do;
lik = cdf('NORMAL',(eta3))-probbnrm(eta1,eta3,r13);
ll_13=log(lik);
end;
if ABIL1=1 and ABIL3=1 then do;
lik = 1-
(cdf('NORMAL',(eta1))-probbnrm(eta1,eta3,r13))-
(cdf('NORMAL',(eta3))-probbnrm(eta1,eta3,r13))-
probbnrm(eta1,eta3,r13);
ll_13=log(lik);
end;
ll=ll_13;
model response123 ~ general(ll);
ods output hessian=hess_13 parameterestimates=parms_13;
where _imputation_=&i;
run;
**ABIL2 and ABIL3**;
proc nlmixed data=mvp.popslda_imp qpoints=5 maxiter=10000 HESSIAN SUBGRADIENT=grad_23
maxfunc=10000 technique=quanew;
parms
int2	2.19
beta_ns2	-1.25
beta_cm2	-0.55
beta_b2	0
int3	1.84
beta_ns3	-0.92
beta_cm3	-0.5
beta_b3	0
r23	0.7304
;
eta2 = int2+beta_ns2*COV+beta_cm2*CGM+beta_b2*BIL/100;
eta3 = int3+beta_ns3*COV+beta_cm3*CGM+beta_b3*BIL/100;
if ABIL2=0 and ABIL3=0 then do;
lik = probbnrm(eta2,eta3,r23);
ll_23 = log(lik);
end;
if ABIL2=0 and ABIL3=1 then do;
lik = cdf('NORMAL',(eta2))-probbnrm(eta2,eta3,r23);
ll_23 = log(lik);
end;
if ABIL2=1 and ABIL3=0 then do;
lik = cdf('NORMAL',(eta3))-probbnrm(eta2,eta3,r23);
ll_23 = log(lik);
end;
if ABIL2=1 and ABIL3=1 then do;
lik = 1-
(cdf('NORMAL',(eta2))-probbnrm(eta2,eta3,r23))-
(cdf('NORMAL',(eta3))-probbnrm(eta2,eta3,r23))-
probbnrm(eta2,eta3,r23);
ll_23 = log(lik);
end;
ll=(ll_23);
model response123 ~ general(ll);
ods output hessian=hess_23 parameterestimates=parms_23;
where _imputation_=&i;
run;
/*drop unnessecary columns in the hessians, gradients, and parameter estimates*/
data H_1;
set hess_12;
drop row Parameter;
run;
data G_1;
set grad_12;
drop observation;
run;
data H_2;
set hess_13;
drop row Parameter;
run;
data G_2;
set grad_13;
drop observation;
run;
data H_3;
set hess_23;
drop row Parameter;
run;
data G_3;
set grad_23;
drop observation;
run;
/*combine the hessians and gradients*/
proc iml;
use H_1; 
read all into H_1; 
close H_1;
use H_2; 
read all into H_2; 
close H_2;
use H_3;
read all into H_3; 
close H_3;
H_comb=block(H_1, H_2, H_3);
create H from H_comb;
append from H_comb;
use G_1; 
read all into G_1; 
close G_1;use G_2; 
read all into G_2; 
close G_2;
use G_3; 
read all into G_3;
close G_3;
G_comb=t(G_1)//t(G_2)//t(G_3);
create G from G_comb;
append from G_comb;
quit; 

**create J and K**;
proc iml;
use G; 
read all into G; 
close G;
nsubjects=ncol(G);
K_1=(G)*t(G);
K=K_1#1/nsubjects;
create K from K;
append from K;
use H; 
read all into H; 
close H;
J=H#1/nsubjects;
create J from J; 
append from J;
quit;

*estimate sigma**;
proc iml;
use J; 
read all into J; 
close J;
use K; 
read all into K;
close K;
use G; 
read all into G; 
close G;
nsubjects=ncol(G);
Sigma=inv(J)*K*inv(J);
Sigma0=Sigma#1/nsubjects;
create Sigma0 from Sigma0; 
append from Sigma0;
quit; 

**parameter estimates;
data theta_est;
set parms_12(keep=estimate) parms_13(keep=estimate) parms_23(keep=estimate);
run;
data theta_est;
set parms_12(keep=parameter estimate) parms_13(keep=parameter estimate) parms_23(keep=parameter estimate);
run;

**create the A matrix**;
proc iml;
A={0.5	0	0	0	0	0	0	0	0	0.5	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	,
0	0.5	0	0	0	0	0	0	0	0	0.5	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	,
0	0	0.5	0	0	0	0	0	0	0	0	0.5	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	,
0	0	0	0.5	0	0	0	0	0	0	0	0	0.5	0	0	0	0	0	0	0	0	0	0	0	0	0	0	,
0	0	0	0	0.5	0	0	0	0	0	0	0	0	0	0	0	0	0	0.5	0	0	0	0	0	0	0	0	,
0	0	0	0	0	0.5	0	0	0	0	0	0	0	0	0	0	0	0	0	0.5	0	0	0	0	0	0	0	,
0	0	0	0	0	0	0.5	0	0	0	0	0	0	0	0	0	0	0	0	0	0.5	0	0	0	0	0	0	,
0	0	0	0	0	0	0	0.5	0	0	0	0	0	0	0	0	0	0	0	0	0	0.5	0	0	0	0	0	,
0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	,
0	0	0	0	0	0	0	0	0	0	0	0	0	0.5	0	0	0	0	0	0	0	0	0.5	0	0	0	0	,
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.5	0	0	0	0	0	0	0	0	0.5	0	0	0	,
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.5	0	0	0	0	0	0	0	0	0.5	0	0	,
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.5	0	0	0	0	0	0	0	0	0.5	0	,
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	,
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	};
use Sigma0; read all into Sigma0; 
close Sigma0;
use Theta_est; 
read all into Theta_est;
close Theta_est;
Theta_est_star=A*Theta_est;
covb=(A*Sigma0*t(A));
create Theta_est_star from Theta_est_star; 
append from Theta_est_star;
create covb from covb;
append from covb;
quit; 


data nlcovb2;
merge parms_seq covb;
_imputation_=&i;
Row=_N_;
run;

data parms_imput2;
merge parms_seq Theta_est_star;
_imputation_=&i;
Row=_N_;
run;

data nlcov;
set nlcov nlcovb2;
run;
data parms_imput;
set parms_imput parms_imput2;
run;
proc datasets nolist lib=work;
delete parms_12 parms_13 parms_23
hess_12 hess_13 hess_23
grad_12 grad_13 grad_23
H_1 G_1 H_2 G_2 H_3 G_3 H G 
K J SIGMA0 theta_est theta_est_star
nlcovb2 parms_imput2 covb;
quit;
%mend pairwise_imput;

%pairwise_imput(i=1);
%pairwise_imput(i=2);
%pairwise_imput(i=3);
%pairwise_imput(i=4);
%pairwise_imput(i=5);
%pairwise_imput(i=6);
%pairwise_imput(i=7);
%pairwise_imput(i=8);
%pairwise_imput(i=9);
%pairwise_imput(i=10);

data parms_imput;
set parms_imput;
rename
COL1=Estimate;
where row is not missing;
drop row;
run;

data cov_imput;
set nlcov;
rename
COL1=int1
COL2=beta_ns1
COL3=beta_cm1
COL4=beta_b1
COL5=int2
COL6=beta_ns2
COL7=beta_cm2
COL8=beta_b2
COL9=r12
COL10=int3
COL11=beta_ns3
COL12=beta_cm3
COL13=beta_b3
COL14=r13
COL15=r23;
where row is not missing;
run;

/*combine the imputed datasets*/
proc mianalyze parms=parms_imput covb=cov_imput;
modeleffects
int1
beta_ns1
beta_cm1
beta_b1
int2
beta_ns2
beta_cm2
beta_b2
int3
beta_ns3
beta_cm3
beta_b3
r12
r13
r23;
ods output ParameterEstimates=result_imput;
run;

proc sort data=result_imput;
by Parm;
run;
