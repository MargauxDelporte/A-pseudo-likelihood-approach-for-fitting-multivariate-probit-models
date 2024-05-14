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

proc nlmixed data=POPS qpoints=5 maxiter=10000 HESSIAN SUBGRADIENT=mvp.grad_complete
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
int3	1.84
beta_ns3	-0.92
beta_cm3	-0.5
beta_b3	0
r12	0.72
r13	0.82
r23	0.7304
;
eta1 = int1+beta_ns1*COV+beta_cm1*CGM+beta_b1*BIL/100;
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
/*sum likelihoods*/
ll=(ll_12)+(ll_23)+(ll_13);
model response123 ~ general(ll);
ods output hessian=mvp.hess_complete parameterestimates=mvp.parms_complete;
run;

*macro for the sandwich estimator;
%macro sandwich(he=,gr=,pa=);
data hessian;
set mvp.&he;
drop row Parameter;
run;
data gradient;
set mvp.&gr;
drop observation;
run;
data parms;
set mvp.&pa;
keep parameter estimate;
run;
proc iml;
use hessian; read all into J; close hessian;
use gradient; read all into G; close gradient;
use parms; read all into P; close parms;
K=t(G)*(G);
Sigma=inv(J)*K*inv(J);
se = sqrt(vecdiag(Sigma)); 
create se from se;append from se;
quit;

data result;
merge parms se;
run;
proc print data=result round;run;
%mend sandwich;
%sandwich(he=hess_complete,gr=grad_complete,pa=parms_complete);
PROC SORT DATA=result;by parameter;run;



/*pseudolikelihood:complete pairs*/
data POPS;
set mvp.popslda ;
response123=-1000000;
WHERE NOT MISSING (BIL) & NOT MISSING (COV) & NOT MISSING (CGM) &
NOT (missing(ABIL1) & missing(abil2) & missing(abil3));
run; 
proc nlmixed data=POPS qpoints=5 maxiter=10000 HESSIAN SUBGRADIENT=mvp.gradient_cp
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
int3	1.84
beta_ns3	-0.92
beta_cm3	-0.5
beta_b3	0
r12	0.72
r13	0.82
r23	0.7304
;
eta1 = int1+beta_ns1*COV+beta_cm1*CGM+beta_b1*BIL/100;
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
/*sum likelihoods*/
ll=(ll_12)+(ll_23)+(ll_13);
model response123 ~ general(ll);
ods output hessian=mvp.hessian_cp parameterestimates=mvp.parms_cp;
run;
%sandwich(he=hessian_cp,gr=gradient_cp,pa=parms_cp);
proc sort data=result;by parameter;run;


/*pseudolikelihood:available cases*/;
data POPS;
set mvp.popslda ;
response123=-1000000;
WHERE NOT MISSING (BIL) & NOT MISSING (COV) & NOT MISSING (CGM)&
NOT (missing(ABIL1) & missing(abil2) & missing(abil3));
run; 

proc nlmixed data=POPS qpoints=5 maxiter=10000 HESSIAN SUBGRADIENT=mvp.grad_av
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
int3	1.84
beta_ns3	-0.92
beta_cm3	-0.5
beta_b3	0
r12	0.72
r13	0.82
r23	0.7304
;
eta1 = int1+beta_ns1*COV+beta_cm1*CGM+beta_b1*BIL/100;
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
/*sum likelihoods*/
ll=(ll_12)+(ll_23)+(ll_13);
model response123 ~ general(ll);
ods output hessian=mvp.hess_av parameterestimates=mvp.parms_av;
run;

%sandwich(he=hess_av,gr=grad_av,pa=parms_av);
proc sort data=result;by parameter;run;



/*Pseudolikelihood after multiple imputation*/
/*imputation step*/
proc mi data=mvp.popslda out=mvp.popslda_imp seed=123 nimpute=10;
class ABIL1 ABIL2 ABIL3 COV CGM;
fcs discrim(ABIL1 ABIL2 ABIL3 COV CGM);
fcs reg(BIL);
var ABIL1 ABIL2 ABIL3 COV CGM BIL;
run;

/*analysis step*/
data mvp.popslda_imp;
set mvp.popslda_imp;
response123=-50000;
run;
proc sort data=mvp.popslda_imp; by _imputation_;run;
proc nlmixed data=mvp.popslda_imp qpoints=5 maxiter=10000 HESSIAN SUBGRADIENT=mvp.grad_imput
technique=newrap cov ecov;
by _imputation_;
parms
int1	2.03
beta_ns1	-1.14
beta_cm1	-0.59
beta_b1	0
int2	2.19
beta_ns2	-1.25
beta_cm2	-0.55
beta_b2	0
int3	1.84
beta_ns3	-0.92
beta_cm3	-0.5
beta_b3	0
r12	0.72
r13	0.82
r23	0.7304
;
eta1 = int1+beta_ns1*COV+beta_cm1*CGM+beta_b1*BIL/100;
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
ll=(ll_12)+(ll_23)+(ll_13);
model response123 ~ general(ll);
ods output ParameterEstimates=mvp.parms_imput
hessian=mvp.hessian_imput CovMatParmEst=nlcovb
AdditionalEstimates=nlparmsa CovMatAddEst=nlcovba;
run;


/*calculate sandwich estimator for each imputed dataset*/
data cov;
run;

%macro sandwich_imput(i=);
data hessian;
set mvp.hessian_imput;
where _imputation_=&i;
drop row Parameter _imputation_;
run;
data gradient;
set mvp.grad_imput;
where _imputation_=&i;
drop observation _imputation_;
run;
proc iml;
use hessian; read all into J; close hessian;
use gradient; read all into G; close gradient;
K=t(G)*(G);
Sigma = inv(J)*K*inv(J);
create Sigma from Sigma;append from Sigma;
quit;
data cov;
set Sigma cov;
run;
data cov;
set cov;
where col1 is not missing;
run;

%mend sandwich_imput;
%sandwich_imput(i=1)
%sandwich_imput(i=2)
%sandwich_imput(i=3)
%sandwich_imput(i=4)
%sandwich_imput(i=5)
%sandwich_imput(i=6)
%sandwich_imput(i=7)
%sandwich_imput(i=8)
%sandwich_imput(i=9)
%sandwich_imput(i=10)

data cov_imput;
merge nlcovb(keep=_imputation_ row parameter) cov;
run;
data cov_imput;
set cov_imput;
rename
COL1=int1
col2=beta_ns1
col3=beta_cm1
col4=beta_b1
col5=int2
col6=beta_ns2
col7=beta_cm2
col8=beta_b2
col9=int3
col10=beta_ns3
col11=beta_cm3
col12=beta_b3
col13=r12
col14=r13
col15=r23;
run;

/*combine the imputed datasets*/
proc mianalyze parms=mvp.parms_imput covb=cov_imput wcov bcov tcov;
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
ods output ParameterEstimates=mvp.result_imput;
run;

proc sort data=mvp.result_imput;by parm;run;
