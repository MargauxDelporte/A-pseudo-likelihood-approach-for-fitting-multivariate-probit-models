libname mvp 'C:\Users\u0118563\OneDrive - KU Leuven\Projecten\MV probit\Data\BIRNH';

****Select complete cases****;
data BIRNH2;
set mvp.BIRNH2 ;
response=100*DRINK+10*FUMKLAS+cholklas;
WHERE NOT MISSING (FUMKLAS) & NOT MISSING (DRINK) & NOT MISSING (gender) & NOT MISSING (STATSOC1) & NOT MISSING (GEWEST1)
& NOT MISSING (BMI2) & NOT MISSING (STATSOC2) &NOT MISSING(AGE) & NOT MISSING(cholklas);
run; 

**variable selection cholklas***;
proc logistic data=birnh2;
class  gender STATSOC1 STATSOC2 GEWEST1;
model cholklas= gender STATSOC1 STATSOC2 BMI2 AGE GEWEST1/link=probit;
run;
proc logistic data=birnh2;
model cholklas= STATSOC1 STATSOC2 BMI2 AGE GEWEST1/link=probit;
run;

*response 1 and 2: FUM-DRINK; 
proc nlmixed data=BIRNH2 qpoints=5 maxiter=10000 maxfunc=10000 technique=newrap HESSIAN SUBGRADIENT=grad_b12;
parms
constant12	0.024201016
r_sex12	-0.08659
r_soc12	0.09297
int_drink1	-2.709749624
int_drink2	-0.231747722
int_drink3	0.156361687
beta_d_gender	-0.548134907
beta_d_statsoc	0.22423817
beta_d_gewest	-0.094606161
int_smoke1	-2.293673314
int_smoke2	-1.659601511
beta_s_gender	-1.238805667
beta_s_bmi	-0.028382134
beta_s_age	-0.014626096
beta_s_statsoc1	-0.157529409
beta_s_statsoc2	-0.254146551
;
y12=(constant12+r_sex12*gender+r_soc12*STATSOC2);
r12=(exp(2 * y12) - 1)/(exp(2 * y12) + 1);
eta_d = beta_d_gender*gender+
beta_d_statsoc*STATSOC1+beta_d_gewest*GEWEST1;
eta_s = beta_s_gender*gender+beta_s_bmi*BMI2+
beta_s_age*AGE+
beta_s_statsoc1*STATSOC1+beta_s_statsoc2*STATSOC2;
/***rij 1***/
/*1*/ if FUMKLAS=1 and DRINK=0 then do;
lik12 = probbnrm(int_smoke1-eta_s,int_drink1-eta_d,r12);
end;
/*2*/if FUMKLAS=2 and DRINK=0 then do;
lik12 = probbnrm(int_smoke2-eta_s,int_drink1-eta_d,r12)-
probbnrm(int_smoke1-eta_s,int_drink1-eta_d,r12);
end;
/*3*/if FUMKLAS=3 and DRINK=0 then do;
lik12 = cdf('NORMAL',(int_drink1-eta_d))-
probbnrm(int_smoke2-eta_s,int_drink1-eta_d,r12);
end;
/***rij 2***/
/*4*/if FUMKLAS=1 and DRINK=1 then do;
lik12 = probbnrm(int_smoke1-eta_s,int_drink2-eta_d,r12)-
probbnrm(int_smoke1-eta_s,int_drink1-eta_d,r12);
end;
/*5*/if FUMKLAS=2 and DRINK=1 then do;
lik12 = probbnrm(int_smoke2-eta_s,int_drink2-eta_d,r12)-
probbnrm(int_smoke2-eta_s,int_drink1-eta_d,r12)-/*4*/
(
probbnrm(int_smoke1-eta_s,int_drink2-eta_d,r12)-
probbnrm(int_smoke1-eta_s,int_drink1-eta_d,r12)
)
;
end;
/*6*/if FUMKLAS=3 and DRINK=1 then do;
lik12 = cdf('NORMAL',(int_drink2-eta_d))- 
cdf('NORMAL',(int_drink1-eta_d))-
probbnrm(int_smoke2-eta_s,int_drink2-eta_d,r12)+
probbnrm(int_smoke2-eta_s,int_drink1-eta_d,r12);
end;
/***rij 3***/
/*7*/if FUMKLAS=1 and DRINK=2 then do;
lik12 = probbnrm(int_smoke1-eta_s,int_drink3-eta_d,r12)-
probbnrm(int_smoke1-eta_s,int_drink2-eta_d,r12);
end;
/*8*/if FUMKLAS=2 and DRINK=2 then do;
lik12 = 
probbnrm(int_smoke2-eta_s,int_drink3-eta_d,r12)-
probbnrm(int_smoke2-eta_s,int_drink2-eta_d,r12)-
probbnrm(int_smoke1-eta_s,int_drink3-eta_d,r12)+
probbnrm(int_smoke1-eta_s,int_drink2-eta_d,r12);
end;
/*9*/if FUMKLAS=3 and DRINK=2 then do;
lik12 = 
cdf('NORMAL',(int_drink3-eta_d))-
cdf('NORMAL',(int_drink2-eta_d))-
probbnrm(int_smoke2-eta_s,int_drink3-eta_d,r12)+
probbnrm(int_smoke2-eta_s,int_drink2-eta_d,r12);
end;
/***rij 4***/
/*10*/if FUMKLAS=1 and DRINK=3 then do;
lik12 = cdf('NORMAL',(int_smoke1-eta_s))-
probbnrm(int_smoke1-eta_s,int_drink3-eta_d,r12);
end;
/*11*/if FUMKLAS=2 and DRINK=3 then do;
lik12 = cdf('NORMAL',(int_smoke2-eta_s))-
probbnrm(int_smoke2-eta_s,int_drink3-eta_d,r12)-
cdf('NORMAL',(int_smoke1-eta_s))+
probbnrm(int_smoke1-eta_s,int_drink3-eta_d,r12);
end;
/*12*/if FUMKLAS=3 and DRINK=3 then do;
lik12 = 1-
cdf('NORMAL',(int_drink3-eta_d))-
cdf('NORMAL',(int_smoke2-eta_s))+
probbnrm(int_smoke2-eta_s,int_drink3-eta_d,r12);
end;
ll=log(lik12);
model response ~ general(ll);
ods output hessian=hessian_b12 parameterestimates=parms_b12;
run;


*response 1 and 3: FUM-CHOLKLAS; 
proc nlmixed data=BIRNH2 qpoints=5 maxiter=10000 maxfunc=10000 technique=newrap HESSIAN SUBGRADIENT=grad_b13;
parms
constant13	-0.03362
r_sex13 -0.08659
r_soc13	0.09297
int_smoke1	-2.293673314
int_smoke2	-1.659601511
beta_s_gender	-1.238805667
beta_s_bmi	-0.028382134
beta_s_age	-0.014626096
beta_s_statsoc1	-0.157529409
beta_s_statsoc2	-0.254146551
int_chol	-2.437889287
beta_c_statsoc1	-0.079429388
beta_c_statsoc2	-0.170497579
beta_c_bmi	-0.035169399
beta_c_age	-0.027805511
beta_c_gewest	0.125759932
;
y13=(constant13+r_sex13*gender+r_soc13*STATSOC2);
r13=(exp(2 * y13) - 1)/(exp(2 * y13) + 1);
eta_s = beta_s_gender*gender+beta_s_bmi*BMI2+
beta_s_age*AGE+
beta_s_statsoc1*STATSOC1+beta_s_statsoc2*STATSOC2;
eta_c = beta_c_bmi*BMI2+
beta_c_age*AGE+beta_c_gewest*GEWEST1+
beta_c_statsoc1*STATSOC1+beta_c_statsoc2*STATSOC2;
/***FUM-CHOLKLAS***/
if FUMKLAS=1 and CHOLKLAS=1 then do;
lik13 = probbnrm(int_smoke1-eta_s,int_chol-eta_c,r13);
end;
if FUMKLAS=2 and CHOLKLAS=1 then do;
lik13 = probbnrm(int_smoke2-eta_s,int_chol-eta_c,r13)-
probbnrm(int_smoke1-eta_s,int_chol-eta_c,r13);;
end;
if FUMKLAS=3 and CHOLKLAS=1 then do;
lik13 = 
cdf('NORMAL',(int_chol-eta_c))-
probbnrm(int_smoke2-eta_s,int_chol-eta_c,r13);
end;
if FUMKLAS=1 and CHOLKLAS=0 then do;
lik13 = cdf('NORMAL',(int_smoke1-eta_s))-
probbnrm(int_smoke1-eta_s,int_chol-eta_c,r13);
end;
if FUMKLAS=2 and CHOLKLAS=0 then do;
lik13 = cdf('NORMAL',(int_smoke2-eta_s))-
cdf('NORMAL',(int_smoke1-eta_s))-
probbnrm(int_smoke2-eta_s,int_chol-eta_c,r13)+
probbnrm(int_smoke1-eta_s,int_chol-eta_c,r13);
end;
if FUMKLAS=3 and CHOLKLAS=0 then do;
lik13 = 1-
cdf('NORMAL',(int_chol-eta_c))-
cdf('NORMAL',(int_smoke2-eta_s))+
probbnrm(int_smoke2-eta_s,int_chol-eta_c,r13);
end;
ll=log(lik13);
model response ~ general(ll);
ods output hessian=hessian_b13 parameterestimates=parms_b13;
run;



*response 2 and 3: DRINK-CHOLKLAS;
proc nlmixed data=BIRNH2 qpoints=5 maxiter=10000 maxfunc=10000 technique=newrap HESSIAN SUBGRADIENT=grad_b23;
parms
constant23	-0.03362
r_sex23	-0.08659
r_soc23	0.09297
int_drink1	-2.709749624
int_drink2	-0.231747722
int_drink3	0.156361687
beta_d_gender	-0.548134907
beta_d_statsoc	0.22423817
beta_d_gewest	-0.094606161
int_chol	-2.437889287
beta_c_statsoc1	-0.079429388
beta_c_statsoc2	-0.170497579
beta_c_bmi	-0.035169399
beta_c_age	-0.027805511
beta_c_gewest	0.125759932
;
y23=(constant23+r_sex23*gender+r_soc23*STATSOC2);
r23=(exp(2 * y23) - 1)/(exp(2 * y23) + 1);
eta_d = beta_d_gender*gender+
beta_d_statsoc*STATSOC1+beta_d_gewest*GEWEST1;
eta_c = beta_c_bmi*BMI2+
beta_c_age*AGE+beta_c_gewest*GEWEST1+
beta_c_statsoc1*STATSOC1+beta_c_statsoc2*STATSOC2;
/***DRINK-CHOLKLAS***/
if DRINK=0 and CHOLKLAS=1 then do;
lik23 = probbnrm(int_drink1-eta_d,int_chol-eta_c,r23);
end;
if DRINK=1 and CHOLKLAS=1 then do;
lik23 = probbnrm(int_drink2-eta_d,int_chol-eta_c,r23)-
probbnrm(int_drink1-eta_d,int_chol-eta_c,r23);
end;
if DRINK=2 and CHOLKLAS=1 then do;
lik23 = probbnrm(int_drink3-eta_d,int_chol-eta_c,r23)-
probbnrm(int_drink2-eta_d,int_chol-eta_c,r23);
end;
if DRINK=3 and CHOLKLAS=1 then do;
lik23 = 
cdf('NORMAL',(int_chol-eta_c))-
probbnrm(int_drink3-eta_d,int_chol-eta_c,r23);
end;
if DRINK=0 and CHOLKLAS=0 then do;
lik23 = cdf('NORMAL',(int_drink1-eta_d))-
probbnrm(int_drink1-eta_d,int_chol-eta_c,r23);
end;
if DRINK=1 and CHOLKLAS=0 then do;
lik23 = cdf('NORMAL',(int_drink2-eta_d))-
cdf('NORMAL',(int_drink1-eta_d))-
probbnrm(int_drink2-eta_d,int_chol-eta_c,r23)+
probbnrm(int_drink1-eta_d,int_chol-eta_c,r23);
end;
if DRINK=2 and CHOLKLAS=0 then do;
lik23 = cdf('NORMAL',(int_drink3-eta_d))-
cdf('NORMAL',(int_drink2-eta_d))-
probbnrm(int_drink3-eta_d,int_chol-eta_c,r23)+
probbnrm(int_drink2-eta_d,int_chol-eta_c,r23);
end;
if DRINK=3 and CHOLKLAS=0 then do;
lik23 = 1-
cdf('NORMAL',(int_chol-eta_c))-
cdf('NORMAL',(int_drink3-eta_d))+
probbnrm(int_drink3-eta_d,int_chol-eta_c,r23);
end;
ll=log(lik23);
model response ~ general(ll);
ods output hessian=hessian_b23 parameterestimates=parms_b23;
run;

/*drop unnessecary columns in the hessians, gradients, and parameter estimates*/
data H_1;
set hessian_b12;
drop row Parameter;
run;
data G_1;
set grad_b12;
drop observation;
run;
data H_2;
set hessian_b13;
drop row Parameter;
run;
data G_2;
set grad_b13;
drop observation;
run;
data H_3;
set hessian_b23;
drop row Parameter;
run;
data G_3;
set grad_b23;
drop observation;
run;
data P_1;
set parms_b12;
estimate_12=estimate; 
keep Parameter estimate_12;
run;
data P_2;
set parms_b13;
estimate_13=estimate; 
keep Parameter estimate_13;
run;
data P_3;
set parms_b23;
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
print nsubjects;
quit;

*estimate sigma**;
proc iml;
use J; 
read all into J; 
close J;
use K; 
read all into K;
close K;
nsubjects=10510;
Sigma=inv(J)*K*inv(J);
Sigma0=Sigma#1/nsubjects;
create Sigma0 from Sigma0; 
append from Sigma0;
var_pm = (vecdiag(Sigma0)); 
create var_pm from var_pm;append from var_pm;
quit; 

**create a list with the variables names in the same ordering as the sigma matrix**;
data parms;
set parms_b12(keep=parameter) parms_b13(keep=parameter) parms_b23(keep=parameter);
run;

**calculate variance of the parameters**;
data result_se;
merge parms var_pm;
run; 

proc sort data=result_se;
by parameter;
run;

proc means data=result_se mean;
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
