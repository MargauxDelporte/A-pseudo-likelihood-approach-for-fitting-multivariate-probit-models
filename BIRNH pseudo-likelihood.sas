libname mvp '';

****Select complete cases****;
data BIRNH2;
set mvp.BIRNH2 ;
response=100*DRINK+10*FUMKLAS+cholklas;
WHERE NOT MISSING (FUMKLAS) & NOT MISSING (DRINK) & NOT MISSING (gender) & NOT MISSING (STATSOC1) & NOT MISSING (GEWEST1)
& NOT MISSING (BMI2) & NOT MISSING (STATSOC2) &NOT MISSING(AGE)& NOT MISSING(cholklas);
run; 


**variable selection cholklas***;
proc logistic data=birnh2;
class  gender STATSOC1 STATSOC2 GEWEST1;
model cholklas= gender STATSOC1 STATSOC2 BMI2 AGE GEWEST1/link=probit;
run;
proc logistic data=birnh2;
model cholklas= STATSOC1 STATSOC2 BMI2 AGE GEWEST1/link=probit;
run;


**correlation depends on parameters via inverse fisher z transform***;
proc nlmixed data=BIRNH2 qpoints=5 maxiter=10000 maxfunc=10000 technique=newrap HESSIAN SUBGRADIENT=mvp.grad_b;
parms
constant12	0.024201016
r_sex12	-0.08659
r_soc12	0.09297
constant13	-0.03362
r_sex13 -0.08659
r_soc13	0.09297
constant23	-0.03362
r_sex23	-0.08659
r_soc23	0.09297
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
int_chol	-2.437889287
beta_c_statsoc1	-0.079429388
beta_c_statsoc2	-0.170497579
beta_c_bmi	-0.035169399
beta_c_age	-0.027805511
beta_c_gewest	0.125759932
;
y12=(constant12+r_sex12*gender+r_soc12*STATSOC2);
r12=(exp(2 * y12) - 1)/(exp(2 * y12) + 1);
y13=(constant13+r_sex13*gender+r_soc13*STATSOC2);
r13=(exp(2 * y13) - 1)/(exp(2 * y13) + 1);
y23=(constant23+r_sex23*gender+r_soc23*STATSOC2);
r23=(exp(2 * y23) - 1)/(exp(2 * y23) + 1);
eta_d = beta_d_gender*gender+
beta_d_statsoc*STATSOC1+beta_d_gewest*GEWEST1;
eta_s = beta_s_gender*gender+beta_s_bmi*BMI2+
beta_s_age*AGE+
beta_s_statsoc1*STATSOC1+beta_s_statsoc2*STATSOC2;
eta_c = beta_c_bmi*BMI2+
beta_c_age*AGE+beta_c_gewest*GEWEST1+
beta_c_statsoc1*STATSOC1+beta_c_statsoc2*STATSOC2;
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
ll=log(lik12)+log(lik13)+log(lik23);
model response ~ general(ll);
ods output hessian=mvp.hessian_b parameterestimates=mvp.parms_b;
run;


*sandwich estimator;
data hessian;
set mvp.hessian_b;
drop row Parameter;
run;
data gradient;
set mvp.grad_b;
drop observation;
run;
data parms;
set mvp.parms_b;
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

