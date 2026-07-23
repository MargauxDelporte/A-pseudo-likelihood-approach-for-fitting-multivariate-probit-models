/* Pairwise fitting approach: fit the pairwise probit model for one response
   pair (1,2), exactly as in %pairanalysis from
   "SIMULATION Pairwise fitting approach binary responses.sas".
   External xlsx import replaced by an inline sample with the same columns
   the macro reads (y1b y2b y3b sex age region). */

data pairdata;
  input y1b y2b y3b sex age region;
  datalines;
0 0 1 1 45 2
1 0 0 0 51 3
1 1 1 1 38 4
0 1 0 0 62 5
1 0 1 1 29 2
0 0 0 1 44 3
1 1 0 0 55 4
0 1 1 1 33 5
1 0 0 0 47 2
0 0 1 1 60 3
1 1 1 0 41 4
0 1 0 1 36 5
1 0 1 0 58 2
0 0 0 1 49 3
1 1 0 0 27 4
0 1 1 0 53 5
1 0 0 1 39 2
0 1 0 0 46 3
1 1 1 1 31 4
0 0 1 0 57 5
;
run;

/*create dummyvariables*/
data pairdata;
set pairdata;
if region=2 then region_2=1;else region_2=0;
if region=3 then region_3=1;else region_3=0;
if region=4 then region_4=1;else region_4=0;
if region=5 then region_5=1;else region_5=0;
run;

/*fit pairwise model for pair 1*/
proc nlmixed data=pairdata qpoints=5 maxiter=10000 maxfunc=10000 technique=newrap HESSIAN SUBGRADIENT=G_12;
parms
r12=0
intercept_1=0
beta_sex_1=0
beta_age_1=0
beta_reg2_1=0
beta_reg3_1=0
beta_reg4_1=0
beta_reg5_1=0
intercept_2=0
beta_sex_2=0
beta_age_2=0
beta_reg2_2=0
beta_reg3_2=0
beta_reg4_2=0
beta_reg5_2=0;

xb_1 = intercept_1+beta_sex_1*sex+beta_age_1*age+
beta_reg2_1*region_2+beta_reg3_1*region_3+beta_reg4_1*region_4+
beta_reg5_1*region_5;
xb_2 = intercept_2+beta_sex_2*sex+beta_age_2*age+
beta_reg2_2*region_2+beta_reg3_2*region_3+beta_reg4_2*region_4+
beta_reg5_2*region_5;

if y1b=0 and y2b=0 then do;
lik_12 = 1-
cdf('NORMAL',(xb_2))-
cdf('NORMAL',(xb_1))+
probbnrm(xb_1,xb_2,r12);
end;
/*2*/
if y1b=0 and y2b=1 then do;
lik_12 = cdf('NORMAL',(xb_2))-
probbnrm(xb_1,xb_2,r12);
end;
/*3*/
if y1b=1 and y2b=0 then do;
lik_12 = cdf('NORMAL',(xb_1))-
probbnrm(xb_1,xb_2,r12);
end;
/*4*/
if y1b=1 and y2b=1 then do;
lik_12 = probbnrm(xb_1,xb_2,r12);
end;
ll_12=log(lik_12);
model y1b ~ general(ll_12);
ods output Hessian=H_12 parameterestimates=parms_12;
run;

proc print data=parms_12;
run;
