/* Pseudo-likelihood (composite) fit of a trivariate probit model.
   Source: %composite_lh_3 from
   "SIMULATION Pseudo-likelihood estimation binary responses.sas".
   External xlsx import replaced by an inline sample with the same
   columns the macro reads (y1b y2b y3b sex age region). */

data simulated_0;
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

%macro composite_lh_3(mydataset,i,j,k);
/*create dummies*/
data &mydataset;
set &mydataset;
if region=2 then region_2=1;else region_2=0;
if region=3 then region_3=1;else region_3=0;
if region=4 then region_4=1;else region_4=0;
if region=5 then region_5=1;else region_5=0;
run;
/*fit pairwise model*/
proc nlmixed data=&mydataset qpoints=5 maxiter=10000 maxfunc=10000 technique=newrap HESSIAN SUBGRADIENT=G;
parms
r&i&j=0
r&i&k=0
r&j&k=0
intercept_&i=0
beta_sex_&i=0
beta_age_&i=0
beta_reg2_&i=0
beta_reg3_&i=0
beta_reg4_&i=0
beta_reg5_&i=0
intercept_&j=0
beta_sex_&j=0
beta_age_&j=0
beta_reg2_&j=0
beta_reg3_&j=0
beta_reg4_&j=0
beta_reg5_&j=0
intercept_&k=0
beta_sex_&k=0
beta_age_&k=0
beta_reg2_&k=0
beta_reg3_&k=0
beta_reg4_&k=0
beta_reg5_&k=0;

xb_&i = intercept_&i+beta_sex_&i*sex+beta_age_&i*age+
beta_reg2_&i*region_2+beta_reg3_&i*region_3+beta_reg4_&i*region_4+
beta_reg5_&i*region_5;
xb_&j = intercept_&j+beta_sex_&j*sex+beta_age_&j*age+
beta_reg2_&j*region_2+beta_reg3_&j*region_3+beta_reg4_&j*region_4+
beta_reg5_&j*region_5;
xb_&k = intercept_&k+beta_sex_&k*sex+beta_age_&k*age+
beta_reg2_&k*region_2+beta_reg3_&k*region_3+beta_reg4_&k*region_4+
beta_reg5_&k*region_5;

/***i en j***/
if y1b=0 and y2b=0 then do;
lik_&i&j = 1-
cdf('NORMAL',(xb_&j))-
cdf('NORMAL',(xb_&i))+
probbnrm(xb_&i,xb_&j,r&i&j);
end;
if y1b=0 and y2b=1 then do;
lik_&i&j = cdf('NORMAL',(xb_&j))-
probbnrm(xb_&i,xb_&j,r&i&j);
end;
if y1b=1 and y2b=0 then do;
lik_&i&j = cdf('NORMAL',(xb_&i))-
probbnrm(xb_&i,xb_&j,r&i&j);
end;
if y1b=1 and y2b=1 then do;
lik_&i&j = probbnrm(xb_&i,xb_&j,r&i&j);
end;
ll_&i&j=log(lik_&i&j);

/***i en k ***/
if y1b=0 and y3b=0 then do;
lik_&i&k = 1-
cdf('NORMAL',(xb_&k))-
cdf('NORMAL',(xb_&i))+
probbnrm(xb_&i,xb_&k,r&i&k);
end;
if y1b=0 and y3b=1 then do;
lik_&i&k = cdf('NORMAL',(xb_&k))-
probbnrm(xb_&i,xb_&k,r&i&k);
end;
if y1b=1 and y3b=0 then do;
lik_&i&k = cdf('NORMAL',(xb_&i))-
probbnrm(xb_&i,xb_&k,r&i&k);
end;
if y1b=1 and y3b=1 then do;
lik_&i&k = probbnrm(xb_&i,xb_&k,r&i&k);
end;
ll_&i&k=log(lik_&i&k);

/***j en k ***/
if y2b=0 and y3b=0 then do;
lik_&j&k = 1-
cdf('NORMAL',(xb_&k))-
cdf('NORMAL',(xb_&j))+
probbnrm(xb_&j,xb_&k,r&j&k);
end;
if y2b=0 and y3b=1 then do;
lik_&j&k = cdf('NORMAL',(xb_&k))-
probbnrm(xb_&j,xb_&k,r&j&k);
end;
if y2b=1 and y3b=0 then do;
lik_&j&k = cdf('NORMAL',(xb_&j))-
probbnrm(xb_&j,xb_&k,r&j&k);
end;
if y2b=1 and y3b=1 then do;
lik_&j&k = probbnrm(xb_&j,xb_&k,r&j&k);
end;
ll_&j&k=log(lik_&j&k);

/*sum log lh*/
ll=ll_&i&j+ll_&j&k+ll_&i&k;
model y1b ~ general(ll);
ods output hessian=H parameterestimates=parms;
run;


proc print data=parms;
run;
%mend composite_lh_3;

%composite_lh_3(simulated_0,1,2,3);
