libname d 'C:\Users\u0118563\OneDrive - KU Leuven\Projecten\MV probit\SAS code\Simulatiestudie\data binair';
libname s 'C:\Users\u0118563\OneDrive - KU Leuven\Projecten\MV probit\SAS code\Simulatiestudie\results_pp_3 bin cl';


/*import macro*/
%macro import_xlsx(filename);
filename xlsx "C:\Users\u0118563\OneDrive - KU Leuven\Projecten\MV probit\Data\Simulatie binary\&filename..xlsx" termstr=LF;
proc import
  datafile=xlsx
  out=d.&filename.
  dbms=xlsx
  replace
;
run;
%mend;

/*import a list with the names of files to import*/
proc import datafile = 'C:\Users\u0118563\OneDrive - KU Leuven\Projecten\MV probit\Data\Simulation study\list.txt'
 out = files
 dbms = dlm
 replace;
run;


/*import files*/
data _null_;
set files;
call execute('%import_xlsx('!!trim(filename)!!');');
run;
ODS _all_ Close;


/**************************************************************************************************************
/*********************************************3 responses pairwise lh*/****************************************************
/**************************************************************************************************************/;
%macro combine(mydataset);
/*drop unnessecary columns in the Hians, Gients, and parameter estimates*/
data H_1;
set H_12;
drop row Parameter;
run;
data G_1;
set G_12;
drop observation;
run;
data H_2;
set H_13;
drop row Parameter;
run;
data G_2;
set G_13;
drop observation;
run;
data H_3;
set H_23;
drop row Parameter;
run;
data G_3;
set G_23;
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

/*combine the Hessians and Gradients*/
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

data s.result_&mydataset;
merge result variances_pm;
by parameter;
run;

proc datasets nolist lib=work;
delete  
G G_1 G_12 G_2 G_23 G_3 G_13
H H_1 H_12 H_2 H_23 H_3 H_13
parameters Parms_12 Parms_13 Parms_23
Result Result_se Sigma0
variances_pm var_pm
P_1 P_2 P_3
estimate_12 estimate_13 estimate_23
J K Sigma se 
parms;
quit;
%mend combine;



%macro pairanalysis(mydataset);
/*create dummyvariables*/
data d.&mydataset;
set d.&mydataset;
if region=2 then region_2=1;else region_2=0;
if region=3 then region_3=1;else region_3=0;
if region=4 then region_4=1;else region_4=0;
if region=5 then region_5=1;else region_5=0;
run;
/*fit pairwise model for pair 1*/
proc nlmixed data=d.&mydataset qpoints=5 maxiter=10000 maxfunc=10000 technique=newrap HESSIAN SUBGRADIENT=G_12;
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

*fit pairwise model for pair 2*/;
proc nlmixed data=d.&mydataset qpoints=5 maxiter=10000 maxfunc=10000 technique=newrap HESSIAN SUBGRADIENT=G_13;
parms
r13=0
intercept_1=0
beta_sex_1=0
beta_age_1=0
beta_reg2_1=0
beta_reg3_1=0
beta_reg4_1=0
beta_reg5_1=0
intercept_3=0
beta_sex_3=0
beta_age_3=0
beta_reg2_3=0
beta_reg3_3=0
beta_reg4_3=0
beta_reg5_3=0;

xb_1 = intercept_1+beta_sex_1*sex+beta_age_1*age+
beta_reg2_1*region_2+beta_reg3_1*region_3+beta_reg4_1*region_4+
beta_reg5_1*region_5;
xb_3 = intercept_3+beta_sex_3*sex+beta_age_3*age+
beta_reg2_3*region_2+beta_reg3_3*region_3+beta_reg4_3*region_4+
beta_reg5_3*region_5;

if y1b=0 and y3b=0 then do;
lik_13 = 1-
cdf('NORMAL',(xb_3))-
cdf('NORMAL',(xb_1))+
probbnrm(xb_1,xb_3,r13);
end;
/*2*/
if y1b=0 and y3b=1 then do;
lik_13 = cdf('NORMAL',(xb_3))-
probbnrm(xb_1,xb_3,r13);
end;
/*3*/
if y1b=1 and y3b=0 then do;
lik_13 = cdf('NORMAL',(xb_1))-
probbnrm(xb_1,xb_3,r13);
end;
/*4*/
if y1b=1 and y3b=1 then do;
lik_13 = probbnrm(xb_1,xb_3,r13);
end;
ll_13=log(lik_13);
model y1b ~ general(ll_13);
ods output Hessian=H_13 parameterestimates=parms_13;
run;

*fit pairwise model for pair 3*/;
proc nlmixed data=d.&mydataset qpoints=5 maxiter=10000 maxfunc=10000 technique=newrap HESSIAN SUBGRADIENT=G_23;
parms
r23=0
intercept_2=0
beta_sex_2=0
beta_age_2=0
beta_reg2_2=0
beta_reg3_2=0
beta_reg4_2=0
beta_reg5_2=0
intercept_3=0
beta_sex_3=0
beta_age_3=0
beta_reg2_3=0
beta_reg3_3=0
beta_reg4_3=0
beta_reg5_3=0;

xb_2 = intercept_2+beta_sex_2*sex+beta_age_2*age+
beta_reg2_2*region_2+beta_reg3_2*region_3+beta_reg4_2*region_4+
beta_reg5_2*region_5;
xb_3 = intercept_3+beta_sex_3*sex+beta_age_3*age+
beta_reg2_3*region_2+beta_reg3_3*region_3+beta_reg4_3*region_4+
beta_reg5_3*region_5;

if y2b=0 and y3b=0 then do;
lik_23 = 1-
cdf('NORMAL',(xb_3))-
cdf('NORMAL',(xb_2))+
probbnrm(xb_2,xb_3,r23);
end;
/*2*/
if y2b=0 and y3b=1 then do;
lik_23 = cdf('NORMAL',(xb_3))-
probbnrm(xb_2,xb_3,r23);
end;
/*3*/
if y2b=1 and y3b=0 then do;
lik_23 = cdf('NORMAL',(xb_2))-
probbnrm(xb_2,xb_3,r23);
end;
/*4*/
if y2b=1 and y3b=1 then do;
lik_23 = probbnrm(xb_2,xb_3,r23);
end;
ll_23=log(lik_23);
model y1b ~ general(ll_23);
ods output Hessian=H_23 parameterestimates=parms_23;
run;
%combine(&mydataset);
%mend pairanalysis;

data s.start00;
run;

/* Start timer */
%let _timer_start = %sysfunc(datetime());
data _null_;
set files;
call execute('%pairanalysis('!!trim(filename)!!');');
run;
/* Stop timer */
data _null_;
  dur = datetime() - &_timer_start;
  put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;

/* TOTAL DURATION:  */





