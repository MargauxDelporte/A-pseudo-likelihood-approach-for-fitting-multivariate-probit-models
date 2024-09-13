libname d 'C:\Users\u0118563\OneDrive - KU Leuven\Projecten\MV probit\SAS code\Simulatiestudie\data binair';
libname s 'C:\Users\u0118563\OneDrive - KU Leuven\Projecten\MV probit\SAS code\Simulatiestudie\results_pl_3 bin cl';


/*import macro*/
%macro import_xlsx(filename);
filename xlsx "C:\Users\u0118563\OneDrive - KU Leuven\Projecten\MV probit\Data\Simulatie binary\&filename..xlsx" termstr=LF;
proc import
  datafile=xlsx
  out=work.&filename.
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
/*********************************************3 responses*/****************************************************
/**************************************************************************************************************/;

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
/*1*/ 
if y1b=0 and y2b=0 then do;
lik_&i&j = 1-
cdf('NORMAL',(xb_&j))-
cdf('NORMAL',(xb_&i))+
probbnrm(xb_&i,xb_&j,r&i&j);
end;
/*2*/
if y1b=0 and y2b=1 then do;
lik_&i&j = cdf('NORMAL',(xb_&j))-
probbnrm(xb_&i,xb_&j,r&i&j);
end;
/*3*/
if y1b=1 and y2b=0 then do;
lik_&i&j = cdf('NORMAL',(xb_&i))-
probbnrm(xb_&i,xb_&j,r&i&j);
end;
/*4*/
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
/*2*/
if y1b=0 and y3b=1 then do;
lik_&i&k = cdf('NORMAL',(xb_&k))-
probbnrm(xb_&i,xb_&k,r&i&k);
end;
/*3*/
if y1b=1 and y3b=0 then do;
lik_&i&k = cdf('NORMAL',(xb_&i))-
probbnrm(xb_&i,xb_&k,r&i&k);
end;
/*4*/
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
/*2*/
if y2b=0 and y3b=1 then do;
lik_&j&k = cdf('NORMAL',(xb_&k))-
probbnrm(xb_&j,xb_&k,r&j&k);
end;
/*3*/
if y2b=1 and y3b=0 then do;
lik_&j&k = cdf('NORMAL',(xb_&j))-
probbnrm(xb_&j,xb_&k,r&j&k);
end;
/*4*/
if y2b=1 and y3b=1 then do;
lik_&j&k = probbnrm(xb_&j,xb_&k,r&j&k);
end;
ll_&j&k=log(lik_&j&k);

/*sum log lh*/
ll=ll_&i&j+ll_&j&k+ll_&i&k;
model y1b ~ general(ll);
ods output hessian=H parameterestimates=parms;
run;

*adapt hessians and gradients for combination;
data hessian;
set H;
drop row Parameter;
run;

data gradient;
set G;
drop observation;
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

data s.result_&mydataset;
merge parms se;
run;

proc datasets nolist lib=work;
 delete  H hessian G gradient  J K Sigma se 
parms;
quit;

%mend composite_lh_3;
/* Start timer */
%let _timer_start = %sysfunc(datetime());
data _null_;
set files;
call execute('%composite_lh_3('!!trim(filename)!!',1,2,3);');
run;
/* Stop timer */
data _null_;
  dur = datetime() - &_timer_start;
  put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;

/* TOTAL DURATION:  */


%composite_lh_3(simulated_0,1,2,3);

data s.start62;
run;

%composite_lh_3(simulated_62,1,2,3);
%composite_lh_3(simulated_63,1,2,3);
%composite_lh_3(simulated_64,1,2,3);
%composite_lh_3(simulated_65,1,2,3);
%composite_lh_3(simulated_66,1,2,3);
%composite_lh_3(simulated_67,1,2,3);
%composite_lh_3(simulated_68,1,2,3);
%composite_lh_3(simulated_69,1,2,3);
%composite_lh_3(simulated_70,1,2,3);
%composite_lh_3(simulated_71,1,2,3);
%composite_lh_3(simulated_72,1,2,3);
%composite_lh_3(simulated_73,1,2,3);
%composite_lh_3(simulated_74,1,2,3);
%composite_lh_3(simulated_75,1,2,3);
%composite_lh_3(simulated_76,1,2,3);
%composite_lh_3(simulated_77,1,2,3);
%composite_lh_3(simulated_78,1,2,3);
%composite_lh_3(simulated_79,1,2,3);
%composite_lh_3(simulated_80,1,2,3);
%composite_lh_3(simulated_81,1,2,3);
%composite_lh_3(simulated_82,1,2,3);
%composite_lh_3(simulated_83,1,2,3);
%composite_lh_3(simulated_84,1,2,3);
%composite_lh_3(simulated_85,1,2,3);
%composite_lh_3(simulated_86,1,2,3);
%composite_lh_3(simulated_87,1,2,3);
%composite_lh_3(simulated_88,1,2,3);
%composite_lh_3(simulated_89,1,2,3);
%composite_lh_3(simulated_90,1,2,3);
%composite_lh_3(simulated_91,1,2,3);
%composite_lh_3(simulated_92,1,2,3);
%composite_lh_3(simulated_93,1,2,3);
%composite_lh_3(simulated_94,1,2,3);
%composite_lh_3(simulated_95,1,2,3);
%composite_lh_3(simulated_96,1,2,3);
%composite_lh_3(simulated_97,1,2,3);
%composite_lh_3(simulated_98,1,2,3);
%composite_lh_3(simulated_99,1,2,3);
