/* Sandwich (robust) variance combination in PROC IML, following the pairwise
   fitting approach: block-diagonal Hessian J from the per-pair Hessians, the
   cross-product K from the stacked subject gradients, and the pseudo-likelihood
   covariance Sigma = inv(J)*K*inv(J) with per-parameter standard errors from
   its diagonal. Source: the PROC IML "create J and K" / "estimate sigma" blocks
   in "SIMULATION Pairwise fitting approach binary responses.sas".
   The per-pair Hessian and subject-gradient blocks (normally the ODS output of
   the three pairwise PROC NLMIXED fits) are supplied here as small numeric
   matrices of the same shape so the combination runs stand-alone. */

data H_1;
  input c1 c2;
  datalines;
6.2 1.1
1.1 4.8
;
run;
data H_2;
  input c1 c2;
  datalines;
5.5 0.7
0.7 5.1
;
run;
data H_3;
  input c1 c2;
  datalines;
7.0 1.4
1.4 6.3
;
run;

/* per-pair subject gradient blocks: one column per subject, one row per parameter */
data G_1;
  input s1 s2 s3 s4 s5;
  datalines;
0.21 -0.13 0.30 -0.05 0.11
-0.18 0.24 -0.09 0.14 -0.07
;
run;
data G_2;
  input s1 s2 s3 s4 s5;
  datalines;
0.16 -0.10 0.22 -0.04 0.08
-0.12 0.19 -0.06 0.10 -0.05
;
run;
data G_3;
  input s1 s2 s3 s4 s5;
  datalines;
0.25 -0.15 0.34 -0.07 0.13
-0.20 0.27 -0.11 0.16 -0.08
;
run;

proc iml;
use H_1; read all into H_1; close H_1;
use H_2; read all into H_2; close H_2;
use H_3; read all into H_3; close H_3;
/* block-diagonal Hessian across the three pairs */
H=block(H_1, H_2, H_3);

use G_1; read all into G_1; close G_1;
use G_2; read all into G_2; close G_2;
use G_3; read all into G_3; close G_3;
/* stack the per-pair subject gradients: parameters down rows, subjects across columns */
G=G_1//G_2//G_3;
nsubjects=ncol(G);

/* J = averaged Hessian, K = averaged gradient cross-product */
J=H#1/nsubjects;
K=(G*t(G))#1/nsubjects;

/* sandwich covariance and standard errors */
Sigma=inv(J)*K*inv(J);
Sigma0=Sigma#1/nsubjects;
var_pm=vecdiag(Sigma0);
se_pm=sqrt(var_pm);
print nsubjects;
print J;
print K;
print se_pm;
create se_pm from se_pm; append from se_pm;
quit;

proc print data=se_pm;
run;
