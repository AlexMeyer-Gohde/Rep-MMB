//**************************************************************************
// A New Comparative Approach to Macroeconomic Modeling and Policy Analysis
//
// Volker Wieland, Tobias Cwik, Gernot J. Mueller, Sebastian Schmidt and
// Maik Wolters
//
// Working Paper, 2009
//**************************************************************************

// Model: NK_CGG99

// Further references:
// Clarida, R., J. Gal?and M. Gertler. 1999. "The Science of Monetary Policy: A New Keynesian Perspective."
// Journal of Economic Literature 37(4), pp. 1661-1707.

// Last edited: 10/08/26 by S. Schmidt

// This file simulates the dynamic response of the model to specific shocks.
// This version uses the extended version of the model (starting on page 1691),
// which adds output persistence and endogenous inflation to the baseline model.
// We use the same parameters as in Rotemberg Woodford (1997) and add a backward looking
// inflation share of 0.48 and a backward looking part of output of 0.44.


var x i pi;



varexo inflation_ demand_;



parameters

 theta sigma phi lambda beta;

theta  =  0.44;
lambda =  0.0244;
phi    =  0.48;
sigma  = -6.25;
beta   =  1/(1+0.035/4);



model(linear);



// Original Model Code:

// Original interest rate rule from table 1 in Clarida et al. (1999)
i = 0.79*i(-1)+(1-0.79)*(2.15*pi(+1)+(0.93/4)*x);
// We devide the original output gap coefficient by 4, since Clarida et al. (2000)
// report that they use annualized quarterly inflation and interest rate data
// to estimate the rule.

x  =  sigma *( i - pi(+1) ) + theta * x(-1)  + (1-theta) * x(+1) + demand_;
pi   =  lambda*x + phi * pi(-1) + (1-phi) * beta *  pi(+1) + inflation_;
end;


shocks;
var demand_;
stderr 0.63;
var inflation_;
stderr 0.4;
end;


//stoch_simul(irf = 0, ar=100, noprint);
stoch_simul (AR=100,IRF=0, noprint,nograph);