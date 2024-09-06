//**************************************************************************
// A New Comparative Approach to Macroeconomic Modeling and Policy Analysis
//
// Volker Wieland, Tobias Cwik, Gernot J. Mueller, Sebastian Schmidt and
// Maik Wolters
//
// Working Paper, 2009
//**************************************************************************

// Model: US_VMDop ("optimistic times")

// Further references:
//      Verona, Martins and Drumond, "(Un)anticipated monetary policy in a DSGE model with a shadow banking system",
//      International Journal of Central Banking 9 (3), 78-124, September 2013

// implemented by Fabio Verona (in September 2012), email: fabio.verona@bof.fi
// Last edited: 12/10/22 by M. Burgert

// Variables

var piU, iU, omegabarSU, RkSXU, nSU, qU,
         cU, wU,rkSU,kbarBU, hU, kbarSU, ReXU,
        FpXU, FwXU,  RcouponXU,uSU,uBU,
          KwXU, KpXU,chiU,KtotU,YU,eps_analystU,eps_couponU,RanalystXU,nBU,stockcapitalU,ntotU,RkBXU,
            btotAU,btotBU,btotSU,levAU,levBU,levSU,CostFinancingSU,SpreadU, lambdanU,PextSU,

    piflexU, iflexU, omegabarSflexU, RkSXflexU, nSflexU, qflexU,
         cflexU, wflexU,rkSflexU,kbarBflexU, hflexU, kbarSflexU, ReXflexU,
        FpXflexU, FwXflexU,  RcouponXflexU,uSflexU,uBflexU,
          KwXflexU, KpXflexU,chiflexU,KtotflexU,YflexU,eps_analystflexU,eps_couponflexU,RanalystXflexU,nBflexU,stockcapitalflexU,ntotflexU,RkBXflexU,
            btotAflexU,btotBflexU,btotSflexU,levAflexU,levBflexU,levSflexU,CostFinancingSflexU,SpreadflexU, lambdanflexU,PextSflexU;



varexo e_xpU;



parameters  KtotUU, ntotUU,KwXUU,KpXUU, RcouponXUU, eps_couponUU,levBUU,btotBUU,btotSUU,levSUU,
            piUU, sUU, rkSUU,rkBUU, iUU, uSUU,uBUU, omegabarSUU, RkSXUU,RkBXUU, nSUU,nBUU, qUU, lambdanUU, cUU, wUU, hUU,
            kbarSUU,kbarBUU, ReXUU,  FpXUU, FwXUU, PextSUU,YUU,eps_analystUU, RanalystXUU, chiUU, eps_biasedUU, RbiasedXUU ,
              gammaSUU,gammaBUU,  sigmaSUU,  sigmaCUU,levAUU,btotAUU,CostFinancingSUU,SpreadUU,stockcapitalUU,gUU;


//PARAMETERS

parameters lambdawUU, sigmaLXUU, betaUU, xiwUU, bUU,lambdafUU,
      psiLXUU, iotaw1UU,  weSUU,weBUU, sigmaaSUU,sigmaaBUU, SdouprXUU,
      alphaUU, deltaUU, xipUU, iota1UU, etagUU, muUU,
      aptilUU, aytilUU, rhotilUU,
      rhoEIS,etaSE, rho_chi, alpha3, alpha4,KwXflexUU,KpXflexUU,FpXflexUU, FwXflexUU %,OMEGA
	  std1_xpU;



std1_xpU=0.25;

alpha3=24;
alpha4=30000;





M_.params(76)=alpha3;
M_.params(75)=0;

pp=load('FIMSteadyStateMMB');
%for i=1:M_.param_nbr;
   % name = deblank(M_.param_names(i,:));
    %if isfield(pp,name);
     %  M_.params(i) = eval(['pp.' name]);
      %// disp([name ' ' num2str(M_.params(i))]);
    %end;
%end;
M_.params(76)=alpha3;
M_.params(75)=0;
load FIMSteadyStateMMB
for jj=1:length(M_.param_names)
set_param_value(M_.param_names{jj},eval(M_.param_names{jj})); 
end;

model;




-nBU + gammaBUU/piU*kbarBU(-1)*qU(-1)*(RkBXU-RcouponXU(-1)) + gammaBUU/piU*(1+RcouponXU(-1))*nBU(-1) + weBUU;
(1+betaUU*RcouponXU)*qU - betaUU*piU(+1)*(  uBU(+1)*rkSU(+1)*( (uSU(+1)*kbarSU/(uBU(+1)*kbarBU))^(1-rhoEIS)) - rkBUU*(exp(sigmaaBUU*(uBU-1))-1)/sigmaaBUU + (1-deltaUU)*qU(+1));
eps_analystU - eps_analystUU - alpha4*(YU(-1)-YUU);
eps_couponU - eps_analystU*(1 + chiU);
chiU - rho_chi*chiU(-1) - (1-rho_chi)*(chiUU + alpha3*(nBU-nBUU));
1 - eps_couponU*(RcouponXU-ReXU)/(1+RcouponXU);
1 - eps_analystU*(RanalystXU-ReXU)/(1+RanalystXU);
-(rkSU/( alphaUU * ((uSU*kbarSU(-1))^(rhoEIS-1)) *((hU / ( ( etaSE*(uSU*kbarSU(-1))^rhoEIS +(1-etaSE)*(uBU*kbarBU(-1))^rhoEIS)^(1/rhoEIS)))^(1-alphaUU)) * (( etaSE*(uSU*kbarSU(-1))^rhoEIS +(1-etaSE)*(uBU*kbarBU(-1))^rhoEIS))^(1/rhoEIS-1))) + ( (wU/(1-alphaUU))^(1-alphaUU/(rhoEIS+alphaUU-alphaUU*rhoEIS))) * ( ( alphaUU/rkSU*((uSU*kbarSU(-1))^(rhoEIS-1)))^(-alphaUU/(rhoEIS+alphaUU-alphaUU*rhoEIS))) *((((( ( etaSE*(uSU*kbarSU(-1))^rhoEIS +(1-etaSE)*(uBU*kbarBU(-1))^rhoEIS)^(1/rhoEIS))^alphaUU)*(hU ^(1-alphaUU)) ))^(alphaUU*(rhoEIS-1)/(rhoEIS+alphaUU-alphaUU*rhoEIS))) * rhoEIS/(rhoEIS+alphaUU-alphaUU*rhoEIS);
rkSU - rkSUU*exp(sigmaaSUU*(uSU-1));
rkSU*( (uSU*kbarSU(-1)/(uBU*kbarBU(-1)))^(1-rhoEIS)) - rkBUU*exp(sigmaaBUU*(uBU-1));
1+ RkBXU -  piU/qU(-1)*(  uBU*rkSU*( (uSU*kbarSU(-1)/(uBU*kbarBU(-1)))^(1-rhoEIS))  - rkBUU*(exp(sigmaaBUU*(uBU-1))-1)/sigmaaBUU  + (1-deltaUU)*qU);
lambdanU*qU*( 1 - SdouprXUU/2*((iU/iU(-1)-1)^2) - SdouprXUU*iU/iU(-1)*(iU/iU(-1)-1)) - lambdanU + betaUU*lambdanU(+1)*qU(+1)*(SdouprXUU * ((iU(+1)/iU)^2)*(iU(+1)/iU-1));
etaSE*kbarSU + (1-etaSE)*kbarBU - (1-deltaUU)*(etaSE*kbarSU(-1) + (1-etaSE)*kbarBU(-1)) - (1 - SdouprXUU/2*(iU/iU(-1)-1)^2)*iU;
1+ RkSXU -  piU/qU(-1)*(  uSU*rkSU - rkSUU*(exp(sigmaaSUU*(uSU-1))-1)/sigmaaSUU  + (1-deltaUU)*qU);
(1- ( normcdf(( (log(omegabarSU(+1))+sigmaSUU^2/2)/sigmaSUU )-sigmaSUU)+omegabarSU(+1)*(1-normcdf(( (log(omegabarSU(+1))+sigmaSUU^2/2)/sigmaSUU ))) ))*(1+RkSXU(+1))/(1+ReXU)  +  ( (1/(omegabarSU(+1)*sigmaSUU))*normpdf(( log(omegabarSU(+1))-sigmaSUU^2/2)/sigmaSUU)+  1 - normcdf(( log(omegabarSU(+1))+sigmaSUU^2/2)/sigmaSUU) - (1/sigmaSUU)*normpdf(( log(omegabarSU(+1))+sigmaSUU^2/2)/sigmaSUU)) / (((1-muUU)/(omegabarSU(+1)*sigmaSUU))* normpdf(( log(omegabarSU(+1))-sigmaSUU^2/2)/sigmaSUU)+ 1- normcdf(( log(omegabarSU(+1))+sigmaSUU^2/2)/sigmaSUU) - (1/sigmaSUU)*normpdf(( log(omegabarSU(+1))+sigmaSUU^2/2)/sigmaSUU))*( (1+RkSXU(+1))/(1+ReXU)*( (1-muUU)*normcdf(( (log(omegabarSU(+1))+sigmaSUU^2/2)/sigmaSUU )-sigmaSUU)+ omegabarSU(+1)*(1-normcdf(( (log(omegabarSU(+1))+sigmaSUU^2/2)/sigmaSUU ))) )-1);
qU(-1)*kbarSU(-1)*(1+RkSXU)*( (1-muUU)*normcdf(( (log(omegabarSU)+sigmaSUU^2/2)/sigmaSUU )-sigmaSUU) + omegabarSU*(1-normcdf(( (log(omegabarSU)+sigmaSUU^2/2)/sigmaSUU ))) )/(nSU(-1)*(1+ReXU(-1))) - qU(-1)*kbarSU(-1)/nSU(-1) + 1;
-nSU + gammaSUU/piU*kbarSU(-1)*qU(-1)*(RkSXU-ReXU(-1)-muUU*(normcdf (( log(omegabarSU)-sigmaSUU^2/2)/sigmaSUU))*(1+RkSXU)) + weSUU + gammaSUU*(1+ReXU(-1))*nSU(-1)/piU;
-lambdanU + betaUU*(1+ReXU)*lambdanU(+1)/piU(+1);
lambdanU - ( (cU-bUU*cU(-1))^(-sigmaCUU) ) + bUU*betaUU*( (cU(+1)-bUU*cU)^(-sigmaCUU) );
YU - gUU - cU - iU  - ( rkSUU*(exp(sigmaaSUU*(uSU-1))-1)/sigmaaSUU )*etaSE*kbarSU(-1) - ( rkBUU*(exp(sigmaaBUU*(uBU-1))-1)/sigmaaBUU )*(1-etaSE)*kbarBU(-1) - muUU* (normcdf (( log(omegabarSU)-sigmaSUU^2/2)/sigmaSUU)) * (1+RkSXU)*qU(-1)*etaSE*kbarSU(-1)/piU;
lambdanU*YU + betaUU * xipUU *( ( (piU^(1-iota1UU))/piU(+1)) ^ (1/(1-lambdafUU)) ) *  FpXU(+1) - FpXU;
lambdanU*YU*lambdafUU*(rkSU/( alphaUU * ((uSU*kbarSU(-1))^(rhoEIS-1)) *((hU / ( ( etaSE*(uSU*kbarSU(-1))^rhoEIS +(1-etaSE)*(uBU*kbarBU(-1))^rhoEIS)^(1/rhoEIS)))^(1-alphaUU)) * (( etaSE*(uSU*kbarSU(-1))^rhoEIS +(1-etaSE)*(uBU*kbarBU(-1))^rhoEIS))^(1/rhoEIS-1))) +  betaUU*xipUU*KpXU(+1)* ( ( (piU^(1-iota1UU))/piU(+1))^(-lambdafUU/(lambdafUU-1))) - KpXU;
hU*(((cU-bUU*cU(-1))^(-sigmaCUU)) - bUU*betaUU*( (cU(+1)-bUU*cU)^(-sigmaCUU)))/lambdawUU +  betaUU * xiwUU *(( (piU^(1-iotaw1UU)))^(1/(1-lambdawUU)))/piU(+1) *( (1/( piU(+1)*wU(+1)/wU )) ^ (lambdawUU/(1-lambdawUU)) ) * FwXU(+1) - FwXU;
( hU^(1+sigmaLXUU)) + betaUU * xiwUU *((  (piU^(1-iotaw1UU)) /( piU(+1)*wU(+1)/wU ) )^(lambdawUU*(1+sigmaLXUU)/(1-lambdawUU))) *KwXU(+1) - KwXU;
KpXU - FpXU*( (1-xipUU*( ((piU(-1)^(1-iota1UU))/piU)^(1/(1-lambdafUU))))/(1-xipUU))^(1-lambdafUU);
KwXU - FwXU*wU/psiLXUU*((1 - xiwUU * ( (  (piU(-1)^(1-iotaw1UU))/( piU*wU/wU(-1) ) ) ^(1/(1-lambdawUU)) ) ) / (1-xiwUU) )^(1-lambdawUU*(1+sigmaLXUU));
ReXU - rhotilUU*ReXU(-1)-(1-rhotilUU)*(ReXUU+aptilUU*(piU(+1)-piUU)+aytilUU*(YU-YUU))+e_xpU/400;
KtotU - ( ( etaSE*(uSU*kbarSU(-1))^rhoEIS +(1-etaSE)*(uBU*kbarBU(-1))^rhoEIS)^(1/rhoEIS));
YU -  ((KtotU^alphaUU)*(hU ^(1-alphaUU)) ) ;
stockcapitalU - etaSE*kbarSU - (1-etaSE)*kbarBU;
ntotU - etaSE*nSU - (1-etaSE)*nBU;
btotBU - qU*kbarBU + nBU;
levBU - qU*kbarBU/nBU;//as in BGG
levSU - qU*kbarSU/nSU;//as in BGG
levAU - etaSE*levSU - (1-etaSE)*levBU;
btotSU - qU*kbarSU + nSU;
btotAU - etaSE*btotSU - (1-etaSE)*btotBU;
CostFinancingSU - omegabarSU(+1)*(1+RkSXU(+1))*qU*kbarSU/(qU*kbarSU-nSU);
SpreadU - CostFinancingSU + RcouponXU;
PextSU + (1+ReXU) - omegabarSU(+1)*(1+RkSXU(+1))*qU*kbarSU/(qU*kbarSU-nSU);

// FLEXIBLE PRICE EQUILIBRIUM EQUATIONS
-nBflexU + gammaBUU/piflexU*kbarBflexU(-1)*qflexU(-1)*(RkBXflexU-RcouponXflexU(-1)) + gammaBUU/piflexU*(1+RcouponXflexU(-1))*nBflexU(-1) + weBUU;
(1+betaUU*RcouponXflexU)*qflexU - betaUU*piflexU(+1)*(  uBflexU(+1)*rkSflexU(+1)*( (uSflexU(+1)*kbarSflexU/(uBflexU(+1)*kbarBflexU))^(1-rhoEIS)) - rkBUU*(exp(sigmaaBUU*(uBflexU-1))-1)/sigmaaBUU + (1-deltaUU)*qflexU(+1));
eps_analystflexU - eps_analystUU - alpha4*(YflexU(-1)-YUU);
eps_couponflexU - eps_analystflexU*(1 + chiflexU);
chiflexU - rho_chi*chiflexU(-1) - (1-rho_chi)*(chiUU + alpha3*(nBflexU-nBUU));
1 - eps_couponflexU*(RcouponXflexU-ReXflexU)/(1+RcouponXflexU);
1 - eps_analystflexU*(RanalystXflexU-ReXflexU)/(1+RanalystXflexU);
-(rkSflexU/( alphaUU * ((uSflexU*kbarSflexU(-1))^(rhoEIS-1)) *((hflexU / ( ( etaSE*(uSflexU*kbarSflexU(-1))^rhoEIS +(1-etaSE)*(uBflexU*kbarBflexU(-1))^rhoEIS)^(1/rhoEIS)))^(1-alphaUU)) * (( etaSE*(uSflexU*kbarSflexU(-1))^rhoEIS +(1-etaSE)*(uBflexU*kbarBflexU(-1))^rhoEIS))^(1/rhoEIS-1))) + ( (wflexU/(1-alphaUU))^(1-alphaUU/(rhoEIS+alphaUU-alphaUU*rhoEIS))) * ( ( alphaUU/rkSflexU*((uSflexU*kbarSflexU(-1))^(rhoEIS-1)))^(-alphaUU/(rhoEIS+alphaUU-alphaUU*rhoEIS))) *((((( ( etaSE*(uSflexU*kbarSflexU(-1))^rhoEIS +(1-etaSE)*(uBflexU*kbarBflexU(-1))^rhoEIS)^(1/rhoEIS))^alphaUU)*(hflexU ^(1-alphaUU)) ))^(alphaUU*(rhoEIS-1)/(rhoEIS+alphaUU-alphaUU*rhoEIS))) * rhoEIS/(rhoEIS+alphaUU-alphaUU*rhoEIS);
rkSflexU - rkSUU*exp(sigmaaSUU*(uSflexU-1));
rkSflexU*( (uSflexU*kbarSflexU(-1)/(uBflexU*kbarBflexU(-1)))^(1-rhoEIS)) - rkBUU*exp(sigmaaBUU*(uBflexU-1));
1+ RkBXflexU -  piflexU/qflexU(-1)*(  uBflexU*rkSflexU*( (uSflexU*kbarSflexU(-1)/(uBflexU*kbarBflexU(-1)))^(1-rhoEIS))  - rkBUU*(exp(sigmaaBUU*(uBflexU-1))-1)/sigmaaBUU  + (1-deltaUU)*qflexU);
lambdanflexU*qflexU*( 1 - SdouprXUU/2*((iflexU/iflexU(-1)-1)^2) - SdouprXUU*iflexU/iflexU(-1)*(iflexU/iflexU(-1)-1)) - lambdanflexU + betaUU*lambdanflexU(+1)*qflexU(+1)*(SdouprXUU * ((iflexU(+1)/iflexU)^2)*(iflexU(+1)/iflexU-1));
etaSE*kbarSflexU + (1-etaSE)*kbarBflexU - (1-deltaUU)*(etaSE*kbarSflexU(-1) + (1-etaSE)*kbarBflexU(-1)) - (1 - SdouprXUU/2*(iflexU/iflexU(-1)-1)^2)*iflexU;
1+ RkSXflexU -  piflexU/qflexU(-1)*(  uSflexU*rkSflexU - rkSUU*(exp(sigmaaSUU*(uSflexU-1))-1)/sigmaaSUU  + (1-deltaUU)*qflexU);
(1- ( normcdf(( (log(omegabarSflexU(+1))+sigmaSUU^2/2)/sigmaSUU )-sigmaSUU)+omegabarSflexU(+1)*(1-normcdf(( (log(omegabarSflexU(+1))+sigmaSUU^2/2)/sigmaSUU ))) ))*(1+RkSXflexU(+1))/(1+ReXflexU)  +  ( (1/(omegabarSflexU(+1)*sigmaSUU))*normpdf(( log(omegabarSflexU(+1))-sigmaSUU^2/2)/sigmaSUU)+  1 - normcdf(( log(omegabarSflexU(+1))+sigmaSUU^2/2)/sigmaSUU) - (1/sigmaSUU)*normpdf(( log(omegabarSflexU(+1))+sigmaSUU^2/2)/sigmaSUU)) / (((1-muUU)/(omegabarSflexU(+1)*sigmaSUU))* normpdf(( log(omegabarSflexU(+1))-sigmaSUU^2/2)/sigmaSUU)+ 1- normcdf(( log(omegabarSflexU(+1))+sigmaSUU^2/2)/sigmaSUU) - (1/sigmaSUU)*normpdf(( log(omegabarSflexU(+1))+sigmaSUU^2/2)/sigmaSUU))*( (1+RkSXflexU(+1))/(1+ReXflexU)*( (1-muUU)*normcdf(( (log(omegabarSflexU(+1))+sigmaSUU^2/2)/sigmaSUU )-sigmaSUU)+ omegabarSflexU(+1)*(1-normcdf(( (log(omegabarSflexU(+1))+sigmaSUU^2/2)/sigmaSUU ))) )-1);
qflexU(-1)*kbarSflexU(-1)*(1+RkSXflexU)*( (1-muUU)*normcdf(( (log(omegabarSflexU)+sigmaSUU^2/2)/sigmaSUU )-sigmaSUU) + omegabarSflexU*(1-normcdf(( (log(omegabarSflexU)+sigmaSUU^2/2)/sigmaSUU ))) )/(nSflexU(-1)*(1+ReXflexU(-1))) - qflexU(-1)*kbarSflexU(-1)/nSflexU(-1) + 1;
-nSflexU + gammaSUU/piflexU*kbarSflexU(-1)*qflexU(-1)*(RkSXflexU-ReXflexU(-1)-muUU*(normcdf (( log(omegabarSflexU)-sigmaSUU^2/2)/sigmaSUU))*(1+RkSXflexU)) + weSUU + gammaSUU*(1+ReXflexU(-1))*nSflexU(-1)/piflexU;
-lambdanflexU + betaUU*(1+ReXflexU)*lambdanflexU(+1)/piflexU(+1);
lambdanflexU - ( (cflexU-bUU*cflexU(-1))^(-sigmaCUU) ) + bUU*betaUU*( (cflexU(+1)-bUU*cflexU)^(-sigmaCUU) );
YflexU - gUU - cflexU - iflexU  - ( rkSUU*(exp(sigmaaSUU*(uSflexU-1))-1)/sigmaaSUU )*etaSE*kbarSflexU(-1) - ( rkBUU*(exp(sigmaaBUU*(uBflexU-1))-1)/sigmaaBUU )*(1-etaSE)*kbarBflexU(-1) - muUU* (normcdf (( log(omegabarSflexU)-sigmaSUU^2/2)/sigmaSUU)) * (1+RkSXflexU)*qflexU(-1)*etaSE*kbarSflexU(-1)/piflexU;
lambdanflexU*YflexU - FpXflexU;
lambdanflexU*YflexU*lambdafUU*(rkSflexU/( alphaUU * ((uSflexU*kbarSflexU(-1))^(rhoEIS-1)) *((hflexU / ( ( etaSE*(uSflexU*kbarSflexU(-1))^rhoEIS +(1-etaSE)*(uBflexU*kbarBflexU(-1))^rhoEIS)^(1/rhoEIS)))^(1-alphaUU)) * (( etaSE*(uSflexU*kbarSflexU(-1))^rhoEIS +(1-etaSE)*(uBflexU*kbarBflexU(-1))^rhoEIS))^(1/rhoEIS-1)))  - KpXflexU;
hflexU*(((cflexU-bUU*cflexU(-1))^(-sigmaCUU)) - bUU*betaUU*( (cflexU(+1)-bUU*cflexU)^(-sigmaCUU)))/lambdawUU - FwXflexU;
( hflexU^(1+sigmaLXUU)) - KwXflexU;
KpXflexU - FpXflexU;
KwXflexU - FwXflexU*wflexU/psiLXUU;
ReXflexU - rhotilUU*ReXflexU(-1)-(1-rhotilUU)*(ReXUU+aptilUU*(piflexU(+1)-piUU)+aytilUU*(YflexU-YUU));
KtotflexU - ( ( etaSE*(uSflexU*kbarSflexU(-1))^rhoEIS +(1-etaSE)*(uBflexU*kbarBflexU(-1))^rhoEIS)^(1/rhoEIS));
YflexU -  ((KtotflexU^alphaUU)*(hflexU ^(1-alphaUU)) ) ;
stockcapitalflexU - etaSE*kbarSflexU - (1-etaSE)*kbarBflexU;
ntotflexU - etaSE*nSflexU - (1-etaSE)*nBflexU;
btotBflexU - qflexU*kbarBflexU + nBflexU;
levBflexU - qflexU*kbarBflexU/nBflexU;//as in BGG
levSflexU - qflexU*kbarSflexU/nSflexU;//as in BGG
levAflexU - etaSE*levSflexU - (1-etaSE)*levBflexU;
btotSflexU - qflexU*kbarSflexU + nSflexU;
btotAflexU - etaSE*btotSflexU - (1-etaSE)*btotBflexU;
CostFinancingSflexU - omegabarSflexU(+1)*(1+RkSXflexU(+1))*qflexU*kbarSflexU/(qflexU*kbarSflexU-nSflexU);
SpreadflexU - CostFinancingSflexU + RcouponXflexU;
PextSflexU + (1+ReXflexU) - omegabarSflexU(+1)*(1+RkSXflexU(+1))*qflexU*kbarSflexU/(qflexU*kbarSflexU-nSflexU);


end;


// Initial values, coinciding with steady-state values
initval;
piU = piUU;
rkSU = rkSUU;
iU = iUU;
omegabarSU = omegabarSUU;
RkSXU = RkSXUU;
nSU = nSUU;
qU = qUU;
cU = cUU;
wU = wUU;
hU = hUU;
kbarSU = kbarSUU;
kbarBU = kbarBUU;
ReXU = ReXUU;
FpXU = FpXUU;
FwXU = FwXUU;
KwXU = KwXUU;
KpXU = KpXUU;
RcouponXU = RcouponXUU;
chiU=0;
KtotU=KtotUU;
YU=YUU;
eps_analystU=eps_analystUU;
eps_couponU=eps_couponUU;
RanalystXU=RanalystXUU;
nBU=nBUU;
stockcapitalU=stockcapitalUU;
ntotU=ntotUU;
uSU=uSUU;
uBU=uBUU;
RkBXU=RkBXUU;
levBU = levBUU;
levSU = levSUU;
btotBU = btotBUU;
btotSU = btotSUU;
levAU= levAUU;
btotAU=btotAUU;
CostFinancingSU = CostFinancingSUU;
SpreadU= SpreadUU;
lambdanU=lambdanUU;
PextSU=PextSUU;

piflexU = piUU;
rkSflexU = rkSUU;
iflexU = iUU;
omegabarSflexU = omegabarSUU;
RkSXflexU = RkSXUU;
nSflexU = nSUU;
qflexU = qUU;
cflexU = cUU;
wflexU = wUU;
hflexU = hUU;
kbarSflexU = kbarSUU;
kbarBflexU = kbarBUU;
ReXflexU = ReXUU;
FpXflexU = FpXflexUU;
FwXflexU = FwXflexUU;
KwXflexU = KwXflexUU;
KpXflexU = KpXflexUU;
RcouponXflexU = RcouponXUU;
chiflexU=0;
KtotflexU=KtotUU;
YflexU=YUU;
eps_analystflexU=eps_analystUU;
eps_couponflexU=eps_couponUU;
RanalystXflexU=RanalystXUU;
nBflexU=nBUU;
stockcapitalflexU=stockcapitalUU;
ntotflexU=ntotUU;
uSflexU=uSUU;
uBflexU=uBUU;
RkBXflexU=RkBXUU;
levBflexU = levBUU;
levSflexU = levSUU;
btotBflexU = btotBUU;
btotSflexU = btotSUU;
levAflexU= levAUU;
btotAflexU=btotAUU;
CostFinancingSflexU = CostFinancingSUU;
SpreadflexU= SpreadUU;
lambdanflexU=lambdanUU;
PextSflexU=PextSUU;



end;


shocks;
var e_xpU; stderr std1_xpU;

end;

//check;
//steady(solve_algo=4);

//stoch_simul(order=1,irf=20) YU cU iU hU ReXU piU nSU btotAU RkXU;

//stoch_simul(order=1,irf=40) YU cU iU YflexU cflexU iflexU ReXU piU piflexU;
stoch_simul (AR=100,IRF=0, noprint,nograph);