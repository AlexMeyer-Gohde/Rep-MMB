//-------------------------------------------------------------------------------------------------------------------------
//      THIS DYNARE FILE REPLICATES
//
//      Christiano, Motto and Rostagno (2010), "Financial Factors in Economic Fluctuations", ECB WP N�1192 - BASELINE MODEL
//
//      ADAPTED for the Macro Model Data Base
//
//      THIS FILE WAS WRITTEN FOR DYNARE 4.2.5 AND MATLAB 7.12 BY:
//
//      Fabio Verona, email: f8fabio@yahoo.it
//
//      August 2012
//
//---------------------------------------------------------------------------------------------------------------------------

// Variables
var   piU, sU, rkU, iU, uU, omegabarU, RkXU, nU, qU, vlU, evU, mbU, RXU, uzcU, lambdazU, mU, RaXU, cU, wU,
      hU, kbarU, ReXU, xU, dmU, RmXU,  pstarU, wstarU, FpXU, FwXU,YU, btotU,
      lambdafU, pitargetU, xbU, muupU, chiiU, gU, muzstarU, gammaU, epsilU,sigmaU, zetacU, zetaiU, tauoU

      piflexU, sflexU, rkflexU, iflexU, uflexU, omegabarflexU, RkXflexU, nflexU, qflexU, vlflexU, evflexU, mbflexU, RXflexU, uzcflexU, lambdazflexU, mflexU, RaXflexU, cflexU, wflexU,
      hflexU, kbarflexU, ReXflexU, xflexU, dmflexU, RmXflexU, FpXflexU, FwXflexU,YflexU, btotflexU;


varexo e_lambdafU, e_pitargetU, e_xbU, e_muupU, e_chiiU, e_gU, e_muzstarU, e_gammaU, e_epsilU, e_sigmaU,  e_zetaiU, e_tauoU e_xpU;// e_zetacU,




parameters btotUU, YUU, piUU, sUU, rkUU, iUU, uUU, omegabarUU, RkXUU, nUU, qUU, vlUU, evUU, mbUU, RXUU, uzcUU, lambdazUU, mUU, RaXUU, cUU, wUU, hUU, kbarUU, ReXUU, xUU, dmUU, RmXUU, pstarUU, wstarUU, FpXUU, FwXUU,  lambdafUU, pitargetUU, xbUU, muupUU, chiiUU, gUU, muzstarUU, gammaUU, epsilUU, sigmaUU, zetacUU, zetaiUU, tauoUU

           FpXflexUU,FwXflexUU,

           rho1_lambdafU, rho1_pitargetU, rho1_xbU, rho1_muupU, rho1_chiiU, rho1_gU, rho1_muzstarU, rho1_gammaU, rho1_epsilU, rho1_sigmaU,  rho1_zetaiU, rho1_tauoU, rho1_zetacU
           std1_lambdafU, std1_pitargetU, std1_xbU, std1_muupU, std1_chiiU, std1_gU, std1_muzstarU, std1_gammaU, std1_epsilU, std1_sigmaU,  std1_zetaiU, std1_tauoU, std1_zetacU

           lambdawUU, sigmaLXUU, thetaUU, vUU, betaUU, xiwUU, HdouprXUU, sigmaqUU, bUU,
           psiLXUU, iotaw1UU, iotamuUU, muUU, weUU, bigthetaUU, sigmaaUU, SdouprXUU,
           upsilUU, alphaUU, psikUU, psilUU, deltaUU, xipUU, iota1UU, ksiUU, chorseUU, etagUU,
           taukUU, taulUU, taucUU, tauUU, aptilUU, aytilUU, rhotilUU, phiUU, pibarUU;




// load the steady state values
pp=load('ss_CMR_baseline_US.mat');

%for i=1:M_.param_nbr;
    %name = deblank(M_.param_names(i,:));
    %if isfield(pp,name);
     %  M_.params(i) = eval(['pp.' name]);
      // disp([name ' ' num2str(M_.params(i))]);
   % end;
%end;

load ss_CMR_baseline_US
for jj=1:length(M_.param_names)
set_param_value(M_.param_names{jj},eval(M_.param_names{jj})); 
end;



model;



// A.1: A measure of marginal cost
-sU+((1-alphaUU)^(alphaUU-1)*alphaUU^(-alphaUU))*((rkU*(1+psikUU*RXU))^alphaUU)*((wU*(1+psilUU*RXU))^(1-alphaUU))/epsilU;

// A.2: Another measure of marginal cost
rkU*(1+psikUU*RXU)/( alphaUU*epsilU*( (upsilUU*muzstarU*( hU*(wstarU ^ (lambdawUU/(lambdawUU-1)) ) )/(uU*kbarU(-1)))^(1-alphaUU) ) )-sU;

// A.3, A.4, and A.5: Conditions associated with Calvo sticky prices
pstarU - ((1-xipUU)*( ( ( (1 - xipUU*(( (pitargetU^iota1UU) * (piU(-1)^(1-iota1UU)) * pibarUU^(1-iota1UU-(1-iota1UU)) )/piU)^(1/(1-lambdafU)) ) / (1-xipUU) ) ) ^ lambdafU ) + xipUU*(((( (pitargetU^iota1UU) * (piU(-1)^(1-iota1UU)) * pibarUU^(1-iota1UU-(1-iota1UU)) )/piU)*pstarU(-1))^( lambdafU / (1-lambdafU) )) )^(1/( lambdafU / (1-lambdafU) ));
lambdazU*( (pstarU^(lambdafU/(lambdafU-1))) * ( epsilU * vlU * ( (uU*kbarU(-1)/(muzstarU*upsilUU))^alphaUU ) * ( ( hU*(wstarU ^ (lambdawUU/(lambdawUU-1)) ) )^(1-alphaUU) ) - phiUU ) ) + ( (( (pitargetU(+1)^iota1UU) * (piU^(1-iota1UU)) * pibarUU^(1-iota1UU-(1-iota1UU)) )/piU(+1)) ^ (1/(1-lambdafU(+1))) ) * betaUU * xipUU * FpXU(+1) - FpXU;
lambdafU*lambdazU*( (pstarU^(lambdafU/(lambdafU-1))) * ( epsilU * vlU * ( (uU*kbarU(-1)/(muzstarU*upsilUU))^alphaUU ) * ( ( hU*(wstarU ^ (lambdawUU/(lambdawUU-1)) ) )^(1-alphaUU) ) - phiUU ) )*sU + betaUU * xipUU *  ( (( (pitargetU(+1)^iota1UU) * (piU^(1-iota1UU)) * pibarUU^(1-iota1UU-(1-iota1UU)) )/piU(+1)) ^ (lambdafU(+1)/(1-lambdafU(+1))) ) * ( FpXU(+1)* (( ( (1 - xipUU*(( (pitargetU(+1)^iota1UU) * (piU^(1-iota1UU)) * pibarUU^(1-iota1UU-(1-iota1UU)) )/piU(+1))^(1/(1-lambdafU(+1))) ) / (1-xipUU) ) )^(1-lambdafU(+1))) ) - ( FpXU* (( ( (1 - xipUU*(( (pitargetU^iota1UU) * (piU(-1)^(1-iota1UU)) * pibarUU^(1-iota1UU-(1-iota1UU)) )/piU)^(1/(1-lambdafU)) ) / (1-xipUU) ) )^(1-lambdafU)) );

// A.6: production function
YU - (pstarU^(lambdafU/(lambdafU-1))) * ( epsilU * vlU * ( (uU*kbarU(-1)/(muzstarU*upsilUU))^alphaUU ) * ( ( hU*(wstarU ^ (lambdawUU/(lambdawUU-1)) ) )^(1-alphaUU) ) - phiUU );

// A.7: supply of capital
lambdazU*qU*( -( ( sqrt(SdouprXUU/2) )*exp(( sqrt(SdouprXUU/2) )*(( zetaiU*iU*muzstarU*upsilUU/iU(-1) )-muzstarUU*upsilUU))-( sqrt(SdouprXUU/2) )*exp(-( sqrt(SdouprXUU/2) )*(( zetaiU*iU*muzstarU*upsilUU/iU(-1) )-muzstarUU*upsilUU)) )*( zetaiU*iU*muzstarU*upsilUU/iU(-1) )+1-( exp(( sqrt(SdouprXUU/2) )*(( zetaiU*iU*muzstarU*upsilUU/iU(-1) )-muzstarUU*upsilUU))+exp(-( sqrt(SdouprXUU/2) )*(( zetaiU*iU*muzstarU*upsilUU/iU(-1) )-muzstarUU*upsilUU))-2 ) )-(lambdazU/muupU)+betaUU*lambdazU(+1)*qU(+1)*( ( ( sqrt(SdouprXUU/2) )*exp(( sqrt(SdouprXUU/2) )*(( zetaiU(+1)*iU(+1)*muzstarU(+1)*upsilUU/iU )-muzstarUU*upsilUU))-( sqrt(SdouprXUU/2) )*exp(-( sqrt(SdouprXUU/2) )*(( zetaiU(+1)*iU(+1)*muzstarU(+1)*upsilUU/iU )-muzstarUU*upsilUU)) )*(( zetaiU(+1)*iU(+1)*muzstarU(+1)*upsilUU/iU )^2) )/(muzstarU(+1)*upsilUU);

// A.8: capital accumulation
kbarU-(1-deltaUU)*(kbarU(-1)/(muzstarU*upsilUU))-(1-( exp(( sqrt(SdouprXUU/2) )*(( zetaiU*iU*muzstarU*upsilUU/iU(-1) )-muzstarUU*upsilUU))+exp(-( sqrt(SdouprXUU/2) )*(( zetaiU*iU*muzstarU*upsilUU/iU(-1) )-muzstarUU*upsilUU))-2 ))*iU;

// A.9: capital utilization
rkU - tauoU*rkUU*exp(sigmaaUU*(uU-1));

// A.10: rate of return on capital
RkXU - ( ( (1-taukUU)*(uU*rkU-tauoU*( rkUU*(exp(sigmaaUU*(uU-1))-1)/sigmaaUU ))+(1-deltaUU)*qU)*piU/(upsilUU*qU(-1)) + taukUU*deltaUU - 1);

// A.11: standard debt contract
(1-( normcdf(( (log(omegabarU(+1))+sigmaU^2/2)/sigmaU )-sigmaU)+omegabarU(+1)*(1-normcdf(( (log(omegabarU(+1))+sigmaU^2/2)/sigmaU ))) ))*( (1+RkXU(+1))/(1+ReXU) )+( (1-( normcdf((log(omegabarU(+1))+sigmaU^2/2)/sigmaU) ))/(1-( normcdf((log(omegabarU(+1))+sigmaU^2/2)/sigmaU) )-muUU*omegabarU(+1)*( normpdf((log(omegabarU(+1))+sigmaU^2/2)/sigmaU)/omegabarU(+1)/sigmaU )) )*(( (1+RkXU(+1))/(1+ReXU) )*( (1-muUU)*normcdf(( (log(omegabarU(+1))+sigmaU^2/2)/sigmaU )-sigmaU)+omegabarU(+1)*(1-normcdf(( (log(omegabarU(+1))+sigmaU^2/2)/sigmaU ))) )-1);

// A.12: zero profit condition
(qU(-1)*kbarU(-1)*(1+RkXU)*( (1-muUU)*normcdf(( (log(omegabarU)+sigmaU(-1)^2/2)/sigmaU(-1) )-sigmaU(-1))+omegabarU*(1-normcdf(( (log(omegabarU)+sigmaU(-1)^2/2)/sigmaU(-1) ))) )/(nU(-1)*(1+ReXU(-1))))-(qU(-1)*kbarU(-1)/nU(-1))+1;

//A.13: law of motion for net worth
-nU+(gammaU/(piU*muzstarU))*(RkXU-ReXU(-1)-( ( normcdf(( (log(omegabarU)+sigmaU(-1)^2/2)/sigmaU(-1) )-sigmaU(-1))+omegabarU*(1-normcdf(( (log(omegabarU)+sigmaU(-1)^2/2)/sigmaU(-1) ))) )-( (1-muUU)*normcdf(( (log(omegabarU)+sigmaU(-1)^2/2)/sigmaU(-1) )-sigmaU(-1))+omegabarU*(1-normcdf(( (log(omegabarU)+sigmaU(-1)^2/2)/sigmaU(-1) ))) ) )*(1+RkXU))*kbarU(-1)*qU(-1)+weUU+gammaU*(1+ReXU(-1))*nU(-1)/(piU*muzstarU);

// A.14: banking services production function
xbU*( ( (evU^(-ksiUU)) ) )*( (1-tauUU)*(mbU(-1)/(piU*muzstarU))*(1-mU)-tauUU*( psilUU*wU*( hU*(wstarU ^ (lambdawUU/(lambdawUU-1)) ) )+psikUU*rkU*kbarU(-1)*uU/(muzstarU*upsilUU) ) )-(( (mbU(-1)/(piU*muzstarU))*(1-mU+chorseUU*dmU(-1)) )+( psilUU*wU*( hU*(wstarU ^ (lambdawUU/(lambdawUU-1)) ) )+psikUU*rkU*kbarU(-1)*uU/(muzstarU*upsilUU) ));

// A.15: ratio of bank excess reserves to their value-added
-evU+((1-tauUU)*(mbU(-1)/(piU*muzstarU))*(1-mU)-tauUU*(psilUU*wU*( hU*(wstarU ^ (lambdawUU/(lambdawUU-1)) ) )+psikUU*rkU*kbarU(-1)*uU/(muzstarU*upsilUU)))/ ((1-vlU)*((uU*kbarU(-1)/(muzstarU*upsilUU))^alphaUU)*(( hU*(wstarU ^ (lambdawUU/(lambdawUU-1)) ) )^(1-alphaUU)));

// A.16: banking efficiency condition
RaXU-(((1-tauUU)*( (1-ksiUU)*xbU*( (evU^(-ksiUU)) ) )-1)/(tauUU*( (1-ksiUU)*xbU*( (evU^(-ksiUU)) ) )+1))*RXU;

// A.17: another banking efficiency condition
lambdazU(+1)/(muzstarU(+1)*piU(+1))*(ReXU-RmXU - chorseUU*RXU(+1)/(( (1-ksiUU)*xbU(+1)*evU(+1)^(-ksiUU) )*tauUU+1));

// A.18: choice of labor
RXU*( (1-alphaUU)*ksiUU*xbU*(( evU^(1-ksiUU) ))*((kbarU(-1)*uU/(muzstarU*upsilUU*( hU*(wstarU ^ (lambdawUU/(lambdawUU-1)) ) )))^alphaUU) )/(1+tauUU*( (1-ksiUU)*xbU*( (evU^(-ksiUU)) ) ))-(1+psilUU*RXU)*wU;

// A.19: marginal utility of consumption
uzcU-(muzstarU*zetacU/(cU*muzstarU-bUU*cU(-1)))+bUU*betaUU*zetacU(+1)/(cU(+1)*muzstarU(+1)-bUU*cU);

// A.20: consumption decision
uzcU -(( zetacU*vUU*( (1+taucUU)*cU*((1/mU)^((1-chiiU)*thetaUU)) * ((1/(1-mU))^((1-chiiU)*(1-thetaUU))) * ((1/dmU(-1))^chiiU)    )^(1-sigmaqUU) )/cU)*((piU*muzstarU/mbU(-1))^(1-sigmaqUU))-(1+taucUU)*lambdazU;

// A.21, A.22, and A.23: Conditions associated with Calvo sticky wages
wstarU - ( ((1-xiwUU)*( ( (1 - xiwUU * ( ( (( (pitargetU^iotaw1UU) * (piU(-1)^(1-iotaw1UU)) * (pibarUU^(1-iotaw1UU-(1-iotaw1UU))) )/( piU*muzstarU*wU/wU(-1) )) * (muzstarUU^(1-iotamuUU)) * (muzstarU^iotamuUU) ) ^ (1/(1-lambdawUU)) ) ) / (1-xiwUU) ) ^ lambdawUU ) + xiwUU*( (( (( (pitargetU^iotaw1UU) * (piU(-1)^(1-iotaw1UU)) * (pibarUU^(1-iotaw1UU-(1-iotaw1UU))) )/( piU*muzstarU*wU/wU(-1) )) * (muzstarUU^(1-iotamuUU)) * (muzstarU^iotamuUU) )*wstarU(-1))^( lambdawUU/(1-lambdawUU) ) ) )^(1/( lambdawUU/(1-lambdawUU) )) );
((wstarU^(lambdawUU/(lambdawUU-1)))*hU*(1-taulUU)*lambdazU/lambdawUU) + betaUU * xiwUU * (muzstarUU^((1-iotamuUU)/(1-lambdawUU))) * ( (( (pitargetU(+1)^iotaw1UU) * (piU^(1-iotaw1UU)) * (pibarUU^(1-iotaw1UU-(1-iotaw1UU))) )^(1/(1-lambdawUU))) / piU(+1) ) * ( (1/( piU(+1)*muzstarU(+1)*wU(+1)/wU )) ^ (lambdawUU/(1-lambdawUU)) ) * (muzstarU(+1)^(iotamuUU/(1-lambdawUU)-1)) * FwXU(+1) - FwXU;
zetacU*((wstarU^(lambdawUU/(lambdawUU-1))*hU)^(1+sigmaLXUU)) + betaUU * xiwUU * (( (( (pitargetU(+1)^iotaw1UU) * (piU^(1-iotaw1UU)) * (pibarUU^(1-iotaw1UU-(1-iotaw1UU))) )/( piU(+1)*muzstarU(+1)*wU(+1)/wU )) * (muzstarUU^(1-iotamuUU)) * (muzstarU(+1)^iotamuUU) )  ^ (lambdawUU*(1+sigmaLXUU)/(1-lambdawUU))) * ( ( ( (1 - xiwUU * ( ( (( (pitargetU(+1)^iotaw1UU) * (piU^(1-iotaw1UU)) * (pibarUU^(1-iotaw1UU-(1-iotaw1UU))) )/( piU(+1)*muzstarU(+1)*wU(+1)/wU )) * (muzstarUU^(1-iotamuUU)) * (muzstarU(+1)^iotamuUU) ) ^ (1/(1-lambdawUU)) ) ) / (1-xiwUU) ) ^ (1-lambdawUU*(1+sigmaLXUU)) ) * wU(+1) * FwXU(+1)/psiLXUU ) - ( ( ( (1 - xiwUU * ( ( (( (pitargetU^iotaw1UU) * (piU(-1)^(1-iotaw1UU)) * (pibarUU^(1-iotaw1UU-(1-iotaw1UU))) )/( piU*muzstarU*wU/wU(-1) )) * (muzstarUU^(1-iotamuUU)) * (muzstarU^iotamuUU) ) ^ (1/(1-lambdawUU)) ) ) / (1-xiwUU) )   ^ (1-lambdawUU*(1+sigmaLXUU)) ) * wU   * FwXU/psiLXUU );

// A.24: choice of Tt
-lambdazU + (betaUU*lambdazU(+1)/(muzstarU(+1)*piU(+1)))*(1+ReXU);

// A.25: choice of Mt
(zetacU*vUU*( (1+taucUU)*cU*((1/mU)^((1-chiiU)*thetaUU)) * ((1/(1-mU))^((1-chiiU)*(1-thetaUU))) * ((1/dmU(-1))^chiiU)    )^(1-sigmaqUU) )*((thetaUU*(1-chiiU)/mU)-(1-chiiU)*(1-thetaUU)/(1-mU))*(((piU*muzstarU)/mbU(-1))^(2-sigmaqUU)) - lambdazU*RaXU -zetacU*( ( sqrt(HdouprXUU/2) )*exp(( sqrt(HdouprXUU/2) )*(( piU(-1)*muzstarU(-1)*mbU(-1)*mU/(mbU(-2)*mU(-1)) )-muzstarUU*piUU))-( sqrt(HdouprXUU/2) )*exp(-( sqrt(HdouprXUU/2) )*(( piU(-1)*muzstarU(-1)*mbU(-1)*mU/(mbU(-2)*mU(-1)) )-muzstarUU*piUU)) )*muzstarU*piU(-1)*piU*muzstarU(-1)/(mbU(-2)*mU(-1)) +betaUU*zetacU(+1)*( ( sqrt(HdouprXUU/2) )*exp(( sqrt(HdouprXUU/2) )*(( piU*muzstarU*mbU*mU(+1)/(mbU(-1)*mU) )-muzstarUU*piUU)) - ( sqrt(HdouprXUU/2) )*exp(-( sqrt(HdouprXUU/2) )*(( piU*muzstarU*mbU*mU(+1)/(mbU(-1)*mU) )-muzstarUU*piUU)) )*((muzstarU*piU)^2)*mbU*mU(+1)/((mbU(-1)*mU)^2);

// A.26: choice of Dm
(betaUU*chiiU(+1)*zetacU(+1)*vUU*(( (1+taucUU)*cU(+1) * ((1/mU(+1))^((1-chiiU(+1))*thetaUU)) * ((1/(1-mU(+1)))^((1-chiiU(+1))*(1-thetaUU))) *(1/dmU)^chiiU(+1) )^(1-sigmaqUU) )*((1/mbU)^(2-sigmaqUU))*((piU(+1)*muzstarU(+1))^(1-sigmaqUU))/dmU) + betaUU*lambdazU(+1)*(1+RmXU)/(piU(+1)*muzstarU(+1)) - lambdazU;

// A.27: choice of Mb
(( zetacU(+1)*vUU*( (1+taucUU)*cU(+1) * ((1/mU(+1))^((1-chiiU(+1))*thetaUU)) * ((1/(1-mU(+1)))^((1-chiiU(+1))*(1-thetaUU))) * (1/dmU)^chiiU(+1) )^(1-sigmaqUU) )*betaUU*(1-thetaUU)*(1-chiiU(+1))*(mbU^(sigmaqUU-2))*((piU(+1)*muzstarU(+1))^(1-sigmaqUU))/(1-mU(+1))  ) +(betaUU*lambdazU(+1)/(piU(+1)*muzstarU(+1)))*(1+RaXU(+1)) - lambdazU;

// A.28: policy rule
log(ReXU)-log(ReXUU) - rhotilUU*(log(ReXU(-1)) -log(ReXUU)) -(1/ReXUU)*(1-rhotilUU)*piUU*( log(pitargetU) -log(piUU) ) - (1/ReXUU)*(1-rhotilUU)*aptilUU*piUU*(log(piU(+1)) - log(pitargetU)) - (1/(4*ReXUU))*(1-rhotilUU)*aytilUU*( cUU/( (cUU+iUU)/(1-etagUU) )*log(cU)+iUU/( (cUU+iUU)/(1-etagUU) )*log(iU)-iUU/( (cUU+iUU)/(1-etagUU) )*log(muupU)+gUU/( (cUU+iUU)/(1-etagUU) )*log(gU)  - (cUU/( (cUU+iUU)/(1-etagUU) )*log(cUU)+iUU/( (cUU+iUU)/(1-etagUU) )*log(iUU)+gUU/( (cUU+iUU)/(1-etagUU) )*log(gUU)) )   - (1/(400*ReXUU))*(e_xpU);

// A.29: law of motion monetary base
mbU - mbU(-1)*(1+xU)/(piU*muzstarU);

// A.30: resource constraint
((pstarU^(lambdafU/(lambdafU-1))) * ( epsilU * vlU * ( (uU*kbarU(-1)/(muzstarU*upsilUU))^alphaUU ) * ( ( hU*(wstarU ^ (lambdawUU/(lambdawUU-1)) ) )^(1-alphaUU) ) - phiUU ) )-gU-cU-iU/muupU-tauoU*( rkUU*(exp(sigmaaUU*(uU-1))-1)/sigmaaUU )*kbarU(-1)/(muzstarU*upsilUU)-( ( ( normcdf(( (log(omegabarU)+sigmaU(-1)^2/2)/sigmaU(-1) )-sigmaU(-1))+omegabarU*(1-normcdf(( (log(omegabarU)+sigmaU(-1)^2/2)/sigmaU(-1) ))) )-( (1-muUU)*normcdf(( (log(omegabarU)+sigmaU(-1)^2/2)/sigmaU(-1) )-sigmaU(-1))+omegabarU*(1-normcdf(( (log(omegabarU)+sigmaU(-1)^2/2)/sigmaU(-1) ))) ) )*(1+RkXU)*qU(-1)*kbarU(-1)/(muzstarU*piU) )-bigthetaUU*(1-gammaU)*(nU-weUU)/gammaU;

// A.32: total bank loans
btotU - psilUU*wU*( hU*(wstarU ^ (lambdawUU/(lambdawUU-1)) ) ) -psikUU*rkU*uU*kbarU(-1)/(muzstarU*upsilUU) - (qU(-1)*kbarU(-1)-nU(-1))/(piU*muzstarU);


// FLEXIBLE PRICE EQUILIBRIUM EQUATIONS
-sflexU+((1-alphaUU)^(alphaUU-1)*alphaUU^(-alphaUU))*((rkflexU*(1+psikUU*RXflexU))^alphaUU)*((wflexU*(1+psilUU*RXflexU))^(1-alphaUU))/epsilU;
rkflexU*(1+psikUU*RXflexU)/( alphaUU*epsilU*( (upsilUU*muzstarU*hflexU/(uflexU*kbarflexU(-1)))^(1-alphaUU) ) )-sflexU;
lambdazflexU*(( epsilU * vlflexU * ( (uflexU*kbarflexU(-1)/(muzstarU*upsilUU))^alphaUU ) * ( hflexU^(1-alphaUU) ) - phiUU ) ) - FpXflexU;
lambdafU*lambdazflexU*(( epsilU * vlflexU * ( (uflexU*kbarflexU(-1)/(muzstarU*upsilUU))^alphaUU ) * (hflexU^(1-alphaUU)) - phiUU ) )*sflexU - FpXflexU;
YflexU - ( epsilU * vlflexU * ( (uflexU*kbarflexU(-1)/(muzstarU*upsilUU))^alphaUU ) * (hflexU^(1-alphaUU) ) - phiUU );
lambdazflexU*qflexU*( -( ( sqrt(SdouprXUU/2) )*exp(( sqrt(SdouprXUU/2) )*(( zetaiU*iflexU*muzstarU*upsilUU/iflexU(-1) )-muzstarUU*upsilUU))-( sqrt(SdouprXUU/2) )*exp(-( sqrt(SdouprXUU/2) )*(( zetaiU*iflexU*muzstarU*upsilUU/iflexU(-1) )-muzstarUU*upsilUU)) )*( zetaiU*iflexU*muzstarU*upsilUU/iflexU(-1) )+1-( exp(( sqrt(SdouprXUU/2) )*(( zetaiU*iflexU*muzstarU*upsilUU/iflexU(-1) )-muzstarUU*upsilUU))+exp(-( sqrt(SdouprXUU/2) )*(( zetaiU*iflexU*muzstarU*upsilUU/iflexU(-1) )-muzstarUU*upsilUU))-2 ) )-(lambdazflexU/muupU)+betaUU*lambdazflexU(+1)*qflexU(+1)*( ( ( sqrt(SdouprXUU/2) )*exp(( sqrt(SdouprXUU/2) )*(( zetaiU(+1)*iflexU(+1)*muzstarU(+1)*upsilUU/iflexU )-muzstarUU*upsilUU))-( sqrt(SdouprXUU/2) )*exp(-( sqrt(SdouprXUU/2) )*(( zetaiU(+1)*iflexU(+1)*muzstarU(+1)*upsilUU/iflexU )-muzstarUU*upsilUU)) )*(( zetaiU(+1)*iflexU(+1)*muzstarU(+1)*upsilUU/iflexU )^2) )/(muzstarU(+1)*upsilUU);
kbarflexU-(1-deltaUU)*(kbarflexU(-1)/(muzstarU*upsilUU))-(1-( exp(( sqrt(SdouprXUU/2) )*(( zetaiU*iflexU*muzstarU*upsilUU/iflexU(-1) )-muzstarUU*upsilUU))+exp(-( sqrt(SdouprXUU/2) )*(( zetaiU*iflexU*muzstarU*upsilUU/iflexU(-1) )-muzstarUU*upsilUU))-2 ))*iflexU;
rkflexU - tauoU*rkUU*exp(sigmaaUU*(uflexU-1));
RkXflexU - ( ( (1-taukUU)*(uflexU*rkflexU-tauoU*( rkUU*(exp(sigmaaUU*(uflexU-1))-1)/sigmaaUU ))+(1-deltaUU)*qflexU)*piflexU/(upsilUU*qflexU(-1)) + taukUU*deltaUU - 1);
(1-( normcdf(( (log(omegabarflexU(+1))+sigmaU^2/2)/sigmaU )-sigmaU)+omegabarflexU(+1)*(1-normcdf(( (log(omegabarflexU(+1))+sigmaU^2/2)/sigmaU ))) ))*( (1+RkXflexU(+1))/(1+ReXflexU) )+( (1-( normcdf((log(omegabarflexU(+1))+sigmaU^2/2)/sigmaU) ))/(1-( normcdf((log(omegabarflexU(+1))+sigmaU^2/2)/sigmaU) )-muUU*omegabarflexU(+1)*( normpdf((log(omegabarflexU(+1))+sigmaU^2/2)/sigmaU)/omegabarflexU(+1)/sigmaU )) )*(( (1+RkXflexU(+1))/(1+ReXflexU) )*( (1-muUU)*normcdf(( (log(omegabarflexU(+1))+sigmaU^2/2)/sigmaU )-sigmaU)+omegabarflexU(+1)*(1-normcdf(( (log(omegabarflexU(+1))+sigmaU^2/2)/sigmaU ))) )-1);
(qflexU(-1)*kbarflexU(-1)*(1+RkXflexU)*( (1-muUU)*normcdf(( (log(omegabarflexU)+sigmaU(-1)^2/2)/sigmaU(-1) )-sigmaU(-1))+omegabarflexU*(1-normcdf(( (log(omegabarflexU)+sigmaU(-1)^2/2)/sigmaU(-1) ))) )/(nflexU(-1)*(1+ReXflexU(-1))))-(qflexU(-1)*kbarflexU(-1)/nflexU(-1))+1;
-nflexU+(gammaU/(piflexU*muzstarU))*(RkXflexU-ReXflexU(-1)-( ( normcdf(( (log(omegabarflexU)+sigmaU(-1)^2/2)/sigmaU(-1) )-sigmaU(-1))+omegabarflexU*(1-normcdf(( (log(omegabarflexU)+sigmaU(-1)^2/2)/sigmaU(-1) ))) )-( (1-muUU)*normcdf(( (log(omegabarflexU)+sigmaU(-1)^2/2)/sigmaU(-1) )-sigmaU(-1))+omegabarflexU*(1-normcdf(( (log(omegabarflexU)+sigmaU(-1)^2/2)/sigmaU(-1) ))) ) )*(1+RkXflexU))*kbarflexU(-1)*qflexU(-1)+weUU+gammaU*(1+ReXflexU(-1))*nflexU(-1)/(piflexU*muzstarU);
xbU*( ( (evflexU^(-ksiUU)) ) )*( (1-tauUU)*(mbflexU(-1)/(piflexU*muzstarU))*(1-mflexU)-tauUU*( psilUU*wflexU*hflexU+psikUU*rkflexU*kbarflexU(-1)*uflexU/(muzstarU*upsilUU) ) )-(( (mbflexU(-1)/(piflexU*muzstarU))*(1-mflexU+chorseUU*dmflexU(-1)) )+( psilUU*wflexU*hflexU+psikUU*rkflexU*kbarflexU(-1)*uflexU/(muzstarU*upsilUU) ));
-evflexU+((1-tauUU)*(mbflexU(-1)/(piflexU*muzstarU))*(1-mflexU)-tauUU*(psilUU*wflexU*hflexU+psikUU*rkflexU*kbarflexU(-1)*uflexU/(muzstarU*upsilUU)))/ ((1-vlflexU)*((uflexU*kbarflexU(-1)/(muzstarU*upsilUU))^alphaUU)*(hflexU^(1-alphaUU)));
RaXflexU-(((1-tauUU)*( (1-ksiUU)*xbU*( (evflexU^(-ksiUU)) ) )-1)/(tauUU*( (1-ksiUU)*xbU*( (evflexU^(-ksiUU)) ) )+1))*RXflexU;
lambdazflexU(+1)/(muzstarU(+1)*piflexU(+1))*(ReXflexU-RmXflexU - chorseUU*RXflexU(+1)/(( (1-ksiUU)*xbU(+1)*evflexU(+1)^(-ksiUU) )*tauUU+1));
RXflexU*( (1-alphaUU)*ksiUU*xbU*(( evflexU^(1-ksiUU) ))*((kbarflexU(-1)*uflexU/(muzstarU*upsilUU*hflexU))^alphaUU) )/(1+tauUU*( (1-ksiUU)*xbU*( (evflexU^(-ksiUU)) ) ))-(1+psilUU*RXflexU)*wflexU;
uzcflexU-(muzstarU*zetacU/(cflexU*muzstarU-bUU*cflexU(-1)))+bUU*betaUU*zetacU(+1)/(cflexU(+1)*muzstarU(+1)-bUU*cflexU);
uzcflexU-(( zetacU*vUU*( (1+taucUU)*cflexU*((1/mflexU)^((1-chiiU)*thetaUU)) * ((1/(1-mflexU))^((1-chiiU)*(1-thetaUU))) * ((1/dmflexU(-1))^chiiU)    )^(1-sigmaqUU) )/cflexU)*((piflexU*muzstarU/mbflexU(-1))^(1-sigmaqUU))-(1+taucUU)*lambdazflexU;
hflexU*(1-taulUU)*lambdazflexU/lambdawUU - FwXflexU;
zetacU*(hflexU^(1+sigmaLXUU)) - wflexU* FwXflexU/psiLXUU;
-lambdazflexU + (betaUU*lambdazflexU(+1)/(muzstarU(+1)*piflexU(+1)))*(1+ReXflexU);
(zetacU*vUU*( (1+taucUU)*cflexU*((1/mflexU)^((1-chiiU)*thetaUU)) * ((1/(1-mflexU))^((1-chiiU)*(1-thetaUU))) * ((1/dmflexU(-1))^chiiU)    )^(1-sigmaqUU) )*((thetaUU*(1-chiiU)/mflexU)-(1-chiiU)*(1-thetaUU)/(1-mflexU))*(((piflexU*muzstarU)/mbflexU(-1))^(2-sigmaqUU)) - lambdazflexU*RaXflexU -zetacU*( ( sqrt(HdouprXUU/2) )*exp(( sqrt(HdouprXUU/2) )*(( piflexU(-1)*muzstarU(-1)*mbflexU(-1)*mflexU/(mbflexU(-2)*mflexU(-1)) )-muzstarUU*piUU))-( sqrt(HdouprXUU/2) )*exp(-( sqrt(HdouprXUU/2) )*(( piflexU(-1)*muzstarU(-1)*mbflexU(-1)*mflexU/(mbflexU(-2)*mflexU(-1)) )-muzstarUU*piUU)) )*muzstarU*piflexU(-1)*piflexU*muzstarU(-1)/(mbflexU(-2)*mflexU(-1)) +betaUU*zetacU(+1)*( ( sqrt(HdouprXUU/2) )*exp(( sqrt(HdouprXUU/2) )*(( piflexU*muzstarU*mbflexU*mflexU(+1)/(mbflexU(-1)*mflexU) )-muzstarUU*piUU)) - ( sqrt(HdouprXUU/2) )*exp(-( sqrt(HdouprXUU/2) )*(( piflexU*muzstarU*mbflexU*mflexU(+1)/(mbflexU(-1)*mflexU) )-muzstarUU*piUU)) )*((muzstarU*piflexU)^2)*mbflexU*mflexU(+1)/((mbflexU(-1)*mflexU)^2);
(betaUU*chiiU(+1)*  zetacU(+1)*vUU*(((1+taucUU)*cflexU(+1) * ((1/mflexU(+1))^((1-chiiU(+1))*thetaUU)) * ((1/(1-mflexU(+1)))^((1-chiiU(+1))*(1-thetaUU))) *(1/dmflexU)^chiiU(+1) )^(1-sigmaqUU) )*((1/mbflexU)^(2-sigmaqUU))*((piflexU(+1)*muzstarU(+1))^(1-sigmaqUU))/dmflexU) + betaUU*lambdazflexU(+1)*(1+RmXflexU)/(piflexU(+1)*muzstarU(+1)) - lambdazflexU;
(( zetacU(+1)*vUU*( (1+taucUU)*cflexU(+1) * ((1/mflexU(+1))^((1-chiiU(+1))*thetaUU)) * ((1/(1-mflexU(+1)))^((1-chiiU(+1))*(1-thetaUU))) * (1/dmflexU)^chiiU(+1) )^(1-sigmaqUU) )*betaUU*(1-thetaUU)*(1-chiiU(+1))*(mbflexU^(sigmaqUU-2))*((piflexU(+1)*muzstarU(+1))^(1-sigmaqUU))/(1-mflexU(+1))  ) +(betaUU*lambdazflexU(+1)/(piflexU(+1)*muzstarU(+1)))*(1+RaXflexU(+1)) - lambdazflexU;
log(ReXflexU)-log(ReXUU) - rhotilUU*(log(ReXflexU(-1)) -log(ReXUU)) -(1/ReXUU)*(1-rhotilUU)*piUU*( log(pitargetU) -log(piUU) ) - (1/ReXUU)*(1-rhotilUU)*aptilUU*piUU*(log(piflexU(+1)) - log(pitargetU)) - (1/(4*ReXUU))*(1-rhotilUU)*aytilUU*( cUU/( (cUU+iUU)/(1-etagUU) )*log(cflexU)+iUU/( (cUU+iUU)/(1-etagUU) )*log(iflexU)-iUU/( (cUU+iUU)/(1-etagUU) )*log(muupU)+gUU/( (cUU+iUU)/(1-etagUU) )*log(gU)  - (cUU/( (cUU+iUU)/(1-etagUU) )*log(cUU)+iUU/( (cUU+iUU)/(1-etagUU) )*log(iUU)+gUU/( (cUU+iUU)/(1-etagUU) )*log(gUU)) );
mbflexU - mbflexU(-1)*(1+xflexU)/(piflexU*muzstarU);
( epsilU * vlflexU * ( (uflexU*kbarflexU(-1)/(muzstarU*upsilUU))^alphaUU ) * (hflexU^(1-alphaUU) ) - phiUU ) -gU-cflexU-iflexU/muupU-tauoU*( rkUU*(exp(sigmaaUU*(uflexU-1))-1)/sigmaaUU )*kbarflexU(-1)/(muzstarU*upsilUU)-( ( ( normcdf(( (log(omegabarflexU)+sigmaU(-1)^2/2)/sigmaU(-1) )-sigmaU(-1))+omegabarflexU*(1-normcdf(( (log(omegabarflexU)+sigmaU(-1)^2/2)/sigmaU(-1) ))) )-( (1-muUU)*normcdf(( (log(omegabarflexU)+sigmaU(-1)^2/2)/sigmaU(-1) )-sigmaU(-1))+omegabarflexU*(1-normcdf(( (log(omegabarflexU)+sigmaU(-1)^2/2)/sigmaU(-1) ))) ) )*(1+RkXflexU)*qflexU(-1)*kbarflexU(-1)/(muzstarU*piflexU) )-bigthetaUU*(1-gammaU)*(nflexU-weUU)/gammaU;
btotflexU - psilUU*wflexU*hflexU -psikUU*rkflexU*uflexU*kbarflexU(-1)/(muzstarU*upsilUU) - (qflexU(-1)*kbarflexU(-1)-nflexU(-1))/(piflexU*muzstarU);



// Shock Equations
// Equations take the form: x = x_ *(1+e_x) + rho_x*(x(-1) - x_),
// where x_ is steady-state value, and e_x is iid shock with stderror std_x

lambdafU = lambdafUU* ( 1 + e_lambdafU) + rho1_lambdafU * ( lambdafU(-1) - lambdafUU );
pitargetU = pitargetUU* ( 1 + e_pitargetU) + rho1_pitargetU * ( pitargetU(-1) - pitargetUU );
xbU = xbUU* ( 1 + e_xbU) + rho1_xbU * ( xbU(-1) - xbUU );
muupU = muupUU* ( 1 + e_muupU) + rho1_muupU * ( muupU(-1) - muupUU );
chiiU = chiiUU* ( 1 + e_chiiU) + rho1_chiiU * ( chiiU(-1) - chiiUU );
gU = gUU* ( 1 + e_gU) + rho1_gU * ( gU(-1) - gUU );
muzstarU = muzstarUU* ( 1 + e_muzstarU) + rho1_muzstarU * ( muzstarU(-1) - muzstarUU );
gammaU = gammaUU* ( 1 + e_gammaU) + rho1_gammaU * ( gammaU(-1) - gammaUU );
epsilU = epsilUU* ( 1 + e_epsilU) + rho1_epsilU * ( epsilU(-1) - epsilUU );
sigmaU = sigmaUU* ( 1 + e_sigmaU) + rho1_sigmaU * ( sigmaU(-1) - sigmaUU );
zetacU = zetacUU;
zetaiU = zetaiUU* ( 1 + e_zetaiU) + rho1_zetaiU * ( zetaiU(-1) - zetaiUU );
tauoU = tauoUU* ( 1 + e_tauoU) + rho1_tauoU * ( tauoU(-1) - tauoUU );

end;

// Initial values, coinciding with steady-state values
initval;

piU = piUU;
sU = sUU;
rkU = rkUU;
iU = iUU;
uU = uUU;
omegabarU = omegabarUU;
RkXU = RkXUU;
nU = nUU;
qU = qUU;
vlU = vlUU;
evU = evUU;
mbU = mbUU;
RXU = RXUU;
uzcU = uzcUU;
lambdazU = lambdazUU;
mU = mUU;
RaXU = RaXUU;
cU = cUU;
wU = wUU;
hU = hUU;
kbarU = kbarUU;
ReXU = ReXUU;
xU = xUU;
dmU = dmUU;
RmXU = RmXUU;
btotU = btotUU;
pstarU = pstarUU;
wstarU = wstarUU;
FpXU = FpXUU;
FwXU = FwXUU;
YU = YUU;

piflexU = piUU;
sflexU = sUU;
rkflexU = rkUU;
iflexU = iUU;
uflexU = uUU;
omegabarflexU = omegabarUU;
RkXflexU = RkXUU;
nflexU = nUU;
qflexU = qUU;
vlflexU = vlUU;
evflexU = evUU;
mbflexU = mbUU;
RXflexU = RXUU;
uzcflexU = uzcUU;
lambdazflexU = lambdazUU;
mflexU = mUU;
RaXflexU = RaXUU;
cflexU = cUU;
wflexU = wUU;
hflexU = hUU;
kbarflexU = kbarUU;
ReXflexU = ReXUU;
xflexU = xUU;
dmflexU = dmUU;
RmXflexU = RmXUU;
btotflexU = btotUU;
FpXflexU = FpXflexUU;
FwXflexU = FwXflexUU;
YflexU = YUU;

lambdafU = lambdafUU;
pitargetU = pitargetUU;
xbU = xbUU;
muupU = muupUU;
chiiU = chiiUU;
gU = gUU;
muzstarU = muzstarUU;
gammaU = gammaUU;
epsilU = epsilUU;
sigmaU = sigmaUU;
zetacU = zetacUU;
zetaiU = zetaiUU;
tauoU = tauoUU;


end;

shocks;
var e_xpU = 0.519^2;
var e_lambdafU; stderr std1_lambdafU;
var e_pitargetU; stderr std1_pitargetU;
var e_xbU; stderr std1_xbU;
var e_muupU; stderr std1_muupU;
var e_chiiU; stderr std1_chiiU;
var e_gU; stderr std1_gU;
var e_muzstarU; stderr std1_muzstarU;
var e_gammaU; stderr std1_gammaU;
var e_epsilU; stderr std1_epsilU;
var e_sigmaU; stderr std1_sigmaU;
//var e_zetacU; stderr std1_zetacU;
var e_zetaiU; stderr std1_zetaiU;
var e_tauoU; stderr std1_tauoU;

//var interest_ = 0;//var interest_; stderr 1;
end;
stoch_simul (AR=100,IRF=0, noprint,nograph);
//check;
//steady(solve_algo=4);

//stoch_simul(order=1,irf=20) YU cU iU hU ReXU piU nU btotU RkXU;

//some_plots
