//DEL NEGRO et al. 2015 Replication with both sticky and flexible price output in Dynare

// DEL NEGRO, M., M. GIANNONI & F. SCHORFHEIDE 2015.: Inflation in the Great Recession and New Keynesian Models,
// American Economic Journal: Macroeconomics
// 2015, 7(1): 168-196.

// This is the DNGS15 version of the SW (2007) model estimated with the same observables as Smets and
%Wouters plus credit spread data.

//------------------------------------------------------------------------------------------------------------------------
//1. Variable declaration
//------------------------------------------------------------------------------------------------------------------------
var c R pi $\pi$ L qk $q^k$ i Rktil $\tilde{R}^k$ rk $r^k$ kbar $\bar{k}$ n y k u mc w wh z ztil mu sigw laf law g b rm pist

//% the flexble-price counterparts
c_f r_f L_f qk_f i_f rk_f y_f k_f u_f kbar_f w_f

//% government spending shock
;



varexo psi_b psi_mu psi_z psi_laf psi_law  psi_sigw psi_pist psi_g psi_rm;



//------------------------------------------------------------------------------------------------------------------------
// 2. Parameter declaration and calibration
//-------------------------------------------------------------------------------------------------------------------------
parameters

alp $\alpha$ zeta_p $\zeta_p$ iota_p $\iota_p$ del $\delta$ Bigphi $\Phi$ s2 h ppsi $\psi$ nu_l $\nu_l$ zeta_w $\zeta_w$ iota_w $\iota_w$ bet $\beta$ psi1 $\psi_1$ psi2 $\psi_2$ psi3 $\psi_3$ sigmac $\sigma_c$ rho $\rho$ epsp epsw
star $g$ rho_g $\rho_g$ rho_b $\rho_b$ rho_mu $\rho_\mu$ rho_z $\rho_z$ rho_laf $\rho_{\lambda_f}$ rho_law $\rho_{\lambda_w}$ rho_rm $\rho_{r^m}$ rho_sigw $\rho_{\sigma_w}$ rho_pist $\rho_{\pi^*}$ eta_gz eta_laf eta_law zstar $\gamma$ rkstar $r^k$ wl_c $\frac{wl}{c}$
cstar $c$ wstar $w$ Lstar $L$ kstar $k$ kbarstar $\bar{k}$ istar $i$ rstar $r$ ystar $y$ zeta_spb $\zeta_{spb}$ gammstar vstar $v$ nstar $n$ zeta_nRk $\zeta_{nRk}$ zeta_nR $\zeta_{nR}$ zeta_nsigw $\zeta_{n\sigma_w}$ zeta_spsigw $\zeta_{sp\sigma_w}$
zeta_nqk $\zeta_{nqk}$ zeta_nn $\zeta_{nn}$ gstar;

alp = 0.2140;
zeta_p = 0.8063;
iota_p = 0.7721;
del = 0.025;
Bigphi = 1.6422;
s2 = 5.2330;
h =  0.2515;
ppsi = 0.7338;
nu_l = 2.1654;
zeta_w = 0.8586;
iota_w = 0.4766;
bet  = 0.9986;
psi1 = 1.2632;
psi2 = 0.0566;
psi3 = 0.2460;
sigmac = 1.6627;
rho = 0.7627;
epsp = 10;
epsw = 10;

//Financial Frictions Parameters
zeta_spb = 0.0500;
gammstar = 0.9900;
vstar = 3.8772;
nstar = 3.8533;
zeta_nRk = 1.3415;
zeta_nR = 0.3424;
zeta_nsigw = 0.0020;
zeta_spsigw = 0.0237;
zeta_nqk = 0.0030;
zeta_nn = 0.9962;

//exogenous processes - level
gstar = 0.1800;

//exogenous processes - autocorrelation

rho_g = 0.9758;
rho_b = 0.9927;
rho_mu = 0.6997;
rho_z = 0.9756;
rho_laf = 0.2600;
rho_law = 0.9770;
rho_rm = 0.0405;

rho_sigw = 0.9899;
rho_pist = 0.9900;

//exogenous processes
eta_gz = 0.8440;
eta_laf = 0.4035;
eta_law = 0.9633;

//Parameters (implicit) -- from steady state

zstar = 0.0034;
rstar = 1.0070;
rkstar = 0.0358;
wstar = 0.6805;
Lstar = 1;
kstar = 5.1777;
kbarstar = 5.1951;
istar = 0.1469;
ystar = 0.8658;
cstar = 0.5631;

wl_c = 0.8056;


//-----------------------------------------------------------------------------------------------------------------------
// 3. The model
//-----------------------------------------------------------------------------------------------------------------------

model(linear);



//consumption Euler equation
c = -(1-h*exp(-zstar))/(sigmac*(1+h*exp(-zstar)))*(R-pi(+1))+b+(h*exp(-zstar))/(1+h*exp(-zstar))*(c(-1)-z)
    + (1/(1+h*exp(-zstar)))*(c(+1)+(1/(1-alp))*(rho_z-1)*ztil)+(sigmac-1)*wl_c/(sigmac*(1+h*exp(-zstar)))*(L-L(+1));

//investment Euler equation
qk = (s2*exp(2*zstar)*(1+bet*exp((1-sigmac)*zstar)))*(i-(1/(1+bet*exp((1-sigmac)*zstar)))*(i(-1)-z)
     -(bet*exp((1-sigmac)*zstar)/(1+bet*exp((1-sigmac)*zstar)))*i(+1)
    -(bet*exp((1-sigmac)*zstar)/(1+bet*exp((1-sigmac)*zstar)))*(1/(1-alp))*(rho_z-1)*ztil-mu);

//evolution of capital
kbar = (1-istar/kbarstar)*(kbar(-1)-z)+(istar/kbarstar)*i+(istar*s2*exp(2*zstar)*(1+bet*exp((1-sigmac)*zstar))/kbarstar)*mu;

//capital utilization
k = u-z+kbar(-1);

//rate of utilization (rental rate of capital)
u = ((1-ppsi)/ppsi)*rk;

//price markup
mc = w+alp*L-alp*k;

//rental rate of capital
k = w-rk+L;

//aggregate production function
y = Bigphi*alp*k+Bigphi*(1-alp)*L+((Bigphi-1)/(1-alp))*ztil;

//ressource constraint
y = gstar*g+cstar/ystar*c+istar/ystar*i+rkstar*kstar/ystar*u-gstar*(1/(1-alp))*ztil;

//Phillips curve
pi = ((1-zeta_p*bet*exp((1-sigmac)*zstar))*(1-zeta_p))/(zeta_p*((Bigphi-1)*epsp+1))*1/(1+iota_p*bet*exp((1-sigmac)*zstar))*mc
    +iota_p*1/(1+iota_p*bet*exp((1-sigmac)*zstar))*pi(-1)+bet*exp((1-sigmac)*zstar)*1/(1+iota_p*bet*exp((1-sigmac)*zstar))*pi(+1)
    +laf;

//evolution of wages
w = ((1-zeta_w*bet*exp((1-sigmac)*zstar))*(1-zeta_w)/(zeta_w*((1.5-1)*epsw+1))*1/(1+bet*exp((1-sigmac)*zstar)))*(wh-w)
    -(1+iota_w*bet*exp((1-sigmac)*zstar))*1/(1+bet*exp((1-sigmac)*zstar))*pi
    +(1/(1+bet*exp((1-sigmac)*zstar)))*(w(-1)-z + iota_w*pi(-1))
    +(bet*exp((1-sigmac)*zstar)*1/(1+bet*exp((1-sigmac)*zstar)))*(w(+1)+ ( 1/(1-alp) )*(rho_z-1)*ztil +pi(+1))+law;

//marginal rate of substitution
wh = (1/(1-h*exp(-zstar)))*(c - h*exp(-zstar)*c(-1) + h*exp(-zstar)*z)+nu_l*L;

//monetary policy rule
R = rho*R(-1)+(1-rho)*psi1*(pi-pist)+(1-rho)*psi2*(y-y_f)+psi3*((y-y_f)-(y(-1)-y_f(-1)))+rm;

//Financial Frictions//
//return to capital
Rktil = pi+rkstar/(rkstar+1-del)*rk+(1-del)/(rkstar+1-del)*qk-qk(-1);
//spreads
Rktil(+1) = R-((sigmac*(1+h*exp(-zstar)))/(1-h*exp(-zstar)))*b+zeta_spb*(qk+kbar-n)+sigw;
//evolution of net worth
n = zeta_nRk*(Rktil-pi)-zeta_nR*(R(-1) - pi)+zeta_nqk*(qk(-1)+kbar(-1))+zeta_nn*n(-1)-zeta_nsigw/zeta_spsigw*sigw(-1)
    -gammstar*vstar/nstar*z;

//exogenous processes//

z = (1/(1-alp))*(rho_z-1)*ztil(-1)+(1/(1-alp))*psi_z;

ztil = rho_z*ztil(-1)+psi_z;

g = rho_g*g(-1)+psi_g+eta_gz*psi_z;

b = rho_b*b(-1)+psi_b;

mu = rho_mu*mu(-1)+psi_mu;

laf = rho_laf*laf(-1)+psi_laf-eta_laf*psi_laf(-1);

law = rho_law*law(-1)+psi_law-eta_law*psi_law(-1);

rm = rho_rm*rm(-1) + psi_rm;

sigw = rho_sigw*sigw(-1)+psi_sigw;

pist = rho_pist*pist(-1)+psi_pist;


// now the flexible-price equations//

c_f = -(1-h*exp(-zstar))/(sigmac*(1+h*exp(-zstar)))*r_f+b+(h*exp(-zstar))/(1+h*exp(-zstar))*(c_f(-1)-z)
    +(1/(1+h*exp(-zstar)))*(c_f(+1)+(1/(1-alp))*(rho_z-1)*ztil)+(sigmac-1)*wl_c/(sigmac*(1+ h*exp(-zstar)))*(L_f-L_f(+1));

qk_f = (s2*exp(2*zstar)*(1+bet*exp((1-sigmac)*zstar)))*(i_f-(1/(1+bet*exp((1-sigmac)*zstar)))*(i_f(-1)-z)
     -(bet*exp((1-sigmac)*zstar)/(1+bet*exp((1-sigmac)*zstar)))*i_f(+1)
    -(bet*exp((1-sigmac)*zstar)/(1+bet*exp((1-sigmac)*zstar)))*(1/(1-alp))*(rho_z-1)* ztil-mu);

kbar_f = (1-istar/kbarstar)*(kbar_f(-1)-z)+(istar/kbarstar)*i_f+(istar*s2*exp(2*zstar)*(1+bet*exp((1-sigmac)*zstar))/kbarstar)*mu;

k_f = u_f-z+kbar_f(-1);

u_f = ((1-ppsi)/ppsi)*rk_f;

w_f = - alp*L_f+alp*k_f;

k_f = w_f-rk_f+L_f;

y_f = Bigphi*alp*k_f+Bigphi*(1-alp)*L_f+((Bigphi-1)/(1-alp))*ztil;

y_f = gstar*g+cstar/ystar*c_f+istar/ystar*i_f+rkstar*kstar/ystar*u_f-gstar*(1/(1-alp))*ztil;

w_f = (1/(1-h*exp(-zstar)))*(c_f-h*exp(-zstar)*c_f(-1)+h*exp(-zstar)*z)+nu_l*L_f;

//assume no financial frictions in the flex-price economy (i.e. arbitrage condition)

qk_f = rkstar/(rkstar+1-del)*rk_f(+1)+(1-del)/(rkstar+1-del)*qk_f(+1)-r_f+(sigmac*(1+h*exp(-zstar)))/(1-h*exp(-zstar))*b;

end;

//--------------------------------------------------------------------------------------------------------------------------
// 4. Steady state
//---------------------------------------------------------------------------------------------------------------------------

steady;

//---------------------------------------------------------------------------------------------------------------------------
// 5. shocks
//---------------------------------------------------------------------------------------------------------------------------

shocks;


var psi_b; stderr 0.0265;
var psi_mu; stderr 0.42853;
var psi_z; stderr 0.4651;
var psi_laf; stderr 0.1853;
var psi_law; stderr 0.2749;
var psi_sigw; stderr 0.0529;
var psi_pist; stderr 0.030;
var psi_rm; stderr 0.2919;
end;

//write_latex_dynamic_model;
//check;
//stoch_simul(irf=20) R pi L y;
stoch_simul (AR=100,IRF=0, noprint,nograph);