//**********************************************************************
//The Implications of Financial Frictions and Imperfect Knowledge in
//the Estimated Model of the US Economy
//
//Yuliya Rychalovska
//
//CERGE-EI Working Paper No. 482
//
//Adaptive Learning version
//Last edited: 13/12/2013 by Sergey Slobodyan
//**********************************************************************
//Model: US_YR13

var ewma epinfma  mc zcap rk k pk
    c inve y lab pinf w r a  b g qs spinf sw kp
    rr er ec nw prem pinf4 eg ms;
	//cl invel labl pinfl pkl rkl wl   ms


varexo  ea eb eqs  epinf ew
em;


parameters


            curvw cgy curvp constelab constepinf constebeta   calfa
            czcap cbeta csadjcost ctou csigma chabb ccs cinvs cfc
            cindw cprobw cindp cprobp csigl clandaw
            crdpi crpi crdy cry crr elast cv clev cff tau
            crhoa crhoas crhob crhog crhols crhoqs crhoms crhopinf crhow cmaw cmap
            ctrendy ctrendc ctrendinv ctrend ctrendw conster cg cgamma clandap cbetabar cr cpie crk cw cikbar cik clk cky ciy ccy crkky cwhlc cwly ro;

// fixed parameters

ctou=       0.025;
clandaw=    1.5;
cg=         0.5787;
curvp=     10.0;
curvw=     10.0;

// estimated parameters

ctrendy =    0.4416;
ctrendc =    0.5026;
ctrendinv =    0.3696;
ctrendw =    0.3974;

constebeta= 0.1443;
constepinf= 0.6605;
constelab=  1.0589;
calfa=      0.1994;
cgy=        0.5787;
csadjcost=  5.9774;
csigma=     1.4874;
chabb=      0.7194;
cprobw=     0.7329;
csigl=      1.5992;
cprobp=     0.5756;
cindw=      0.2020;
cindp=      0.4184;
czcap=      0.1356;
cfc=        1.5886;
crpi=       1.6692;
crr=        0.9209;
cry=        0.1372;
crdy=       0.1004;

crhoa=      0.9828;
crhob=      0.3881;
crhog=      0.9249;
crhols=     0.9928;
crhoqs=     0.4265;
crhopinf=   0.3206;
crhow=      0.5702;
cmap =      0.5020;
cmaw  =     0.4520;

tau=0.025;
cff=1.0129;
clev=3.4525;
cv=0.9447;
elast=0.0197;
ro = 0.9881;


// derived parameters from steady state : see stst_f19.m
ctrend=(ctrendy+ctrendc+ctrendinv+ctrendw)/4;
cpie=     1+constepinf/100;
cgamma=   1+ctrend/100 ;
cbeta=    1/(1+constebeta/100);
conster=  (cr-1)*100;
clandap=  cfc;
cbetabar= cbeta*cgamma^(-csigma);
cr=       cpie/(cbeta*cgamma^(-csigma));
crk=      (cbeta^(-1))*(cgamma^csigma) - (1-ctou);
cw =      (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));
cikbar=   (1-(1-ctou)/cgamma);
cik=      (1-(1-ctou)/cgamma)*cgamma;
clk=      ((1-calfa)/calfa)*(crk/cw);
cky=      cfc*(clk)^(calfa-1);
ciy=      cik*cky;
ccy=      1-cg-cik*cky;
crkky=    crk*cky;
cwhlc=    (1/clandaw)*(1-calfa)/calfa*crk*cky/ccy;
cwly=     1-crk*cky;
cff=crk+1-ctou;


model(linear);

//Annual inflation:
            pinf4= 0.25 * (4*pinf + 4*pinf(-1) + 4*pinf(-2) + 4*pinf(-3));
// Original Model Code:
//stst_f19_traver;


            mc =  calfa*rk+(1-calfa)*(w) - 1*a - 0*(1-calfa)*a ;
            zcap =  (1/(czcap/(1-czcap)))* rk ;
            rk =  w+lab-k ;
            k =  kp(-1)+zcap ;
            inve = (1/(1+cbetabar*cgamma))* (  inve(-1) + cbetabar*cgamma*inve(1)+(1/(cgamma^2*csadjcost))*pk ) +qs ;
            pk = -r+pinf(1)-b +0*(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b + ((cff-1+ctou)/cff)*rk(1) +  ((1-ctou)/cff)*pk(1) ;
            c = (chabb/cgamma)/(1+chabb/cgamma)*c(-1) + (1/(1+chabb/cgamma))*c(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(lab-lab(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf(+1) + b) +0*b ;
            y = ccy*c+ciy*inve+g  +  1*((cff-1+ctou)*(ciy/cik))*zcap + (ciy/cik)*cff*(1-(cr/cff))*(1-1/clev)*(rr + pk(-1) +k);


            y = cfc*( calfa*k+(1-calfa)*lab +a );
            pinf =  (1/(1+cbetabar*cgamma*cindp)) * ( cbetabar*cgamma*pinf(1) +cindp*pinf(-1)
               +((1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1)*(mc)  )  + spinf ;
            w =  (1/(1+cbetabar*cgamma))*w(-1)
               +(cbetabar*cgamma/(1+cbetabar*cgamma))*w(1)
               +(cindw/(1+cbetabar*cgamma))*pinf(-1)
               -(1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*pinf
               +(cbetabar*cgamma)/(1+cbetabar*cgamma)*pinf(1)
               +(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1))*
                   (csigl*lab + (1/(1-chabb/cgamma))*c - ((chabb/cgamma)/(1-chabb/cgamma))*c(-1) -w)
               +1*sw ;
            r = crpi*(1-crr)*pinf +cry*(1-crr)*(y-cfc*a) +crdy*(y-y(-1)-cfc*(a-a(-1))) +crr*r(-1) +ms  ;
            kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs ;

            a     = crhoa*a(-1)   + ea;
            b     = crhob*b(-1)   + eb;
            g     = crhog*(g(-1)) + eg + cgy*ea;
            qs    = crhoqs*qs(-1)+ eqs;
            ms    =  em;
            spinf = crhopinf*spinf(-1) + epinfma - cmap*epinfma(-1);
	          epinfma=epinf;
	      sw = crhow*sw(-1) + ewma - cmaw*ewma(-1) ;
	          ewma=ew;

rr + pk(-1)= ((cff-1+tau)/cff)*rk + (1-tau)/cff*pk ;
nw/(cv*cff) = clev*(rr-ec(-1))+ ec(-1) + nw(-1);
er + pk = ((cff-1+tau)/cff)*rk(1) + (1-tau)/cff*pk(1)  ; //+ ETA_Q;
ec-(r-pinf(1)+b)=-elast*(nw -pk -k);
ec=er;
prem=ec-(r-pinf(1)+b);

            //cl = c(-1);
            //invel = inve(-1);
            //labl = lab(-1);
            //pinfl = pinf(-1);
            //pkl = pk(-1);
            //rkl = rk(-1);
            //wl = w(-1);

           end;

shocks;
var ea;
stderr 0.4853;
var eb;
stderr 1.7918;

var eqs;
stderr 0.3959;
var em; stderr 0.2005;
var epinf;
stderr 0.1852;
var ew;
stderr 0.2277;
end;

stoch_simul (AR=100,IRF=0, noprint,nograph);
//stoch_simul(irf = 0, ar=100, periods=10000);

