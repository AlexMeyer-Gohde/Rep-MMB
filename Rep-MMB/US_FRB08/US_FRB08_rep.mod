//**************************************************************************
// A New Comparative Approach to Macroeconomic Modeling and Policy Analysis
//
// Volker Wieland, Tobias Cwik, Gernot J. Mueller, Sebastian Schmidt and 
// Maik Wolters
//
// Working Paper, 2009
//**************************************************************************

// Model: US_FRB08

// Last edited: 09/07/23 by M. Wolters with M. Burgert

//This model file contains the linearized Version of the FRB-US Model. 
//The model equations and parameters are of February 2008 and kindly provided by Thomas Laubach.
//This version of the model assumes rational expectations.


var xgap2 lur picnia picxfe pigdp drffe ec eco qec ecd qecd echsh ech eh qeh 
    ecnia kcd kh epdc qepdc epdo qepdo epd eps qeps ki kpdc kpdo kps hks ks 
    rpd rtpdc rtpdo rtps tapdda vpdc vpdo vps ex emo ceng emp em fcbnr fnicnr 
    fnilnr fynicnr fynilnr rfynic rfynil xgdp xb xg xp efs hmfpt mfpt veoa 
    xgpot hgx empt hgemp xgvt hgxgv xhhit hgxhhi xgap xgdpt hggdpt lhp qlhp 
    lww lep leo leh lfpr left lest leppot lurnat lprdt jccanr qynidnr rrmet 
    tryh wdnfcnr xbnr xgdpnr xgnr ycsnr yh yhibnr yhinr yhl yhp yhpntnr yht 
    ynicpnr ynidnr yniish yniinr ynilnr yninr ypnr pcniar pxncr dlpl plr qpcniar 
    qpxncr qplr qpxpr qpxgr ppdcr poilr pcengr pcer pcfr uces ucfs pcxfer pmor 
    pgdpr plminr pxbr pxgr hgpmor pimor dpxnc dlpxp egfi egfl egfo egsi egsl 
    egso egt eg gfdbtnr gfintnr gfrecnr gfexpnr gfdefnr gftrd gsdbtnr gsintnr 
    gsrecnr gsexpnr gssrpnr gstrd rgfint tfpnr trfci trfp trfpt trsci trsib 
    trsp trspt trssi tspnr rfftay rffe rtbe rg5e rg5p rg10e rg10p rbaae rbaap 
    rcbe rcbp rcar rme req reqp wpsnr rcgain wpo fxgap fgdp fgdpt fpi10 fpi10t  
    fpic frs10 frl10 fpxr fpxrr rstar zeco zecd zgapc2 zeh zvpdc zvpdo zvps 
    zxbc zxbo zxbs zpib5 zpi5 zpi10 zpic30 zlhp zyh zyhp zyht zynid zlurc zlurl 
    zlurnc zpc zpl zpnc zdivgr zrff5 zrff10 zrff30 zgap05 zgap10 zgap30 ecol1 
    ecdl1 echshl1 ehl1 ehl2 qehf1 epdcl1 epdcl2 epdol1 epdol2 qepdol1 epsl1 
    epsl2 qepsl1 kil1 kil2 kil3 kil4 vpdcf1 vpdof1 vpsf1 cengl1 cengl2 rfyncl1 
    rfynll1 xbf1 xbl1 xgl1 efsl1 efsl2 efsl3 efsl4 hmfptl1 hgxl1 xgdptl1 xgdptl2 
    hggdpl1 hggdpl2 xgap2l1 lhpl1  qlhpl1 lurl1 leppol1 leppol2 lurnal1 wdnfcl1 
    wdnfcl2 xgdpnl1 yhibnl1 ynidnl1 ynilnl1 yninrl1 picnil1 picnil2 picnil3 
    dpxncl1 dpxncl2 dlpll1 dlpll2 poilrl1 pcengl1 pcerl1 pcfrl1 pcfrl2 pcfrl3 
    ucesl1 ucfsl1 picxfl1 picxfl2 picxfl3 pxbrl1 dlpxpl1 dlpxpl2 dlpxpl3 egfil1 
    egfil2 egfll1 egfll2 egfol1 egfol2 egsil1 egsil2 egsll1 egsll2 egsol1 egsol2 
    gfdbtl1 gsdbtl1 rgfinl1 trfpl1 trfptl1 trssil1 rtbel1  rg5el1  rmel1 rmel2 
    fxgapl1 fgdpl1 fgdpl2 fgdpl3 fgdpl4 fgdpl5 fpi10l1 fpi10l2 fpi10l3 frl10l1 
    fpxrrl1 lzeco lzecd lzgapc2 lzeh zehf2 lzvpdc zvpdcf2 lzvpdo zvpdof2 lzvpol1 
    lzvps zvpsf2 lzvpsl1 lzxbc  zxbcf2 lzxbo zxbof2 lzxbol1 lzxbs zxbsf2 lzxbsl1 
    lzpib5 lzlhp lzynid lzlurc lzlurl lzlurnc lzpc lzpl lzpnc one;



varexo 
ceng_ 
dlpl_ 
ecd_ 
ech_ 
eco_ 
egsi_ 
egsl_ 
egso_ 
eh_ 
emo_ 
epdc_ 
epdo_ 
eps_ 
ex_ 
%fiscal_ 
fpi10_      
fpxrr_ 
frl10_ 
fxgap_ 
gftrd_ 
gstrd_ 
hmfpt_ 
%interest_ 
ki_ 
lhp_ 
mfpt_  
pcengr_ 
pcer_ 
pcfr_ 
pcniar_ 
pmor_ 
poilr_ 
pxncr_ 
rbaap_ 
rcar_ 
rcbp_ 
rcgain_
reqp_ 
rfynic_ 
rfynil_ 
rg10p_ 
rg5p_  
rgfint_ 
rme_ 
rtbe_ 
trfci_ 
trfp_ 
trsci_ 
trsib_ 
trsp_ 
trssi_  
uces_ 
ucfs_ 
wdnfcnr_ 
yhibnr_ 
ynidnr_    
yniish_
egfi_ egfl_;            


parameters 

           fsector datet taxon qsector gsector dglprd rstar1  trfpt1  trfpt2  
           trspt1  trspt2  tayr1   tayp0   tayp1  tayp2 tayp3 tayx0 tayx1  
           zeroadj rhozero lagadj mca mco lurb1 lurb2 lurb3 lurb4 ecs1 ecs2 
           ecos1 ecoa1 ecoa2 ecod1 qecs1 qecs2 qecs3 qecs4 qecs5 qecb1 qecb2
           qecb3 qecb4 qecb5 ecda1 ecda2 ecdb2 qecds1 qecds2 qecda3 echs1 echa1 
           echa2 echa3 eha1 eha2 eha3 ehb2 ehb3 qehs1 qehs2 qeha2 ecnias1 ecnias2
           ecnias3 kcd1 kh1 epdcb1 epdcb2a epdcb2b epdcb3 qepdcs1 epdob1 epdob2a 
           epdob2b epdob3 qepdos1 epds1 epsb1 epsb2a epsb2b epsb3 qepss1 kib1 
           kib2 kib3 kib4 kpdcs1 kpdos1 kpss1 hkss1 hkss2 hkss3 rtpdcs1 rtpdcs2 
           rtpdos1 rtpdos2 rtpss1 tapdda1 exa1 exd0 exd1 emoa1 emod0 emod1 emod2 
           cenga1 cengb1a cengb1b cengb2 cengb3a cengb3b emps1 emps2 ems1 ems2 
           fcbns1 fcbns2 fcbns2a fcbns2b fcbns3 fcbns4 fcbns5 fnicns1 fnicns2 
           fnicns3 fnilns1 fnilns2 fnilns3 fynicns1 fynilns1 rfynica1 rfynica2 
           rfynica3 rfynila1 rfynila2 rfynila3 rfynila4 rfynila5 rfynila6 rfynila7 
           rfynila8 xgdps1 xgdps2 xgdps3 xgdps4 xgdps5 xgdps6 xgdps7 xgdps8 
           xgdps9 xgdps10 egfs1 egfs2 egfs3 egss1 egss2 egss3 xbs1 xbs2 xbs3 
           xgvs1 xgs1 xps1 xps2 xps3 xps4 xps5 xps6 xps7 xps8 xps9 efss1 efss2 
           efss3 hmfptrho veoaa1 empta1 hgempa1 xgvta1 hgxgva1 xhhita1 hgxhhia1 
           xgdpts1 xgdpts2 xgdpts3 xgdpts4 xgdpts5 lhpb2 lhpb3 lhpb4 lhpb5 lhpb6 
           lwwa2 lwwa3 leors1 leors2 leors3 leors4 lehs1 lehs2 lfpra1 lfpra2 
           lefta1 lesta1 leppots1 leppots2 leppots3 jccans1 jccans2 jccans3 
           jccans4 jccans5 qynidb2 qynids1 rrmeta1 tryhs1 tryhs2 tryhs3 tryhs4 
           yhptns1 yhptns2 yhptns3 wdnfca1 wdnfcb1 wdnfcb2 wdnfcb3 wdnfcs1 
           wdnfcs2 xbns1 xbns2 xbns3 xgdpns1 xgdpns2 xgdpns3 xgdpns4 xgdpns5 
           xgdpns6 xgns1 ycsns1 ycsns2 ycsns3 ycsns4 jccacs1 jccacs2 jccacs3 
           yhs1 yhs2 yhs3 yhibna1 yhibna3 yhibna4 yhibna6 yhibna7 yhibna8 yhins1 
           yhins2 yhins3 yhins4 yhls1 yhls2 yhps1 yhps2 yhpntns1 yhpntns2 yhpntns3 
           yhpntns4 yhpntns5 yhpntns6 yhpntns7 yhpntns8 yhts1 ynicpns1 ynicpns2 
           ynicpns3 ynicpns4 ynicpns5 ynicpns6 ynicpns7 ynidna1 ynidna2 yniina1 
           yniina2 yniina3 yniina4 yniina5 yniins1 yniins2 yniins3 yniins4 yniins5 
           yniins6 yniins7 yniins8 ynilns1 ynilns2 ynilns3 ynins1 ynins2 ynins3 
           ynins4 ypns1 ypns2 ypns3 pcniaa1 pcniaa21 pcniaa22 pcniaa23 pcniaa24 
           pcniaa4 pcniaa5 pcniaa6 pcniaa8 pcnias1 pcnias2 pxnca1 pxnca21 pxnca22 
           pxnca23 pxnca4 pxnca5 pxnca6 pxnca8 pipla1 pipla21 pipla22 pipla23 
           pipla4 pipla6 pipla7 dlpls1 dlpls2 tsiens1 qpcnias1 uqpca2 qpxnca1 
           qpxnca2 qpxps1 poilra1 poilrb1 pcengrb1 pcengrb2 pcengrb3 pcengrb4 
           pcera1 pcera2 pcerb1 pcerb2 pcerb3 pcfra1 pcfrb11 pcfrb12 pcfrb13 
           ucesa1 ucesa2 ucesa3 ucesb1 ucesb2 ucesb3 ucfsa1 ucfsa2 ucfsb1 ucfsb3 
           pmoa1 pmoa2 pmoa3 hgpmora1 pimors1 egfia1 egfia2a egfia2b egfia3 
           egfia4a egfia4b egfla1 egfla2a egfla2b egfla3 egfla4a egfla4b egfoa1 
           egfoa2a egfoa2b egfoa3 egfoa4a egfoa4b egsia1 egsia2a egsia2b egsia3 
           egsia4a egsia4b egsla1 egsla2a egsla2b egsla3 egsla4a egsla4b egsoa1 
           egsoa2a egsoa2b egsoa3 egsoa4a egsoa4b egta1 egs1 egs2 gfdbtns2 
           gfdbtns3 gfdbtns4 gfdbtns5 gfdbtns1 gfintns1 gfrecns1 gfrecns2 gfrecns3 
           gfrecns4 trfcis1 gfexpns1 gfexpns2 gfexpns3 gfexpns4 gfexpns5 gfexpns6 
           gfptns1 gfdefns1 gftrda1 gftrda2 gsdbtns2 gsdbtns3 gsdbtns4 gsdbtns5 
           gsdbtns1 gsintns1 gsrecns1 gsrecns2 gsrecns3 gsrecns4 gsrecns5 trscis1 
           trsibs1 trssis1 gsexpns2 gsexpns3 gsexpns4 gsexpns5 gsexpns1 gsptns1 
           gssrpns1 gstrda1 gstrda2 rgfinta1 rgfinta2 rgfintb1 rgfintb2 rgfintb3 
           tfpns1 tfpns2 tfpns3 tfpns4 trfcia1 trfcia4 trfcia5 trfpa21 trfpa22 
           trfpa3 trscia1 trscia31 trscia32 trscia4 trsiba1 trsiba3 trspa1 
           trspa21 trspa22 trspa3 trspa4 trssia11 trssia12 trssia3 tspns1 tspns2 
           tspns3 tspns4 rtbeb11 rtbeb12 rtbeb21 rtbeb22 rg5pa1 rg5pa2 rg10pa1 
           rg10pa2 rbaapa1 rbaapa2 rcbpa1 rcbpa2 rcara3 rcara41 rcara42 rmea1 
           rmea21 rmea22 rmea3 reqpb1 wpsns1 rcgaina1 rcgaina2 rcgaina5 wpos2 
           wpos3 wpos4 wpos5 wpos6 wpos7 wpos1 ydns1 ydns2 ydns3 fxgapb1a fxgapb1b  
           fxgapb2 fxgapb3 fxgapb4 fgdpta1 fpi10b2 fpi10b3a fpi10b3b fpi10tb1 
           fpi10tb2 fpicb1 fpicb2 frs10b1 frs10b2 frs10b3 frs10b4 frl10b1 frl10b2 
           frl10b3 frl10b4 fpxrb1 fpxrb2 fpxrs1 fpxrs2 fprho fpxrrb1 fpxrrb2 
           zfecoa1 zfecoa2 zfecoa6 zfecoa7 zbeco1 zbeco2 zbeco3 zbeco4 zbeco5 
           zbeco6 zbeco7 zbeco8 zbeco9 zbeco10 zbeco11 zbeco12 zbeco13 zbeco14 
           zbeco15 zbeco16 zbeco17 zbeco18 zbeco19 zbeco20 zbeco21 zbeco22 
           zbeco23 zbeco24 zbeco25 zbeco26 zbeco27 zbeco28 zbeco29 zbeco30 
           zbeco31 zfecda1 zfecda2 zfecda6 zfecda7 zbecda1 zbecda2 zbecda3 
           zbecda4 zbecda5 zbecda6 zbecda7 zbecda8 zbecda9 zbecda10 zbecda11 
           zbecda12 zbecda13 zbecda14 zbecda15 zbecda16 zbecda17 zbecda18 
           zbecda19 zbecda20 zbecda21 zbecda22 zbecda23 zbecda24 zbecda25 
           zbecda26 zbecda27 zbecda28 zbecda29 zbecda30 zbecda31 zfgapc2a1 
           zfgapc2a6 zfgapc2a7 zbgapc2a1 zbgapc2a2 zbgapc2a3 zbgapc2a4 zbgapc2a5 
           zbgapc2a6 zbgapc2a7 zbgapc2a8 zbgapc2a9 zbgapc2a10 zbgapc2a11 zbgapc2a12 
           zbgapc2a13 zbgapc2a14 zfeha1 zfeha2 zfeha3 zfeha6 zfeha7 zfeha8 
           zbeha1 zbeha2 zbeha3 zbeha4 zbeha5 zbeha6 zbeha7 zbeha8 zbeha9 
           zbeha10 zbeha11 zbeha12 zbeha13 zbeha14 zbeha15 zbeha16 zbeha17 
           zbeha18 zbeha19 zbeha20 zbeha21 zbeha22 zbeha23 zbeha24 zbeha25 
           zbeha26 zbeha27 zbeha28 zbeha29 zbeha30 zbeha31 zfvpdc1 zfvpdc2 
           zfvpdc3 zfvpdc6 zfvpdc7 zfvpdc8 zbvpdca1 zbvpdca2 zbvpdca3 zbvpdca4 
           zbvpdca5 zbvpdca6 zbvpdca7 zbvpdca8 zbvpdca9 zbvpdca10 zbvpdca11 
           zbvpdca12 zbvpdca13 zbvpdca14 zbvpdca15 zbvpdca16 zbvpdca17 zbvpdca18 
           zbvpdca19 zbvpdca20 zbvpdca21 zbvpdca22 zfvpdo1 zfvpdo2 zfvpdo3 
           zfvpdo6 zfvpdo7 zfvpdo8 zbvpdoa1 zbvpdoa2 zbvpdoa3 zbvpdoa4 zbvpdoa5 
           zbvpdoa6 zbvpdoa7 zbvpdoa8 zbvpdoa9 zbvpdoa10 zbvpdoa11 zbvpdoa12 
           zbvpdoa13 zbvpdoa14 zbvpdoa15 zbvpdoa16 zbvpdoa17 zbvpdoa18 zbvpdoa19 
           zbvpdoa20 zbvpdoa21 zbvpdoa22 zfvps1 zfvps2 zfvps3 zfvps6 zfvps7 
           zfvps8 zbvpsa1 zbvpsa2 zbvpsa3 zbvpsa4 zbvpsa5 zbvpsa6 zbvpsa7 
           zbvpsa8 zbvpsa9 zbvpsa10 zbvpsa11 zbvpsa12 zbvpsa13 zbvpsa14 zbvpsa15 
           zbvpsa16 zbvpsa17 zbvpsa18 zbvpsa19 zbvpsa20 zbvpsa21 zbvpsa22 zfxbc1 
           zfxbc2 zfxbc3 zfxbc6 zfxbc7 zfxbc8 zbxbca1 zbxbca2 zbxbca3 zbxbca4 
           zbxbca5 zbxbca6 zbxbca7 zbxbca8 zbxbca9 zbxbca10 zbxbca11 zbxbca12 
           zbxbca13 zbxbca14 zbxbca15 zbxbca16 zbxbca17 zbxbca18 zbxbca19 zbxbca20 
           zbxbca21 zbxbca22 zbxbca23 zfxbo1 zfxbo2 zfxbo3 zfxbo6 zfxbo7 zfxbo8 
           zbxboa1 zbxboa2 zbxboa3 zbxboa4 zbxboa5 zbxboa6 zbxboa7 zbxboa8 
           zbxboa9 zbxboa10 zbxboa11 zbxboa12 zbxboa13 zbxboa14 zbxboa15 zbxboa16 
           zbxboa17 zbxboa18 zbxboa19 zfxbs1 zfxbs2 zfxbs3 zfxbs6 zfxbs7 zfxbs8 
           zbxbsa1 zbxbsa2 zbxbsa3 zbxbsa4 zbxbsa5 zbxbsa6 zbxbsa7 zbxbsa8 zbxbsa9 
           zbxbsa10 zbxbsa11 zbxbsa12 zbxbsa13 zbxbsa14 zbxbsa15 zbxbsa16 zbxbsa17 
           zbxbsa18 zbxbsa19 zfpib5a1 zfpib5a6 zbpib5a1 zbpib5a2 zbpib5a3 zbpib5a4 
           zbpib5a5 zbpib5a6 zbpib5a7 zbpib5a8 zbpib5a9 zbpib5a10 zbpib5a11
           zbpib5a12 zbpib5a13 zbpib5a14 zbpib5a15 zbpib5a16 zbpib5a17 zbpib5a18 
           zfpi5a1 zfpi5a6 zbpi5a1 zbpi5a2 zbpi5a3 zbpi5a4 zbpi5a5 zbpi5a6 
           zbpi5a7 zbpi5a8 zbpi5a9 zbpi5a10 zbpi5a11 zbpi5a12 zbpi5a13 zbpi5a14 
           zfpi10a1 zfpi10a6 zbpi10a1 zbpi10a2 zbpi10a3 zbpi10a4 zbpi10a5 zbpi10a6 
           zbpi10a7 zbpi10a8 zbpi10a9 zbpi10a10 zbpi10a11 zbpi10a12 zbpi10a13 
           zbpi10a14 zfpic30a1 zfpic30a6 zbpic30a1 zbpic30a2 zbpic30a3 zbpic30a4 
           zbpic30a5 zbpic30a6 zbpic30a7 zbpic30a8 zbpic30a9 zbpic30a10 zbpic30a11 
           zbpic30a12 zbpic30a13 zbpic30a14 zflhpa1 zflhpa2 zflhpa6 zflhpa7 
           zblhpa1 zblhpa2 zblhpa3 zblhpa4 zblhpa5 zblhpa6 zblhpa7 zblhpa8 
           zblhpa9 zblhpa10 zblhpa11 zblhpa12 zblhpa13 zblhpa14 zblhpa15 zblhpa16 
           zfynida1 zfynida2 zfynida6 zfynida7 zbynida1 zbynida2 zbynida3 zbynida4 
           zbynida5 zbynida6 zbynida7 zbynida8 zbynida9 zbynida10 zbynida11 
           zbynida12 zbynida13 zbynida14 zbynida15 zbynida16 zbynida17 zbynida18 
           zbynida19 zflurca1 zflurca6 zblurca1 zblurca2 zblurca3 zblurca4 
           zblurca5 zblurca6 zblurca7 zblurca8 zblurca9 zblurca10 zblurca11 
           zblurca12 zblurca13 zblurca14 zblurca15 zblurca16 zblurca17 zblurca18 
           zblurca19 zblurca20 zblurca21 zblurca22 zblurca23 zblurca24 zblurca25 
           zblurca26 zblurca27 zblurca28 zblurca29 zblurca30 zblurca31 zblurca32 
           zblurca33 zflurla1 zflurla6 zblurla1 zblurla2 zblurla3 zblurla4 
           zblurla5 zblurla6 zblurla7 zblurla8 zblurla9 zblurla10 zblurla11 
           zblurla12 zblurla13 zblurla14 zblurla15 zblurla16 zblurla17 zblurla18 
           zblurla19 zblurla20 zblurla21 zblurla22 zblurla23 zblurla24 zblurla25 
           zblurla26 zblurla27 zblurla28 zblurla29 zblurla30 zblurla31 zblurla32 
           zblurla33 zflurnca1 zflurnca6 zblurnca1 zblurnca2 zblurnca3 zblurnca4 
           zblurnca5 zblurnca6 zblurnca7 zblurnca8 zblurnca9 zblurnca10 zblurnca11 
           zblurnca12 zblurnca13 zblurnca14 zblurnca15 zblurnca16 zblurnca17 
           zblurnca18 zblurnca19 zblurnca20 zblurnca21 zblurnca22 zblurnca23 
           zblurnca24 zblurnca25 zblurnca26 zblurnca27 zblurnca28 zblurnca29 
           zblurnca30 zblurnca31 zblurnca32 zblurnca33 zfpca1 zfpca6 zbpca1 
           zbpca2 zbpca3 zbpca4 zbpca5 zbpca6 zbpca7 zbpca8 zbpca9 zbpca10 
           zbpca11 zbpca12 zbpca13 zbpca14 zbpca15 zbpca16 zbpca17 zbpca18 
           zbpca19 zbpca20 zbpca21 zbpca22 zbpca23 zbpca24 zbpca25 zbpca26
           zbpca27 zbpca28 zbpca29 zbpca30 zbpca31 zbpca32 zbpca33 zfpla1 
           zfpla6 zbpla1 zbpla2 zbpla3 zbpla4 zbpla5 zbpla6 zbpla7 zbpla8 
           zbpla9 zbpla10 zbpla11 zbpla12 zbpla13 zbpla14 zbpla15 zbpla16 
           zbpla17 zbpla18 zbpla19 zbpla20 zbpla21 zbpla22 zbpla23 zbpla24 
           zbpla25 zbpla26 zbpla27 zbpla28 zbpla29 zbpla30 zbpla31 zbpla32
           zbpla33 zfpnca1 zfpnca6 zbpnca1 zbpnca2 zbpnca3 zbpnca4 zbpnca5 
           zbpnca6 zbpnca7 zbpnca8 zbpnca9 zbpnca10 zbpnca11 zbpnca12 zbpnca13 
           zbpnca14 zbpnca15 zbpnca16 zbpnca17 zbpnca18 zbpnca19 zbpnca20 zbpnca21 
           zbpnca22 zbpnca23 zbpnca24 zbpnca25 zbpnca26 zbpnca27 zbpnca28 zbpnca29 
           bpnca30 zbpnca31 zbpnca32 zbpnca33 zfdivgra1 zfdivgra6 zbdivgra1 
           zbdivgra2 zbdivgra3 zbdivgra4 zbdivgra5 zbdivgra6 zbdivgra7 zbdivgra8 
           zbdivgra9 zbdivgra10 zbdivgra11 zbdivgra12 zbdivgra13 zbdivgra14 
           zbdivgra15 zbdivgra16 zbdivgra17 zbdivgra18 zbdivgra19 zfrff5a1 
           zfrff5a6 zbrff5a1 zbrff5a2 zbrff5a3 zbrff5a4 zbrff5a5 zbrff5a6 
           zbrff5a7 zbrff5a8 zbrff5a9 zbrff5a10 zbrff5a11 zbrff5a12 zbrff5a13 
           zbrff5a14 zfrff10a1 zfrff10a6 zbrff10a1 zbrff10a2 zbrff10a3 zbrff10a4 
           zbrff10a5 zbrff10a6 zbrff10a7 zbrff10a8 zbrff10a9 zbrff10a10 zbrff10a11 
           zbrff10a12 zbrff10a13 zbrff10a14 zfrff30a1 zfrff30a6 zbrff30a1 
           zbrff30a2 zbrff30a3 zbrff30a4 zbrff30a5 zbrff30a6 zbrff30a7 zbrff30a8 
           zbrff30a9 zbrff30a10 zbrff30a11 zbrff30a12 zbrff30a13 zbrff30a14 
           zfgap05a1 zfgap05a6 zbgap05a1 zbgap05a2 zbgap05a3 zbgap05a4 zbgap05a5 
           zbgap05a6 zbgap05a7 zbgap05a8 zbgap05a9 zbgap05a10 zbgap05a11 zbgap05a12 
           zbgap05a13 zbgap05a14 zfgap10a1 zfgap10a6 zbgap10a1 zbgap10a2 zbgap10a3
           zbgap10a4 zbgap10a5 zbgap10a6 zbgap10a7 zbgap10a8 zbgap10a9 zbgap10a10 
           zbgap10a11 zbgap10a12 zbgap10a13 zbgap10a14 zfgap30a1 zfgap30a6 
           zbgap30a1 zbgap30a2 zbgap30a3 zbgap30a4 zbgap30a5 zbgap30a6 zbgap30a7 
           zbgap30a8 zbgap30a9 zbgap30a10 zbgap30a11 zbgap30a12 zbgap30a13 
           zbgap30a14 zfyha1 zfyha2 zbyha1 zbyha2 zbyha3 zbyha4 zbyha5 zbyha6 
           zbyha7 zbyha8 zbyha9 zbyha10 zbyha11 zbyha12 zbyha13 zbyha14 zbyha15 
           zbyha16 zbyha17 zbyha18 zfyhta1 zfyhta2 zbyhta1 zbyhta2 zbyhta3 
           zbyhta4 zbyhta5 zbyhta6 zbyhta7 zbyhta8 zbyhta9 zbyhta10 zbyhta11 
           zbyhta12 zbyhta13 zbyhta14 zbyhta15 zbyhta16 zbyhta17 zbyhta18 zbyhta19 
           zbyhta20 zbyhta21 zbyhta22 zfyhpa1 zfyhpa2 zbyhpa1 zbyhpa2 zbyhpa3 
           zbyhpa4 zbyhpa5 zbyhpa6 zbyhpa7 zbyhpa8 zbyhpa9 zbyhpa10 zbyhpa11 
           zbyhpa12 zbyhpa13 zbyhpa14 zbyhpa15 zbyhpa16 zbyhpa17 zbyhpa18 
           zbyhpa19 zbyhpa20 zbyhpa21 zbyhpa22 hprdtwa1 zyhsta1 zyhtsta1 zyhpsta1 
           yhshrs1 yhtshrs1 yhpshrs1 ptra1 ptra2 ptra3 ptra4 ptra5 ptra6 ptra7 
           ptra8 ptra9 ptra10 ptra11 ptra12 ptra13 ptra14 rrtra1;


fsector =  1;
datet   =  0; // datet=1: Model-Consistent solution based on time t expectations;
              // datet=0: Same timing as in VAR expectations
taxon   =  1; // taxon=1: Tax non-neutralities in place; taxon=0: non-neutralities removed (eh,epd)
qsector =  1; // Switch for Stochastic Simulations: Ad hoc forecasting equations (usually exog.)
gsector =  1; // Switch for Stochastic Simualtions: Ad hoc government forecasting equations (usually exog.)
dglprd = 1; // Switch for treatment of productivity growth in government sector
rstar1  =  0.95; //

trfpt1  =  0.05; // tax rate reaction functions to debt/GDP
trfpt2  =  0.50; // tax rate reaction functions to debt/GDP
trspt1  =  0.05; // tax rate reaction functions to debt/GDP
trspt2  =  0.50; // tax rate reaction functions to debt/GDP

tayr1   =   0;
tayp0   =       0.125;
tayp1   =       tayp0;
tayp2   =       tayp0;
tayp3   =       tayp0;

tayx0   =       1.0;
tayx1   =       0;

zeroadj = 0; //.25
rhozero = 0;
lagadj=0;

mca=0;
mco=0;

// LUR
lurb1 =     0.87997794000000000;
lurb2 =     0.31117806000000000;
lurb3 =    -0.06414897700000000;
lurb4 =    -0.15320231000000000;

// EC
ecs1 =     0.72472800000000000;
ecs2 =     0.15066900000000000;

// ECO
ecos1 =     0.79700700000000000;
ecoa1 =     0.09658366291970000;
ecoa2 =     0.16567146260700000;
ecod1 =     0.54061698114000000;

// QEC
qecs1 =     1.07403100000000000;
qecs2 =     0.15977500000000000;
qecs3 =     0.29365800000000000;
qecs4 =     5.77367900000000000;
qecs5 =     0.32485900000000000;
qecb1 =     0.82223284447400000;
qecb2 =     0.02744726492880000;
qecb3 =     0.17776715552600000;
qecb4 =    -0.52831383004700010;
qecb5 =     0.03785160449050000;

// ECD
ecda1 =     0.17934374348000000;
ecda2 =    -0.14488268355800000;
ecdb2 =     1.96861385021000000;

// QECD
qecds1 =    15.53093800000000000;
qecds2 =     0.03651200000000000;
qecda3 =    -0.65797457059600000;

// ECH
echs1 =    10.58128800000000000;

// ECHSH
echa1 =    -0.04189322518140000;
echa2 =     0.60465997056300000;
echa3 =     0.00251399113239000;

// EH
eha1 =     0.09320978458430000;
eha2 =     0.25327096841500000;
eha3 =     0.23275999948600000;
ehb2 =    -0.01858143360210000;
ehb3 =    -0.01025202488966000;

// QEH
qehs1 =    86.88568700000000000;
qehs2 =     0.18190400000000000;
qeha2 =    -0.27575823838800000;

// ECNIA
ecnias1 =     0.72430200000000000;
ecnias2 =     0.12455800000000000;
ecnias3 =     0.15114000000000000;

// KCD
kcd1 =     0.92922200000000000;

// KH
kh1 =     0.98935200000000000;

// EPDC
epdcb1 =     0.05000000000000000;
epdcb2a =     0.27186586842900000;
epdcb2b =     0.35821673881200000;
epdcb3 =     0.30220198578600000;

// QEPDC
qepdcs1 =     3.07231400000000000;

// EPDO
epdob1 =     0.05000000000000000;
epdob2a =    -0.01651866085250000;
epdob2b =     0.17719174122800000;
epdob3 =     0.37500148337200000;

// QEPDO
qepdos1 =     6.02899300000000000;

// EPD
epds1 =     0.59788100000000000;

// EPS
epsb1 =     0.02500000000000000;
epsb2a =     0.24062440268500000;
epsb2b =     0.16250920002900000;
epsb3 =     0.55426950201200000;

// QEPS
qepss1 =    16.04211000000000000;

// KI
kib1 =     0.52534821535152040;
kib2 =     0.24359511952804890;
kib3 =     0.23105666512043030;
kib4 =    -0.12431414744102960;

// KPDC
kpdcs1 =     0.90090700000000000;

// KPDO
kpdos1 =     0.96090800000000000;

// KPS
kpss1 =     0.98925600000000000;

// HKS
hkss1 =     0.19408000000000000;
hkss2 =     0.35962700000000000;
hkss3 =     0.38140300000000000;

// RTPDC
rtpdcs1 =     0.02166600000000000;
rtpdcs2 =    -0.46801500000000000;

// RTPDO
rtpdos1 =     0.04512000000000000;
rtpdos2 =    -0.46801500000000000;

// RTPS
rtpss1 =     0.11369800000000000;

// TAPDDA
tapdda1 =    -0.01785700000000000;

// EX
exa1 =    -0.10494112578145360;
exd0 =     2.13717932268070500;
exd1 =     0.90951689811986520;

// EMO
emoa1 =    -0.26712712016987720;
emod0 =     1.41485636953648900;
emod1 =     1.29215661261832800;
emod2 =     0.42680041941100160;

// CENG
cenga1 =    -0.24542150515617510;
cengb1a =     0.70240755106276310;
cengb1b =     0.73071030490526220;
cengb2 =     0.02260721843178481;
cengb3a =    -0.24629590222339240;
cengb3b =    -0.18682195374457920;

// EMP
emps1 =     2.60137500000000000;
emps2 =    -1.60188500000000000;

// EM
ems1 =     0.91410800000000000;
ems2 =     0.08860500000000000;

// FCBN
fcbns1 =    -3.49160700000000000;
fcbns2 =     4.35484300000000000;
fcbns2a =     0.92200900000000000;
fcbns2b =     0.07799100000000000;
fcbns3 =    -1.10004000000000000;
fcbns4 =     0.96544300000000010;
fcbns5 =     0.24817500000000000;

// FNICN
fnicns1 =     0.97281700000000000;
fnicns2 =    -0.01477000000000000;
fnicns3 =     1.00123900000000000;

// FNILN
fnilns1 =     0.97091500000000000;
fnilns2 =    -0.01232800000000000;
fnilns3 =     0.99896500000000000;

// FYNICN
fynicns1 =     0.17468300000000000;

// FYNILN
fynilns1 =     0.23974900000000000;

// RFYNIC
rfynica1 =    -0.15244691438846210;
rfynica2 =     0.15556971510600210;
rfynica3 =     0.65214382926511630;

// RFYNIL
rfynila1 =    -0.24380622977402910;
rfynila2 =     0.09318731044804289;
rfynila3 =     0.08245051399975228;
rfynila4 =     0.03826362037625722;
rfynila5 =     0.10591076073997820;
rfynila6 =     0.13293345589883290;
rfynila7 =     0.26092418096450080;
rfynila8 =     0.06140865827026788;

// XGDP
xgdps1 =     0.68801900000000000;
xgdps2 =     0.04864300000000000;
xgdps3 =     0.08252700000000000;
xgdps4 =     0.02889300000000000;
xgdps5 =     0.59812600000000000;
xgdps6 =    -0.59402600000000010;
xgdps7 =     0.10535100000000000;
xgdps8 =    -0.13953600000000000;
xgdps9 =     0.06467400000000000;
xgdps10 =     0.11733800000000000;
egfs1 =     0.46659400000000000;
egfs2 =     0.13513400000000000;
egfs3 =     0.39827200000000000;
egss1 =     0.20762000000000000;
egss2 =     0.19392600000000000;
egss3 =     0.59845400000000000;

// XB
xbs1 =     1.29070600000000000;
xbs2 =    -0.14600600000000000;
xbs3 =    -0.14469900000000000;
xgvs1 =     0.26812800000000000;

// XG
xgs1 =     0.98571400000000000;

// XP
xps1 =     0.66185300000000000;
xps2 =     0.04674400000000000;
xps3 =     0.07942500000000000;
xps4 =     0.02780500000000000;
xps5 =     0.10138500000000000;
xps6 =     0.00839300000000000;
xps7 =     0.02905600000000000;
xps8 =     0.02189800000000000;
xps9 =     0.02344000000000000;

// EFS
efss1 =     1.00531400000000000;
efss2 =    -0.77599600000000000;
efss3 =     0.77067000000000000;

// HMFPT
hmfptrho =     0.93999999999999990;

// VEOA
veoaa1 =     0.98800000000000000;

// EMPT
empta1 =     0.10000000000000000;

// HGEMP
hgempa1 =     0.90000000000000000;

// XGVT
xgvta1 =     0.10000000000000000;

// HGXGV
hgxgva1 =     0.90000000000000000;

// XHHIT
xhhita1 =     0.10000000000000000;

// HGXHHI
hgxhhia1 =     0.90000000000000000;

// XGDPT
xgdpts1 =     0.77733000000000000;
xgdpts2 =    -0.01111500000000000;
xgdpts3 =     0.00859400000000000;
xgdpts4 =     0.11310500000000000;
xgdpts5 =     0.11208700000000000;

// LHP
lhpb2 =     0.21565800000000000;
lhpb3 =     0.17399100000000000;
lhpb4 =     0.60923900000000000;
lhpb5 =     0.39076100000000000;
lhpb6 =    -0.06798899999999999;

// LWW
lwwa2 =     0.45995660070444500;
lwwa3 =     0.35017544283190200;

// LEOR
leors1 =    20.04675900000000000;
leors2 =    16.26665500000000000;
leors3 =     2.50796700000000000;
leors4 =     0.27213700000000000;

// LEH
lehs1 =     1.50191200000000000;
lehs2 =     0.01052900000000000;

// LFPR
lfpra1 =     0.25350100000000000;
lfpra2 =    -0.00035700000000000;

// LEFT
lefta1 =     0.90000000000000000;

// LEST
lesta1 =     0.90000000000000000;

// LEPPOT
leppots1 =    -0.01306000000000000;
leppots2 =    -0.01695000000000000;
leppots3 =    -0.15566600000000000;

// JCCAN
jccans1 =     0.13602400000000000;
jccans2 =     0.19787800000000000;
jccans3 =     0.33572100000000000;
jccans4 =     0.14835300000000000;
jccans5 =     0.18202400000000000;

// QYNIDN
qynidb2 =     1.00000000000000000;
qynids1 =    -1.36633900000000000;

// RRMET
rrmeta1 =     0.90480000000000000;

// TRYH
tryhs1 =     0.11919800000000000;
tryhs2 =     0.03075600000000000;
tryhs3 =    -0.10398900000000000;
tryhs4 =    -0.04596500000000000;
yhptns1 =     0.33344000000000000;
yhptns2 =     0.48582300000000000;
yhptns3 =     0.18073700000000000;

// WDNFCN
wdnfca1 =    -0.04153184313843793;
wdnfcb1 =     0.16178999656411170;
wdnfcb2 =     0.31537468081475100;
wdnfcb3 =     0.00198634617075808;
wdnfcs1 =    -2.86404600000000000;
wdnfcs2 =     1.86404600000000000;

// XBN
xbns1 =     1.29070600000000000;
xbns2 =    -0.14600600000000000;
xbns3 =    -0.14469900000000000;

// XGDPN
xgdpns1 =     1.03951100000000000;
xgdpns2 =     0.59812600000000000;
xgdpns3 =    -0.59402600000000010;
xgdpns4 =    -0.13953600000000000;
xgdpns5 =     0.02574300000000000;
xgdpns6 =     0.07019100000000000;

// XGN
xgns1 =     0.98571400000000000;

// YCSN
ycsns1 =     0.69987200000000000;
ycsns2 =    -0.97099600000000000;
ycsns3 =    -0.41619000000000000;
ycsns4 =     0.71631800000000000;
jccacs1 =     1.43690200000000000;
jccacs2 =    -0.24135900000000000;
jccacs3 =    -0.19554299999999990;

// YH
yhs1 =     0.57921900000000000;
yhs2 =     0.14751600000000000;
yhs3 =     0.27326500000000000;

// YHIBN
yhibna1 =     1.00000000000003600;
yhibna3 =     0.05781590372088752;
yhibna4 =     0.43976145398781750;
yhibna6 =     0.00203964594842794;
yhibna7 =     0.04794771504698590;
yhibna8 =     0.00478661146982463;

// YHIN
yhins1 =     0.59332200000000000;
yhins2 =     0.22813600000000000;
yhins3 =     0.00004300000000000;
yhins4 =     0.17849900000000000;

// YHL
yhls1 =    -1.17677700000000000;
yhls2 =    -1.14156000000000000;

// YHP
yhps1 =     0.93453500000000000;
yhps2 =     0.06546500000000000;

// YHPNTN
yhpntns1 =     3.52589700000000000;
yhpntns2 =    -4.51690400000000000;
yhpntns3 =   -17.73142700000000000;
yhpntns4 =    12.82002000000000000;
yhpntns5 =    -8.36559600000000000;
yhpntns6 =    -1.05293800000000000;
yhpntns7 =    -1.81789200000000000;
yhpntns8 =    -0.64552400000000000;

// YHT
yhts1 =     9.24293800000000100;

// YNICPN
ynicpns1 =     8.81205800000000000;
ynicpns2 =    -5.74014400000000000;
ynicpns3 =    -0.64495700000000000;
ynicpns4 =    -0.73607500000000000;
ynicpns5 =    -0.73030100000000000;
ynicpns6 =    -6.86777600000000000;
ynicpns7 =     0.03941900000000000;

// YNIDN
ynidna1 =     0.04854176156750000;
ynidna2 =     0.44400097074800000;

// YNIIN
yniina1 =     0.79747139382107790;
yniina2 =     0.06504200397968549;
yniina3 =     0.20252860617892220;
yniina4 =     0.93528053487107640;
yniina5 =     0.78692469720201120;
yniins1 =     4.83864600000000000;
yniins2 =     2.86404600000000000;
yniins3 =     0.18308900000000000;
yniins4 =     0.07925900000000000;
yniins5 =    -0.00379700000000000;
yniins6 =     0.19852800000000000;
yniins7 =     0.13360100000000000;
yniins8 =     0.01975400000000000;

// YNILN
ynilns1 =     0.83353100000000000;
ynilns2 =     0.04469000000000000;
ynilns3 =     0.12177900000000000;

// YNIN
ynins1 =     1.13274500000000000;
ynins2 =     0.03914700000000000;
ynins3 =    -0.03427800000000000;
ynins4 =    -0.13782500000000000;

// YPN
ypns1 =     0.60335900000000000;
ypns2 =     0.13070700000000000;
ypns3 =     0.26593400000000000;

// PCNIA
pcniaa1 =     0.05345963250200000;
pcniaa21 =     0.19089296942400000;
pcniaa22 =     0.13020188243000000;
pcniaa23 =     0.08779723489500001;
pcniaa24 =     0.18975504477200000;
pcniaa4 =    -0.01439537957800000;
pcniaa5 =     1.00698383632000000;
pcniaa6 =     0.38907365960000000;
pcniaa8 =     0.85406145598700000;
pcnias1 =     0.04936200000000000;
pcnias2 =     0.13882200000000000;

// PXNC
pxnca1 =     0.05609013937800000;
pxnca21 =     0.33560112810700000;
pxnca22 =     0.23601350656000000;
pxnca23 =     0.07411259138400000;
pxnca4 =    -0.01406663313800000;
pxnca5 =     0.78678934778300000;
pxnca6 =     1.35803877798000000;
pxnca8 =     0.63294484245200000;

// PIPL
pipla1 =     0.06495028437800000;
pipla21 =     0.20098403513700000;
pipla22 =     0.24002339229900000;
pipla23 =     0.15628289361100000;
pipla4 =    -0.01434662855500000;
pipla6 =     0.82504538883400000;
pipla7 =     0.15307031699000000;
dlpls1 =     0.07273100000000000;
dlpls2 =     0.20016100000000000;
tsiens1 =     8.06539800000000000;

// QPCNIA
qpcnias1 =     0.18656000000000000;
uqpca2 =    -0.53567007410000000;

// QPXNC
qpxnca1 =     2.98507462687000000;
qpxnca2 =    -1.98507462687000000;

// QPXP
qpxps1 =     0.74786800000000000;

// POILR
poilra1 =    -0.19931742478769890;
poilrb1 =     0.36606340553734860;

// PCENGR
pcengrb1 =    -0.05265757387586706;
pcengrb2 =    -0.02591203086750279;
pcengrb3 =     0.02958683282242650;
pcengrb4 =     0.75614666793982940;

// PCER
pcera1 =    -0.11054933532390130;
pcera2 =     0.04306306304073697;
pcerb1 =     0.13453388478857130;
pcerb2 =     0.20729166053040690;
pcerb3 =     0.04371777025197865;

// PCFR
pcfra1 =    -0.16688185405526880;
pcfrb11 =     0.41520648000181230;
pcfrb12 =    -0.00824831453140233;
pcfrb13 =     0.40917924877889530;

// UCES
ucesa1 =    -0.48384684805145880;
ucesa2 =     0.41995735880580240;
ucesa3 =     0.24277664936872510;
ucesb1 =    -0.04354515853187443;
ucesb2 =     0.76646901217701860;
ucesb3 =     0.26847039509614130;

// UCFS
ucfsa1 =    -0.10152120411058390;
ucfsa2 =     0.06943223210796481;
ucfsb1 =     0.00563710195204836;
ucfsb3 =     0.37434306786105860;

// PMO
pmoa1 =     0.43912200000000000;
pmoa2 =     0.30691300000000000;
pmoa3 =     0.69308700000000000;

// HGPMOR
hgpmora1 =     0.90000000000000000;

// PIMOR
pimors1 =     0.12345600000000000;

// EGFI
egfia1 =    -0.18394330238045290;
egfia2a =    -0.24680170860429420;
egfia2b =    -0.17820833976183210;
egfia3 =     2.01129676230384100;
egfia4a =     0.01046625988566389;
egfia4b =    -0.00824285257827215;

// EGFL
egfla1 =    -0.07986102297274635;
egfla2a =     0.25758086510406960;
egfla2b =    -0.10115236614004040;
egfla3 =     1.33960527247441100;
egfla4a =    -0.00049941403603505;
egfla4b =     0.00070796792063541;

// EGFO
egfoa1 =    -0.17144831350528240;
egfoa2a =    -0.23549871893128250;
egfoa2b =    -0.10695124365892500;
egfoa3 =     1.78020075416061700;
egfoa4a =    -0.00038516765044011;
egfoa4b =     0.00117488981088158;

// EGSI
egsia1 =    -0.20316532027428990;
egsia2a =     0.03363569843997749;
egsia2b =    -0.10711876241400460;
egsia3 =     1.48498254405308900;
egsia4a =     0.00767889395960698;
egsia4b =    -0.00538194005479192;

// EGSL
egsla1 =    -0.12505734167481360;
egsla2a =     0.15110849292099590;
egsla2b =    -0.02057494855531095;
egsla3 =     0.66081178813101190;
egsla4a =    -0.00070396094438263;
egsla4b =     0.00112583020388717;

// EGSO
egsoa1 =    -0.08128878632774075;
egsoa2a =     0.52242478077076290;
egsoa2b =     0.18029049832237940;
egsoa3 =     0.30958319493669360;
egsoa4a =    -0.00136019325294800;
egsoa4b =     0.00139634568679049;

// EGT
egta1 =    -0.10000000000000000;

// EG
egs1 =     0.35502900000000000;
egs2 =     0.64497100000000000;

// GFDBTN
gfdbtns2 =     0.12866400000000000;
gfdbtns3 =    -0.12194800000000000;
gfdbtns4 =     0.00568600000000000;
gfdbtns5 =    -0.00581600000000000;
gfdbtns1 =     0.99341400000000010;

// GFINTN
gfintns1 =    16.04946300000000000;

// GFRECN
gfrecns1 =     0.46260700000000000;
gfrecns2 =     0.11455000000000000;
gfrecns3 =     0.04683000000000000;
gfrecns4 =     0.37601300000000000;
trfcis1 =     4.49050800000000000;

// GFEXPN
gfexpns1 =     0.13101000000000000;
gfexpns2 =     0.15369300000000000;
gfexpns3 =     0.42844700000000000;
gfexpns4 =     0.12735300000000000;
gfexpns5 =     0.02239300000000000;
gfexpns6 =     0.13710400000000010;
gfptns1 =    11.75350600000000000;

// GFDEFN
gfdefns1 =     9.10059300000000000;

// GFTRD
gftrda1 =     0.63905428450608790;
gftrda2 =    -0.00021061937132746;

// GSDBTN
gsdbtns2 =     0.22216600000000000;
gsdbtns3 =    -0.22637900000000000;
gsdbtns4 =     0.04333500000000000;
gsdbtns5 =    -0.02175100000000000;
gsdbtns1 =     0.98262900000000000;

// GSINTN
gsintns1 =     8.01736400000000100;

// GSRECN
gsrecns1 =     0.18761800000000000;
gsrecns2 =     0.03139000000000000;
gsrecns3 =     0.54262000000000000;
gsrecns4 =     0.01253400000000000;
gsrecns5 =     0.22583799999999990;
trscis1 =    25.31811800000000000;
trsibs1 =    10.65960900000000000;
trssis1 =   404.76472300000000000;

// GSEXPN
gsexpns2 =     0.20596300000000000;
gsexpns3 =     0.20223100000000000;
gsexpns4 =     0.00865000000000000;
gsexpns5 =     0.00288500000000000;
gsexpns1 =     0.58027100000000000;
gsptns1 =    40.29360300000000000;

// GSSRPN
gssrpns1 =    20.70475300000000000;

// GSTRD
gstrda1 =     0.74926798279551620;
gstrda2 =    -0.00003881189014207;

// RGFINT
rgfinta1 =     0.04424422209453022;
rgfinta2 =     0.06204177367387913;
rgfintb1 =    -0.17133869041964870;
rgfintb2 =     0.19987202525810740;
rgfintb3 =     0.04423541135296243;

// TFPN
tfpns1 =     1.14578300000000000;
tfpns2 =    -0.11459800000000000;
tfpns3 =    -0.03118500000000000;
tfpns4 =     8.63565000000000000;

// TRFCI
trfcia1 =     0.71407249944833140;
trfcia4 =     0.00172820870546111;
trfcia5 =     0.00382589778662357;

// TRFP
trfpa21 =     0.56817223354206310;
trfpa22 =     0.34942685988922830;
trfpa3 =     0.00034427536674426;

// TRSCI
trscia1 =     0.75946029100417440;
trscia31 =    -0.00062346124814265;
trscia32 =     0.00055480299350093;
trscia4 =     0.12357749586377120;

// TRSIB
trsiba1 =     0.91371152109125410;
trsiba3 =    -0.00004236235627602;

// TRSP
trspa1 =     0.62131574157998890;
trspa21 =     0.91979333369980030;
trspa22 =    -0.54110907528479400;
trspa3 =     0.00000670961673159;
trspa4 =     0.02523637603295423;

// TRSSI
trssia11 =     1.12142509967081700;
trssia12 =    -0.18529742056667950;
trssia3 =    -0.00000279390504096;

// TSPN
tspns1 =     1.14578300000000000;
tspns2 =    -0.11459800000000000;
tspns3 =    -0.03118500000000000;
tspns4 =    33.04425500000000000;

// RTBE
rtbeb11 =     0.75934130984098370;
rtbeb12 =     0.12609272581749860;
rtbeb21 =     0.78269426723110300;
rtbeb22 =    -0.66812830288956260;

// RG5P
rg5pa1 =    -0.32277000000000000;
rg5pa2 =     0.86383000000000000;

// RG10P
rg10pa1 =    -0.55387000000000000;
rg10pa2 =     0.85364000000000000;

// RBAAP
rbaapa1 =    -0.61186000000000000;
rbaapa2 =     0.90088000000000000;

// RCBP
rcbpa1 =    -0.94450000000000000;
rcbpa2 =     0.82646000000000000;

// RCAR
rcara3 =     0.62933284517933360;
rcara41 =     0.22492622732911550;
rcara42 =     0.14574092749152240;

// RME
rmea1 =     0.76720028583715290;
rmea21 =     0.97258818579414850;
rmea22 =    -0.73978847163100300;
rmea3 =    -0.06825682782284648;

// REQP
reqpb1 =     0.93427108611037020;

// WPSN
wpsns1 =     0.25176300000000000;

// RCGAIN
rcgaina1 =     0.41802403831786060;
rcgaina2 =     0.58197596168213880;
rcgaina5 =     0.16255919025238090;

// WPO
wpos2 =     0.06870700000000000;
wpos3 =    -0.06437400000000000;
wpos4 =    -0.00177800000000000;
wpos5 =     0.00802300000000000;
wpos6 =    -0.00632100000000000;
wpos7 =     0.00250200000000000;
wpos1 =     0.99574300000000000;
ydns1 =     1.14745500000000000;
ydns2 =    -0.11804900000000000;
ydns3 =    -0.03039800000000000;

// FXGAP
fxgapb1a =     1.18791873322402700;
fxgapb1b =    -0.36028137863804400;
fxgapb2 =    -0.07043152293092064;
fxgapb3 =     0.03625296278160418;
fxgapb4 =     0.00000000000000000;

// FGDPT
fgdpta1 =    -0.10000000000000000;

// FPI10
fpi10b2 =     0.32062602000000000;
fpi10b3a =     0.38869524000000000;
fpi10b3b =     5.14967970000000000;

// FPI10T
fpi10tb1 =     0.95000000000000000;
fpi10tb2 =     0.05000000000000000;

// FPIC
fpicb1 =     0.65858227904952550;
fpicb2 =     0.34141772095047420;

// FRS10
frs10b1 =     1.00000000000000000;
frs10b2 =     1.00000000000000000;
frs10b3 =     0.50000000000000000;
frs10b4 =     1.00000000000000000;

// FRL10
frl10b1 =    -0.08230000000000000;
frl10b2 =     0.16091000000000000;
frl10b3 =     0.40743000000000000;
frl10b4 =     0.21512000000000000;

// FPXR
fpxrb1 =     0.06000000000000003;
fpxrb2 =     0.19999999999999850;
fpxrs1 =     0.61697700000000000;
fpxrs2 =    -0.73609500000000010;

// FPXRR
fprho =     0.00000000000000000;
fpxrrb1 =     0.18293566655131470;
fpxrrb2 =     0.07323158412381960;

// ZECO
zfecoa1 =     0.09423738003239999;
zfecoa2 =    -0.01468293965926700;
zfecoa6 =     1.04663514379776800;
zfecoa7 =    -0.15580802070493000;

// ZBECO
zbeco1 =    -0.00020185600841500;
zbeco2 =    -0.00010770711388800;
zbeco3 =    -0.00001338884816020;
zbeco4 =     0.00002705572067910;
zbeco5 =     0.00000837898829290;
zbeco6 =     0.00015927187044300;
zbeco7 =    -0.00004385510561060;
zbeco8 =     0.00006708057376320;
zbeco9 =    -0.00009837724391610;
zbeco10 =     0.00003150933560740;
zbeco11 =    -0.00005291089803830;
zbeco12 =    -0.00008327579982630;
zbeco13 =     0.00029589624978500;
zbeco14 =    -0.00019087632688800;
zbeco15 =    -0.00007591010362540;
zbeco16 =     0.00008634911633850;
zbeco17 =     0.00012927932410200;
zbeco18 =    -0.00002898975570000;
zbeco19 =     0.00001771373411670;
zbeco20 =    -0.00004066650916700;
zbeco21 =     0.00001188876353250;
zbeco22 =    -0.00002891856273070;
zbeco23 =    -0.00005175269643740;
zbeco24 =    -0.00006478953217690;
zbeco25 =    -0.00000840906651595;
zbeco26 =     0.00000987830694352;
zbeco27 =     0.73121282946699990;
zbeco28 =     0.00492542542754000;
zbeco29 =    -0.00743357386590000;
zbeco30 =    -0.00236655199096000;
zbeco31 =     0.00268006848485000;

// ZECD
zfecda1 =     0.17942469395300000;
zfecda2 =     0.02512800035771100;
zfecda6 =     0.66125807091722000;
zfecda7 =     0.14004761442866300;

// ZBECD
zbecda1 =    -0.00078349984932000;
zbecda2 =    -0.00063621853633000;
zbecda3 =    -0.00053723271603800;
zbecda4 =    -0.00040248227711100;
zbecda5 =    -0.00161294241054000;
zbecda6 =     0.00207246868793000;
zbecda7 =    -0.00045775208328700;
zbecda8 =     0.00047444469884100;
zbecda9 =     0.00037582027728600;
zbecda10 =    -0.00030605622070600;
zbecda11 =    -0.00056865794095800;
zbecda12 =     0.00019947691554300;
zbecda13 =     0.00235943337880000;
zbecda14 =    -0.00047621889295100;
zbecda15 =     0.00024731432280800;
zbecda16 =    -0.00038864759451100;
zbecda17 =    -0.00014640626294500;
zbecda18 =     0.00018358690423600;
zbecda19 =     0.00000588721219744;
zbecda20 =     0.00003652261066340;
zbecda21 =     0.00003282019505060;
zbecda22 =     0.00004841590303630;
zbecda23 =    -0.00003782420118300;
zbecda24 =     0.00031726155638900;
zbecda25 =     0.00007685142936440;
zbecda26 =     0.00043304041057000;
zbecda27 =     1.03036609913000000;
zbecda28 =     0.02070207959270000;
zbecda29 =     0.00530545512357000;
zbecda30 =    -0.01408637493270000;
zbecda31 =    -0.01357568977010000;

// ZGAPC2
zfgapc2a1 =     0.03565066659701600;
zfgapc2a6 =     0.66125807091722000;
zfgapc2a7 =     0.14004761442866300;

// ZBGPC2
zbgapc2a1 =    -0.03417384897760000;
zbgapc2a2 =    -0.00942618734464000;
zbgapc2a3 =    -0.00483030060919000;
zbgapc2a4 =     0.00548473735054000;
zbgapc2a5 =    -0.05035021597730000;
zbgapc2a6 =    -0.01075451297140000;
zbgapc2a7 =     0.02347953924390000;
zbgapc2a8 =    -0.00199327599458000;
zbgapc2a9 =     0.15067844681500000;
zbgapc2a10 =    -0.00561735717748000;
zbgapc2a11 =    -0.02044715177540000;
zbgapc2a12 =    -0.00532923566468000;
zbgapc2a13 =     0.04294559958080000;
zbgapc2a14 =     0.03961846569940000;

// ZEH
zfeha1 =     0.09373761037809999;
zfeha2 =    -0.02234683784639400;
zfeha3 =    -0.02054297025152100;
zfeha6 =     1.13596398491524200;
zfeha7 =    -0.01924379752798400;
zfeha8 =    -0.21915397852215800;

// ZBEH
zbeha1 =    -0.00054758352828400;
zbeha2 =    -0.00037554155272200;
zbeha3 =    -0.00051413247866400;
zbeha4 =    -0.00006383837989620;
zbeha5 =     0.00008430516195090;
zbeha6 =     0.00104649751949000;
zbeha7 =     0.00007374531712390;
zbeha8 =     0.00029776290633000;
zbeha9 =     0.00063141959616600;
zbeha10 =    -0.00125748752670000;
zbeha11 =    -0.00082355569590900;
zbeha12 =     0.00030693367968300;
zbeha13 =     0.00150109593957000;
zbeha14 =    -0.00302677981831000;
zbeha15 =     0.00060448182531500;
zbeha16 =    -0.00032930355300200;
zbeha17 =     0.00038590219172800;
zbeha18 =     0.00061968249898300;
zbeha19 =    -0.00023298429664600;
zbeha20 =     0.00002779540716980;
zbeha21 =    -0.00004690360733010;
zbeha22 =    -0.00000945130159196;
zbeha23 =     0.00005731769689090;
zbeha24 =     0.00039714132723800;
zbeha25 =    -0.00016475596383200;
zbeha26 =     0.00050229398876300;
zbeha27 =     0.55862755912600000;
zbeha28 =    -0.01298838003050000;
zbeha29 =    -0.01856956220230000;
zbeha30 =    -0.02157250102640000;
zbeha31 =    -0.00962731742328000;

// ZVPDC
zfvpdc1 =     0.05000000000000000;
zfvpdc2 =    -0.01271096764600600;
zfvpdc3 =    -0.01685753644179700;
zfvpdc6 =     1.19742855106042000;
zfvpdc7 =     0.08293137591583299;
zfvpdc8 =    -0.33715072883594400;

// ZBVPDC
zbvpdca1 =    -0.00069669830347700;
zbvpdca2 =    -0.00096710947136800;
zbvpdca3 =    -0.00058734955470300;
zbvpdca4 =    -0.00003462061678850;
zbvpdca5 =     0.00124708033597000;
zbvpdca6 =     0.00014618413722700;
zbvpdca7 =    -0.00076851848711400;
zbvpdca8 =     0.00014656503899500;
zbvpdca9 =    -0.00077131102507900;
zbvpdca10 =     0.00228577794634000;
zbvpdca11 =    -0.00084420120581500;
zbvpdca12 =     0.00024415220311900;
zbvpdca13 =     0.00082914961447300;
zbvpdca14 =    -0.00021698547799700;
zbvpdca15 =    -0.00000000000000189;
zbvpdca16 =     0.00000000000000060;
zbvpdca17 =    -0.00000000000000420;
zbvpdca18 =    -0.00000000000000125;
zbvpdca19 =     0.00519324083494000;
zbvpdca20 =    -0.01038606273960000;
zbvpdca21 =    -0.02125303036360000;
zbvpdca22 =    -0.00694250517399000;

// ZVPDO
zfvpdo1 =     0.05000000000000000;
zfvpdo2 =     0.00096340104241200;
zfvpdo3 =    -0.00833857246549300;
zfvpdo6 =     0.91481171236454990;
zfvpdo7 =     0.18603947015811200;
zfvpdo8 =    -0.16677144930986400;

// ZBVPDO
zbvpdoa1 =     0.00047639431931000;
zbvpdoa2 =    -0.00002628001535660;
zbvpdoa3 =     0.00019692311161100;
zbvpdoa4 =     0.00008106193209880;
zbvpdoa5 =     0.00024014956273000;
zbvpdoa6 =     0.00031638032545900;
zbvpdoa7 =    -0.00012171375451200;
zbvpdoa8 =     0.00018884090212900;
zbvpdoa9 =    -0.00062365703580600;
zbvpdoa10 =    -0.00072809934766300;
zbvpdoa11 =    -0.00034209638831400;
zbvpdoa12 =     0.00000670597110055;
zbvpdoa13 =    -0.00002180133402510;
zbvpdoa14 =     0.00020058700945100;
zbvpdoa15 =     0.00000000000000023;
zbvpdoa16 =     0.00000000000000035;
zbvpdoa17 =    -0.00000000000000054;
zbvpdoa18 =     0.00000000000000067;
zbvpdoa19 =     0.00863159343403000;
zbvpdoa20 =     0.00444051246810000;
zbvpdoa21 =     0.00640345062465000;
zbvpdoa22 =     0.00993772001366000;

// ZVPS
zfvps1 =     0.02500000000000000;
zfvps2 =    -0.00569935499061300;
zfvps3 =    -0.00382380897484200;
zfvps6 =     1.19131191463130000;
zfvps7 =    -0.07502184063082200;
zfvps8 =    -0.15295235899369500;

// ZBVPS
zbvpsa1 =     0.00025407183453800;
zbvpsa2 =    -0.00000563212048438;
zbvpsa3 =     0.00016405037043300;
zbvpsa4 =     0.00014998009969600;
zbvpsa5 =    -0.00001568896704760;
zbvpsa6 =     0.00012426274692200;
zbvpsa7 =    -0.00017385735012100;
zbvpsa8 =     0.00000036150249237;
zbvpsa9 =     0.00006492206775420;
zbvpsa10 =    -0.00056247018418300;
zbvpsa11 =    -0.00031190380977300;
zbvpsa12 =    -0.00009254579012430;
zbvpsa13 =     0.00006294934434650;
zbvpsa14 =     0.00002089066202520;
zbvpsa15 =     0.00000000000000063;
zbvpsa16 =     0.00000000000000069;
zbvpsa17 =    -0.00000000000000125;
zbvpsa18 =     0.00000000000000018;
zbvpsa19 =    -0.00065517948190500;
zbvpsa20 =    -0.00854451453699000;
zbvpsa21 =    -0.00449314225851000;
zbvpsa22 =    -0.00525514423712000;

// ZXBC
zfxbc1 =     0.05000000000000000;
zfxbc2 =    -0.01271096764600600;
zfxbc3 =    -0.01685753644179700;
zfxbc6 =     1.19742855106042000;
zfxbc7 =     0.08293137591583299;
zfxbc8 =    -0.33715072883594400;

// ZBXBC
zbxbca1 =     0.00001168690696110;
zbxbca2 =     0.00000027489525101;
zbxbca3 =    -0.00007041002248690;
zbxbca4 =     0.00006203428244370;
zbxbca5 =    -0.00071061701021800;
zbxbca6 =     0.00003520401254230;
zbxbca7 =     0.00049607791838900;
zbxbca8 =     0.00022437821662300;
zbxbca9 =    -0.00004504313733650;
zbxbca10 =    -0.00000358606216894;
zbxbca11 =    -0.00045527882330600;
zbxbca12 =     0.00022081937939800;
zbxbca13 =     0.00013507772778500;
zbxbca14 =     0.00002263081160320;
zbxbca15 =     0.16637475207700000;
zbxbca16 =     0.11377234483300000;
zbxbca17 =     0.05463218213920000;
zbxbca18 =     0.02498842320650000;
zbxbca19 =     0.00000000000000003;
zbxbca20 =     0.00000000000000006;
zbxbca21 =    -0.00000000000000001;
zbxbca22 =    -0.00000000000000001;
zbxbca23 =     0.00000000000000009;

// ZXBO
zfxbo1 =     0.05000000000000000;
zfxbo2 =     0.00096340104241200;
zfxbo3 =    -0.00833857246549300;
zfxbo6 =     0.91481171236454990;
zfxbo7 =     0.18603947015811200;
zfxbo8 =    -0.16677144930986400;

// ZBXBO
zbxboa1 =     0.00013358221184500;
zbxboa2 =     0.00025449248606000;
zbxboa3 =     0.00024216223862500;
zbxboa4 =     0.00019922980751100;
zbxboa5 =    -0.00120695652995000;
zbxboa6 =     0.00021539884025500;
zbxboa7 =     0.00062732626277400;
zbxboa8 =     0.00046728377112400;
zbxboa9 =    -0.00010305234420600;
zbxboa10 =    -0.00082946674404100;
zbxboa11 =    -0.00073267166716700;
zbxboa12 =     0.00040062722861600;
zbxboa13 =     0.00020087801825300;
zbxboa14 =     0.00006406613031330;
zbxboa15 =     0.26517731226400000;
zbxboa16 =     0.20737912110100000;
zbxboa17 =     0.11484594641300000;
zbxboa18 =     0.05920911545860000;
zbxboa19 =     0.00000000000000017;

// ZXBS
zfxbs1 =     0.02500000000000000;
zfxbs2 =    -0.00569935499061300;
zfxbs3 =    -0.00382380897484200;
zfxbs6 =     1.19131191463130000;
zfxbs7 =    -0.07502184063082200;
zfxbs8 =    -0.15295235899369500;

// ZBXBS
zbxbsa1 =     0.00006466146382360;
zbxbsa2 =     0.00018206700053800;
zbxbsa3 =     0.00018025347459000;
zbxbsa4 =     0.00015930032265100;
zbxbsa5 =    -0.00097330518389900;
zbxbsa6 =     0.00016714871311500;
zbxbsa7 =     0.00050141505961000;
zbxbsa8 =     0.00037568165587000;
zbxbsa9 =    -0.00007094024469490;
zbxbsa10 =    -0.00058628226160300;
zbxbsa11 =    -0.00059228248496000;
zbxbsa12 =     0.00031145580158700;
zbxbsa13 =     0.00014761584440900;
zbxbsa14 =     0.00004150818983730;
zbxbsa15 =     0.20551375996400000;
zbxbsa16 =     0.16267287731300000;
zbxbsa17 =     0.09069643445880000;
zbxbsa18 =     0.04720769571600000;
zbxbsa19 =     0.00000000000000013;

// ZPIB5
zfpib5a1 =     0.05882353000000000;
zfpib5a6 =     0.94117647000000000;

// ZBPIB5
zbpib5a1 =     0.05567364513686002;
zbpib5a2 =     0.03213388983239129;
zbpib5a3 =     0.04831377185678321;
zbpib5a4 =     0.00136204448210009;
zbpib5a5 =    -0.24800594648467000;
zbpib5a6 =    -0.04939153761848081;
zbpib5a7 =    -0.05430346140256770;
zbpib5a8 =     0.02028681910282160;
zbpib5a9 =     0.33141412640241090;
zbpib5a10 =     0.71070127703966430;
zbpib5a11 =     0.12892933910636320;
zbpib5a12 =     0.06387504744496386;
zbpib5a13 =     0.03301212283629473;
zbpib5a14 =     0.02443749879731572;
zbpib5a15 =     0.06277961198499149;
zbpib5a16 =     0.03917746207061756;
zbpib5a17 =     0.01543767554087694;
zbpib5a18 =     0.03442062205617091;

// ZPI5
zfpi5a1 =     0.05882353000000000;
zfpi5a6 =     0.94117647000000000;

// ZBPI5
zbpi5a1 =     0.05391077624399176;
zbpi5a2 =     0.03901712531502367;
zbpi5a3 =     0.03927185388265116;
zbpi5a4 =     0.01575352521491724;
zbpi5a5 =    -0.21106105124113900;
zbpi5a6 =    -0.02988969416361464;
zbpi5a7 =    -0.05156341078808968;
zbpi5a8 =     0.01860931728420345;
zbpi5a9 =     0.27390483890847440;
zbpi5a10 =     0.85204671934369540;
zbpi5a11 =     0.10802754829794630;
zbpi5a12 =     0.04842487667477174;
zbpi5a13 =     0.03256529521752238;
zbpi5a14 =     0.03168493732897470;

// ZPI10
zfpi10a1 =     0.03448275800000000;
zfpi10a6 =     0.96551724200000000;

// ZBPI10
zbpi10a1 =     0.02680390337344685;
zbpi10a2 =     0.02213469369401856;
zbpi10a3 =     0.02319157150217065;
zbpi10a4 =     0.01066838693985439;
zbpi10a5 =    -0.11964153609454870;
zbpi10a6 =    -0.01172045837545768;
zbpi10a7 =    -0.03358213149248879;
zbpi10a8 =     0.01294319161637936;
zbpi10a9 =     0.15200093434590220;
zbpi10a10 =     0.91720144449084520;
zbpi10a11 =     0.04302632217545288;
zbpi10a12 =     0.03080385595946715;
zbpi10a13 =     0.01768194919194954;
zbpi10a14 =     0.01632950937387312;

// ZPIC30
zfpic30a1 =     0.02439024400000000;
zfpic30a6 =     0.97560975600000000;

// ZBPC30
zbpic30a1 =     0.03451669387236234;
zbpic30a2 =     0.01241131681941162;
zbpic30a3 =     0.01293539821964619;
zbpic30a4 =     0.00587661239283005;
zbpic30a5 =    -0.06714286824325572;
zbpic30a6 =    -0.00691976376837769;
zbpic30a7 =    -0.01860664153930819;
zbpic30a8 =     0.00711045840305668;
zbpic30a9 =     0.08555881514681835;
zbpic30a10 =     0.93425997869697150;
zbpic30a11 =     0.02522702110104343;
zbpic30a12 =     0.01713160832385359;
zbpic30a13 =     0.00999442576057152;
zbpic30a14 =     0.00925754106250087;

// ZLHP
zflhpa1 =     0.21565800000000000;
zflhpa2 =    -0.03603665805531100;
zflhpa6 =     0.93916634000000000;
zflhpa7 =    -0.16710095640000000;

// ZBLHP
zblhpa1 =    -0.00012164308818930;
zblhpa2 =    -0.00003948324929119;
zblhpa3 =    -0.00001862502664693;
zblhpa4 =     0.00000267400011575;
zblhpa5 =    -0.00017544590287942;
zblhpa6 =    -0.00001691596159103;
zblhpa7 =     0.00008612447214424;
zblhpa8 =     0.00002411618017216;
zblhpa9 =     0.00008212121215405;
zblhpa10 =     0.00017707736401169;
zblhpa11 =     0.00006555019879728;
zblhpa12 =    -0.00002785201098344;
zblhpa13 =    -0.00006353388753394;
zblhpa14 =    -0.00000596397134507;
zblhpa15 =     0.07878941473545197;
zblhpa16 =     0.70924982618563740;

// ZYNID
zfynida1 =     0.04854176156750000;
zfynida2 =    -0.02069910672318000;
zfynida6 =     1.36755002499689000;
zfynida7 =    -0.42641853230637900;

// ZBYNID
zbynida1 =     0.00030898234145384;
zbynida2 =     0.00036985484591112;
zbynida3 =     0.00032821158181724;
zbynida4 =     0.00023012815985315;
zbynida5 =    -0.00022480329593001;
zbynida6 =    -0.00012162539209289;
zbynida7 =    -0.00014716245522356;
zbynida8 =     0.00005494825807706;
zbynida9 =     0.00043864288516904;
zbynida10 =    -0.00123717692903458;
zbynida11 =    -0.00155683951903894;
zbynida12 =     0.00044103731598766;
zbynida13 =    -0.00026691548775063;
zbynida14 =     0.00019681764825921;
zbynida15 =    -0.01060485115865744;
zbynida16 =    -0.00888229828213106;
zbynida17 =    -0.00157516768991789;
zbynida18 =    -0.00262486325554603;
zbynida19 =     0.49665067337795680;

// ZLURC
zflurca1 =     0.00326076272320900;
zflurca6 =     0.93900514143851230;

// ZBLURC
zblurca1 =     0.28393169218200000;
zblurca2 =     0.18010662842600000;
zblurca3 =     0.06835938868020000;
zblurca4 =    -0.01019496229590000;
zblurca5 =     0.01606736651920000;
zblurca6 =     0.00768810718250000;
zblurca7 =     1.62961937992000000;
zblurca8 =     0.98583759951700000;
zblurca9 =     1.29581757262000000;
zblurca10 =     0.39243565933500000;
zblurca11 =     0.11904743932600000;
zblurca12 =     0.05509979329950000;
zblurca13 =     0.01294018746180000;
zblurca14 =     0.00644262635261000;
zblurca15 =    -0.00720796498812000;
zblurca16 =    -0.00326012093182000;
zblurca17 =     0.00412920933315000;
zblurca18 =    -5.03694462581000000;
zblurca19 =    -0.15777950338800000;
zblurca20 =     0.02205444549860000;
zblurca21 =     3.84601287561000000;
zblurca22 =    -4.37836583490000000;
zblurca23 =     0.04727475224700000;
zblurca24 =    -0.02860123603570000;
zblurca25 =     0.00052595425340400;
zblurca26 =     0.00162857864099000;
zblurca27 =     1.41987635642000000;
zblurca28 =    -0.01482248175270000;
zblurca29 =    -0.45481657396100000;
zblurca30 =    -0.02501565206850000;
zblurca31 =    -0.92513776404300000;
zblurca32 =    -0.00612657855148401;
zblurca33 =     0.00000016093328871;

// ZLURL
zflurla1 =     0.00468679810814500;
zflurla6 =     0.92784022190159510;

// ZBLURL
zblurla1 =     0.29540669976300000;
zblurla2 =     0.18912542666600000;
zblurla3 =     0.07225684868720000;
zblurla4 =    -0.01197834798230000;
zblurla5 =     0.01861238557840000;
zblurla6 =     0.00904673259490000;
zblurla7 =     2.14036280782000000;
zblurla8 =     1.28755297763000000;
zblurla9 =     1.83134136939000000;
zblurla10 =     0.52470702606200000;
zblurla11 =     0.11804178447000000;
zblurla12 =     0.05487522293530000;
zblurla13 =     0.01290125674610000;
zblurla14 =     0.00747747121274000;
zblurla15 =    -0.00836745377148000;
zblurla16 =    -0.00379655839378000;
zblurla17 =     0.00478140860494000;
zblurla18 =    -6.54244397096000000;
zblurla19 =    -0.18025023489100000;
zblurla20 =     0.02310529341190000;
zblurla21 =     3.34919704459000000;
zblurla22 =    -3.90594746287000000;
zblurla23 =     0.07116999316470000;
zblurla24 =    -0.04298395704400000;
zblurla25 =     0.00077540858289400;
zblurla26 =     0.00244545908018000;
zblurla27 =     2.10579874240000000;
zblurla28 =    -0.01775891357070000;
zblurla29 =    -0.68271299680300000;
zblurla30 =    -0.03755730111230000;
zblurla31 =    -1.36767345775000000;
zblurla32 =    -0.00709238863012249;
zblurla33 =     0.00000001255221325;

// ZLURNC
zflurnca1 =     0.00350350921202600;
zflurnca6 =     0.93753787651665910;

// ZBLRNC
zblurnca1 =     0.26075316409400000;
zblurnca2 =     0.16689887871500000;
zblurnca3 =     0.06375388441250000;
zblurnca4 =    -0.01053724806550000;
zblurnca5 =     0.01638442456540000;
zblurnca6 =     0.00796287062009000;
zblurnca7 =     1.86997843555000000;
zblurnca8 =     1.13543591714000000;
zblurnca9 =     1.58335356728000000;
zblurnca10 =     0.45692015011400000;
zblurnca11 =     0.10431184624500000;
zblurnca12 =     0.04848801849250000;
zblurnca13 =     0.01139867638940000;
zblurnca14 =     0.00658159145810000;
zblurnca15 =    -0.00736554135807000;
zblurnca16 =    -0.00334137925594000;
zblurnca17 =     0.00420920360855000;
zblurnca18 =    -5.71527099232000000;
zblurnca19 =    -0.15873356618500000;
zblurnca20 =     0.02039387605290000;
zblurnca21 =     2.97123843662000000;
zblurnca22 =    -3.46261042228000000;
zblurnca23 =     0.06051627989170000;
zblurnca24 =    -0.03677865156520000;
zblurnca25 =     0.00068704150139400;
zblurnca26 =     0.00209325042595000;
zblurnca27 =     1.81803689880000000;
zblurnca28 =    -0.02181839087810000;
zblurnca29 =    -0.58474335490600000;
zblurnca30 =    -0.03214836935590000;
zblurnca31 =    -1.17924296297000000;
zblurnca32 =    -0.00624367342072864;
zblurnca33 =     0.00000002007067833;

// ZPC
zfpca1 =     0.02448046144615500;
zfpca6 =     0.93900514143851230;

// ZBPC
zbpca1 =     0.03390311811020000;
zbpca2 =     0.02250001702610000;
zbpca3 =     0.00884745547591000;
zbpca4 =    -0.00226956866028000;
zbpca5 =     0.00284904560626000;
zbpca6 =     0.00144758023470000;
zbpca7 =    -0.00146357682230000;
zbpca8 =     0.00088330306785200;
zbpca9 =    -0.00485408795065000;
zbpca10 =     0.00360734692946000;
zbpca11 =     0.01222117203420000;
zbpca12 =     0.00575642883151000;
zbpca13 =     0.00137131989204000;
zbpca14 =     0.00132408561551000;
zbpca15 =    -0.00134474309815000;
zbpca16 =    -0.00068075003887100;
zbpca17 =     0.00076719703702100;
zbpca18 =     0.31648417187000000;
zbpca19 =     0.02573479728440000;
zbpca20 =     0.00251587485526000;
zbpca21 =     0.33515201221900000;
zbpca22 =    -0.40039856641000000;
zbpca23 =    -0.00050450954969700;
zbpca24 =     0.00033745185396700;
zbpca25 =    -0.00001014182311530;
zbpca26 =    -0.00001952585926900;
zbpca27 =    -0.02103518529650000;
zbpca28 =     0.00131417340137000;
zbpca29 =     0.00550912383069000;
zbpca30 =     0.00029887435156300;
zbpca31 =     0.01391625058210000;
zbpca32 =    -0.00113776738029370;
zbpca33 =     0.00000003829701120;

// ZPL
zfpla1 =     0.02905944107132800;
zfpla6 =     0.92784022190159510;

// ZBPL
zbpla1 =     0.01594066716320000;
zbpla2 =     0.01017207691320000;
zbpla3 =     0.00387665420964000;
zbpla4 =     0.00627414353060000;
zbpla5 =    -0.00894525117045000;
zbpla6 =    -0.00457191227565000;
zbpla7 =     0.01610456683080000;
zbpla8 =     0.01338988229670000;
zbpla9 =     0.00376196405976000;
zbpla10 =     0.01015620027030000;
zbpla11 =     0.02527629273130000;
zbpla12 =     0.01221257691340000;
zbpla13 =     0.00292567659960000;
zbpla14 =    -0.00386814047146000;
zbpla15 =     0.00407329995437000;
zbpla16 =     0.00197501639189000;
zbpla17 =    -0.00234066693159000;
zbpla18 =     0.29628555251900000;
zbpla19 =     0.00194223789102000;
zbpla20 =    -0.01012078609630000;
zbpla21 =     0.58661265123700000;
zbpla22 =    -0.21388421318700000;
zbpla23 =    -0.00052807249083500;
zbpla24 =     0.00035194574324400;
zbpla25 =    -0.00001043504380550;
zbpla26 =    -0.00002033844805070;
zbpla27 =    -0.02170401500220000;
zbpla28 =     0.00132450254056000;
zbpla29 =     0.00573888606939000;
zbpla30 =     0.00031190358223100;
zbpla31 =     0.01433874180440000;
zbpla32 =     0.00347010297157318;
zbpla33 =     0.00000000240219906;

// ZPNC
zfpnca1 =     0.02212862975317600;
zfpnca6 =     0.93753787651665910;

// ZBPNC
zbpnca1 =     0.03712734358180000;
zbpnca2 =     0.02490802001500000;
zbpnca3 =     0.00985692756190000;
zbpnca4 =    -0.00276823974025000;
zbpnca5 =     0.00342660519733000;
zbpnca6 =     0.00180736316642000;
zbpnca7 =    -0.00212615652951000;
zbpnca8 =     0.00032520517486000;
zbpnca9 =    -0.00535082183368000;
zbpnca10 =     0.00317719869026000;
zbpnca11 =     0.01256301347000000;
zbpnca12 =     0.00593789319748000;
zbpnca13 =     0.00139450968329000;
zbpnca14 =     0.00154389822344000;
zbpnca15 =    -0.00161207931325000;
zbpnca16 =    -0.00080119175073900;
zbpnca17 =     0.00091631124331100;
zbpnca18 =     0.26391838288400000;
zbpnca19 =     0.02330096034520000;
zbpnca20 =     0.00271041999057000;
zbpnca21 =     0.28460956002700000;
zbpnca22 =    -0.35648166176500000;
zbpnca23 =    -0.00050404463747800;
zbpnca24 =     0.00033527684788600;
zbpnca25 =    -0.00000987073436105;
zbpnca26 =    -0.00001936308588840;
zbpnca27 =    -0.02054200255490000;
zbpnca28 =     0.00123748010925000;
zbpnca29 =     0.00546277904903000;
zbpnca30 =     0.00029660122282000;
zbpnca31 =     0.01357057928830000;
zbpnca32 =    -0.00135747338304021;
zbpnca33 =    -0.00000003126870382;

// ZDIVGR
zfdivgra1 =     0.02439024300000000;
zfdivgra6 =     0.97560975700000000;

// ZBDVGR
zbdivgra1 =     0.00621201294074858;
zbdivgra2 =    -0.00019787038664765;
zbdivgra3 =    -0.01089541245364467;
zbdivgra4 =    -0.00004989707758885;
zbdivgra5 =     0.04647696371965138;
zbdivgra6 =    -0.01296848700462050;
zbdivgra7 =    -0.06536180677042758;
zbdivgra8 =    -0.02332466874883909;
zbdivgra9 =     0.05517799880354345;
zbdivgra10 =     0.00493116697828100;
zbdivgra11 =    -0.17298153825903380;
zbdivgra12 =     0.00393719877631504;
zbdivgra13 =    -0.00489670023933186;
zbdivgra14 =     0.02200680363042373;
zbdivgra15 =     0.01911217425225434;
zbdivgra16 =    -0.00358256292639525;
zbdivgra17 =    -0.00408853516721142;
zbdivgra18 =     0.00014323147627970;
zbdivgra19 =     0.98841569236441280;

// ZRFF5
zfrff5a1 =     0.05882353000000000;
zfrff5a6 =     0.94117647000000000;

// ZBRFF5
zbrff5a1 =     0.01738406617687738;
zbrff5a2 =    -0.01351066357405079;
zbrff5a3 =    -0.02841001959931608;
zbrff5a4 =    -0.04003646394372262;
zbrff5a5 =     0.20776639730071160;
zbrff5a6 =    -0.06937750790574720;
zbrff5a7 =     0.15077774160941050;
zbrff5a8 =    -0.05046426011133809;
zbrff5a9 =     0.76129762910680140;
zbrff5a10 =     0.06457308094049631;
zbrff5a11 =     0.31523118713586270;
zbrff5a12 =    -0.12151695018981950;
zbrff5a13 =    -0.03042930815673905;
zbrff5a14 =     0.00421089389491502;

// ZRFF10
zfrff10a1 =     0.03448275800000000;
zfrff10a6 =     0.96551724200000000;

// ZBRF10
zbrff10a1 =    -0.00098386813091711;
zbrff10a2 =    -0.01154947002098994;
zbrff10a3 =    -0.01856067673930048;
zbrff10a4 =    -0.02258952591560004;
zbrff10a5 =     0.12748318726368140;
zbrff10a6 =    -0.03445563439517701;
zbrff10a7 =     0.08994970038820871;
zbrff10a8 =    -0.02763882351885405;
zbrff10a9 =     0.84466157026166690;
zbrff10a10 =     0.05368354080736326;
zbrff10a11 =     0.17631476563255800;
zbrff10a12 =    -0.07238344095970957;
zbrff10a13 =    -0.02203885485429868;
zbrff10a14 =    -0.00122594837008183;

// ZRFF30
zfrff30a1 =     0.02439024400000000;
zfrff30a6 =     0.97560975600000000;

// ZBRF30
zbrff30a1 =     0.00023181962342189;
zbrff30a2 =    -0.00638200035294065;
zbrff30a3 =    -0.01049331700928818;
zbrff30a4 =    -0.01296291883159351;
zbrff30a5 =     0.07229552305430142;
zbrff30a6 =    -0.02004695882116786;
zbrff30a7 =     0.05113301186255675;
zbrff30a8 =    -0.01593371671633315;
zbrff30a9 =     0.91255214062042690;
zbrff30a10 =     0.02960641657066043;
zbrff30a11 =     0.10082747132960990;
zbrff30a12 =    -0.04122459154708843;
zbrff30a13 =    -0.01224726450184477;
zbrff30a14 =    -0.00043267589033968;

// ZGAP05
zfgap05a1 =     0.05882353000000000;
zfgap05a6 =     0.94117647000000000;

// ZBGP05
zbgap05a1 =    -0.32608959566429200;
zbgap05a2 =    -0.07805403033857826;
zbgap05a3 =    -0.02097929467812309;
zbgap05a4 =     0.05260004208055552;
zbgap05a5 =    -0.42160743179976860;
zbgap05a6 =     0.06443033493220028;
zbgap05a7 =     0.02514503685340025;
zbgap05a8 =     0.10530588117600510;
zbgap05a9 =     0.22672617883792790;
zbgap05a10 =     0.37252287860063580;
zbgap05a11 =     0.50579493794377250;
zbgap05a12 =     0.05709287894937924;
zbgap05a13 =    -0.09109345447553752;
zbgap05a14 =    -0.06657514710154722;

// ZGAP10
zfgap10a1 =     0.03448275800000000;
zfgap10a6 =     0.96551724200000000;

// ZBGP10
zbgap10a1 =    -0.16236874115912860;
zbgap10a2 =    -0.04131276973974999;
zbgap10a3 =    -0.01401788495372624;
zbgap10a4 =     0.02229013220658072;
zbgap10a5 =    -0.18126026332729640;
zbgap10a6 =     0.02970451549982506;
zbgap10a7 =     0.02537548745190779;
zbgap10a8 =     0.04723152051622648;
zbgap10a9 =     0.07894873985910435;
zbgap10a10 =     0.19540926364617450;
zbgap10a11 =     0.26320801431592170;
zbgap10a12 =     0.01550603962355681;
zbgap10a13 =    -0.04939835794385838;
zbgap10a14 =    -0.03425621234444853;

// ZGAP30
zfgap30a1 =     0.02439024300000000;
zfgap30a6 =     0.97560975700000000;

// ZBGP30
zbgap30a1 =    -0.09410103481690320;
zbgap30a2 =    -0.02379999939484767;
zbgap30a3 =    -0.00789824937590626;
zbgap30a4 =     0.01317879893602168;
zbgap30a5 =    -0.10709361486502510;
zbgap30a6 =     0.01736786549391348;
zbgap30a7 =     0.01380949450164327;
zbgap30a8 =     0.02771970349707603;
zbgap30a9 =     0.04819655137216727;
zbgap30a10 =     0.11262048465175780;
zbgap30a11 =     0.15170199552027650;
zbgap30a12 =     0.00988507129522261;
zbgap30a13 =    -0.02832234490372199;
zbgap30a14 =    -0.01976033119310949;

// ZYH
zfyha1 =     0.75000000000000000;
zfyha2 =     1.00752819544453400;

// ZBYH
zbyha1 =    -0.00058624990254500;
zbyha2 =     0.00034868237334100;
zbyha3 =     0.00015952289554200;
zbyha4 =     0.00047593359786500;
zbyha5 =    -0.00130805476445000;
zbyha6 =     0.00020752206145500;
zbyha7 =     0.00076808633468800;
zbyha8 =     0.00029997657217800;
zbyha9 =     0.00418891644276000;
zbyha10 =     0.00038792791642900;
zbyha11 =    -0.00075573669929600;
zbyha12 =    -0.00086362218550500;
zbyha13 =    -0.00039788896424500;
zbyha14 =     0.00003246979615440;
zbyha15 =     0.00248270505543000;
zbyha16 =     0.00055539672932300;
zbyha17 =    -0.00010467535389300;
zbyha18 =     0.00007878181517560;

// ZYHT
zfyhta1 =     0.75000000000000000;
zfyhta2 =     1.00752819544453400;

// ZBYHT
zbyhta1 =    -0.00111067670859000;
zbyhta2 =    -0.00036433283800000;
zbyhta3 =    -0.00028405434747100;
zbyhta4 =     0.00008732944024670;
zbyhta5 =    -0.00131264748228000;
zbyhta6 =     0.00104584195368000;
zbyhta7 =     0.00027953519293100;
zbyhta8 =     0.00024045067540200;
zbyhta9 =     0.00577510431414000;
zbyhta10 =     0.00160457304817000;
zbyhta11 =     0.00037341312478200;
zbyhta12 =    -0.00003424624729560;
zbyhta13 =     0.00167173445376000;
zbyhta14 =    -0.00025318033967200;
zbyhta15 =     0.00178443995452000;
zbyhta16 =     0.00033853103066900;
zbyhta17 =     0.00028446833904500;
zbyhta18 =     0.00019448266182100;
zbyhta19 =     0.00235601784297000;
zbyhta20 =     0.00059398660282600;
zbyhta21 =    -0.00024832909144700;
zbyhta22 =     0.00029149089482300;

// ZYHP
zfyhpa1 =     0.75000000000000000;
zfyhpa2 =     1.00752819544453400;

// ZBYHP
zbyhpa1 =     0.00078343366647300;
zbyhpa2 =     0.00138459071210000;
zbyhpa3 =     0.00074961824536100;
zbyhpa4 =     0.00074529971719200;
zbyhpa5 =    -0.00105229888020000;
zbyhpa6 =    -0.00055673890405200;
zbyhpa7 =     0.00113523987890000;
zbyhpa8 =    -0.00014068228924700;
zbyhpa9 =     0.00259014439019000;
zbyhpa10 =     0.00035834441381100;
zbyhpa11 =    -0.00089078028706200;
zbyhpa12 =    -0.00077093463407200;
zbyhpa13 =    -0.00366294234115000;
zbyhpa14 =     0.00061448019461600;
zbyhpa15 =     0.00246508562522000;
zbyhpa16 =     0.00045041845212200;
zbyhpa17 =    -0.00014012587929800;
zbyhpa18 =    -0.00003108856079700;
zbyhpa19 =     0.00203187061582000;
zbyhpa20 =    -0.00012337090958700;
zbyhpa21 =     0.00011601875463800;
zbyhpa22 =    -0.00008048707953920;

// HPRDTW
hprdtwa1 =     0.96303478353000000;

// ZYHST
zyhsta1 =     0.05000000000000000;

// ZYHTST
zyhtsta1 =     0.05000000000000000;

// ZYHPST
zyhpsta1 =     0.05000000000000000;

// YHSHR
yhshrs1 =     0.74226400000000000;

// YHTSHR
yhtshrs1 =     0.10953300000000000;

// YHPSHR
yhpshrs1 =     0.20280300000000000;

// PTR
ptra1 =     0.05117863942941380;
ptra2 =     0.00940211718360560;
ptra3 =     0.04990577450556040;
ptra4 =    -0.01048653111857970;
ptra5 =     0.01165361904939040;
ptra6 =    -0.06407319326233719;
ptra7 =     0.02534706938518230;
ptra8 =     0.01648869473929550;
ptra9 =     0.01058381008846820;
ptra10 =     0.90000000000000000;
ptra11 =     0.07374166005854491;
ptra12 =    -0.06702131852844340;
ptra13 =    -0.01238841393000660;
ptra14 =     0.01932690073569770;

// RRTR
rrtra1 =     0.97000000000000000;


 
model(linear);



// Original Model Code:


rfftay =  (1-tayr1)*(rstar+.25*(picxfe+picxfe(-1)+picxfl1(-1)+picxfl2(-1)))+ tayr1*rffe(-1) + tayp0*(picxfe    ) + tayp1*(picxfe(-1) ) + tayp2*(picxfl1(-1) ) + tayp3*(picxfl2(-1) ) + tayx0*xgap2+ tayx1*xgap2(-1);

xgap2   = 100*(xgdp - xgdpt);
lur     = lurnat + (lurb1 + lurb2) * (lur(-1) - lurnat(-1)) - lurb2 * (lurl1(-1) - lurnal1(-1)) + (lurb3 + lurb4) * xgap2 - lurb4 * xgap2(-1);
picnia  = 400*(pcniar - pcniar(-1) + dlpxp);
picxfe  = 400*(pcxfer - pcxfer(-1) + dlpxp);
pigdp   = 400*(pgdpr - pgdpr(-1) + dlpxp);
drffe   = rffe - rffe(-1);
ec      = ecs1*eco + ecs2*ech + (1-ecs1-ecs2)*kcd;
eco     = eco(-1) + ecoa1*(qec(-1)-eco(-1)) + ecoa2*(eco(-1)-ecol1(-1)) + datet*zeco + (1-datet)*lzeco(-1) + ecod1*(ecos1*(yhl-yhl(-1))+(1-ecos1)*(yht-yht(-1))-(datet*zeco+ (1-datet)*lzeco(-1)) ) + eco_;
qec     = qecs1 * (qecb1 + qecb2) * zyh +  qecs2 * (qecb3 - qecb2) * zyht +  qecs3 * qecb4 * zyhp  +  qecs4 * qecb5 * (qecs5 * (wpsnr-pcniar)+(1-qecs5) * wpo);
ecd     = ecd(-1) +  ecda1 * (qecd(-1)-ecd(-1)) +  ecda2 * (ecd(-1)-ecdl1(-1)) +  datet * zecd + (1-datet) * lzecd(-1) +  0.01 * ecdb2 * (datet * zgapc2 + (1-datet) * lzgapc2(-1)) +  ecd_;
qecd    = qec + qecds1 * hggdpt/400 + qecda3 * qecds2 * (rcar-zpi5);
echsh   = (1 + echa1) * echshl1 + echa2 * (echshl1 - echshl1(-1)) + echa3 * 0.01 * rrmet + ech_;
ech     =   echs1 * echsh + kh(-1);
eh      = eh(-1)   +  eha1*(qeh(-1)-eh(-1)) +  eha2*(eh(-1)-ehl1(-1)) +  eha3*(ehl1(-1)-ehl2(-1))+datet*zeh+ (1-datet)*lzeh(-1) + ehb2*(rme(-1)-rmel1(-1)) +  ehb3 * (rmel1(-1)-rmel2(-1))+  eh_;
qeh     = qec + qehs1 * hggdpt/400 + qeha2 * (qehs2 * ((1-taxon * .28) * rme-zpi10) - pcniar);
ecnia   = ecnias1 * eco + ecnias2 * ecd + ecnias3 * ech;
kcd     = kcd1*kcd(-1) + (1-kcd1)*ecd;
kh      = kh1*kh(-1) + (1-kh1)*eh;
epdc    = epdc(-1)+ epdcb1 * (qepdc(-1)-epdc(-1)) + epdcb2a * (epdc(-1)-epdcl1(-1))+ epdcb2b * (epdcl1(-1)-epdcl2(-1)) + epdcb3 * ( datet * (zxbc+zvpdc) + (1-datet) * (lzxbc(-1)+lzvpdc(-1)) )+ (1 - epdcb3) * (xb(-1)-xbl1(-1)) + epdc_;
qepdc   = xb + vpdc + qepdcs1 * 4 * hgx/400;
epdo    = epdo(-1) +  epdob1*(qepdol1(-1)-epdol1(-1)) +  epdob2a*(epdo(-1)-epdol1(-1)) +  epdob2b*(epdol1(-1)-epdol2(-1))+ epdob3*(datet*(zxbo(-1)+zvpdo(-1)) + (1-datet)*(lzxbol1(-1)+lzvpol1(-1))) + (1 - epdob3) * (xb(-1)-xbl1(-1)) + epdo_;
qepdo   = xb + vpdo + qepdos1 * 4 * hgx/400;
epd     = epds1 * epdo + (1-epds1) * epdc;
eps     = eps(-1)+ epsb1 * (qepsl1(-1)-epsl1(-1))+ epsb2a * (eps(-1)-epsl1(-1))+ epsb2b * (epsl1(-1)-epsl2(-1)) + epsb3 * (datet * (zxbs(-1)+zvps(-1))+ (1-datet) * (lzxbsl1(-1)+lzvpsl1(-1)) )+ (1 - epsb3) * (xb(-1)-xbl1(-1))+ eps_;
qeps    =  xb + vps + qepss1 * 4 * hgx/400;
ki      = ki(-1) +  kib1 * (ki(-1)-kil1(-1)) +  kib2 * (efs(-1)-efsl1(-1)) +  kib3 * (efsl1(-1)-efsl2(-1)) +  kib4 * (kil3(-1)-kil4(-1)-(efsl3(-1)-efsl4(-1)))+  ki_;
kpdc    = kpdcs1 * kpdc(-1) + (1-kpdcs1) * epdc;
kpdo    = kpdos1 * kpdo(-1) + (1-kpdos1) * epdo;
kps     = kpss1 * kps(-1) + (1-kpss1) * eps;
hks     = 400 * (hkss1 * (kpdc-kpdc(-1)) + hkss2 * (kpdo-kpdo(-1)) + hkss3 * (kps-kps(-1)) + (1 - hkss1 - hkss2 - hkss3) * (ki-ki(-1)) );
ks      = ks(-1) + hks/400;
rpd     = 0.5 * (1-0.35) * (rg5e+rbaae-rg10e) - 0.5 * (datet * zpib5+(1-datet) * lzpib5(-1))  + 0.5 * req;
rtpdc   = rtpdcs1 * rpd + rtpdcs2 * tapdda  + 0.5 * (ppdcr+ppdcr(-1)) - pxbr - dlpxp;
rtpdo   =  rtpdos1 * rpd + rtpdos2 * tapdda  - pxbr - dlpxp;
rtps    =  rtpss1 * rpd - pxbr - dlpxp;
tapdda  =  tapdda1 * (rpd + datet * zpib5 + (1-datet) * lzpib5(-1));
vpdc    =  - rtpdc;
vpdo    =  - rtpdo;
vps     =  - rtps;
ex      = ex(-1)+ exa1 * (ex(-1)-fgdp(-1))+ exa1 * 0.94039 * (fpxr(-1)-pcniar(-1))+ exd0 * (fgdp-fgdpl1(-1))/2 + exd1 * (fgdpl1(-1)-fgdpl5(-1))/4 + ex_;
emo     = emo(-1) + emoa1*(emo(-1)-.33 * yh(-1) - 0.67 * xgdp(-1))+ emoa1 * 0.743170 * (pmor(-1)-pgdpr(-1)) - emoa1 * 0.443955 * (pmor(-1)+fpxr(-1)-pcniar(-1)) + emod0 * (xgap2-xgap2(-1))/100 + emod1 * (xgap2(-1)-xgap2l1(-1))/100+ emod2 * xgap2l1(-1)/100+ emo_;
ceng    =  ceng(-1) + cenga1 * (ceng(-1)-xg(-1)-veoa(-1))+ cengb1a * (xg-xg(-1)) + cengb1b * (xg(-1)-xgl1(-1)) + cengb2 * (veoa-veoa(-1)) + cengb3a * (ceng(-1)-cengl1(-1))+ cengb3b * (cengl1(-1)-cengl2(-1))+ ceng_;
emp     =  emps1 * ceng + emps2 * xgpot;
em      =  ems1 * emo + ems2 * emp;
fcbnr   = fcbns1 * (pxbr+ex) + fcbns2 * (fcbns2a * (pmor+emo)  + fcbns2b * (poilr + pxbr + emp)) + fcbns3 * fynicnr + fcbns4 * fynilnr + fcbns5 * (pxgr + xgpot);
fnicnr  = fnicns1 * fnicns3 *fnicnr(-1) + 0.54 * fnicns1 * fpic/400 - 0.67 * fnicns1 * (fpxr - fpxr(-1) + fpic/400 - dlpxp)- 0.87 * fnicns2 * fcbnr;
fnilnr  = fnilns1 * fnilns3 * fnilnr(-1) + 0.32 * fnilns1 * pigdp/400  - 0.06 * fnilns1 * (fpxr - fpxr(-1) + fpic/400 - dlpxp) - 1.87 * fnilns2 * fcbnr;
fynicnr = fynicns1 * rfynic + fnicnr(-1);
fynilnr = fynilns1 * rfynil + fnilnr(-1);
rfynic  = rfynic(-1) + rfynica1 * (rfynic(-1) - rfynil(-1)) + rfynica2 * (rfynic(-1) - rfyncl1(-1))+ rfynica3 * (rfynil - rfynil(-1))+ rfynic_;
rfynil  = (1 + rfynila1) * rfynil(-1) + rfynila2 * rg10e(-1)  + rfynila3 * rtbe(-1) + rfynila4 * reqp(-1) + rfynila5 * (rfynil(-1) - rfynll1(-1)) + rfynila6 * (rg10e - rg10e(-1))+ rfynila7 * (rtbe - rtbe(-1)) + rfynila8 * (reqp - reqp(-1)) + rfynil_;
xgdp    = xgdps1*ecnia + xgdps2*eh + xgdps3*epd + xgdps4*eps + xgdps5*ki + xgdps6*ki(-1) + xgdps7*ex + xgdps8*em + xgdps9*(egfs1*egfo + egfs2*egfi + egfs3*egfl) + xgdps10*(egss1*egso + egss2*egsi + egss3*egsl);
xb      =  xbs1 * xgdp + xbs2 * (xgvs1 * egfl + (1 - xgvs1) * egsl) + xbs3 * xgpot;
xg      =  xgs1 * xb + (1 - xgs1) * emp;
xp      =  xps1 * ecnia + xps2 * eh + xps3 * epd + xps4 * eps + xps5 * ex  + xps6 * egfi + xps7 * egfo  + xps8 * egsi   + xps9 * egso;
efs     = efss1 * xb +  efss2 * ki + efss3 * ki(-1);
hmfpt   = hmfptrho*hmfpt(-1) + hmfpt_;
mfpt    =  mfpt(-1) + hmfpt/400 + mfpt_;
veoa    =  veoaa1*veoa(-1)- (1-veoaa1)*pcengr(-1);
xgpot   =  (.7*leppot + .265*ks + .035*veoa + mfpt)/(1-.035);
hgx     =  (.7*400*(leppot - leppot(-1)) + .265*hks  + .035 * 400 * (veoa - veoa(-1)) + hmfpt)/(1-.035);
empt    =  empt(-1) + empta1*(emp(-1)-empt(-1))+ hgemp/400;
hgemp   =  hgempa1 * hgemp(-1)+ (1 - hgempa1) * 400 * (emp - emp(-1));
xgvt    =  xgvt(-1)+ xgvta1 * (xgvs1 * egfl(-1) + (1-xgvs1) * egsl(-1)- xgvt(-1))+ hgxgv/400;
hgxgv   =  hgxgva1 * hgxgv(-1)  + (1 - hgxgva1) * 400 * (xgvs1 * (egfl - egfl(-1)) + (1- xgvs1) * (egsl - egsl(-1)));
xhhit   = xhhit(-1)+ xhhita1 * (xgpot(-1) - xhhit(-1)) + hgxhhi/400;
hgxhhi  =  hgxhhia1 * hgxhhi(-1)+ (1 - hgxhhia1) * 400 * (xgpot - xgpot(-1));
xgap    =  100*(xg-xgpot);
xgdpt   =  xgdpts1*xgpot  + xgdpts2*empt+ xgdpts3*xgpot + xgdpts4*xgvt+ xgdpts5*xhhit;
hggdpt  =  xgdpts1*hgx  + xgdpts2*hgemp   + xgdpts3*hgx      + xgdpts4*hgxgv             + xgdpts5*hgxhhi;
lhp     = lhp(-1) + lhpb2*(qlhp(-1)-lhp(-1)- (mfpt-mfpt(-1))/0.965) + lhpb3*(lhp(-1)-lhpl1(-1)) + lhpb4*(datet*zlhp+ (1-datet)*lzlhp(-1)) + lhpb5*(xg- xgl1 - (hgxl1/400 - leppol1 + leppol2) - (hmfpt - hmfptl1)/(400 * 0.965)) + lhpb6*(xgl1 - xgl1(-1) - (hgxl1(-1)/400 - leppol2 + leppol2(-1))-(hmfptl1-hmfptl1(-1))/(400 * 0.965)) + lhp_;
qlhp    =  xg - lprdt;
lww     =  (1 - lwwa2) * lww(-1) + lwwa3 * (lhp - 2*lhpl1 + lhpl1(-1));
lep     =  lhp - lww;
leo     =  leors1 * leh - leors2 * lep  - leors3 * (egsl - dglprd * lprdt)        - leors4 * (egfl - dglprd * lprdt);
leh     =  lehs1 * lfpr - lehs2 * lur;
lfpr    =  (1-lfpra1) * lfpr(-1) + lfpra2 * (lur(-1)-lurnat(-1));
left    =  lefta1 * left(-1) + (1-lefta1)*(egfl - dglprd * lprdt);
lest    = lesta1 * lest(-1)+ (1-lesta1)*(egsl - dglprd * lprdt);
leppot  =  leppots1 * lurnat + leppots2 * left   + leppots3 * lest;
lurnat  =  0*lurnat(-1);
lprdt   =  xgpot - leppot;
jccanr  =  jccans1*kh(-1) + jccans2*kpdc(-1)+ jccans3*kpdo(-1)+ jccans4*kps(-1)+ jccans5*(xgdpt+pgdpr);
qynidnr =  qynidb2 * (ynicpnr + qynids1 * (trfci + trsci));
rrmet   =  rrmeta1*rrmet(-1) + (1-rrmeta1)*(rme- zpi10);
tryh    =  tryhs1 * tfpnr + tryhs2 * tspnr+ tryhs3 * (yhl + pcniar)   + tryhs4 * (yhptns1*xbnr + yhptns2*yhinr + yhptns3*ynidnr);
wdnfcnr =  wdnfcnr(-1)+ wdnfca1 * (wdnfcnr(-1) + wdnfcs1 * yninr(-1) + wdnfcs2 * ynilnr(-1)) + wdnfcb1 * (wdnfcnr(-1)-wdnfcl1(-1))+ wdnfcb2 * (wdnfcl1(-1)-wdnfcl2(-1)) + wdnfcb3 * xgap2+ wdnfcnr_;
xbnr    =  xbns1 * xgdpnr + xbns2 * (plr  - dglprd * lprdt + xgvs1 * egfl + (1 - xgvs1) * egsl)+ xbns3 * (pxbr + xgpot);
xgdpnr  =  xgdpns1 * xp + xgdpns2 * (pxbr + ki) + xgdpns3 * (pxbr + ki(-1)) + xgdpns4 * (ems1 * (pmor + emo) + ems2 * (poilr + pxbr + emps1*ceng + emps2*xgpot))      + xgdpns5 * (plr - dglprd * lprdt + egfl)        + xgdpns6 * (plr - dglprd * lprdt + egsl);
xgnr    =  xgns1*xbnr + (1-xgns1)*(emps1*ceng + emps2*xgpot + poilr + pxbr);
ycsnr   =  ycsns1 * ynicpnr + ycsns2 * (trfci + trsci)+ ycsns3 * ynidnr+ ycsns4 * (jccacs1 * jccanr + jccacs2 * (pgdpr + xgdpt)+jccacs3 * kh(-1));
yh      =  yhs1*yhl + yhs2*yht  + yhs3*yhp;
yhibnr  =  yhibnr(-1)+ yhibna1 * ((picxfe + picxfe(-1) + picxfl1(-1) + picxfl2(-1))/1600 - (dlpxp + dlpxp(-1) + dlpxpl1(-1) + dlpxpl2(-1))/4 )+ yhibna3 * (ecnia(-1) + pcniar(-1) - yhibnr(-1))+ yhibna4 * (yhibnr(-1)-yhibnl1(-1)  - (picxfe(-1) + picxfl1(-1) + picxfl2(-1) + picxfl3(-1))/1600- (dlpxp(-1) + dlpxpl1(-1) + dlpxpl2(-1) + dlpxpl3(-1))/4 ) + yhibna6 * rcar(-1) + yhibna7 * (pcxfer(-1)-pcniar(-1)+ecd(-1)-ecnia(-1)) + yhibna8 * (rffe-rffe(-1))+ yhibnr_;
yhinr   =  yhins1*yniinr + yhins2*gfintnr + yhins3*gsintnr + yhins4*yhibnr;
yhl     =  ynilnr - pcniar + yhls1*tryh + yhls2 * trssi;
yhp     =    yhps1 * (yhptns1*xbnr + yhptns2*yhinr + yhptns3*ynidnr)+ yhps2 * (yhpntnr - tryh) - pcniar;
yhpntnr =   yhpntns1 * (pcxfer+kcd(-1)) + yhpntns2 * yhibnr   + yhpntns3 * (trfci + trsci) + yhpntns4 * ynicpnr   + yhpntns5 * ynidnr + yhpntns6 * zpi10+ yhpntns7 * gfdbtnr  + yhpntns8 * gsdbtnr;
yht     =  xgdpt + yhts1 * (gftrd + gstrd);
ynicpnr =  ynicpns1 * yninr+ ynicpns2 * ynilnr+ ynicpns3 * yniinr+ ynicpns4 * xbnr+ ynicpns5 * (ecnia + pcniar)+ ynicpns6 * trsib+ ynicpns7 * (pgdpr + xgdpt);
ynidnr  =  ynidnr(-1) + pxbr - pxbr(-1)+ ynidna1 * (qynidnr(-1)-ynidnr(-1))+ ynidna2 * (ynidnr(-1)-ynidnl1(-1)  - (pxbr(-1) - pxbrl1(-1))) + datet * zynid + (1-datet) * lzynid(-1)  + ynidnr_;
yniish  =  yniina1 * yniish(-1) - (yniina2*yniins3 + yniina3*yniins4 + yniina5*yniins5) * (yniins2*yninr(-1) + (1-yniins2)*ynilnr(-1)) + yniina2*yniins3*(kh(-1) + yniins6*rrmet) + yniina3*yniins4*(wdnfcnr(-1) + yniins7*rbaae) + yniina4*yniins4*(wdnfcnr(-1) - wdnfcl1(-1) + yniins7*(rbaae-rbaae(-1))- yniins2*(yninr(-1)-yninrl1(-1)) - (1-yniins2)*(ynilnr(-1)-ynilnl1(-1)))  + yniina5*(yniins8*fnicnr(-1)-(yniins8-yniins5)*fnicnr(-1)) + yniish_;
yniinr  =  yniins1 * yniish + yniins2 * yninr(-1) + (1 - yniins2) * ynilnr(-1);
ynilnr  =  plr + ynilns1 * lhp + ynilns2 * (egfl - dglprd * lprdt)+ ynilns3 * (egsl - dglprd * lprdt);
yninr   =  ynins1 * xgdpnr + ynins2 * fynicnr + ynins3 * fynilnr  + ynins4 * jccanr;
ypnr    =  ypns1 * (yhl + pcniar) + ypns2 * (yht + pcniar) + ypns3 * (yhptns1*xbnr + yhptns2*yhinr + yhptns3*ynidnr);
picnia/400 =   pcniaa1 * (qpcniar(-1)-pcniar(-1)) +  pcniaa21 * picnia(-1)/400 +  pcniaa22 * picnil1(-1)/400+  pcniaa23 * picnil2(-1)/400 +  pcniaa24 * picnil3(-1)/400 +   datet * lzpc +  (1-datet) * lzpc(-1) +  pcniaa4 * (datet * lzlurc +  (1-datet) * lzlurc(-1))+ pcniaa5 * pcnias1 * (pcer-pcer(-1))+  pcniaa6 * pimor +  pcniaa8 * pcnias2 * (pcfr-pcfr(-1)) +  pcniar_;
dpxnc   =   pxnca1 * (qpxncr(-1)-pxncr(-1)) +  pxnca21 * dpxnc(-1)+  pxnca22 * dpxncl1(-1) +  pxnca23 * dpxncl2(-1) +   datet * lzpnc  +  (1-datet) * lzpnc(-1)+  pxnca4 * (datet * lzlurnc +  (1-datet) * lzlurnc(-1)) +  pxnca5 * pcnias1 * (pcer-pcer(-1))+pxnca6 * pimor+  pxnca8 * pcnias2 * (pcfr-pcfr(-1)) +  pxncr_;
dlpl    =  pipla1 * (qplr(-1)-plr(-1)) + pipla21 * dlpl(-1)+ pipla22 * dlpll1(-1) + pipla23 * dlpll2(-1)+ datet * zpl + (1-datet) * lzpl(-1)+ pipla4 * (datet * zlurl + (1-datet) * lzlurl(-1)) + pipla6 * dlpls1 * (ynilnr - ynilnr(-1) + tsiens1 * (trssi - trssi(-1)) - (lhp - lhp(-1)) - dlpl) + pipla7 * dlpls2 * (plminr - plminr(-1)) +dlpl_;
plr    =  plr(-1) + dlpl - dlpxp;
qpcniar  = qpxpr + qpcnias1 * uqpca2 * (pmor - pcniar);
qpxncr = pxncr + qpxnca1 * qpxpr + qpxnca2 * (qpcniar-pcniar);
qplr   =  plr + pxgr - qpxgr;
qpxpr  =  qpxps1 * (qpxgr - pxgr);
qpxgr  =  plr - lprdt;
ppdcr  =  0*ppdcr(-1);
poilr  =  qsector * ((1 + poilra1) *poilr(-1) + poilrb1*(poilr(-1) - poilrl1(-1)))+ poilr_;
pcengr =  (1 + pcengrb2) *pcengr(-1) + pcengrb1 * (pcengr(-1) - pcengl1(-1))+ pcengrb3 * poilr(-1)+ pcengrb4 * (poilr-poilr(-1)) + pcengr_;
pcer   =  (1 + pcera1) * pcer(-1)+ pcera2 * (pcengr(-1) + pxbr(-1) - pcxfer(-1)) + pcerb1 * (pcer(-1) - pcerl1(-1))+ pcerb2 * (pcengr + pxbr - pcxfer - (pcengr(-1) + pxbr(-1) - pcxfer(-1)))+ pcerb3 * (poilr + pxbr - pcxfer - (poilr(-1) + pxbr(-1) - pcxfer(-1))) + pcer_;
pcfr   =   (1 + pcfra1) * pcfr(-1) + pcfrb11 * (pcfr(-1) - pcfrl1(-1)) + pcfrb12 * (pcfrl1(-1) - pcfrl2(-1)) + pcfrb13 * (pcfrl2(-1) - pcfrl3(-1)) + pcfr_;
uces   =  (1 + ucesa1) * uces(-1) + ucesa2 * pcer(-1) + ucesa3 * (ceng(-1) - xg(-1)) + ucesb1 * (uces(-1) - ucesl1(-1))+ ucesb2 * (pcer - pcer(-1)) + ucesb3 * (ceng - xg - (ceng(-1) - xg(-1))) + uces_;
ucfs   =  (1 + ucfsa1) * ucfs(-1)+ ucfsa2 * pcfr(-1)+ ucfsb1 * (ucfs(-1)-ucfsl1(-1))+ ucfsb3 * (pcfr-pcfr(-1))+ ucfs_;
pcxfer =  pcniar - pcnias1 * pcer - pcnias2 * pcfr;
pmor   =  pmor(-1) - dlpxp + pmoa1 * (0.64 * (pcniar(-1)-fpxr(-1)) + 0.36 * pxbr(-1) - pmor(-1))+ pmoa2 * (pcniar - pcniar(-1) - (fpxr - fpxr(-1)))+ pmoa3 * (pxbr - pxbr(-1) + dlpxp) + pmor_;
pgdpr  =  xgdpnr - xgdp;
plminr  = 0*one;
pxbr   =  xbnr - xb;
pxgr   =  xgnr - xg;
hgpmor  =  hgpmora1 * hgpmor(-1)+ (1-hgpmora1) * 400 * (pmor+dlpxp-pmor(-1)-dlpxp(-1));
pimor =  pimors1 * (pmor - pmor(-1) - hgpmor/400);
dpxnc   = pxncr - pxncr(-1) + dlpxp;
dlpxp   =  xps1 * picnia/400+ (1 - xps1) * dpxnc;

// Fiscal policy shock is implemented to affect all three components of federal government spending:
egfi    = gsector*(egfi(-1) +egfia1*(egfi(-1)-egt(-1)) +egfia2a*(egfi(-1)-egfil1(-1)) +egfia2b*(egfil1(-2)-egfil2(-1)) + egfia3*(egt-egt(-1))+egfia4a*xgap2+egfia4b*xgap2(-1)) + egfi_;
egfl    = gsector*(egfl(-1) +egfla1*(egfl(-1)-egt(-1)) +egfla2a*(egfl(-1)-egfll1(-1)) +egfla2b*(egfll1(-2)-egfll2(-1)) + egfla3*(egt-egt(-1))+egfla4a*xgap2+egfla4b*xgap2(-1)) + egfl_;
egfo    = gsector*(egfo(-1) +egfoa1*(egfo(-1)-egt(-1)) +egfoa2a*(egfo(-1)-egfol1(-1)) +egfoa2b*(egfol1(-2)-egfol2(-1)) + egfoa3*(egt-egt(-1))+egfoa4a*xgap2+egfoa4b*xgap2(-1)) ;%+ fispol

egsi    = gsector*(egsi(-1) +egsia1*(egsi(-1)-egt(-1)) +egsia2a*(egsi(-1)-egsil1(-1)) +egsia2b*(egsil1(-2)-egsil2(-1)) + egsia3*(egt-egt(-1))+egsia4a*xgap2+egsia4b*xgap2(-1)) + egsi_;
egsl    = gsector*(egsl(-1) +egsla1*(egsl(-1)-egt(-1)) +egsla2a*(egsl(-1)-egsll1(-1)) +egsla2b*(egsll1(-2)-egsll2(-1)) + egsla3*(egt-egt(-1))+egsla4a*xgap2+egsla4b*xgap2(-1)) + egsl_;
egso    = gsector*(egso(-1) +egsoa1*(egso(-1)-egt(-1)) +egsoa2a*(egso(-1)-egsol1(-1)) +egsoa2b*(egsol1(-2)-egsol2(-1)) + egsoa3*(egt-egt(-1))+egsoa4a*xgap2+egsoa4b*xgap2(-1)) + egso_;
egt     =  egt(-1) + egta1 * (egt(-1) - xgdpt(-1))+ (hggdpt + hggdpl1 + hggdpl2 + hggdpl2(-1))/1600;
eg      =  egs1*(egfs1*egfo + egfs2*egfi + egfs3*egfl)+ egs2*(egss1*egso + egss2*egsi + egss3*egsl);
gfdbtnr =  gfdbtns1 * gfdbtnr(-1)+ gfdbtns2 * gfexpnr + gfdbtns3 * gfrecnr + gfdbtns4 * egfi + gfdbtns5 * (pgdpr + xgdpt);
gfintnr = gfdbtnr(-1) + gfintns1 * rgfint;
gfrecnr =  gfrecns1 * tfpnr+ gfrecns2 * (trfcis1 * trfci + ynicpnr) + gfrecns3 * (ecnia + pcniar)+ gfrecns4 * ynilnr;
gfexpnr =  gfexpns1 * (plr - dglprd * lprdt + egfl) + gfexpns2 * egfo + gfexpns3 * (xgdpt + pcniar + gfptns1*gftrd) + gfexpns4 * gfintnr  + gfexpns5 * (pgdpr + xgdpt)  + gfexpns6 * (pgdpr + xgdpt);
gfdefnr =  gfdefns1 * gfexpnr + (1 - gfdefns1) * gfrecnr;
gftrd   =  gftrda1 *gftrd(-1) + gftrda2 * xgap2  + gftrd_;
gsdbtnr =  gsdbtns1 * gsdbtnr(-1) + gsdbtns2 * gsexpnr+ gsdbtns3 * gsrecnr   + gsdbtns4 * egsi+ gsdbtns5 * (pgdpr + xgdpt);
gsintnr =  gsintns1 * (gsdbtnr(-1) + gfintns1*rgfint)  + (1 - gsintns1) * xbnr;
gsrecnr =  gsrecns1 * tspnr + gsrecns2 * (trscis1*trsci + ynicpnr) + gsrecns3 * (trsibs1*trsib + ecnia + pcniar)+ gsrecns4 * (trssis1*trssi + ynilnr) + gsrecns5 * (pgdpr + xgdpt);
gsexpnr =  gsexpns1 * (plr - dglprd * lprdt + egsl) + gsexpns2 * egso + gsexpns3 * (xgdpt + gsptns1*gstrd + pcniar) + gsexpns4 * gsintnr + gsexpns5 * (pgdpr + xgpot);
gssrpnr =  gssrpns1 * gsrecnr + (1 - gssrpns1) * gsexpnr;
gstrd   =  gstrda1 * gstrd(-1) + gstrda2 * xgap2 + gstrd_;
rgfint  =  rgfint(-1) + rgfinta1 * (rg5e(-1)/100 - rgfint(-1)) + rgfinta2 * (rtbe(-1)/100 - rg5e(-1)/100) + rgfintb1 * (rgfint(-1) - rgfinl1(-1))+ rgfintb2 * (rg5e(-1) - rg5el1(-1))/100 + rgfintb3 * (rtbe - rtbe(-1))/100  + rgfint_;
tfpnr   =  tfpns4 * trfp+ tfpns1 * ypnr + tfpns2 * (xgdpt + gfptns1*gftrd + pcniar)+ tfpns3 * (xgdpt + gsptns1*gstrd + pcniar);
trfci   =  trfcia1 * trfci(-1) + trfcia4 * xgap2 + trfcia5 * picnia  + trfci_;
trfp    =  trfpt   + trfpa21 * (trfp(-1)-trfpt(-1))+ trfpa22 * (trfpl1(-1)-trfptl1(-1)) + trfpa3 * xgap2(-1)  + trfp_;
trfpt   =  trfpt(-1)+ trfpt1 * (gfdbtnr(-1)-xgdpnr(-1)) + trfpt2 * (gfdbtnr(-1)-xgdpnr(-1)- (gfdbtl1(-1)-xgdpnl1(-1)));
trsci   =  trscia1 *trsci(-1)  + trscia31 * xgap2+ trscia32 * xgap2(-1)+ trscia4 * (trfci-trfci(-1))+ trsci_;
trsib   =  trsiba1 * trsib(-1) + trsiba3 * xgap2  + trsib_;
trsp    =  trspa1 * trsp(-1)  + trspa21 * trspt + trspa22 * trspt(-1) +trspa3 * xgap2(-1) + trspa4 * (trfp-trfp(-1))+trsp_;
trspt   =  trspt(-1) + trspt1 * (gsdbtnr(-1)-xgdpnr(-1)) + trspt2 * (gsdbtnr(-1)-xgdpnr(-1)- (gsdbtl1(-1)-xgdpnl1(-1)));
trssi   =  trssia11 * trssi(-1)+ trssia12 * trssil1(-1)+ trssia3 * xgap2+ trssi_;
tspnr   =  tspns4 * trsp + tspns1 * ypnr + tspns2 * (xgdpt + gfptns1*gftrd + pcniar) + tspns3 * (xgdpt + gsptns1*gstrd + pcniar);
rffe    =  rfftay;
rtbe    =  rtbeb11 * rtbe(-1)+ rtbeb12 * rtbel1(-1)+ rtbeb21 * rffe + rtbeb22 * rffe(-1) + rtbe_;
rg5e    = zrff5 + rg5p;
rg5p    =  rg5pa1 * zgap05+ rg5pa2 * rg5p(-1) - rg5pa1 * rg5pa2 * zgap05(-1)  + rg5p_ ;
rg10e   = zrff10 + rg10p;
rg10p   =  rg10pa1 * zgap10  + rg10pa2 * rg10p(-1)- rg10pa1 * rg10pa2 * zgap10(-1) + rg10p_ ;
rbaae   = zrff10 + rbaap;
rbaap   =  rbaapa1 * zgap10 + rbaapa2 * rbaap(-1)- rbaapa1 * rbaapa2 * zgap10(-1) + rbaap_ ;
rcbe    = zrff30 + rcbp;
rcbp    =  rcbpa1 * zgap30+ rcbpa2 * rcbp(-1) - rcbpa1 * rcbpa2 * zgap30(-1)+ rcbp_ ;
rcar    =  rcara3 * rcar(-1)+ rcara41 * rg5e+ rcara42 * rg5e(-1) + rcar_;
rme     =  rmea1 * rme(-1) + rmea21 * rbaae + rmea22 * rbaae(-1) + rmea3 * (rbaae-rtbe) + rme_;
req     =  rcbe - zpic30 + reqp;
reqp    =  reqpb1 * reqp(-1) + reqp_;
wpsnr   =  ynidnr - wpsns1 * (req - zdivgr);
rcgain  = rcgaina1 * rcgain(-1) + rcgaina2 * (picxfe+picxfe(-1)+picxfl1(-1)+picxfl2(-1))/4 + rcgaina5 * xgap2 + rcgain_;
wpo     =  wpos1 * wpo(-1) + wpos2 * (ydns1*ypnr+ydns2*tfpnr+ydns3*tspnr-pcniar)+ wpos3 * ecnia + wpos4 * (yhibnr-pcniar)+ wpos5 * (ecd+pcxfer-pcniar)+ wpos6 * (kcd(-1)+pcxfer-pcniar) + wpos7 * (rcgain-picnia);
fxgap   =  fsector*( fxgapb1a*fxgap(-1)+ fxgapb1b*fxgapl1(-1) + fxgapb2*(frl10(-1)-fpi10t(-1))+ fxgapb3*xgap2(-1)+ fxgapb4*rstar(-1)+ fxgap_ );
fgdp    = fgdpt + fxgap/100;
fgdpt   =  fgdpt(-1) + fgdpta1 * (fgdpt(-1) - xgdpt(-1))   + (hggdpt + hggdpl1 + hggdpl2 + hggdpl2(-1))/1600;
fpi10   =  fsector*( (fpi10(-1)+fpi10l1(-1)+fpi10l2(-1)+fpi10l3(-1))/4 + fpi10b2*fxgap(-1)+ fpi10b3a*(poilr-poilr(-1))+ fpi10b3b*(poilr(-1)-poilrl1(-1))+ fpi10_);
fpi10t  = fsector*(fpi10tb1*fpi10t(-1) + fpi10tb2*fpi10);
fpic    =  fsector*(fpicb1*fpi10 + fpicb2*fpic(-1));
frs10   =  fsector*( frs10b1*rstar(-1)+ frs10b2*(fpi10+fpi10(-1)+fpi10l2+fpi10l2(-1))/4 + frs10b3*(fpi10+fpi10(-1)+fpi10l2+fpi10l2(-1))/4+ frs10b4*fxgap);
frl10   =  fsector * (frl10(-1)+ frl10b1 * (frl10(-1) - frs10(-1))+ frl10b2 * (frl10(-1) - frl10l1(-1))  + frl10b3 * (frs10 - frs10(-1))+ frl10b4 * (fxgap - fxgap(-1))+ frl10_);
fpxr    =  fpxrr + fpxrb1 * (rg10e-zpi10-(frl10-fpi10t)) + fpxrb2 * (fpxrs1 * (fnicnr - xgdpnr) + fpxrs2 * (fnilnr - xgdpnr));
fpxrr   =  (fprho - fpxrrb1) * fpxrr(-1)+ fpxrrb2 * (fpxrr(-1) - fpxrrl1(-1))+ fpxrr_;
rstar   =  rstar1*rstar(-1) + (1-rstar1)*(rffe-picnia);
zeco    =  zfecoa1*(qec-qec(-1))+ zfecoa2*(qec(+1)-qec) + zfecoa6*zeco(+1) + zfecoa7*lzeco(+1);
zecd    =   zfecda1*(qecd-qecd(-1)) + zfecda2*(qecd(+1)-qecd)    + zfecda6*zecd(+1)+ zfecda7*lzecd(+1);
zgapc2  =  zfgapc2a1*xgap2+ zfgapc2a6*zgapc2(+1)+ zfgapc2a7*lzgapc2(+1);
zeh     =    zfeha1*(qeh-qeh(-1)) + zfeha2*(qeh(+1)-qeh)  + zfeha3*(qehf1(+1)-qeh(+1))+ zfeha6*zeh(+1)+ zfeha7*lzeh(+1)+zfeha8*zehf2(+1);
zvpdc   = zfvpdc1*(vpdc-vpdc(-1))+ zfvpdc2*(vpdc(+1)-vpdc)+ zfvpdc3*(vpdcf1(+1)-vpdc(+1))+ zfvpdc6*zvpdc(+1)  + zfvpdc7*lzvpdc(+1)+ zfvpdc8*zvpdcf2(+1);
zvpdo   = zfvpdo1*(vpdo-vpdo(-1))+ zfvpdo2*(vpdo(+1)-vpdo)+ zfvpdo3*(vpdof1(+1)-vpdo(+1))+ zfvpdo6*zvpdo(+1)  + zfvpdo7*lzvpdo(+1)+ zfvpdo8*zvpdof2(+1);
zvps    = zfvps1*(vps-vps(-1))+ zfvps2*(vps(+1)-vps)+ zfvps3*(vpsf1(+1)-vps(+1))+ zfvps6*zvps(+1)  + zfvps7*lzvps(+1)+ zfvps8*zvpsf2(+1);
zxbc    = zfxbc1*(xb-xb(-1))+ zfxbc2*(xb(+1)-xb)+ zfxbc3*(xbf1(+1)-xb(+1))+ zfxbc6*zxbc(+1)  + zfxbc7*lzxbc(+1)+ zfxbc8*zxbcf2(+1);
zxbo    = zfxbo1*(xb-xb(-1))+ zfxbo2*(xb(+1)-xb)+ zfxbo3*(xbf1(+1)-xb(+1))+ zfxbo6*zxbo(+1)  + zfxbo7*lzxbo(+1)+ zfxbo8*zxbof2(+1);
zxbs    = zfxbs1*(xb-xb(-1))+ zfxbs2*(xb(+1)-xb)+ zfxbs3*(xbf1(+1)-xb(+1))+ zfxbs6*zxbs(+1)  + zfxbs7*lzxbs(+1)+ zfxbs8*zxbsf2(+1);
zpib5   = zfpib5a1*400*(pxbr - pxbr(-1) + dlpxp)+ zfpib5a6*zpib5(+1);
zpi5  =   zfpi5a1*picnia + zfpi5a6*zpi5(+1);
zpi10  =  zfpi10a1*picnia  + zfpi10a6*zpi10(+1);
zpic30 =  zfpic30a1*picnia + zfpic30a6*zpic30(+1);
zlhp   =  zflhpa1*(qlhp-qlhp(-1))+ zflhpa2*(qlhp(+1)-qlhp)  + zflhpa6*zlhp(+1) + zflhpa7*lzlhp(+1);
zyh    =  (1-(zfyha1^.25)*(zfyha2))*yh+ (zfyha1^.25)*zyh(+1);
zyhp   =  (1-(zfyhpa1^.25)*(zfyhpa2))*yhp+ (zfyhpa1^.25)*zyhp(+1);
zyht   =  (1-(zfyhta1^.25)*(zfyhta2))*yht+ (zfyhta1^.25)*zyht(+1);
zynid  =  zfynida1*(qynidnr-qynidnr(-1)-(pxgr-pxgr(-1))) + zfynida2*(qynidnr(+1)-qynidnr-(pxgr(+1)-pxgr))+ zfynida6*zynid(+1) + zfynida7*lzynid(+1);
zlurc  =  zflurca1*(lur-lurnat) + zflurca6*zlurc(+1);
zlurl  =  zflurla1*(lur-lurnat) + zflurla6*zlurl(+1);
zlurnc =  zflurnca1*(lur-lurnat)+ zflurnca6*zlurnc(+1);
zpc   =   zfpca1*(qpcniar-qpcniar(-1)+dlpxp)+ zfpca6*zpc(+1);
zpl    =  zfpla1*(qplr-qplr(-1)+dlpxp) + zfpla6*zpl(+1);
zpnc   =  zfpnca1*(qpxncr-qpxncr(-1)+dlpxp)+ zfpnca6*zpnc(+1);
zdivgr =  zfdivgra1*400*(ynidnr - pxgr - (ynidnr(-1) - pxgr(-1)))+ zfdivgra6*zdivgr(+1);
zrff5  =  zfrff5a1*rffe + zfrff5a6*zrff5(+1);
zrff10  =  zfrff10a1*rffe + zfrff10a6*zrff10(+1);
zrff30  =  zfrff30a1*rffe  + zfrff30a6*zrff30(+1);
zgap05  =  zfgap05a1*xgap + zfgap05a6*zgap05(+1);
zgap10 =  zfgap10a1*xgap + zfgap10a6*zgap10(+1);
zgap30 =  zfgap30a1*xgap + zfgap30a6*zgap30(+1);

ecol1 = eco(-1);
ecdl1   =ecd(-1);
echshl1    = echsh(-1);
ehl1    =eh(-1);
ehl2   = ehl1(-1);
qehf1 = qeh(+1);
epdcl1  =epdc(-1);
epdcl2  =  epdcl1(-1);
epdol1  =epdo(-1);
epdol2  =  epdol1(-1);
qepdol1 = qepdo(-1);
epsl1   =eps(-1);
epsl2  =  epsl1(-1);
qepsl1  =qeps(-1);
kil1    =ki(-1);
kil2    =kil1(-1);
kil3    =kil2(-1);
kil4  =  kil3(-1);
vpdcf1  =  vpdc(+1);
vpdof1  =  vpdo(+1);
vpsf1  =  vps(+1);
cengl1 = ceng(-1);
cengl2  = cengl1(-1);
rfyncl1 = rfynic(-1);
rfynll1 = rfynil(-1);
xbf1  =  xb(+1);
xbl1 = xb(-1);
xgl1 = xg(-1);
efsl1   =efs(-1);
efsl2   = efsl1(-1);
efsl3   = efsl2(-1);
efsl4   = efsl3(-1);
hmfptl1 = hmfpt(-1);
hgxl1 = hgx(-1);
xgdptl1 = xgdpt(-1);
xgdptl2  = xgdptl1(-1);
hggdpl1 = hggdpt(-1);
hggdpl2  = hggdpl1(-1);
xgap2l1  =xgap2(-1);
lhpl1 = lhp(-1);
qlhpl1   =qlhp(-1);
lurl1   =lur(-1);
leppol1   =leppot(-1);
leppol2   =leppol1(-1);
lurnal1   =lurnat(-1);
wdnfcl1    =wdnfcnr(-1);
wdnfcl2    =wdnfcl1(-1);
xgdpnl1 = xgdpnr(-1);
yhibnl1 =yhibnr(-1);
ynidnl1 = ynidnr(-1);
ynilnl1 = ynilnr(-1);
yninrl1 = yninr(-1);
picnil1  = picnia(-1);
picnil2  = picnil1(-1);
picnil3  = picnil2(-1);
dpxncl1 = dpxnc(-1);
dpxncl2 = dpxncl1(-1);
dlpll1 = dlpl(-1);
dlpll2 = dlpll1(-1);
poilrl1  =poilr(-1);
pcengl1=pcengr(-1);
pcerl1  = pcer(-1);
pcfrl1  = pcfr(-1);
pcfrl2  = pcfrl1(-1);
pcfrl3  = pcfrl2(-1);
ucesl1 = uces(-1);
ucfsl1 = ucfs(-1);
picxfl1  = picxfe(-1);
picxfl2  = picxfl1(-1);
picxfl3 = picxfl2(-1);
pxbrl1   =pxbr(-1);
dlpxpl1  = dlpxp(-1);
dlpxpl2    = dlpxpl1(-1);
dlpxpl3    = dlpxpl2(-1);
egfil1   = egfi(-1);
egfil2   = egfil1(-1);
egfll1   = egfl(-1);
egfll2   = egfll1(-1);
egfol1   = egfo(-1);
egfol2   = egfol1(-1);
egsil1   = egsi(-1);
egsil2   = egsil1(-1);
egsll1   = egsl(-1);
egsll2   = egsll1(-1);
egsol1   = egso(-1);
egsol2   = egsol1(-1);
gfdbtl1 = gfdbtnr(-1);
gsdbtl1 = gsdbtnr(-1);
rgfinl1 = rgfint(-1);
trfpl1 = trfp(-1);
trfptl1 = trfpt(-1);
trssil1 = trssi(-1);
rtbel1 = rtbe(-1);
rg5el1 = rg5e(-1);
rmel1  =  rme(-1);
rmel2   =rmel1(-1);
fxgapl1 = fxgap(-1);
fgdpl1  =fgdp(-1);
fgdpl2  = fgdpl1(-1);
fgdpl3  =fgdpl2(-1);
fgdpl4  = fgdpl3(-1);
fgdpl5  =fgdpl4(-1);
fpi10l1  =fpi10(-1);
fpi10l2 = fpi10l1(-1);
fpi10l3  =fpi10l2(-1);
frl10l1  =frl10(-1);
fpxrrl1  =fpxrr(-1);
lzeco  =   zeco(+1);
lzecd  =  zecd(+1);
lzgapc2 = zgapc2(+1);
lzeh    = zeh(+1);
zehf2 = lzeh(+1);
lzvpdc  =  zvpdc(+1);
zvpdcf2  = lzvpdc(+1);
lzvpdo  =  zvpdo(+1);
zvpdof2  = lzvpdo(+1);
lzvpol1 = lzvpdo(-1);
lzvps  =  zvps(+1);
zvpsf2  = lzvps(+1);
lzvpsl1  =lzvps(-1);
lzxbc   =  zxbc(+1);
zxbcf2   = lzxbc(+1);
lzxbo   =  zxbo(+1);
zxbof2   = lzxbo(+1);
lzxbol1 = lzxbo(-1);
lzxbs   =  zxbs(+1);
zxbsf2   = lzxbs(+1);
lzxbsl1 =lzxbs(-1);
lzpib5 =  zpib5(+1);
lzlhp   = zlhp(+1);
lzynid  =  zynid(+1);
lzlurc =  zlurc(+1);
lzlurl =  zlurl(+1);
lzlurnc =  zlurnc(+1);
lzpc  =  zpc(+1);
lzpl   =  zpl(+1);
lzpnc  =  zpnc(+1);
one    =   one(-1);
end;

shocks;
var ceng_ ; stderr 0.01;
var dlpl_ ; stderr 0.01;
var ecd_ ; stderr 0.01;
var ech_ ; stderr 0.01;
var eco_ ; stderr 0.01;
var egsi_ ; stderr 0.01;
var egsl_ ; stderr 0.01;
var egso_ ; stderr 0.01;
var eh_ ; stderr 0.01;
var emo_ ; stderr 0.01;
var epdc_ ; stderr 0.01;
var epdo_ ; stderr 0.01;
var eps_ ; stderr 0.01;
var ex_ ; stderr 0.01;

var fpi10_      ; stderr 0.01;
var fpxrr_ ; stderr 0.01;
var frl10_ ; stderr 0.01;
var fxgap_ ; stderr 0.01;
var gftrd_ ; stderr 0.01;
var gstrd_ ; stderr 0.01;
var hmfpt_ ; stderr 0.01;

var ki_ ; stderr 0.01;
var lhp_ ; stderr 0.01;
var mfpt_  ; stderr 0.01;
var pcengr_ ; stderr 0.01;
var pcer_ ; stderr 0.01;
var pcfr_ ; stderr 0.01;
var pcniar_ ; stderr 0.01;
var pmor_ ; stderr 0.01;
var poilr_ ; stderr 0.01;
var pxncr_ ; stderr 0.01;

end;

%load('Sigma_US_FRB08.mat');

%M_.Sigma_e = Sigma_e_;
%options_.Schur_vec_tol = 1e-6;
stoch_simul (AR=100,IRF=0, noprint,nograph);