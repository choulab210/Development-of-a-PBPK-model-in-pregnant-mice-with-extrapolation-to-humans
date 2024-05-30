GMicePBPK.code <- '
$PROB
# Gestational F-53B PBPK model for female C57BL/6J mice
- Author    : Wei-Chun Chou, Jing Zhang
- Date      : Aug, 2023
- Structure : GI tract, Plasma, Liver, Fat, Mammary gland, Kidney, Filtrate, PTC, rest of Fetuss, fetal liver, fetal brain, amniotic fluid  
- Note.1    : Initial physiological parameters and optimized parameters ares matched with the values in non-pregnnat PBPK models
- Note.2    : Growth equation of physioloical parameters values were taken from O’Flaherty et al. (1992), O’Flaherty et al. (1995)

$PARAM @annotated
// #+ mice parameters
BW                  : 0.025   : kg,                  Bodyweight (measrument data if available)     ; Value obtained from TK study
Htc                 : 0.48    : Unitless,            Hematocrit for mice                           ; Value obtained from Hejtmancik et al., 2002
QCC                 : 16.5    : L/h/kg^0.75,         Cardiac output                                ; Value obtained from Brown et al., 1997
QLC                 : 0.161   : Unitless,            Fraction blood flow to liver                  ; Value obtained from Brown et al., 1997
QKC                 : 0.091	  : Unitless,            Fraction blood flow to kidney                 ; Value obtained from Brown et al., 1997
QMC                 : 0.002   : Unitless,            Fraction blood flow to Mammary gland          ; Value obtained from Hanwell and Linzell, 1973
QFC                 : 0.07    : Unitless,            Fraction blood flow to Fat                    ; Value obtained from Brown et al., 1997
VLC                 : 0.055   : Unitless,            Fraction of liver volume                      ; Value obtained from Brown et al., 1997
VKC                 : 0.017   : Unitless,            Fraction of kidney volume                     ; Value obtained from Brown et al., 1997
VMC                 : 0.01    : Unitless,            Fraction of Mammary gland tissue              ; Value obtained from Hanwell and Linzell, 1973
VFC                 : 0.068   : Unitless,            Fraction of fat tissue                        ; Value obtained from Brown et al., 1997
VPlasC              : 0.049   : L/kg BW,             Fraction of plasma volume                     ; Value obtained from Brown et al., 1997
VFilC               : 0.0017  : L/kg BW,             Fraction of filtrate (10% of Kidney volume)   ; Value obtained from Worley et al., 2017
VPTCC               : 1.35e-4 : L/g kidney,          Volume of proximal tubule cells               ; Value calculated from Hsu et al., 2014 (60 million PTC cells/gram kidney, 1 PTC = 2250 um3)
protein             : 2.0e-6  : mg protein/PTCs,     Amount of protein in proximal tubule cells    ; Value obtained from Addis et al., 1936
PL                  : 2.10    : Unitless,            Liver-to-plasma partition coefficient         ; Value obtained from TK study
PK                  : 0.80    : Unitless,            Kidney-to-plasma partition coefficient        ; Value obtained from Loccisano et al., 2012
PM                  : 0.16    : Unitless,            Mammary gland-to-plasma partition coefficient ; Value obtained from Loccisano et al., 2012
PF                  : 0.13    : Unitless,            Fat-to-plasma partition coefficient           ; Value obtained from Loccisano et al., 2012
PRest               : 0.22    : Unitless,            Rest of body partition coefficient            ; Value was optimized by Chou and Lin, 2019 
PPla                : 0.59    : Unitless,            placenta-to-plasmapartition coefficient       ; VValue obtained from TK study
MW                  : 570.67  : g/mol,               F-53B molecular mass                          
Free                : 0.009   : Unitless,            Free fraction of F-53B in maternal plasma     ; Value obtained from TK study
KbileC              : 3.90e-4 : 1/(h*BW^-0.25),      Biliary elimination rate                      ; Value was optimized by Chou and Lin, 2019 
KurineC             : 0.015   : 1/(h*BW^-0.25),      Urinary elimination rate                      ; Value obtained from TK study
K0C                 : 1.000   : 1/(h*BW^-0.25),      Rate of absorption of F-53B in stomach        ; Value was optimized by Chou and Lin, 2019
KabsC               : 2.43    : 1/(h*BW^-0.25),      Rate of absorption of F-53B in small intestines; Value obtained from TK study
KunabsC             : 5.40e-4 : 1/(h*BW^-0.25),      Rate of unabsorbed dose to appear in feces    ; Value obtained from TK study
GFRC                : 59      : L/hr/kg kiney,       Glomerular filtration rate (female)           ; Value obtained from Corley, 2005
GEC                 : 0.54    : 1/(h*BW^0.25),       Gastric emptying constant                     ; Value obtained from Yang et al., 2013
Vmax_baso_invitro   : 393.45  : pmol/mg Protein/min, Vmax of basolateral transporter               ; Value obtained from Worley and Fisher, 2015; assumed to be same as PFOA
Km_baso             : 27.2    : mg/L,                Km of basolateral transporter                 ; Value calculated from Worley et al., 2007
Vmax_apical_invitro : 4185    : pmol/mg protein/min, Vmax of apical transporter                    ; Value optimzed from Chou and Lin, 2019 
Km_apical           : 52.3    : mg/L,                Km of apical transporter                      ; Value optimzed from Chou and Lin, 2019
RAFbaso             : 3.99    : Unitless             Relative activity factor                      ; Value optimzed from Chou and Lin, 2019 
RAFapi              : 2.81    : Unitless             Relative acitivty factor                      ; Value optimzed from Chou and Lin, 2019 
Kdif                : 4.6e-5  : L/h,                 Diffusion rate from PTCs to kidney serum      ; Value optimzed from Chou and Lin, 2019 
KeffluxC            : 5.60    : 1/(h*BW^-0.25),      Rate of clearance of F-53B from PTCs into blood; Value optimzed from Chou and Lin, 2019  
Ktrans1C            : 1.27    : L/h/kg^0.75,         Mother-to-fetus placental transfer rate       ; Value from Loccisano et al., 2013; assumed to be same as PFOA
Ktrans2C            : 1       : L/h/kg^0.75,         Fetus-to-mother placental transfer rate       ; Value from Loccisano et al., 2013; assumed to be same as PFOA
Ktrans3C            : 0.23    : L/h/kg^0.75,         Fetus-to-amniotic fluid transfer rate         ; Value from Loccisano et al., 2012; assumed to be same as PFOA 
Ktrans4C            : 0.001   : L/h/kg^0.75,         Amniotic fluid-to-fetus transfer rate         ; Value from Loccisano et al., 2012; assumed to be same as PFOA 

// #+ Fetal parameters
N                   : 8       : number,              Number of fetus                               ; Value from TK study
VPlasC_Fet          : 0.049   : Unitless,            Fraction volume of fetal plasma               ; Value was assumed to be same as the mother
Free_Fet            : 0.009   : Unitless,            Free fraction of F-53B in plasma for Fetus     ; Value was assumed to be same as the mother
PL_Fet              : 2.10    : Unitless,            Liver-to-plasma partition coefficient in Fetus; Value was assumed to be same as the mother
PRest_Fet           : 0.22    : Unitless,            Rest of body-to-plasma Partition coefficient  ; Value was assumed to be same as the mother
QLC_Fet             : 0.161   : Unitless,            Blood flows to fetal liver                    ; Value from Itskovitz, 1987 and Loccisano et al. (2012) (Sheep); 
PB_Fet              : 1.55    : Unitless,            Brain-to-plasma partition coefficient in Fetus; Value was assumed to TK study 
QBC_Fet             : 0.1055  : Unitless,            Blood flows to fetal brain                    ; Carter, A. M. and W. Gu (1988). 

$MAIN
// #+ Time varabiles: Gestational day and age
// #+ GD             : Day, Gestational days
// #+ GA             : Week, Gestational age

double GD            = TIME/ 24;          
double GA            = TIME/ 168;          

// #+ Placenta volume (L); all equations from O’Flaherty et al.(1992), O.Clarke et al.(1993), adjusted from rat to mice

double VPla          = 1e-7;

if (GD <= 6)         {VPla = 0;} 
else if (GD <= 9.25) {VPla = (N*(8*(GD - 6)))/(1.0e6);} 
else if (GD > 9.25)  {VPla = ((N*32*exp(-0.23*(GD-9.25)))  + (N*40*(exp(0.28*(GD-9.25)) - 1))) /(1.0e6);}

// #+ Blood flow to mice placenta (L/h); equations from O’Flaherty et al. (1992), O.Clarke et al. (1993), adjusted from rat to mice
// #+ Qdec1          : L/d, Plasma flow to the placenta during GD6 - 9.25; yolk sac placenta predominates
// #+ Qdec2          : L/d, Plasma flow to the placenta during GD9.25 -11; yolk sac placenta disappears
// #+ Qcap           : L/d, Chorioallantoic placenta appears
// #+ QPla1          : L/h, 25% of incr. in maternal blood flow assoc. with
// #+ QPla           : L/h, Plasma flow to placenta

double Qdec1         = 0;
double Qdec2         = 0;
double Qcap          = 0;

if (GD <= 5.5)       {QPla1 = 1e-7;} 
else if (GD <= 9.25) {Qdec1 = 0.55*(GD - 5.5);}                                 
else if (GD <= 11)   {Qdec2 = (2.2*exp(-0.23*(GD-9.25)));}                         
else if (GD > 11)    {Qcap  = pow((0.1207*(GD-11)), 4.36);}   

double Qdec          = Qdec1 + Qdec2;
double QPla1         = (N*(0.02*Qdec + Qcap))/24;    
double QPla          = QPla1*(1-Htc); 


// #+ Changes of tissue volumes during pregnancy
// #+ VM             : L, Volume of mammary gland at GD0
// #+ VM_P           : L, Volume of mammary gland during pregnancy; Equation from Yoon et al. (2009), Loccisano et al. (2012) and Lin et al. (2013), adjusted from rat to mice
// #+ VF             : L, Volume of fat at GD0
// #+ VF_P           : L, Volume of fat gland during pregnancy; Equation from Yoon et al. (2009), Loccisano et al. (2012) and Lin et al. (2013), adjusted from rat to mice
// #+ VL             : L, Volume of liver
// #+ VK             : L, Volume of kidney
// #+ VKb            : L, Volume of kidney serum
// #+ VPlas          : L, Volume of plasma
// #+ VFil           : L, Volume of flitrate

double VM            = VMC*BW;  
double VM_P          = VM*(1 + 0.27*GD*(30.0 / 300)); 
double VF            = VFC*BW; 
double VF_P          = VF*(1 + 0.0165*GD*(30.0 / 300)); 
double VL            = VLC*BW;
double VK            = VKC*BW;
double VPlas         = VPlasC*BW;

// #+ Volumes of kidney and related compartmnet
double MK            = VKC*BW*1000;                                                          
double ML            = VLC*BW*1000;                     
double VPTC          = MK*VPTCC;                       
double VKb           = VK*0.16;          
double VFil          = VFilC*BW; 

// #+ Growth equations of tissues for the fetus druing pregnancy
// #+ VFet_1,        : kg, Fetus volume for one fetus; Equlation from O’Flaherty et al. (1995)
// #+ VFet,          : kg, Fetus volume for whole litter; Equlation from O’Flaherty et al. (1995)
// #+ VAmX,          : L,  Amniotic fluid volume for one fetus; Equlation from Loccisano et al. (2012), adjusted from rat to mice
// #+ VAm,           : L,  Amniotic fluid volume for whole litter; Equlation from Loccisano et al. (2012), adjusted from rat to mice
// #+ VPlas_Fet,     : L,  Plasma volume for fetus; Equlation from Loccisano et al. (2012)
// #+ VL_Fet,        : L,  Liver volume for fetus; Equlation from PK study data
// #+ VRest_Fet,     : L,  Rest of body volume for fetus; Equlation from Loccisano et al. (2012) 
// #+ VB_Fet,        : L,  Brain volume for fetus; Equlation from Sikov, M. R. and J. M. Thomas (1970). (2012) 
// #+ VBal_Fet,      : L,  Check the balance for fetus volume

double VFet_1        = 0;

if (GD <= 8.6)       {VFet_1 = pow((0.12*GD),4.53)/(1.0e6);} 
else if (GD <= 15.8) {VFet_1 = pow((0.12*8.6),4.53)/(1.0e6) + pow((1.2*(GD-8.6)),2.6)/(1.0e6);} 
else if (GD > 15.8)  {VFet_1 = (273.38+((GD-15.8)/(19-15.8))*((1.25*1.0e3)-273.38))/(1.0e6);}

double VFet          = VFet_1*N;
double VAmX          = GD >= 8 ? ((-4E-6)*pow(GD, 3) + 0.0002*pow(GD, 2) - 0.0023*GD + 0.0099)*(30.0 / 300) : 1E-7;
double VAm           = VAmX*N;
double VPlas_Fet     = VPlasC_Fet*VFet;                      
double VL_Fet        = (0.406/(1 + exp((14.716 - GD)/0.907))/1000);
double VRest_Fet     = 0.93*VFet - VPlas_Fet - VL_Fet - VB_Fet;                
double VB_Fet        = (4.191*exp(-exp(2.554-0.06726*GD))/1000)*N*(30.0 / 300);
double VBal_Fet      = 0.93*VFet - (VRest_Fet + VPlas_Fet + VL_Fet + VB_Fet);    

// #+ Fetal blood flows for fetus druing pregnancy
// #+ QC_Fet,        : L/h, Plasma flow to the fetus
// #+ QL_Fet,        : L/h, Plasma flow to the fetal liver
// #+ QRest_Fet,     : L/h, Plasma flow to the fetal rest of body
// #+ QFetBal,       : L/h, Check the balance of fetal blood flow

double QC_Fet        = (QPla/(1 + 20000*exp(-0.55*GD)));     
double QL_Fet        = QC_Fet*QLC_Fet;
double QB_Fet        = QC_Fet*QBC_Fet;
double QRest_Fet     = QC_Fet - QL_Fet- QB_Fet; 
double QFetBal       = QC_Fet - (QRest_Fet + QL_Fet + QB_Fet);

// #+ Changes of maternal body weight during pregnancy
// #+ BW_P           : Mice BW during pregnancy
// #+ BWinc          : BW increased during pregnancy
// #+ VRest          : Volume of rest of body in Dam
// #+ VBal           : Volume Balance check 

double BW_P          = BW + (VF_P - VF) + (VM_P - VM) + VPla + VFet + VAm; 
double BWinc         = (VF_P - VF) + (VM_P - VM) + VPla + VFet + VAm; 
double VRest         = (0.93*BW_P) - (VL + VK + VM_P + VF_P + VPlas + VPla + VFet + VAm); 
double VBal          = (0.93*BW_P) - VRest -(VL + VK + VM_P + VF_P + VPlas + VPla + VFet + VAm); 

// #+ Equations for dam physiological parameters; growth equation obtained from Loccisano et al. (2012)
// #+ QC_P           : L/h, Cardiac output (adjusted for plasma) 
// #+ QC             : L/h, Cardiac output at GD0 (adjusted for plasma)
// #+ QF             : L/h, Plasma flow to fat at GD0
// #+ QF_P           : L/h, Plasma flow to fat during pregnacy
// #+ QM             : L/h, Plasma flow to mammary gland at GD0
// #+ QM_P           : L/h, Plasma flow to mammary gland during pregnacy
// #+ QL             : L/h, Plasma flow to liver
// #+ QK             : L/h, Plasma flow to kidney
// #+ QRest          : L/h, Plasma flow to the rest of body
// #+ QBal           : L/h, Blood flow balance check

double QC_P          = QC+((N*(Qdec+Qcap))/24)*(1-Htc);                        
double QC            = QCC*pow(BW,(0.75))*(1-Htc);                               
double QF            = QFC*QC;                                      
double QF_P          = QF*(VF_P/VF);                              
double QM            = QMC*QC;                                      
double QM_P          = QM*(VM_P/VM);                              
double QL            = QLC*QC;                                      
double QK            = QKC*QC;
double QRest         = QC_P - (QK + QL + QM_P + QF_P + QPla);                
double QBal          = QC_P - (QK + QL + QRest + QPla + QM_P + QF_P);

// #+ Pregnant mice Kinetics parameters
// #+ GFR            : L/h, Glomerular Filtration Rate (GFR); GFR was assumed 50% of renal plasma flow
// #+ GE             : 1/h, Gastric emptying rate
// #+ K0             : 1/h, Rate of uptake from the stomach into the liver
// #+ Kbile          : 1/h, Biliary elimination rate constant, liver to feces storage
// #+ Kurine         : 1/h, Urinary elimination rate constant
// #+ Kabs           : 1/h, Rate constant of absorption of F-53B from small intestine to liver
// #+ Kunabs         : 1/h, Rate constant of unabsorbed dose to appear in feces
// #+ Kinetics parmaters for kidney and its subcompartment
// #+ PTC            : cells/kg BW, Number of PTC (cells/kg BW) (based on 60 million PTC/gram kidney)
// #+ MPTC           : g, mass of the proximal tubule cells (assuming density 1 kg/L) 
// #+ Vmax_basoC     : mg/h/kg BW^0.75, Vmax of basolateral transporters (average Oat1 and Oat3)
// #+ Vmax_baso      : 
// #+ Vmax_apicalC   : mg/h/kg BW^0.75, Vmax of apical transporters in in vitro studies (Oatp1a1)
// #+ Vmax_apical    : 
// #+ Kefflux        : 1/h, Efflux clearance rate from PTC to blood 
// #+ Ktrans_1       : L/h, Rate constant for placental transfer; dam to fetus
// #+ Ktrans_2       : L/h, Rate constant for placental transfer; fetus to dam
// #+ Ktrans_3       : L/h, Amniotic fluid transfer rate; fetus to fluid
// #+ Ktrans_4       : L/h, Amniotic fluid transfer rate; fluid to fetus

double GFR           = GFRC*(MK/1000);
double GE            = GEC*pow(BW_P,(-0.25));                
double K0            = K0C*pow(BW_P,(-0.25));               
double Kbile         = KbileC*pow(BW_P,(-0.25));                   
double Kurine        = KurineC*pow(BW_P,(-0.25));                 
double Kabs          = KabsC*pow(BW_P,(-0.25));                    
double Kunabs        = KunabsC*pow(BW_P,(-0.25));                
double PTC           = VKC*6e7*1000;                                       
double MPTC          = VPTC*1000;                                         
double Vmax_basoC    = (Vmax_baso_invitro*RAFbaso*PTC*protein*60*(MW/1e12)*1000);       
double Vmax_baso     = Vmax_basoC*pow(BW_P,0.75);                          
double Vmax_apicalC  = (Vmax_apical_invitro*RAFapi*PTC*protein*60*(MW/1e12)*1000);   
double Vmax_apical   = Vmax_apicalC*pow(BW_P,0.75);                      
double Kefflux       = KeffluxC*pow(BW_P,(-0.25));                        
double Ktrans_1      = Ktrans1C*(pow(VFet_1,0.75)*N);                       
double Ktrans_2      = Ktrans2C*(pow(VFet_1,0.75)*N);                      
double Ktrans_3      = Ktrans3C*(pow(VFet_1,0.75)*N);                       
double Ktrans_4      = Ktrans4C*(pow(VFet_1,0.75)*N);                      
  
// #+ Mass balance adjusted factor; avoding the negative occur at time = 0
double KDOSE        = (TIME==0)?0:1;  

$INIT @annotated
// #+ Set up the initial concentration; 
ADOSE              :0.002: mg, Amount of input dose; assumed a virtual compartment for validating the model mass balance
APlas_free         : 0   : mg, Amount of free F-53B in the plasma compartment
APTC               : 0   : mg, Amount of F-53B in the proximal tubule cells subcompartment
AFil               : 0   : mg, Amount of F-53B in the filtrate subcompartment
Aurine             : 0   : mg, Amount of F-53B in the urine virtual compartment
AKb                : 0   : mg, Amount of F-53B in the kidney blood compartment
ARest              : 0   : mg, Amount of F-53B in the rest of body compartment
Afeces             : 0   : mg, Amount of F-53B in the feces virtual compartment 
AL                 : 0   : mg, Amount of F-53B in the liver virtual compartment
AM                 : 0   : mg, Amount of F-53B in the mammary gland
AF                 : 0   : mg, Amount of F-53B in the fat compartment
A_baso             : 0   : mg, Amount of F-53B in the baso subcompartment 
A_apical           : 0   : mg, Amount of F-53B in the apical subcompartment 
Adif               : 0   : mg, Amount of F-53B in the diffusion virtual compartment 
Aefflux            : 0   : mg, Amount of F-53B in the efflux virtual compartment; simulation of F-53B from PTC to blood
Atrans_1           : 0   : mg, Amount of F-53B in the placental transfer from dam to Fetus  
Atrans_2           : 0   : mg, Amount of F-53B in the placental transfer from Fetus to dam 
Atrans_3           : 0   : mg, Amount of F-53B in the Amniotic fluid transfer from Amniotic fluid to fetus  
Atrans_4           : 0   : mg, Amount of F-53B in the Amniotic fluid transfer from fetus to Amniotic fluid 
APla               : 0   : mg, Amount of F-53B in the placenta compartment
ASI                : 0   : mg, Amount of F-53B in the small intestine compartment
AST                : 0   : mg, Amount of F-53B in the stomach compartment
AabsST             : 0   : mg, Amount of absorbed F-53B in the stomach compartment
AabsSI             : 0   : mg, Amount of absorbed F-53B in the small intestine compartment
AAm                : 0   : mg, Amount of F-53B in the Amniotic fluid compartment
AL_Fet             : 0   : mg, Amount of F-53B in the Fetal liver compartment
AB_Fet             : 0   : mg, Amount of F-53B in the Fetal brain compartment
ARest_Fet          : 0   : mg, Amount of F-53B in the Fetal rest of Fetus body compartment
APlas_Fet_free     : 0   : mg, Amount of F-53B in the Fetal plasma compartment
AUC_CPlas          : 0   : mg/L*hr, Area under curve of F-53B in maternal plasma
AUC_CPlas_Fet      : 0   : mg/L*hr, Area under curve of F-53B in fetal plasma
AUC_CL             : 0   : mg/L*hr, Area under curve of F-53B in maternal liver
AUC_CPla           : 0   : mg/L*hr, Area under curve of F-53B in placenta
AUC_CL_Fet         : 0   : mg/L*hr, Area under curve of F-53B in fetal liver
AUC_CB_Fet         : 0   : mg/L*hr, Area under curve of F-53B in fetal brain

$ODE 
// #+ Concentrations in the tissues and in the venous plasma leaving each of the tissues (Unit: mg/L) 
// #+ Concentrations for dam; CX indicate the concentration in the tissue (e.g., CL); CVX represent the concentration of chemical in tissue leaving tissue 
// #+ CPlas_free   : mg/L, Free F-53B concentration in the plasma
// #+ CL           : mg/L, Concentration of F-53B in the liver compartment
// #+ CKb          : mg/L, Concetraitons of F-53B in venous plasma leaving kidney
// #+ CK           : mg/L, Concetraitons of F-53B in Kidney compartment
// #+ CM           : mg/L, Concentration of F-53B in the mammary gland compartment
// #+ CF           : mg/L, Concentration of F-53B in the fat compartment
// #+ CAm          : mg/L, Concentration of F-53B in the amniotic fluid compartment
// #+ CRest        : mg/L, Concentration of F-53B in the rest of the body
// #+ CPTC         : mg/L, Concetraitons of F-53B in proximal tubule cells (PTC)
// #+ CFil         : mg/L, Concetraitons of F-53B in filtrate (Fil)
// #+ CPla         : mg/L, Concetraitons of F-53B in placenta

double CPlas_free  = APlas_free/ VPlas;    
double CPlas       = CPlas_free / Free;
double CL          = AL/VL;                                             
double CKb         = AKb/VKb;                                            
double CK          = CVK*PK;                                            
double CM          = AM/VM;                                              
double CF          = AF/VF;                                               
double CAm         = AAm/(VAm + 1E-7);
double CRest       = ARest/VRest;                                     
double CPTC        = APTC/VPTC;                                       
double CFil        = AFil/VFil;                                       
double CPla        = APla/VPla;
double CVL         = CL/PL;                                            
double CVK         = CKb;                                              
double CVM         = CM/PM;                                            
double CVF         = CF/PF;                                            
double CVRest      = CRest/PRest;                                   
double CVPla       = CPla/PPla;

// #+ Concentrations for fetus; 
// #+ CPlas_Fet_free   : mg/L, Free F-53B concentration in the plasma
// #+ CRest_Fet        : mg/L, Free F-53B concentration in the rest of body
// #+ CL_Fet           : mg/L, Free F-53B concentration in the liver
// #+ CB_Fet           : mg/L, Free F-53B concentration in the brain

double CPlas_Fet_free  = APlas_Fet_free /(VPlas_Fet + 1e-7); 
double CPlas_Fet       = CPlas_Fet_free/Free_Fet;
double CRest_Fet       = ARest_Fet/(VRest_Fet+ 1e-7);
double CVRest_Fet      = CRest_Fet/ PRest_Fet;
double CL_Fet          = AL_Fet/ (VL_Fet+ 1e-7);
double CVL_Fet         = CL_Fet/ PL_Fet;
double CB_Fet          = AB_Fet/ (VB_Fet+ 1e-7);
double CVB_Fet         = CB_Fet/ PB_Fet;

// #+ Equation for estimation the rates of compartments in the PBPK model for the pregnant mice
// #+ RA_baso      : mg/h, Rate of basolateral transporters
// #+ RA_apical    : mg/h, Rate of apical transporter
// #+ Rdif         : mg/h, Rate of diffusion from blood into the PTC
// #+ RAefflux     : mg/h, Rate of efflux clearance rate from PTC to blood
// #+ RCI          : mg/h, Rate of clerance (CL) to via glomerular filtration (GFR)
// #+ RPTC         : mg/h, Rate of change in PTC
// #+ RFil         : mg/h, Rate of change in Fil
// #+ RKb          : mg/h, Rate of change in Kidney-serum compartment
// #+ RST          : mg/h, Rate of change in stomach compartment
// #+ RSI          : mg/h, Rate of change in small intestines
// #+ RabsST       : mg/h, Rate of absorption in Stomach
// #+ RabsSI       : mg/h, Rate of absorption in small intestines
// #+ RL           : mg/h, Rate of change in liver compartment
// #+ RF           : mg/h, Rate of chnage in fat comaprtment
// #+ RM           : mg/h, Rate of change in mammary tissues
// #+ RRest        : mg/h, Rate of change in rest of body
// #+ RPla         : mg/h, Rate of change in placenta compartment
// #+ RPlas_free   : mg/h, Rate of free F-53B change in the maternal plasma
// #+ RPlas        : mg/h, Rate of total F-53B change in the maternal plasma
// #+ Rurine       : mg/h, Rate of change in urine
// #+ Rtrans_1     : mg/h, Rate of change in placenta tranfer from mother to fetus
// #+ Rtrans_2     : mg/h, Rate of change in placenta tranfer from fetus to mother
// #+ Rtrans_3     : mg/h, Rate of change in amniotic fluid tranfer from fetus to fluid
// #+ Rtrans_4     : mg/h, Rate of change in amniotic fluid tranfer from fluid to fetus
// #+ RAm          : mg/h, Rate of change in amniotic fluid compartment 
// #+ Rbile        : mg/h, Rate of change in bile compartment 
// #+ RB_Fet       : mg/h, Rate of change in fetal brain
// #+ RL_Fet       : mg/h, Rate of change in fetal liver
// #+ RRest_Fet    : mg/h, Rate of change in fetal rest of body
// #+ RPlas_Fet    : mg/h, Rate of free F-53B change in the fetal plasma

double RA_baso        = (Vmax_baso*CKb)/(Km_baso + CKb);                 
double RA_apical      = (Vmax_apical*CFil)/(Km_apical + CFil);        
double Rdif           = Kdif*(CKb - CPTC);                               
double RAefflux       = Kefflux*APTC;                                
double RCI            = CPlas*GFR*Free;                             
double RPTC           = Rdif + RA_apical + RA_baso - RAefflux;           
double RFil           = RCI - RA_apical - AFil*Kurine;          
double RKb            = QK*(CPlas - CVK)*Free - CPlas*GFR*Free - Rdif - RA_baso;  
double RST            = -K0*AST - GE*AST;
double RabsST         = K0*AST;
double RSI            = GE*AST - Kabs*ASI - Kunabs*ASI;
double RabsSI         = Kabs*ASI;
double RL             = QL*(CPlas - CVL)*Free - Kbile*AL + Kabs*ASI + K0*AST;   
double RF             = QF_P*(CPlas - CVF)*Free;            
double RM             = QM_P*(CPlas - CVM)*Free;                            
double RRest          = QRest*(CPlas - CVRest)*Free;                         
double Rurine         = Kurine*AFil;                                   
double Rfeces         = Kbile*AL + Kunabs*ASI;                         
double Rtrans_1       = Ktrans_1*CVPla*Free;                          
double Rtrans_2       = Ktrans_2*CPlas_Fet*Free;                           
double Rtrans_3       = Ktrans_3*CVRest_Fet*Free_Fet;                              
double Rtrans_4       = Ktrans_4*CAm;
double RAm            = Rtrans_3 - Rtrans_4;                               
double RPla           = QPla*(CPlas - CVPla)*Free + Rtrans_2 - Rtrans_1;          
double RPlas_free     = (QRest*CVRest*Free) + (QK*CVK*Free) + (QL*CVL*Free) + (QM_P*CVM*Free) + (QF_P*CVF*Free) + (QPla*CVPla*Free) - 
                        (QC_P*CPlas*Free) + RAefflux;  
double RL_Fet         = QL_Fet*(CPlas_Fet - CVL_Fet)*Free_Fet;     
double RB_Fet         = QB_Fet*(CPlas_Fet - CVB_Fet)*Free_Fet;     
double RRest_Fet      = QRest_Fet*(CPlas_Fet - CVRest_Fet)*Free_Fet - Rtrans_3 + Rtrans_4;     
double RPlas_Fet_free = (QRest_Fet*CVRest_Fet*Free_Fet) + (QL_Fet*CVL_Fet*Free_Fet)+ (QB_Fet*CVB_Fet*Free_Fet) - (QC_Fet*CPlas_Fet*Free_Fet) + Rtrans_1 - Rtrans_2;

// #+ ODE equation for compartments in the pregnant mice
dxdt_A_baso           = RA_baso;                                           
dxdt_A_apical         = RA_apical;                                       
dxdt_Adif             = Rdif;                                              
dxdt_Aefflux          = RAefflux;                                      
dxdt_APTC             = RPTC;                                            
dxdt_AFil             = RFil;                                            
dxdt_AKb              = RKb;                                                
dxdt_AST              = RST;
dxdt_ASI              = RSI;
dxdt_AabsST           = RabsST;
dxdt_AabsSI           = RabsSI;
dxdt_AL               = RL;                                                  
dxdt_AF               = RF;                                            
dxdt_AM               = RM;                                            
dxdt_ARest            = RRest;                                            
dxdt_Aurine           = Rurine;                                          
dxdt_Afeces           = Rfeces;                                          
dxdt_Atrans_1         = Rtrans_1;                                
dxdt_Atrans_2         = Rtrans_2;                                
dxdt_Atrans_3         = Rtrans_3;                                 
dxdt_Atrans_4         = Rtrans_4;
dxdt_APla             = RPla;                                        
dxdt_APlas_free       = RPlas_free;                                  
dxdt_AAm              = RAm;                                          
dxdt_AUC_CPlas        = CPlas;
dxdt_AUC_CPlas_Fet    = CPlas_Fet;
dxdt_AUC_CL           = CL;
dxdt_AUC_CPla         = CPla;
dxdt_AUC_CL_Fet       = CL_Fet;
dxdt_AUC_CB_Fet       = CB_Fet;


// #+ ODE equation for fetus compartment
dxdt_ARest_Fet        = RRest_Fet;                                 
dxdt_AL_Fet           = RL_Fet;  
dxdt_AB_Fet           = RB_Fet; 
dxdt_APlas_Fet_free   = RPlas_Fet_free;                         

// #+ Virtural compartment; input dose
dxdt_ADOSE         = 0;

// #+ Mass Balance check (pregnant mice)
double ATissue     = AF + AM + ARest + APlas_free + AKb + AL + APla + AFil + APTC + AST + ASI;
double ALoss       = Aurine + Atrans_1 - Atrans_2 + Afeces;
double ATotal      = ATissue + ALoss;
double Mbal        = ADOSE - ATotal*KDOSE; 

// #+ Mass Balance check (Fetus)
double ATissueF    = APlas_Fet_free + ARest_Fet + AL_Fet + AB_Fet;
double ALossF      = Atrans_2 + Atrans_3 - Atrans_4;
double ATotalF     = ATissueF + ALossF;
double DoseF       = Atrans_1;
double MbalF       = DoseF - ATotalF; 

$TABLE
capture Plasma        = CPlas;
capture Liver         = CL;
capture Plasma_Fet    = CPlas_Fet;
capture Placenta      = CPla;
capture AmnioticFluid = CAm;
capture Brain_Fet     = CB_Fet;
capture Liver_Fet     = CL_Fet;
capture Kidney        = CK;
capture Fat           = CF;
capture Rest          = CRest;
capture Bal           = Mbal; 
capture Bal_Fet       = MbalF;
'