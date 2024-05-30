PreGHumanPBPK.code <- '
$PROB
# Pre-pregnant F-53B PBPK model 
- Author    : Wei-Chun Chou, Jing Zhang
- Date      : Aug, 2023
- Strucutre : GI tract, Plasma, Liver, Fat, Mammary gland, Kidney, Filtrate, PTC
- Note1     : Pysiological parameters were obtained from Brown, 1997; Yoon et al., 2011; Loccisano et al., 2013; Corley, 2005 
- Note2     : Chemical-specific parameters were converted according to "k human= k mice*(bw human/bw mice)^-0.25"

$PARAM @annotated
BW                  : 54        : kg,                  Body weight                                   ; Value from Yoon et al. (2011)
Htc                 : 0.44      : Unitless             Hematocrit for human                          ; ICRP Publication 89 (2003)
QCC                 : 16.4      : L/h/kg^0.75,         Cardiac output                                ; Value obtained from yoon et al. (2011); Brown 1997; Forsyth 1968
QLC                 : 0.25      : Unitless,            Fraction blood flow to liver (%QCC)           ; Value obtained from Brown 1997; Fisher 2000
QKC                 : 0.141	    : Unitless,            Fraction blood flow to kidney (%QCC)          ; Value obtained from Brown 1997, Forsyth 1968
QMC                 : 0.027     : Unitless,            Fraction blood flow to Mammary gland (%QCC)   ; Value obtained from yoon et al., 2011
QFC                 : 0.052     : Unitless,            Fraction blood flow to Fat (%QCC)             ; Value obtained from Brown 1997; yoon et al., 2011
VLC                 : 0.026     : Unitless,            Fraction of liver tissue volume (%BW)         ; Value obtained from Brown 1997
VKC                 : 0.004     : Unitless,            Fraction of kidney tissue volume (%BW)        ; Value obtained from Brown 1997
VMC                 : 0.0062    : Unitless,            Fraction of Mammary gland volume (%BW)        ; Value obtained from yoon et al., 2011
VFC                 : 0.214     : Unitless,            Fraction of fat tissue (%BW)                  ; Value was assumed to be 1.5 times of male parameter (0.0723) based on ICRP, 1975.
VPlasC              : 0.0428    : L/kg BW,             Fraction of plasma volume (%BW)               ; Value obtained from Davies and Morries, 1993
VFilC               : 8.4e-4    : L/kg BW,             Fraction of filtrate (10% of Kidney volume)   ; Value obatedin from Worley et al., 2017
VPTCC               : 1.35e-4   : L/kg kidney,         Volume of proximal tubule cells               ; Value calculated from Hsu et al., 2014 (60 million PTC cells/gram kidney, 1 PTC = 2250 um3)
protein             : 2.0e-6    : mg protein/PTCs,     Amount of protein in proximal tubule cells    ; Value obtaedin from Addis et al., 1936
PL                  : 2.1       : Unitless,            Liver-to-plasma partition coefficient         ; Value obtained from mice PBPK model
PK                  : 1.26      : Unitless,            Kidney-to-plasma partition coefficient        ; Value obtained from TK study
PM                  : 0.16      : Unitless,            Mammary gland-to-plasma partition coefficient ; Value from Loccisano et al., 2013
PF                  : 0.27      : Unitless,            Fat-to-plasma partition coefficient           ; Value obtained from TK study
PRest               : 0.43      : Unitless,            Rest of body-to-plasma partition coefficient  ; Value from mice PBPK model
MW                  : 570.67    : g/mol,               F-53B molecular mass                          ; Value from Worley et al. (2017)
Free                : 0.067     : Unitless,            Free fraction of F-53B in plasma              ; Value obtained from mice PBPK model
KbileC              : 1.35e-6   : 1/(h*BW^-0.25),      Biliary elimination rate                      ; Value optimized from Chou and Lin, 2019                        
KurineC             : 0.005     : 1/(h*BW^-0.25),      Urinary elimination rate                      ; Value optimized from Chou and Lin, 2019                        
K0C                 : 0.827     : 1/(h*BW^-0.25),      Rate of absorption of F-53B in stomach        ; Value from Worley et al. (2017); Assumed to be same as PFOA
KabsC               : 0.356     : 1/(h*BW^-0.25),      Rate of absorption of F-53B in small intestines; Value from Worley et al. (2017); Assumed to be same as PFOA 
KunabsC             : 7.921e-5  : 1/(h*BW^-0.25),      Rate of unabsorbed dose to appear in feces    ; Value from Worley et al. (2017); Assumed to be same as PFOA  
GFRC                : 27.28     : L/hr/kg kiney,       Glomerular filtration rate (female)           ; Value from Corley, 2005
GEC                 : 3.510     : 1/(h*BW^0.25),       Gastric emptying constant                     ; Value from Yang et al., 2014
Vmax_baso_invitro   : 479       : pmol/mg Protein/min, Vmax of basolateral transporter               ; Value optimized from Chou and Lin, 2019
Km_baso             : 64.40     : mg/L,                Km of basolateral transporter                  ; Value calculated from Worley et al. (2007)
Vmax_apical_invitro : 51803     : pmol/mg protein/min, Vmax of apical transporter                     ; Value optimized from Chou and Lin, 2019 
Km_apical           : 20.1      : mg/L,                Km of apical transporter                       ; Value optimized from Chou and Lin, 2019
RAFbaso             : 1         : Unitless             Relative activity factor for baso             ; Value calculated from Worley et al. (2017); Assumed to be same as PFOA 
RAFapi              : 0.001     : Unitless             Relative activity factor for api              ; Value optimized from Chou and Lin, 2019 
Kdif                : 1.46e-6   : L/h,                 Diffusion rate from PTCs to kidney serum      ; Value from Worley et al. (2017); Assumed to be same as PFOA
KeffluxC            : 0.821     : 1/(h*BW^-0.25),      Rate of clearance of F-53B from PTCs into blood; Value optimized from Chou and Lin, 2019

$MAIN
// #+ Time varabiles: Gestational day and age
// #+ YEAR          

double YEAR           = TIME/24*365;                                     

// #+ Physiological parameters: Blood flows
// #+ QC            : L/h,       Cardiac output 
// #+ QK            : L/h,       Blood flows of kidney
// #+ QL            : L/h,       Blood flows of liver  
// #+ QM            : L/h,       Blood flows of mammary gland 
// #+ QF            : L/h,       Blood flows of fat tissue 
// #+ QRest         : L/h,       Blood flows of rest of body  
// #+ QBal          : L/h,       Mass balance check for total blood flows 

double QC           = QCC*pow(BW, 0.75)*(1-Htc);                  
double QK           = QKC*QC;                                     
double QL           = QLC*QC;                                     
double QM           = QMC*QC;                                     
double QF           = QFC*QC;                                     
double QRest        = QC - QK - QL - QM - QF;                          
double QBal         = QC - (QK + QL + QRest + QM + QF);

// Mass balance adjusted factor; avoding the negative occur at time = 0
double KDOSE        = (TIME==0)?0:1;

// #+ Physiological parameters: tissue volume
// #+ VL            : L,         Volume of liver
// #+ VK            : L,         Volume of Kidney
// #+ VKb           : L,         Volume of blood in the kidney  
// #+ VFil          : L,         Volume of filtrate
// #+ VPTC          : L,         Volume of proximal tubule cells (PTCs)
// #+ MK            : g,         Kidney weight                                 ; based on density of kidney = 1.0 g/mL
// #+ ML            : g,         Liver weight                                  ; based on density of liver = 1.05 g/mL, Overmeyer 1987
// #+ VBal          : L,         Mass balance check for total tissue volume

double VL           = VLC*BW;                                     
double VK           = VKC*BW;                                     
double VKb          = VK*0.16;                                   
double VM           = VMC*BW;                                     
double VF           = VFC*BW;                                     
double VPlas        = VPlasC*BW;                               
double VFil         = VFilC*BW;                                  
double VPTC         = VK*VPTCC;                                  
double VRest        = (0.93*BW) - VL - VK - VM - VF - VPlas;  
double MK           = VKC*BW*1000;                                                                            
double ML           = VLC*1.05*BW*1000;                                       
double VBal         = (0.93*BW) - VRest - (VL + VK + VM + VF + VPlas);  

// #+ Kinetic parameters
// #+ Kbile         : 1/h, Biliary elimination, liver to feces storage
// #+ Kurine        : 1/h, Urinary elimination
// #+ Kabs          : 1/h, Rate of absorption of F-53B from small intestine to liver
// #+ Kunabs        : 1/h, Rate of unabsorbed dose to appear in feces
// #+ PTC           : cells/kg BW, Number of PTC (cells/kg BW) (based on 60 million PTC/gram kidney)
// #+ MPTC          : g, mass of the proximal tubule cells (assuming density 1 kg/L)
// #+ Vmax_basoC    : mg/h/kg BW^0.75, Vmax of basolateral transporters (average Oat1 and Oat3)
// #+ Vmax_apicalC  : mg/h/kg BW^0.75, Vmax of apical transporters in in vitro studies (Oatp1a1)
// #+ Vmax_baso     : mg/h,
// #+ Vmax_apical   : mg/h,
// #+ Kefflux       : 1/h, Efflux clearance rate from PTC to blood  
// #+ GFR           : L/h, Glomerular filtration rate
// #+ GE            : 1/h, Gasric emptying rate
// #+ K0            : 1/h, Rate of uptake from the stomach into the liver

double Kbile        = KbileC*pow(BW,(-0.25));                  
double Kurine       = KurineC*pow(BW,(-0.25));                 
double Kabs         = KabsC*pow(BW,(-0.25));                    
double Kunabs       = KunabsC*pow(BW,(-0.25));                
double PTC          = VKC*6e7*1000;                              
double MPTC         = VPTC*1000;                                
double Vmax_basoC   = (Vmax_baso_invitro*RAFbaso*PTC*protein*60*(MW/1e12)*1000);         
double Vmax_apicalC = (Vmax_apical_invitro*RAFapi*PTC*protein*60*(MW/1e12)*1000);      
double Vmax_baso    = Vmax_basoC*pow(BW,0.75);             
double Vmax_apical  = Vmax_apicalC*pow(BW,0.75);         
double Kefflux      = KeffluxC*pow(BW,(-0.25));              
double GFR          = GFRC*(MK/1000);                            
double GE           = GEC*pow(BW,(-0.25));                        
double K0           = K0C*pow(BW,(-0.25));                        


$CMT ADOSE A_baso A_apical Adif Aefflux ACI APlas_free AUCCPlas_free APTC AUCCPTC AFil AUCFil Aurine AKb AUCKb ARest AUCCRest AST
AabsST ASI AabsSI Afeces AL AUCCL AF AUCCF AM AUCCM

$ODE
// #+ Concentrations in the tissues and in the venous palsma leaving each of the tissues (Unit: mg/L) 
// #+ CPlas_free  : mg/L, Free F-53B concentration in the plasma
// #+ CPlas       : mg/L, Concentration of total F-53B in the plasma
// #+ CL          : mg/L, Concentration of F-53B in the liver compartment
// #+ CKb         : mg/L, Concentration of F-53B in venous plasma leaving kidney
// #+ CK          : mg/L, Concentration of F-53B in Kidney compartment
// #+ CM          : mg/L, Concentration of F-53B in mammary gland compartment
// #+ CF          : mg/L, Concentration of F-53B in fat compartment
// #+ CRest       : mg/L, Concentration of F-53B in rest of the body
// #+ CPTC        : mg/L, Concentration of F-53B in proximal tubule cells (PTC)
// #+ CFil        : mg/L, Concentration of F-53B in Filtrate (fil)
// #+ CVL         : mg/L, Concentration in the venous blood leaving the liver
// #+ CVK         : mg/L, Concentration in the venous blood leaving the kidney
// #+ CVM         : mg/L, Concentration in the venous blood leaving the mammary gland
// #+ CVF         : mg/L, Concentration in the venous blood leaving the fat
// #+ CVRest      : mg/L, Concentration in the venous blood leaving the rest of body

double CPlas_free = APlas_free/VPlas;                      
double CPlas      = CPlas_free/Free;                               
double CL         = AL/VL;                                      
double CKb        = AKb/VKb;                                   
double CK         = CKb*PK;                                     
double CM         = AM/VM;                                      
double CF         = AF/VF;                                      
double CRest      = ARest/VRest;                             
double CPTC       = APTC/VPTC;                                       
double CFil       = AFil/VFil;                                       
double CVL        = CL/PL;                                     
double CVK        = CKb;                                       
double CVM        = CM/PM;                                     
double CVF        = CF/PF;                                     
double CVRest     = CRest/PRest;                    

// #+ Equation for estimating the rate of compartment in mother PBPK
// #+ Virtural kidney sub-compartment 
// #+ RA_baso     : mg/h, Rate of basolateral transporters
// #+ RA_apical   : mg/h, Rate of apical transporter
// #+ Rdif        : mg/h, Rate of diffusion from into the PTC
// #+ RAefflux    : mg/h, Rate of efflux clearance rate from PTC to blood
// #+ RCI         : mg/h, Rate of clearance(CL) to via glomerular filtration (GFR)
// #+ RPTC        : mg/h, Rate of change in PTC
// #+ RFil        : mg/h, Rate of change in Fil
// #+ RKb         : mg/h, Rate of change in Kidney-serum compartment
// #+ RST         : mg/h, Rate of change in stomach compartment
// #+ RSI         : mg/h, Rate of change in small intestines
// #+ RabsST      : mg/h, Rate of change of absorption in Stomach
// #+ RabsSI      : mg/h, Rate of change of absorption in small intestines
// #+ RL          : mg/h, Rate of change in liver compartment
// #+ RF          : mg/h, Rate of change in fat compartment
// #+ RM          : mg/h, Rate of change in mammary tissues
// #+ RRest       : mg/h, Rate of change in rest of body
// #+ Rurine      : mg/h, Rate of change in urine
// #+ Rbile       : mg/h, Rate of change in bile compartment 
// #+ RPlas_free  : mg/h, Rate of change in the plasma

double RA_baso    = (Vmax_baso*CKb)/(Km_baso + CKb);                                            
double RA_apical  = (Vmax_apical*CFil)/(Km_apical + CFil);                                   
double Rdif       = Kdif*(CKb - CPTC);                                                           
double RAefflux   = Kefflux*APTC;                                                             
double RCI        = CPlas*GFR*Free;                                                                   
double RPTC       = Rdif + RA_apical + RA_baso - RAefflux;                                        
double RFil       = RCI - RA_apical - AFil*Kurine;                                        
double RKb        = QK*(CPlas - CVK)*Free - RCI - Rdif - RA_baso;                               
double RST        = - K0*AST - GE*AST;                                                             
double RSI        = GE*AST - Kabs*ASI - Kunabs*ASI;                                                
double RabsST     = K0*AST;                                                                     
double RabsSI     = Kabs*ASI;
double RL         = QL*(CPlas - CVL)*Free - Kbile*AL + Kabs*ASI + K0*AST;                                
double RF         = QF*(CPlas - CVF)*Free;                                                               
double RM         = QM*(CPlas - CVM)*Free;   
double RRest      = QRest*(CPlas - CVRest)*Free;                                                      
double Rurine     = Kurine*AFil;                                                               
double Rbile      = Kbile*AL;                                                                    
double Rfeces     = Rbile + Kunabs*ASI; 
double RPlas_free = (QRest*CVRest*Free) + (QK*CVK*Free) + (QL*CVL*Free) + (QM*CVM*Free) + (QF*CVF*Free) - (QC*CPlas*Free) + RAefflux;  

// #+ ODE equation 
dxdt_A_baso       = RA_baso;                                                                      
dxdt_A_apical     = RA_apical;                                                                  
dxdt_Adif         = Rdif;                                                                           
dxdt_Aefflux      = RAefflux;                                                                    
dxdt_ACI          = RCI;                                                                             
dxdt_APTC         = RPTC;                                                                           
dxdt_AFil         = RFil;                                                                           
dxdt_AKb          = RKb;                                                                             
dxdt_AST          = RST; 
dxdt_AabsST       = RabsST;                                                                       
dxdt_ASI          = RSI;                                                                            
dxdt_AabsSI       = RabsSI;                                                                       
dxdt_AL           = RL;                                                                               
dxdt_AF           = RF;                                                                               
dxdt_AM           = RM; 
dxdt_ARest        = RRest;                                                                         
dxdt_Aurine       = Rurine;                                          
dxdt_Afeces       = Rfeces;                                          
dxdt_APlas_free   = RPlas_free;                                                               

// #+ Virtural compartment for estimating input dose
dxdt_ADOSE        = 0;

// #+ AUC equation 
dxdt_AUCCPlas_free = CPlas_free;                                                                  
dxdt_AUCCPTC      = CPTC;                                                                       
dxdt_AUCFil       = CFil;                                                                      
dxdt_AUCKb        = CK;                                                                            
dxdt_AUCCRest     = CRest;                                                                     
dxdt_AUCCL        = CL;                                                                            
dxdt_AUCCM        = CM;                                                                            
dxdt_AUCCF        = CF;                                                                            

// #+  Mass Balance check
double ATissue    = APlas_free + ARest + AKb + AFil + APTC + AL + AM + AF + AST + ASI;
double ALoss      = Aurine + Afeces;
double ATotal     = ATissue + ALoss;
double Mbal       = ADOSE - ATotal;

$TABLE
capture Plasma    = CPlas*1000;
capture Liver     = CL;
capture Kidney    = CK;
capture Bal       = Mbal;
capture AT        = ATissue;
'

GHumanPBPK.code <- '
$PROB
# Gestational F-53B PBPK model for human
- Author : Wei-Chun Chou, Jing Zhang
- Date   : Aug, 2023
- Structure: GI tract, Plasma, Liver, Fat, Mammary gland, Kidney, Filtrate, PTC, rest of fetuses, fetal liver, amniotic fluid  
- Note1  : Initial physiological parameters and optimized parameters is matched with the values in non-pregnantPBPK models
- Note2  : Growth equation of physiological parameters values and fetus-specific parameters was taken from Loccisano et al. (2013); Yoon et al. (2011); Kapraun et al., 2019 

$PARAM @annotated
// #+  Maternal parameters
BW                  : 60        : kg,                  Body weight                                   ; Value from Yoon et al. (2011)
QCC                 : 16.4      : L/h/kg^0.75,         Cardiac output                                ; Value obtained from yoon et al. (2011); Brown 1997; Forsyth 1968
QLC                 : 0.25      : Unitless,            Fraction blood flow to liver (%QCC)           ; Value obtained from Brown 1997; Fisher 2000
QKC                 : 0.141	    : Unitless,            Fraction blood flow to kidney (%QCC)          ; Value obtained from Brown 1997, Forsyth 1968
QMC                 : 0.027     : Unitless,            Fraction blood flow to Mammary gland (%QCC)   ; Value obtained from yoon et al., 2011
QFC                 : 0.052     : Unitless,            Fraction blood flow to Fat (%QCC)             ; Value obtained from Brown 1997; yoon et al., 2011
VLC                 : 0.026     : Unitless,            Fraction of liver tissue volume (%BW)         ; Value obtained from Brown 1997
VKC                 : 0.004     : Unitless,            Fraction of kidney tissue volume (%BW)        ; Value obtained from Brown 1997
VMC                 : 0.0062    : Unitless,            Fraction of Mammary gland volume (%BW)        ; Value obtained from yoon et al., 2011
VFC                 : 0.214     : Unitless,            Fraction of fat tissue (%BW)                  ; Value was assumed to be 1.5 times of male parameter (0.0723) based on ICRP, 1975.
VPlasC              : 0.0428    : L/kg BW,             Fraction of plasma volume (%BW)               ; Value obtained from Davies 1993
VFilC               : 8.4e-4    : L/kg BW,             Fraction of filtrate (10% of Kidney volume)   ; Value obatedin from Worley et al., 2017
VPTCC               : 1.35e-4   : L/kg kidney,         Volume of proximal tubule cells               ; Value calculated from Hsu et al., 2014 (60 million PTC cells/gram kidney, 1 PTC = 2250 um3)
protein             : 2.0e-6    : mg protein/PTCs,     Amount of protein in proximal tubule cells    ; Value obtaedin from Addis et al., 1936
PL                  : 2.1       : Unitless,            Liver-to-plasma partition coefficient         ; Value obtained from TK study
PK                  : 1.26      : Unitless,            Kidney-to-plasma partition coefficient        ; Value obtained from TK study
PM                  : 0.16      : Unitless,            Mammary gland-to-plasma partition coefficient ; Value from Loccisano et al., 2013
PF                  : 0.27      : Unitless,            Fat-to-plasma partition coefficient           ; Value obtained from TK study
PRest               : 0.43      : Unitless,            Rest of body-to-plasma partition coefficient  ; Value obtained from TK study
PPla                : 0.07      : Unitless,            Placenta-to-plasma partition coefficient      ; Value obtained from TK study
MW                  : 570.67    : g/mol,               F-53B molecular mass                           ; Value from Worley et al. (2017)
Free                : 0.067     : Unitless,            Free fraction of F-53B in plasma               ; Value optimized from Chou and Lin, 2019
KbileC              : 1.31e-6   : 1/(h*BW^-0.25),      Biliary elimination rate                      ; Value optimized from Chou and Lin, 2019                        
KurineC             : 0.005     : 1/(h*BW^-0.25),      Urinary elimination rate                      ; Value optimized from Chou and Lin, 2019                        
K0C                 : 0.806     : 1/(h*BW^-0.25),      Rate of absorption of F-53B in stomach         ; Value from Worley et al. (2017); Assumed to be same as PFOA
KabsC               : 0.347     : 1/(h*BW^-0.25),      Rate of absorption of F-53B in small intestines; Value from Worley et al. (2017); Assumed to be same as PFOA 
KunabsC             : 7.715e-5  : 1/(h*BW^-0.25),      Rate of unabsorbed dose to appear in feces    ; Value from Worley et al. (2017); Assumed to be same as PFOA  
GFRC                : 27.28     : L/hr/kg kiney,       Glomerular filtration rate (female)           ; Value from Corley, 2005
GEC                 : 3.510     : 1/(h*BW^0.25),       Gastric emptying constant                     ; Value from Yang et al., 2014
Vmax_baso_invitro   : 479       : pmol/mg Protein/min, Vmax of basolateral transporter               ; Value optimized from Chou and Lin, 2019
Km_baso             : 64.4      : mg/L,                Km of basolateral transporter                  ; Value calculated from Worley et al. (2007)
Vmax_apical_invitro : 51803     : pmol/mg protein/min, Vmax of apical transporter                     ; Value optimized from Chou and Lin, 2019 
Km_apical           : 20.1      : mg/L,                Km of apical transporter                       ; Value optimized from Chou and Lin, 2019
RAFbaso             : 1         : Unitless             Relative activity factor for baso             ; Value calculated from Worley et al. (2017); Assumed to be same as PFOA 
RAFapi              : 0.001     : Unitless             Relative activity factor for api              ; Value optimized from Chou and Lin, 2019 
Kdif                : 1.42e-6   : L/h,                 Diffusion rate from PTCs to kidney serum      ; Value from Worley et al. (2017); Assumed to be same as PFOA
KeffluxC            : 0.800     : 1/(h*BW^-0.25),      Rate of clearance of F-53B from PTCs into blood; Value optimized from Chou and Lin, 2019
Ktrans1C            : 0.200     : L/h/kg^0.75,         Mother-to-fetus placental transfer rate       ; Value from Loccisano et al., 2013
Ktrans2C            : 0.084     : L/h/kg^0.75,         Fetus-to-mother placental transfer rate       ; Value from Loccisano et al., 2013
Ktrans3C            : 0.033     : L/h/kg^0.75,         Fetus-to-amniotic fluid transfer rate         ; Value from Loccisano et al., 2013 
Ktrans4C            : 1.429e-4  : L/h/kg^0.75,         Amniotic fluid-to-fetus transfer rate         ; Value from Loccisano et al., 2013 
VPlasC_Fet          : 0.0428    : Unitless,            Fraction of fetal plasma volume (%BW)         ; Value was assumed to be same with mother; Loccisano et al. (2013)
Free_Fet            : 0.070     : Unitless,            Free fraction of F-53B in plasma at t = 0     ; Value was assumed to be same with mother
PRest_Fet           : 0.22      : Unitless,            Rest of body-to-plasma partition coefficient  ; Value was assumed to be same with mother
PL_Fet              : 2.10      : Unitless,            Liver-to-plasma partition coefficient in Fetus; Value was assumed to be same with mother
PB_Fet              : 1.55      : Unitless,            Brain-to-plasma partition coefficient in Fetus; Value was assumed to be same with mother; Sharma et al., (2017)

$MAIN
// #+ Time varabiles: Gestational day and age
// #+ GD            : Day,       Gestational days
// #+ GA            : Week,      Gestational age

double GD           = TIME/24;                                     
double GA           = GD/7;                                        

// #+ Growth equations for maternal or fetal physiological parameters; Equation from Kapraun et al. (2019): https://doi.org/10.1371/journal.pone.0215906 
// #+ VPlas         : L,         Volume of plasma                              ; Equation from Kapraun et al., 2019                                         
// #+ Htc           : Unitless,  Hematocrit                                    ; Equation from Kapraun et al., 2019
// #+ GFR           : L/h,       Glomerular filtration rate                    ; Equation from Kapraun et al., 2019; orinral unit is mL/min; using (60/1000) convert unit (mL/min) to L/h
// #+ VAm           : L,         Growth equations for amniotic fluid volume (L); Equation from Kapraun et al., 2019 
// #+ VPla          : L,         Volume of placenta                            ; Equation from Kapraun et al., 2019
// #+ VF_0          : L,         Volume of fat tissue at GD0
// #+ VF_P          : L,         Volume of fat tissue during pregnnacy         ; Equation from Kapraun et al., 2019; 0.95 kg/L is the mean density (Martin et al., 1994)
// #+ VM_0          : L,         Volume of mammary gland at GD0
// #+ VM_P          : L,         Volume of mammary gland during pregnnacy      ; Equation from Loccisano et al. (2013)
// #+ VL            : L,         Volume of liver                               
// #+ VK            : L,         Volume of Kidney
// #+ VKb           : L,         Volume of blood in the kidney  
// #+ VFil          : L,         Volume of filtrate
// #+ VPTC          : L,         Volume of proximal tubule cells (PTC)
// #+ MK            : g,         Kidney weight                                 ; based on density of kidney = 1.0 g/mL
// #+ ML            : g,         Liver weight                                  ; based on density of liver = 1.05 g/mL, Overmeyer 1987
                                         
double VPlas        = (1.2406/(1 + exp(-0.31338*(GA - 17.813)))) + 2.4958;
double Htc          = (39.192 - 0.10562*GA - (7.1045E-4)*pow(GA, 2))/100;
double GFR          = (113.73 + 3.5784*GA - 0.067272*pow(GA, 2))*(0.06); 
double VAm          = (822.34/(1 + exp(-0.26988*(GA - 20.150))))/1000;
double VPla         = GA < 2? 0 : (-1.7646*GA + 0.91775*pow(GA, 2) - 0.011543*pow(GA, 3))/1000;
double VF_P         = (1/0.95)*(17.067 + 0.14937*GA); 
double VF_0         = BW*VFC;
double VM_0         = BW*VMC;
double VM_P         = BW*(((VMC + (0.0065*exp(-7.444868*exp(-0.000678*(GD*24)))))));
double VL           = VLC*BW;
double VK           = VKC*BW;
double MK           = VKC*BW*1000;                                                                            
double ML           = VLC*BW*1.05*1000;                                       
double VKb          = VK*0.16;                            
double VFil         = VFilC*BW;                                      
double VPTC         = VK*VPTCC;  


// #+ Changes of maternal body weight during pregnnacy
// #+ BW_P          : kg,        Maternal BW in pregnant women; Equation from Loccisano et al. (2013)
// #+ BWinc         : kg,        BW increased during prengncy
// #+ VRest         : L,         Volume of rest of body
// #+ VBal          : L,         Mass balance check for total tissue volume

double BW_P         = BW + (VF_P - VF_0) + (VM_P - VM_0) + VPla + VFet + VAm; 
double BWinc        = (VF_P - VF_0) + (VM_P - VM_0) + VPla + VFet + VAm;  
double VRest        = (0.93*BW_P) - (VL + VK + VM_P + VF_P + VPlas + VPla + VFet + VAm); 
double VBal         = (0.93*BW_P) - VRest -(VL + VK + VM_P + VF_P + VPlas + VPla + VFet + VAm);  
                   
// #+ Growth equations for blood flows; Euqations from  Loccisano et al. (2013) or Kapraun et al., 2019
// #+ QC            : L/h,       Cardiac output during pregnancy 
// #+ QC_0          : L/h,       Cardiac output at GD0           
// #+ QF_P          : L/h,       Blood flows of liver during pregnancy   
// #+ QF            : L/h,       Blood flows of fat tissue for non-pregnant women
// #+ QK_P          : L/h,       Blood flows of kidney during prengnacy
// #+ QK            : L/h,       Blood flows of kidney for non-pregnant women
// #+ QL_P          : L/h,       Blood flows of liver during pregnancy
// #+ QL            : L/h,       Blood flows of liver for non-pregnant women 
// #+ QPla          : L/h,       Blood flows of placenta; Equation from Kapraun et al., 2019
// #+ QM            : L/h,       Blood flows of mammary gland for non-pregnant women
// #+ QM_P          : L/h,       Blood flows of mammary gland during pregnnacy
// #+ QRest         : L/h,       Blood flows of rest of body for non-pregnant women 
// #+ QBal          : L/h,       Mass balance check for total blood flows 

double QC_0         = QCC*pow(BW,0.75)*(1-Htc);
double QC           = QC_0 + 3.2512*GA + 0.15947*pow(GA, 2) - 0.0047059*pow(GA, 3);
double QF_P         = ((0.01)*(8.5 + (-0.0175)*GA))*QC;
double QF           = QFC*QC_0;
double QK_P         = (0.01)*(17 + (-0.01)*GA)*QC;
double QK           = QKC*QC_0;
double QL_P         = (0.01)*(27 + (-0.175)*GA)*QC;
double QL           = QLC*QC_0;
double QPla         = GA < 3.6 ? 0 : (0.00022)*(GA - 3.6)*(0.4 + 0.29*GA)*QC;
double QM           = QMC*QC_0;
double QM_P         = QM*(VM_P/VM_0);
double QRest        = QC - (QK_P + QL_P + QM_P + QF_P + QPla);
double QBal         = QC - (QK_P + QL_P + QRest + QPla + QM_P + QF_P);

// #+ Kinetic parameters
// #+ Kbile         : 1/h, Biliary elimination, liver to feces storage
// #+ Kurine        : 1/h, Urinary elimination
// #+ Kabs          : 1/h, Rate of absorption of F-53B from small intestine to liver
// #+ Kunabs        : 1/h, Rate of unabsorbed dose to appear in feces
// #+ PTC           : cells/kg BW, Number of PTC (cells/kg BW) (based on 60 million PTC/gram kidney)
// #+ MPTC          : g, mass of the proximal tubule cells (assuming density 1 kg/L)
// #+ Vmax_basoC    : mg/h/kg BW^0.75, Vmax of basolateral transporters (average Oat1 and Oat3)
// #+ Vmax_apicalC  : mg/h/kg BW^0.75, Vmax of apical transporters in in vitro studies (Oatp1a1)
// #+ Vmax_baso     : mg/h,
// #+ Vmax_apical   : mg/h,
// #+ Kefflux       : 1/h, Efflux clearance rate from PTC to blood  
// #+ Ktrans_1      : L/h, Rate constant for placental transfer; Mother to fetus
// #+ Ktrans_2      : L/h, Rate constant for placental transfer; fetus to Mother
// #+ Ktrans_3      : L/h, Amniotic fluid transfer rate; fetus to fluid
// #+ Ktrans_4      : L/h, Amniotic fluid transfer rate; fluid to fetus
// #+ GE            : 1/h, Gasric emptying time
// #+ K0            : 1/h, Rate of uptake from the stomach into the liver

double Kbile        = KbileC*pow(BW_P,(-0.25));                   
double Kurine       = KurineC*pow(BW_P,(-0.25));                
double Kabs         = KabsC*pow(BW_P,(-0.25));                     
double Kunabs       = KunabsC*pow(BW_P,(-0.25));                
double PTC          = VKC*6e7*1000;                                       
double MPTC         = VPTC*1000;                                         
double Vmax_basoC   = (Vmax_baso_invitro*RAFbaso*PTC*protein*60*(MW/1e12)*1000);       
double Vmax_apicalC = (Vmax_apical_invitro*RAFapi*PTC*protein*60*(MW/1e12)*1000);  
double Vmax_baso    = Vmax_basoC*pow(BW_P,0.75);                          
double Vmax_apical  = Vmax_apicalC*pow(BW_P,0.75);                      
double Kefflux      = KeffluxC*pow(BW_P,(-0.25));                         
double Ktrans_1     = Ktrans1C*(pow(VFet,0.75));                      
double Ktrans_2     = Ktrans2C*(pow(VFet,0.75));                     
double Ktrans_3     = Ktrans3C*(pow(VFet,0.75));                      
double Ktrans_4     = Ktrans4C*(pow(VFet,0.75));                     
double GE           = GEC*pow(BW_P,(-0.25));                        
double K0           = K0C*pow(BW_P,(-0.25));                        

// #+ Fetus physiological parameters
// #+ Htc_Fet       : Unitness,  Fetal Hematocrit; Equation from Kapraun et al., 2019; (0.5) fetal hematocrit (same as newborn); Sisson, et al 1959
// #+ VFet          : L,         Volume or mass (kg) of a human fetus; Equation obtained from Kapraun et al., 2019
// #+ VPlas_Fet     : L,         Plasma volume of fetal 
// #+ VL_Fet        : L,         Volume of fetal liver; Equation obtained from Kapraun et al., 2019
// #+ VRest_Fet     : L,         Volume of fetal rest of body
// #+ VFetBal       : L,         Mass balance for fetal volume
// #+ QCC_Fet       : L/h/kg,    Cardiac output of Fetus; value obtained from Mendes et al. (2015); Table 3     
// #+ QC_Fet        : L/h,       Cardiac output of Fetus; equation obatained from Loccisano et al. (2013)
// #+ QL_Fet        : L/h,       Blood flows of fetal Liver; Equation modifed from Kapraun et al., 2019 (blood flow in the intra-abdominal umbilical vein "Qpla_f" was assumed as 0)
// #+ QRest_Fet     : L/h,       Blood flows of fetal rest of body
// #+ QBal_Fet      : L/h,       Mass balance check for fetal blood flows

double Htc_Fet      = (4.5061*GA - 0.18487*pow(GA,2) + 0.0026766*pow(GA,3))/100;
double VFet         = (0.0018282*exp((15.12691)*(1-exp(-0.077577*GA))))/1000;
double VPlas_Fet    = VPlasC_Fet*VFet;
double VL_Fet       = (0.0075*exp(10.68*(1-exp((-0.062)*GA))))/1050;
double VB_Fet       = ((0.01574*exp(10.91*(1-exp((-0.065)*GA))))/1040)/1000;
double VRest_Fet    = (0.93*VFet) - VPlas_Fet - VL_Fet - VB_Fet;
double VFetBal      = (0.93*VFet) - (VRest_Fet + VPlas_Fet + VL_Fet + VB_Fet);

double QCC_Fet      = 54;
double QC_Fet       = QCC_Fet*VPlas_Fet*(1 - Htc_Fet); 
double QL_Fet       = (6.5/54)*(1 - 26.5/75)*QC_Fet;
double QB_Fet       = (14.3/75)*QC_Fet; 
double QRest_Fet    = QC_Fet - QL_Fet - QB_Fet;
double QBal_Fet     = QC_Fet - (QRest_Fet + QL_Fet + QB_Fet);
            
// Mass balance adjusted factor; avoding the negative occur at time = 0
double KDOSE        = (TIME==0)?0:1;

$INIT @annotated
// #+ Set up the initial concentration; the initialconcentration was obtained from the steady-state concentration from non-pregnant model
ADOSE             : 0   : mg, Amount of input dose; virtual compartment
APlas_free        : 0   : mg, Amount of F-53B in the plasma 
APTC              : 0   : mg, Amount of F-53B in the proximal tubule cells
AFil              : 0   : mg, Amount of F-53B in the filtrate
Aurine            : 0   : mg, Amount of F-53B in urine
AKb               : 0   : mg, Amount of F-53B in the kidney blood 
ARest             : 0   : mg, Amount of F-53B in rest of body
Afeces            : 0   : mg, Amount of F-53B in the feces  
AL                : 0   : mg, Amount of F-53B in the liver
AM                : 0   : mg, Amount of F-53B in the mammary gland
AF                : 0   : mg, Amount of F-53B in fat compartment
A_baso            : 0   : mg, Amount of F-53B in the apical subcompartment 
A_apical          : 0   : mg, Amount of F-53B in the apical subcompartment 
Adif              : 0   : mg, Amount of of F-53B in the diffused into the proximal tubule cells (PTC)
Aefflux           : 0   : mg, Amount of F-53B in the efflux virtual compartment; simulation of F-53B from PTC to blood
Atrans_1          : 0   : mg, Amount of F-53B in the placental transfer from mother to fetus  
Atrans_2          : 0   : mg, Amount of F-53B in the placental transfer from fetus to mother 
Atrans_3          : 0   : mg, Amount of F-53B in the Amniotic fluid transfer from Amniotic fluid to fetus  
Atrans_4          : 0   : mg, Amount of F-53B in the Amniotic fluid transfer from fetus to Amniotic fluid 
APla              : 0   : mg, Amount of F-53B in the placenta compartment
ASI               : 0   : mg, Amount of F-53B in the small intestine compartment
AST               : 0   : mg, Amount of F-53B in the stomach compartment
AabsST            : 0   : mg, Amount of absorbed F-53B in the stomach compartment
AabsSI            : 0   : mg, Amount of absorbed F-53B in the small intestine compartment
AAm               : 0   : mg, Amount of F-53B in the Amniotic fluid compartment
AL_Fet            : 0   : mg, Amount of F-53B in the Fetal liver compartment
AB_Fet            : 0   : mg, Amount of F-53B in the Fetal brain compartment
ARest_Fet         : 0   : mg, Amount of F-53B in the rest of body compartment
APlas_Fet_free    : 0   : mg, Amount of F-53B in the plasma compartment
AUC_CPlas         : 0   : mg/L*hr, Area under curve of F-53B in plasma
AUC_CPlas_Fet     : 0   : mg/L*hr, Area under curve of F-53B in fetal plasma


$ODE
// #+ Concentrations in the tissues and in the venous palsma leaving each of the tissues (Unit: mg/L) 
// #+ Concentrations for mother; CX indicate the cocentration in the tissue (e.g., CL); CVX represent the concentration of chemical in tissue leaving tissue 
// #+ CPlas_free  : mg/L, Free F-53B concentration in the plasma
// #+ CPlas       : mg/L, Concentration of total F-53B in the plasma
// #+ CL          : mg/L, Concentration of F-53B in the liver compartment
// #+ CKb         : mg/L, Concentration of F-53B in venous plasma leaving kidney
// #+ CK          : mg/L, Concentration of F-53B in Kidney compartment
// #+ CM          : mg/L, Concentration of F-53B in the mammary gland compartment
// #+ CF          : mg/L, Concentration of F-53B in the fat compartment
// #+ CAm         : mg/L, Concentration of F-53B in the Amniotic fluid compartment
// #+ CRest       : mg/L, Concentration of F-53B in the rest of the body
// #+ CPTC        : mg/L, Concentration of F-53B in proximal tubule cells (PTC)
// #+ CFil        : mg/L, Concentration of F-53B in Filtrate (fil)
// #+ CPla        : mg/L, Concentration of F-53B in placenta
// #+ CVL         : mg/L, Concentration in the venous blood leaving liver
// #+ CVK         : mg/L, Concentration in the venous blood leaving kidney
// #+ CVM         : mg/L, Concentration in the venous blood leaving mammary gland
// #+ CVF         : mg/L, Concentration in the venous blood leaving fat
// #+ CVRest      : mg/L, Concentration in the venous blood leaving rest of body
// #+ CVPla       : mg/L, Concentration in the venous blood leaving placenta

double CPlas_free = APlas_free/VPlas;                             
double CPlas      = CPlas_free/Free;                                      
double CL         = AL/VL;                                             
double CKb        = AKb/VKb;                                          
double CK         = CKb*PK;                                            
double CM         = AM/VM_P;                                             
double CF         = AF/VF_P;                                            
double CAm        = AAm/(VAm + 1E-7);
double CRest      = ARest/VRest;                                    
double CPTC       = APTC/VPTC;                                       
double CFil       = AFil/VFil;                                       
double CPla       = APla/(VPla+ 1E-7);
double CVL        = CL/PL;                                            
double CVK        = CKb;                                              
double CVM        = CM/PM;                                            
double CVF        = CF/PF;                                            
double CVRest     = CRest/PRest;                                   
double CVPla      = CPla/PPla;

// #+ Concentrations for fetus; Fet represent the fetus
// #+ CPlas_Fet_free  : mg/L, Free F-53B concentration in the fetal plasma
// #+ CPlas_Fet       : mg/L, Total F-53B concentration in the fetal plasma
// #+ CRest_Fet       : mg/L, F-53B concentration in the fetal rest of body
// #+ CL_Fet          : mg/L, F-53B concentration in the fetal liver
// #+ CVL_Fet         : mg/L, F-53B concentration in the venous blood leaving fetal liver
// #+ CB_Fet          : mg/L, F-53B concentration in the fetal brain
// #+ CVB_Fet         : mg/L, F-53B concentration in the venous blood leaving fetal brain

double CPlas_Fet_free  = APlas_Fet_free/(VPlas_Fet + 1E-7);
double CPlas_Fet  = CPlas_Fet_free/Free_Fet;
double CRest_Fet  = ARest_Fet/(VRest_Fet + 1E-7);
double CVRest_Fet = ARest_Fet/((VRest_Fet + 1E-7)*PRest_Fet);
double CL_Fet     = AL_Fet/ (VL_Fet + 1E-7);
double CVL_Fet    = CL_Fet/ PL_Fet;
double CB_Fet     = AB_Fet/ (VB_Fet + 1E-7);
double CVB_Fet    = CB_Fet/ PB_Fet;

// #+ Equation for estimating the rate of compartment in mother PBPK
// #+ Virtural kidney sub-compartment 
// #+ RA_baso     : mg/h, Rate of basolateral transporters
// #+ RA_apical   : mg/h, Rate of apical transporter
// #+ Rdif        : mg/h, Rate of diffusion from into the PTC
// #+ RAefflux    : mg/h, Rate of efflux clearance rate from PTC to blood
// #+ RCI         : mg/h, Rate of clearance(CL) to via glomerular filtration (GFR)
// #+ RPTC        : mg/h, Rate of change in PTC
// #+ RFil        : mg/h, Rate of change in Fil
// #+ RKb         : mg/h, Rate of change in Kidney-serum compartment
// #+ RST         : mg/h, Rate of change in stomach compartment
// #+ RSI         : mg/h, Rate of change in small intestines
// #+ RabsST      : mg/h, Rate of change of absorption in Stomach
// #+ RabsSI      : mg/h, Rate of change of absorption in small intestines
// #+ RL          : mg/h, Rate of change in liver compartment
// #+ RF          : mg/h, Rate of change in fat compartment
// #+ RM          : mg/h, Rate of change in mammary tissues
// #+ RRest       : mg/h, Rate of change in rest of body
// #+ RPla        : mg/h, Rate of change in placenta compartment
// #+ Rurine      : mg/h, Rate of change in urine
// #+ Rtrans_1    : mg/h, Rate of change in placenta tranfer from mother to fetus
// #+ Rtrans_2    : mg/h, Rate of change in placenta tranfer from fetus to mother
// #+ Rtrans_3    : mg/h, Rate of change in amniotic fluid tranfer from fetus to fluid
// #+ Rtrans_4    : mg/h, Rate of change in amniotic fluid tranfer from fluid to fetus
// #+ RAm         : mg/h, Rate of change in amniotic fluid compartment 
// #+ Rbile       : mg/h, Rate of change in bile compartment 
// #+ RL_Fet      : mg/h, Rate of change in fetal liver
// #+ RL_Fet      : mg/h, Rate of change in fetal brain
// #+ RRest_Fet   : mg/h, Rate of change in fetal rest of body
// #+ RPlas_Fet_free : mg/h, Rate of free F-53B change in the fetal plasma

double RA_baso    = (Vmax_baso*CKb)/(Km_baso + CKb);                
double RA_apical  = (Vmax_apical*CFil)/(Km_apical + CFil);        
double Rdif       = Kdif*(CKb - CPTC);                               
double RAefflux   = Kefflux*APTC;                                
double RCI        = CPlas*GFR*Free;                              
double RPTC       = Rdif + RA_apical + RA_baso - RAefflux;           
double RFil       = RCI - RA_apical - AFil*Kurine;           
double RKb        = QK_P*(CPlas - CVK)*Free - CPlas*GFR*Free - Rdif - RA_baso;  
double RST        = -K0*AST - GE*AST;
double RSI        = GE*AST - Kabs*ASI - Kunabs*ASI;
double RabsST     = K0*AST;
double RabsSI     = Kabs*ASI;
double RL         = QL_P*(CPlas - CVL)*Free - Kbile*AL + Kabs*ASI + K0*AST;   
double RF         = QF_P*(CPlas - CVF)*Free;                            
double RM         = QM_P*(CPlas - CVM)*Free;                            
double RRest      = QRest*(CPlas - CVRest)*Free;                         
double Rfeces     = Kbile*AL + Kunabs*ASI;                       
double Rurine     = Kurine*AFil;                                   
double Rtrans_1   = Ktrans_1*CVPla*Free;                          
double Rtrans_2   = Ktrans_2*CPlas_Fet*Free_Fet;                           
double Rtrans_3   = Ktrans_3*CVRest_Fet*Free_Fet;                              
double Rtrans_4   = Ktrans_4*CAm;
double RPla       = QPla*(CPlas - CVPla)*Free + Rtrans_2 - Rtrans_1;     
double RPlas_free = (QRest*CVRest*Free) + (QK_P*CVK*Free) + (QL_P*CVL*Free) + (QM_P*CVM*Free) + (QF_P*CVF*Free) + (QPla*CVPla*Free) - (QC*CPlas*Free) + RAefflux;  
double RAm        = Rtrans_3 - Rtrans_4;                               
double RL_Fet     = QL_Fet*(CPlas_Fet - CVL_Fet)*Free_Fet;
double RB_Fet     = QB_Fet*(CPlas_Fet - CVB_Fet)*Free_Fet;
double RRest_Fet  = QRest_Fet*(CPlas_Fet - CVRest_Fet)*Free_Fet - Rtrans_3 + Rtrans_4;     
double RPlas_Fet  = (QRest_Fet*CVRest_Fet*Free_Fet) + (QL_Fet*CVL_Fet*Free_Fet) + (QB_Fet*CVB_Fet*Free_Fet) - (QC_Fet*CPlas_Fet*Free_Fet) + Rtrans_1 - Rtrans_2;
double Abile      = Kbile*AL;

// #+ ODE equation for mother compartment
dxdt_A_baso       = RA_baso;                                          
dxdt_A_apical     = RA_apical;                                       
dxdt_Adif         = Rdif;                                              
dxdt_Aefflux      = RAefflux;                                       
dxdt_APTC         = RPTC;                                           
dxdt_AFil         = RFil;                                            
dxdt_AKb          = RKb;                                                
dxdt_AST          = RST;
dxdt_AabsST       = RabsST;
dxdt_ASI          = RSI;
dxdt_AabsSI       = RabsSI;
dxdt_AL           = RL;                                                  
dxdt_AF           = RF;                                            
dxdt_AM           = RM;                                            
dxdt_ARest        = RRest;                                            
dxdt_Aurine       = Rurine;                                          
dxdt_Afeces       = Rfeces;                                          
dxdt_Atrans_1     = Rtrans_1;                                
dxdt_Atrans_2     = Rtrans_2;
dxdt_APla         = RPla;                                        
dxdt_APlas_free   = RPlas_free;                                  
dxdt_Atrans_3     = Rtrans_3;                                 
dxdt_Atrans_4     = Rtrans_4;                                  
dxdt_AAm          = RAm;                                         
dxdt_AL_Fet       = RL_Fet;
dxdt_AB_Fet       = RB_Fet;
dxdt_ARest_Fet    = RRest_Fet;                                   
dxdt_APlas_Fet_free = RPlas_Fet;                         
dxdt_AUC_CPlas     = CPlas;
dxdt_AUC_CPlas_Fet = CPlas_Fet;

// #+ Virtural compartment for estmating input dose
dxdt_ADOSE        = 0;

// #+ Mass Balance check (Mother)
double ATissue    = AF + AM + ARest + APlas_free + AKb + AL + APla + AFil + APTC + AST + ASI;
double ALoss      = Aurine + Atrans_1 - Atrans_2 + Afeces;
double ATotal     = ATissue + ALoss;
double Mbal       = ADOSE - ATotal*KDOSE; 

// #+ Mass Balance check (Fetus)
double ATissueF   = APlas_Fet_free + ARest_Fet + AL_Fet +AB_Fet;
double ALossF     = Atrans_2 + Atrans_3 - Atrans_4;
double ATotalF    = ATissueF + ALossF;
double DoseF      = Atrans_1;
double MbalF      = DoseF - ATotalF; 

// Embryo/fetus concentration
double CFtotal    = ATissueF/(VFet + 1.0e-7)/1000; 
            
$TABLE
capture Plasma    = CPlas;
capture Liver     = CL;
capture Kidney    = CK;
capture Placenta  = CPla;
capture CordB     = CPlas_Fet;
capture Fbrain    = CB_Fet;
capture Bal       = Mbal;
capture Bal_Fet   = MbalF;
'
