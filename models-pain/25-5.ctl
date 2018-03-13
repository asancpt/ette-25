;APPENDIX 25.5 NONMEM CONTROL FILE TO ESTIMATE
;PAIN AND REMEDICATION MODEL PARAMETERS
$PROB Analgesic Pain Model for nonrandomly censored data
$INPUT ISIM ID=L1 TIME MDV PRLF=DV QUIT PTIM DOSE CL VC Q VP KA ALAG FBIO
; ISIM ID TIME MDV AMT PRLF QUIT PTIM DOSE
; ID = subject ID number
; TIME = time of dose or observation
; MDV = event ID (0=obs, 1=dose)
; AMT = dose amount (MDV=1) or "." (MDV=0)
; PRLF = Ordinal pain relief score (0=no relief thru 4=full relief)
; TQT = Time to remedication
; QUIT = Indicator of remedication (0=stay in study; 1=quit study)
; PTIM = Time of previous pain observation
; DOSE = nominal dose amount
$DATA ../data/analgesicSim.csv IGNORE=#
$THETA
; Pain Model Parameters
0.4 ; KE0 [1/h]
-2 ; BT1
-2 ; BT2
-2 ; BT3
-2 ; BT4
100 ; EC50 [ng/mL]
20 ; EMAX
0 FIX ; ALPH
2 ; GAMM
5 ; AA
;
; Remedication Model Parameters
0.2 ; LAM0
0.002 ; LAM1
0.002 ; LAM2
0.001 FIX ; LAM3
0.0001 FIX ; LAM4
$OMEGA
; Pain and Remedication Model Parameters
1 ; Pain
0.0 FIX ; Remedication
$ABBREVIATED DERIV2=NOCOMMON
$PRED
; Specify Pain Model Parameters
KE0 = THETA(1)
BT1 = THETA(2)
BT2 = THETA(3)
BT3 = THETA(4)
BT4 = THETA(5)
EC50 = THETA(6)
EMAX = THETA(7)
ALPH = THETA(8)
GAMM = THETA(9)
AA = THETA(10)
ZPAN = ETA(1)
;
; Specify Remedication Model Parameters
LAM0 = THETA(11)
LAM1 = THETA(12)
LAM2 = THETA(13)
LAM3 = THETA(14)
LAM4 = THETA(15)
ZRMD = ETA(2)
; Calculate concentration in effect compartment
K20 = CL/VC
K23 = Q/VC
K32 = Q/VP
BET1 = K23+K32+K20
BET2 = SQRT(BET1**2 - 4*K32*K20)
BETA = 0.5*(BET1 - BET2)
ALFA = K32*K20/BETA
BSL = KE0*KA*DOSE*FBIO/VC
M1 = ALFA - KA
M2 = -M1
Q1 = BETA - KA
Q2 = -Q1
R1 = KE0 - KA
R2 = -R1
S1 = BETA - ALFA
S2 = -S1
T1 = KE0 - ALFA
T2 = -T1
U1 = KE0 - BETA
U2 = -U1
Z1 = K32 - KA
Z2 = K32 - ALFA
Z3 = K32 - BETA
Z4 = K32 - KE0
TIM2 = TIME-ALAG
E1 = EXP(-KA*TIM2)
E2 = EXP(-ALFA*TIM2)
E3 = EXP(-BETA*TIM2)
E4 = EXP(-KE0*TIM2)
CE1 = Z1*E1/(M1*Q1*R1)
CE2 = Z2*E2/(M2*S1*T1)
CE3 = Z3*E3/(Q2*S2*U1)
CE4 = Z4*E4/(R2*T2*U2)
IF (TIME .LE. ALAG) THEN
CE = 0.0
ELSE
CE = BSL*(CE1 + CE2 + CE3 + CE4)
ENDIF
; Specify placebo effect for each cumulative probability
PEFF = EXP(-ALPH*TIME) - EXP(-GAMM*TIME)
PEFF1 = BT1 + AA*PEFF
PEFF2 = PEFF1 + BT2
PEFF3 = PEFF2 + BT3
PEFF4 = PEFF3 + BT4
; Specify drug effect
DEFF = EMAX * CE/(EC50 + CE)
; Logits for cummulative probabilities
LGT1 = PEFF1 + DEFF + ETA(1)
LGT2 = PEFF2 + DEFF + ETA(1)
LGT3 = PEFF3 + DEFF + ETA(1)
LGT4 = PEFF4 + DEFF + ETA(1)
; Exponentiate logit
C1 = EXP(LGT1)
C2 = EXP(LGT2)
C3 = EXP(LGT3)
C4 = EXP(LGT4)
; Calculate cumulative probability of response j
P1 = C1/(1+C1) ; P(Y<=1|X)
P2 = C2/(1+C2) ; P(Y<=2|X)
P3 = C3/(1+C3) ; P(Y<=3|X)
P4 = C4/(1+C4) ; P(Y==4|X)
; Likelihood (Yj), by pain relief (j = 0 to 4)
Y0 = 1 - P1
Y1 = P1 - P2
Y2 = P2 - P3
Y3 = P3 - P4
Y4 = P4
IND0=0
IND1=0
IND2=0
IND3=0
IND4=0
IF (PRLF .EQ. 0) IND0=1
IF (PRLF .EQ. 1) IND1=1
IF (PRLF .EQ. 2) IND2=1
IF (PRLF .EQ. 3) IND3=1
IF (PRLF .EQ. 4) IND4=1
; Calculate likelihood of pain score
YP = Y0*IND0 + Y1*IND1 + Y2*IND2 + Y3*IND3 + Y4*IND4
; Calculate likelihood of remedication
LAMM = LAM0
IF(PRLF .EQ. 1) LAMM = LAM1
IF(PRLF .EQ. 2) LAMM = LAM2
IF(PRLF .EQ. 3) LAMM = LAM3
IF(PRLF .EQ. 4) LAMM = LAM4
LAMM = LAMM + ETA(2)
; Prob that subject has not remedicated upto time=TIME
YR0 = EXP(-LAMM*TIME)
;
; Prob that subject will remedicate at time=TIME, and that subject has
; not remedicated at time=PTIM
ETIM = TIME - PTIM
IF(MDV .EQ. 0) THEN
YR10 = 1 - EXP(-LAMM*ETIM)
ELSE
YR10 = 1
ENDIF
YR11 = EXP(-LAMM*PTIM)
YR1 = YR10*YR11
;
; Prob for remedication model
YR = YR0*(1-QUIT) + YR1*QUIT
; Overall likelihood
IF (MDV .EQ. 0) THEN
Y = YP*YR
ELSE
Y = 0
ENDIF
$ESTIMATION SIG=3 MAX=9999 PRINT=1 METHOD=COND LAPLACE LIKE NOABORT
$COV PRINT=ER
$TABLE NOPRINT ONEHEADER FILE=analgesicEst3.tab
ISIM ID TIME PTIM MDV DOSE QUIT
P1 P2 P3 P4 Y0 Y1 Y2 Y3 Y4 YR0 YR1 YP YR Y

