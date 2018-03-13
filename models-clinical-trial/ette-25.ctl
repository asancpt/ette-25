$PROB Analgesic Pain Model for nonrandomly censored data
$INPUT ID=L1 TIME MDV PRLF=DV QUIT PTIM DOSE
;
; ID = subject ID number
; TIME = time of dose or observation
; EVID = event ID (0=obs, 1=dose)
; PRLF = Ordinal pain relief score (0=no relief thru 4=full relief)
; TQT = Time to remedication
; QUIT = Indicator of remedication (0=stay in study; 1=quit study)
; PTIM = Time of previous pain observation
; DOSE = nominal dose amount

$DATA ../data/analgesicTemplate.csv IGNORE=#

$THETA

  ; PK Model Parameters
  2.0 FIX ; CL/F [L/h]
  10.0 FIX ; VC/F [L]
  1.0 FIX ; Q/F [L/h]
  20.0 FIX ; VP/F [L]
  2.0 FIX ; KA [1/h]
  0.1 FIX ; ALAG [h]
  1 FIX ; FBIO [1]
  ;
  ; Pain Model Parameters
  0.5 ; KE0 [1/h]
  -2.5 ; BT1
  -2 ; BT2
  -1.5 ; BT3
  -1 ; BT4
  40 ; EC50 [ng/mL]
  10 ; EMAX
  0 FIX ; ALPH
  1 ; GAMM
  3 ; AA
  ;
  ; Remedication Model Parameters
  0.5 ; LAM0
  0.005 ; LAM1
  0.005 ; LAM2
  0.001 ; LAM3
  0.0001 ; LAM4

$OMEGA
  ; PK Model Parameters
  0.09 ; CL/F
  0.09 ; VC/F
  0 FIX ; Q/F
  0 FIX ; VP/F
  0.49 ; KA
  1.0 ; ALAG
  0.49 ; FBIO
  ;
  ; Pain and Remedication Model Parameters
  1 ; Pain
  0.0 FIX ; Remediation

$ABBREVIATED DERIV2=NOCOMMON

$PRED
  ; Specify PK Model Parameters
  CL = THETA(1)*EXP(ETA(1))
  VC = THETA(2)*EXP(ETA(2))
  Q = THETA(3)*EXP(ETA(3))
  VP = THETA(4)*EXP(ETA(4))
  KA = THETA(5)*EXP(ETA(5))
  ALAG = THETA(6)*EXP(ETA(6))
  FBIO = THETA(7)*EXP(ETA(7))
  ; Specify Pain Model Parameters
  KE0 = THETA(8)
  BT1 = THETA(9)
  BT2 = THETA(10)
  BT3 = THETA(11)
  BT4 = THETA(12)
  EC50 = THETA(13)
  EMAX = THETA(14)
  ALPH = THETA(15)
  GAMM = THETA(16)
  AA = THETA(17)
  ZPAN = ETA(8)
  ;
  ; Specify Remedication Model Parameters
  LAM0 = THETA(18)
  LAM1 = THETA(19)
  LAM2 = THETA(20)
  LAM3 = THETA(21)
  LAM4 = THETA(22)
  ZRMD = ETA(9)
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
  LGT1 = PEFF1 + DEFF + ETA(8)
  LGT2 = PEFF2 + DEFF + ETA(8)
  LGT3 = PEFF3 + DEFF + ETA(8)
  LGT4 = PEFF4 + DEFF + ETA(8)
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
  ; If PREDPP is being called for simulation, then . . .
  ; . . . generate uniform random number (source #2), and
  ; . . . call this random number UNIF1
  ; . . . use UNIF1 to assign the level of pain relief
  ; . . . NOTE: PRLF is named PRLS for simulation
  IF(ICALL .EQ. 4) THEN
  CALL RANDOM(2, R)
  UNIF1 = R
  PRLF = 0
  IF (P1 .GT. UNIF1) PRLF=1
  IF (P2 .GT. UNIF1) PRLF=2
  IF (P3 .GT. UNIF1) PRLF=3
  IF (P4 .GT. UNIF1) PRLF=4
  PRLS = PRLF
  ENDIF
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
  YP = Y0*IND0+Y1*IND1+Y2*IND2+Y3*IND3+Y4*IND4
  ; Calculate likelihood of remedication
  LAMM = LAM0
  IF(PRLF .EQ. 1) LAMM = LAM1
  IF(PRLF .EQ. 2) LAMM = LAM2
  IF(PRLF .EQ. 3) LAMM = LAM3
  IF(PRLF .EQ. 4) LAMM = LAM4
  LAMM = LAMM + ETA(9)
  ; Probability that subject has not remedicated upto time=TIME
  YR0 = EXP(-LAMM*TIME)
  ;
  ; Probability that subject will remedicate at time=TIME, given that
  ; the subject has not remedicated upto time=PTIM
  ETIM = TIME - PTIM
  YR10 = 1 - EXP(-LAMM*ETIM)
  YR11 = EXP(-LAMM*PTIM)
  YR1 = YR10*YR11
  ; If PREDPP is being called for simulation, then . . .
  ; . . . generate uniform random number (source #2), and
  ; . . . call this random number UNIF2
  ; . . . use UNIF2 to determine whether subject quit
  IF(ICALL .EQ. 4) THEN
  CALL RANDOM(2, R)
  UNIF2 = R
  QUIT = 0
  IF (YR1 .GT. UNIF2) QUIT=1
  ENDIF
  ; Get simulation iteration number
  ISIM = 0
  IF (ICALL .EQ. 4) ISIM = IREP
  ;
  ; Likelihood for remedication model
  YR = YR0*(1-QUIT) + YR1*QUIT
  Y = YP*YR

$SIM (55555) (54321 UNIFORM) ONLY SUB=1

$TABLE NOPRINT ONEHEADER NOAPPEND FILE=analgesicSim.tab
ISIM ID TIME MDV PRLS QUIT PTIM DOSE CL VC Q VP KA ALAG FBIO

