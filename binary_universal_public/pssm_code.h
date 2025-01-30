o = pssm(s, scansite_human_BRCT_BRCA1_pssm, 15)/2.258e+05;
if (o > 0) {
  o = log(o);
  o = 2.31706127099737e-07+(0.210947016967266-2.31706127099737e-07)/(1+exp(3.04234*(-5.94869-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tBRCT\tany_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0699854212223765);
}

o = pssm(s, scansite_human_KIN_AMPK_group_pssm, 15)/4.536e+04;
if (o > 0) {
  o = log(o);
  o = 0.00188008701167695+(0.106172644182624-0.00188008701167695)/(1+exp(0.598562*(-3.67134-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tAMPK_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0282842712474619);
}

o = pssm(s, scansite_human_KIN_CDK5_pssm, 15)/1.034e+05;
if (o > 0) {
  o = log(o);
  o = 2.75740855529468e-08+(0.360744375400965-2.75740855529468e-08)/(1+exp(0.97807*(-3.18552-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tCDK2_CDK3_CDK1_CDK5_group\t%.6f\t%.6f\n", name, *n, c, o, 0.04);
}

o = pssm(s, scansite_human_KIN_CK1_group_pssm, 15)/1.399e+05;
if (o > 0) {
  o = log(o);
  o = 9.95701610076016e-07+(0.145180988444129-9.95701610076016e-07)/(1+exp(0.993866*(-4.0268-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tCK1_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0529150262212918);
}

o = pssm(s, scansite_human_KIN_DMPK1_pssm, 15)/2.688e+03;
if (o > 0) {
  o = log(o);
  o = 0.00918460458511681+(0.291871921182266-0.00918460458511681)/(1+exp(0.63519*(-1.97074-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tDMPK_group\t%.6f\t%.6f\n", name, *n, c, o, 0.04);
}

o = pssm(s, scansite_human_KIN_EGFR_pssm, 15)/9.997e+05;
if (o > 0) {
  o = log(o);
  o = 0.0196920514889329+(0.0813793809332912-0.0196920514889329)/(1+exp(3.73875*(-3.21577-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tEGFR_group\t%.6f\t%.6f\n", name, *n, c, o, 0.04);
}

o = pssm(s, scansite_human_KIN_InsR_pssm, 15)/2.190e+04;
if (o > 0) {
  o = log(o);
  o = 0.0237693672837333+(0.102735436709296-0.0237693672837333)/(1+exp(3.09673*(-2.53617-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tInsR_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0346410161513775);
}

o = pssm(s, scansite_human_KIN_Lck_pssm, 15)/1.877e+05;
if (o > 0) {
  o = log(o);
  o = 0.0237856170384623+(0.151954202989566-0.0237856170384623)/(1+exp(0.999999*(-3-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tSrc_group\t%.6f\t%.6f\n", name, *n, c, o, 0.066332495807108);
}

o = pssm(s, scansite_human_KIN_MAPK14_pssm, 15)/1.485e+04;
if (o > 0) {
  o = log(o);
  o = 1.19929673242453e-08+(0.148220256486223-1.19929673242453e-08)/(1+exp(43.6605*(-2.82847-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tp38_group\t%.6f\t%.6f\n", name, *n, c, o, 0.04);
}

o = pssm(s, scansite_human_KIN_MAPK3_pssm, 15)/2.383e+05;
if (o > 0) {
  o = log(o);
  o = 0.000500417034004853+(0.331151825906694-0.000500417034004853)/(1+exp(0.659862*(-2.14883-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tMAPK3_MAPK1_MAPK7_NLK_group\t%.6f\t%.6f\n", name, *n, c, o, 0.04);
}

o = pssm(s, scansite_human_KIN_PKD1_pssm, 15)/1.814e+06;
if (o > 0) {
  o = log(o);
  o = 3.91054175771969e-08+(0.116863312249919-3.91054175771969e-08)/(1+exp(1.64269*(-4.77718-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tPKD_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0346410161513775);
}

o = pssm(s, scansite_human_PTB_SHC1_pssm, 15)/2.583e+04;
if (o > 0) {
  o = log(o);
  o = 1.44098589574965e-08+(0.577823050045167-1.44098589574965e-08)/(1+exp(1.87041*(-1.73069-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tPTB\tSHC1_SHC2_SHC3_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0494871659305393);
}

o = pssm(s, scansite_human_SH2_PLCG1_1_pssm, 15)/4.975e+04;
if (o > 0) {
  o = log(o);
  o = 0.00804220390320895+(0.131040568715424-0.00804220390320895)/(1+exp(1.1319*(-2.14053-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tSH2\tPLCG_1_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0404061017820884);
}

o = pssm(s, scansite_human_SH2_PLCG1_2_pssm, 15)/5.961e+04;
if (o > 0) {
  o = log(o);
  o = 0.0285162648018724+(0.179889183290924-0.0285162648018724)/(1+exp(0.999999*(-2-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tSH2\tPLCG_2_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0404061017820884);
}

o = pssm(s, scansite_LiLabOPAL_human_SH2_INPP5D_pssm, 15)/1.387e-03;
if (o > 0) {
  o = log(o);
  o = 0.0122117256875409+(0.0823386607323488-0.0122117256875409)/(1+exp(1.27666*(-3.09424-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_LiLabOPAL\thuman\tSH2\tINPP5D\t%.6f\t%.6f\n", name, *n, c, o, 0.0285714285714286);
}

o = pssm(s, scansite_LiLabOPAL_human_SH2_SH3BP2_pssm, 15)/2.135e-05;
if (o > 0) {
  o = log(o);
  o = 0.017911233061943+(0.116550893235419-0.017911233061943)/(1+exp(3.84889*(-2.44766-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_LiLabOPAL\thuman\tSH2\tSH3BP2\t%.6f\t%.6f\n", name, *n, c, o, 0.0285714285714286);
}

o = pssm(s, scansite_LiLabOPAL_human_SH2_SHC2_pssm, 15)/1.528e-01;
if (o > 0) {
  o = log(o);
  o = 0.0207364427656338+(0.101906853305672-0.0207364427656338)/(1+exp(1.31642*(-2.20547-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_LiLabOPAL\thuman\tSH2\tSHC2_SHC3_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0404061017820884);
}

o = pssm(s, scansite_LiLabOPAL_human_SH2_SHE_pssm, 15)/2.430e-03;
if (o > 0) {
  o = log(o);
  o = 0.0116298990401533+(0.133810106332926-0.0116298990401533)/(1+exp(0.608993*(-1.51017-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_LiLabOPAL\thuman\tSH2\tSHE\t%.6f\t%.6f\n", name, *n, c, o, 0.0285714285714286);
}

o = pssm(s, scansite_LiLabOPAL_human_SH2_Syk_1_pssm, 15)/4.793e-03;
if (o > 0) {
  o = log(o);
  o = 0.00656765550958794+(0.10175596578118-0.00656765550958794)/(1+exp(1.527*(-3.1099-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_LiLabOPAL\thuman\tSH2\tSyk_1\t%.6f\t%.6f\n", name, *n, c, o, 0.0285714285714286);
}

o = pssm(s, scansite_LiLabOPAL_human_SH2_TNS4_pssm, 15)/1.657e-03;
if (o > 0) {
  o = log(o);
  o = 0.00236464040993422+(0.137621910237852-0.00236464040993422)/(1+exp(0.75074*(-2.58628-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_LiLabOPAL\thuman\tSH2\tTNS4\t%.6f\t%.6f\n", name, *n, c, o, 0.0285714285714286);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_CLK2_pssm, 15)/1.413e+03;
if (o > 0) {
  o = log(o);
  o = 0.00503259201635941+(0.421209858103062-0.00503259201635941)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tCLK_group\t%.6f\t%.6f\n", name, *n, c, o, 0.04);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_DAPK3_pssm, 15)/1.074e+03;
if (o > 0) {
  o = log(o);
  o = 9.25307706661765e-11+(0.112107567281039-9.25307706661765e-11)/(1+exp(99.4089*(-13.4476-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tDAPK_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0346410161513775);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_ICK_pssm, 15)/1.384e+03;
if (o > 0) {
  o = log(o);
  o = 0.00436129308904892+(0.386594515151589-0.00436129308904892)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tRCK_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0346410161513775);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_LKB1_pssm, 15)/3.879e+02;
if (o > 0) {
  o = log(o);
  o = 0.000122769270670808+(0.110169491525424-0.000122769270670808)/(1+exp(4.61123*(-2.49789-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tLKB1\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_MST1_pssm, 15)/1.351e+03;
if (o > 0) {
  o = log(o);
  o = 0.00356383307469715+(0.339755861506917-0.00356383307469715)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tMST_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0282842712474619);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_MST4_pssm, 15)/1.178e+03;
if (o > 0) {
  o = log(o);
  o = 0.00436129308904892+(0.386594515151589-0.00436129308904892)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tYSK_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0346410161513775);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_PAK2_pssm, 15)/5.623e+02;
if (o > 0) {
  o = log(o);
  o = 0.010462108110282+(0.190793697897193-0.010462108110282)/(1+exp(1.95923*(-3.42444-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tPAK_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0489897948556636);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_Pim1_pssm, 15)/2.950e+03;
if (o > 0) {
  o = log(o);
  o = 0.00356383307469715+(0.339755861506917-0.00356383307469715)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tPim3_Pim1_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0282842712474619);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_Pim2_pssm, 15)/1.935e+03;
if (o > 0) {
  o = log(o);
  o = 0.0025226437265188+(0.266792809839167-0.0025226437265188)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tPim2\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_SLK_pssm, 15)/5.803e+02;
if (o > 0) {
  o = log(o);
  o = 0.00356383307469715+(0.339755861506917-0.00356383307469715)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tSLK_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0282842712474619);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_TGFbR2_pssm, 15)/6.494e+02;
if (o > 0) {
  o = log(o);
  o = 0.00436129308904892+(0.386594515151589-0.00436129308904892)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tACTR2_ACTR2B_TGFbR2_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0346410161513775);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_TLK1_pssm, 15)/1.861e+03;
if (o > 0) {
  o = log(o);
  o = 0.00356383307469715+(0.339755861506917-0.00356383307469715)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tTLK_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0282842712474619);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_TNIK_pssm, 15)/1.359e+03;
if (o > 0) {
  o = log(o);
  o = 0.00503259201635941+(0.421209858103062-0.00503259201635941)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tMSN_group\t%.6f\t%.6f\n", name, *n, c, o, 0.04);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_p70S6K_pssm, 15)/1.975e+03;
if (o > 0) {
  o = log(o);
  o = 0.00412223670003529+(0.0660059353402013-0.00412223670003529)/(1+exp(1.06901*(-3.53051-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tp70S6K_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0282842712474619);
}

o = pssm(s, scansite_TurkLabOPL1_human_KIN_EphA3_pssm, 15)/3.702e+02;
if (o > 0) {
  o = log(o);
  o = 0.048195566335724+(0.144279729622026-0.048195566335724)/(1+exp(8.36245*(-1.70993-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL1\thuman\tKIN\tEph_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0748331477354788);
}

o = pssm(s, scansite_YaffeLabOPL_human_KIN_NEK2_pssm, 15)/7.903e+02;
if (o > 0) {
  o = log(o);
  o = 0.0296541137013766+(0.144003522831755-0.0296541137013766)/(1+exp(2.00528*(-1.93177-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_YaffeLabOPL\thuman\tKIN\tNEK1_NEK5_NEK3_NEK4_NEK11_NEK2_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0489897948556636);
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YAR018C_pssm, 15)/3.189e+02;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYAR018C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YAR019C_pssm, 15)/4.717e+02;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYAR019C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YCL024W_pssm, 15)/3.061e+02;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYCL024W\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YDL025C_pssm, 15)/1.277e+03;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYDL025C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YDL028C_pssm, 15)/1.452e+02;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYDL028C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YDR122W_pssm, 15)/4.473e+02;
if (o > 0) {
  o = log(o);
  o = 0.012019760540367+(0.349551678433319-0.012019760540367)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.0282842712474619) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYDR122W_YLR096W_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0282842712474619);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YDR247W_pssm, 15)/6.496e+02;
if (o > 0) {
  o = log(o);
  o = 0.012019760540367+(0.349551678433319-0.012019760540367)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.0282842712474619) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYDR247W_YPL026C_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0282842712474619);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YDR283C_pssm, 15)/1.163e+02;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYDR283C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YDR477W_pssm, 15)/9.748e+02;
if (o > 0) {
  o = log(o);
  o = 2.49769204574875e-09+(0.0802516202821197-2.49769204574875e-09)/(1+exp(1.53381*(-3.9888-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYDR477W\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YDR507C_pssm, 15)/1.146e+03;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYDR507C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YFL029C_pssm, 15)/3.397e+02;
if (o > 0) {
  o = log(o);
  o = 0.00400184893465361+(0.0502484815019326-0.00400184893465361)/(1+exp(4.6856*(-3.07458-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYFL029C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YFL033C_pssm, 15)/1.099e+02;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYFL033C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YFR014C_pssm, 15)/7.929e+02;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYFR014C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YGL179C_pssm, 15)/4.697e+02;
if (o > 0) {
  o = log(o);
  o = 0.012019760540367+(0.349551678433319-0.012019760540367)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.0282842712474619) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYGL179C_YER129W_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0282842712474619);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YGL180W_pssm, 15)/3.347e+02;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYGL180W\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YHR082C_pssm, 15)/1.195e+03;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYHR082C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YHR135C_pssm, 15)/2.363e+00;
if (o > 0) {
  o = log(o);
  o = 0.0205376114650051+(0.0883511500815992-0.0205376114650051)/(1+exp(4.63185*(-3.4214-o)));
  if (o > 0.0346410161513775) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYPL204W_YER123W_YHR135C_YNL154C_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0346410161513775);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YJL141C_pssm, 15)/1.423e+03;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYJL141C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YJL164C_pssm, 15)/1.738e+03;
if (o > 0) {
  o = log(o);
  o = 2.32510810941229e-06+(0.28675732310358-2.32510810941229e-06)/(1+exp(0.997692*(-2.52641-o)));
  if (o > 0.0346410161513775) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYPL203W_YJL164C_YKL166C_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0346410161513775);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YJR059W_pssm, 15)/4.571e+02;
if (o > 0) {
  o = log(o);
  o = 0.012019760540367+(0.349551678433319-0.012019760540367)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.0282842712474619) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYKL198C_YJR059W_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0282842712474619);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YKL101W_pssm, 15)/4.898e+02;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYKL101W\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YKL116C_pssm, 15)/3.849e+02;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYKL116C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YKL126W_pssm, 15)/1.284e+03;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYKL126W\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YKL171W_pssm, 15)/6.788e+02;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYKL171W\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YLR113W_pssm, 15)/3.263e+02;
if (o > 0) {
  o = log(o);
  o = 2.72343193141226e-06+(0.146411354102443-2.72343193141226e-06)/(1+exp(1.64015*(-3.02234-o)));
  if (o > 0.0489897948556636) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYPR054W_YKL161C_YHR030C_YLR113W_YGR040W_YBL016W_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0489897948556636);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YLR248W_pssm, 15)/8.764e+02;
if (o > 0) {
  o = log(o);
  o = 0.012019760540367+(0.349551678433319-0.012019760540367)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.0282842712474619) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYGL158W_YLR248W_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0282842712474619);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YMR104C_pssm, 15)/8.174e+02;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYMR104C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YMR216C_pssm, 15)/2.858e+02;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYMR216C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YNL020C_pssm, 15)/5.033e+02;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYNL020C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YNL298W_pssm, 15)/1.252e+03;
if (o > 0) {
  o = log(o);
  o = 0.0091258504323059+(0.0821667060674633-0.0091258504323059)/(1+exp(0.51152*(-4.8582-o)));
  if (o > 0.0346410161513775) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYHL007C_YNL298W_YOL113W_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0346410161513775);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YNR047W_pssm, 15)/5.700e+02;
if (o > 0) {
  o = log(o);
  o = 0.012019760540367+(0.349551678433319-0.012019760540367)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.0282842712474619) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYCR091W_YNR047W_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0282842712474619);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YOL016C_pssm, 15)/5.933e+02;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYOL016C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YOL100W_pssm, 15)/5.727e+02;
if (o > 0) {
  o = log(o);
  o = 0.0146814798279921+(0.396928911214138-0.0146814798279921)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.0346410161513775) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYDR466W_YDR490C_YOL100W_group\t%.6f\t%.6f\n", name, *n, c, o, 0.0346410161513775);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YOR233W_pssm, 15)/2.292e+03;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYOR233W\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YOR351C_pssm, 15)/9.962e+02;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYOR351C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YPL141C_pssm, 15)/6.564e+02;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYPL141C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YPL204W_pssm, 15)/5.283e+01;
if (o > 0) {
  o = log(o);
  o = 0.00852928157631675+(0.27536231884058-0.00852928157631675)/(1+exp(0.716035*(-1.40086-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYPL204W\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

o = pssm(s, scansite_TurkLabOPL_yeast_KIN_YPL209C_pssm, 15)/1.161e+03;
if (o > 0) {
  o = log(o);
  o = 4.77735835211401e-09+(0.105691056910569-4.77735835211401e-09)/(1+exp(1.93103*(-3.55124-o)));
  if (o > 0.02) {
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\tyeast\tKIN\tYPL209C\t%.6f\t%.6f\n", name, *n, c, o, 0.02);
  }
}

