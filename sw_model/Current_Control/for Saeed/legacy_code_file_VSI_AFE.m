def = legacy_code('initialize');
def.SourceFiles = {'epwm1_isr.c'};
def.HeaderFiles = {'math.h'};
def.SFunctionName = 'ex_sfun_ctrl_VSI_AFE';
def.OutputFcnSpec = 'void epwm1_isr(double y1[1], double y2[1], double y3[1], double y4[1], double y5[1], double y6[1], double y7[1], double y8[1], double y9[1], double y10[1], double u1, double u2, double u3, double u4, double u5, double u6, double u7, double u8, double u9, double u10, double u11, double u12, double u13, double u14, double u15, double u16, double u17, double u18, double u19, double u20, double u21, double u22, double u23, double u24)';
def.SampleTime = 1/20e3;
legacy_code('sfcn_cmex_generate',def);
legacy_code('compile',def);
legacy_code('slblock_generate',def)