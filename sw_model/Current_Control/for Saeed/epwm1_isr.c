#define PWM_Mode 1											//set PWM mode
//PWM_Mode=1, center aligned SVM
//PWM_Mode=2, Discontinious PWM
//------------------------------------------------constant for calculation----------------------------------------------------------
const double One_3				= 0.333333333333333333333333333;			//	1/3
const double PIover6 			= 0.52359877559829887307710723054658;		//	pi/6
const double PIover3 			= 1.0471975511965977461542144610932;		//	pi/3
const double PI			 		= 3.1415926535897932384626433832795;		//	pi
const double Two_PI				= 6.283185307179586476925286766559;			//	pi
const double SIN_PIover3	 	= 0.86602540378443864676372317075294;		//	sin(pi/3)
const double TwooverSQRT3		= 1.1547005383792515290182975610039;		//	2/sqrt(3)
const double SQRT3				= 1.7320508075688772935274463415059;		//	sqrt(3)
const double One_SQRT3			= 0.577350269;								// 
const double SQRT2				= 1.4142135623730950488016887242097;		//	sqrt(2)
const double SQRT_2over3		= 0.81649658092773;							//	sqrt(2/3)
const double SQRT_3over2		= 1.224744871;								//  sqrt(3/2)
const double SQRT_1over2		= 0.707106781;								//  sqrt(1/2)
const double SQRT3_2			= 0.86602540;								//	SQRT(3)/2
const double SIN_PI_6			= 0.5;										//	sin(pi/6)
const double COS_PI_6			= 0.86602540378444;							//	cos(pi/6)
// const double fline              = 60;
//----------------------------------------------------------------------------------------------------------------------------------
void inital_para(void);									// set all parameters to zero when stop
// current and voltage protection in DSP
const double I_peak_limit 		= 40;
const double V_peak_limit 		= 320;
const double V_dc_limit			= 500;
double f_line                   = 400;					//line frequency (Hz) for ac bus
//--------------------------------------variables for VSI Current and Voltage Loop Calculation------------------------------------------------------
const double Lout_VSI			= 1000e-6;
// const double Cout_VSI			= 30.8e-6;
const double R_Lout_VSI			= 110e-3;
// const double R_load_VSI			= 15;
double Vdc_VSI					= 270;
// power control
double Pref_VSI                 = 0;
double Qref_VSI                 = 0;
// power controller
double Kppq_VSI					= 0;
double Kipq_VSI					= 0;

// current loop:
double idref_VSI				= 0;
double idslop_VSI				= 0;
double iqref_VSI				= 0;
double iqslop_VSI				= 0;
// const double vdref_VSI 			= 99.6;
// double vdslop_VSI	  			= 0;
// const double vqref_VSI 			= 0;
          
const double Tsw_VSI			= 0.00005;
const double bandwidth_i_VSI 	= 1000;	
// const double bandwidth_v_VSI 	= 100;
double Kpi_VSI					= 0;
double Kii_VSI					= 0;
// double Kpv_VSI					= 0;
// double Kiv_VSI					= 0;
double K_decoup_VSI				= 0;
double delta_id_antiwind_VSI 	= 0;
double delta_iq_antiwind_VSI 	= 0;
// double delta_vd_antiwind_VSI 	= 0;
// double delta_vq_antiwind_VSI 	= 0;
//variables for Voltage Generation (Open loop or Close loop)
// int phase_counter_VSI 			= 0;

double Vd_out_VSI		  		= 0;					//reference voltage D value in DQ axis
double Vq_out_VSI		  		= 0;					//reference voltage Q value in DQ axis
double out_vector_alpha_VSI 	= 0;
double out_vector_beta_VSI		= 0;
double In_ud_VSI		  		= 0;
double In_uq_VSI 	 			= 0;
double In_idref_VSI 			= 0;
double In_iqref_VSI 			= 0;
void protection_VSI(void);
void I_loop_VSI(void);
void V_loop_VSI(void);
void PQ_loop_VSI(void);
void O_loop_VSI(void);
double M_Index_VSI 				= 0;

double Ia_VSI, Ib_VSI, Ic_VSI, Va_VSI, Vb_VSI, Vc_VSI, Vab_VSI, Vbc_VSI, Vca_VSI, Vno_VSI, Vpo_VSI, Per_i_1, Per_i_2;
double Ialpha_VSI, Ibeta_VSI, Id_VSI, Iq_VSI, Valpha_VSI, Vbeta_VSI, Vd_VSI, Vq_VSI; 
double P_VSI, Q_VSI;

// variables for VSI PLL
double DQ_theta_VSI                 = 0;					//reference voltage angle in DQ axis
double SIN_VP_theta_VSI             = 0;
double COS_VP_theta_VSI             = 0;
double SIN_VP_theta_VSI_ID          = 0;                    //angle for current and duty cycle dq transfermation
double COS_VP_theta_VSI_ID          = 0;                    //angle for current and duty cycle dq transfermation
double FW_VP_VSI[2]                 = {0,0};
// double phase_angle_PLL_VSI			= 0;
int PLL_Voltage_OK_VSI              = 0;
double counter_PLL_VSI              = 0;
const double counter_PLL_Max_VSI    = 600;
double Kp_PLL_VSI                   = 0.8921;//1;//
double Ki_PLL_VSI                   = 39.64;//2;//
double N                            = 0.000;// (0.015 is unstable 0.005 is stable for 1 and 2 pair)  // (0.0065 is unstable 0.005 is stable for 1.784 and 158.6 pair)
void PLL_VSI(void);

//--------------------------------------variables for AFE Current and Voltage Loop Calculation------------------------------------------------------
const double Lout_AFE			= 468e-6;
const double R_Lout_AFE			= 85e-3;
const double C_dc_AFE			= 100e-6;
const double R_load_AFE			= 96;
double Vdc_AFE					= 0;
// current loop:
double idref_AFE				= 0;
double idslop_AFE				= 0;
double iqref_AFE				= 0;
double iqslop_AFE				= 0;
const double vdcref_AFE 		= 270;
double vdcslop_AFE	  			= 134;
const double Tsw_AFE			= 0.00005;
const double bandwidth_i_AFE 	= 1000;
const double bandwidth_v_AFE 	= 100;
double Kpi_AFE					= 0;
double Kii_AFE					= 0;
double Kpv_AFE					= 0;
double Kiv_AFE					= 0;
double K_decoup_AFE				= 0;
double delta_id_antiwind_AFE 	= 0;
double delta_iq_antiwind_AFE 	= 0;
double delta_vdc_antiwind_AFE	= 0;
double Y                        = 0;
//variables for Voltage Generation (Open loop or Close loop)

// variables for AFE PLL
double DQ_theta_AFE                 = 0;					//reference voltage angle in DQ axis
double SIN_VP_theta_AFE             = 0;
double COS_VP_theta_AFE             = 0;
double FW_VP_AFE[2]                 = {0,0};
// double phase_angle_PLL_VSI			= 0;
int PLL_Voltage_OK_AFE              = 0;
double counter_PLL_AFE              = 0;
const double counter_PLL_Max_AFE    = 60000;
double Kp_PLL_AFE                   = 0.2;
double Ki_PLL_AFE                   = 1.0;
void PLL_AFE(void);

// double VLL_DQ_theta_AFE 		= 0;					//reference line to line voltage angle in DQ axis
// double SIN_VLL_theta_AFE	 	= 0;					//sin of reference line to line voltage angle in DQ axis
// double COS_VLL_theta_AFE		= 0;					//cos of reference line to line voltage angle in DQ axis
// double SIN_VPN_theta_AFE	 	= 0;					//sin of reference phase to neutral voltage angle in DQ axis
// double COS_VPN_theta_AFE		= 0;					//cos of reference phase to neutral voltage angle in DQ axis
double Vd_out_AFE		  		= 0;					//reference voltage D value in DQ axis
double Vq_out_AFE		  		= 0;					//reference voltage Q value in DQ axis
double Vd_out_AFE_final			= 0;
double Vq_out_AFE_final			= 0;
double out_vector_alpha_AFE 	= 0;
double out_vector_beta_AFE		= 0;
double Vout_max_AFE				= 0;
double Vout_vector_AFE			= 0;
double Vout_vector_AFE_final	= 0;
double In_ud_AFE		  		= 0;
double In_uq_AFE 	 			= 0;
double In_idref_AFE 			= 0;
const double One_Vpd     		= 0.010040160642570281124497991967871;		//1/Vsd: 1/99.6 when input phase to neutral voltage is 57.5
void protection_AFE(void);
void PLL(void);
void I_loop_AFE(void);
void V_loop_AFE(void);
double M_Index_AFE 				= 0;

double Ia_AFE, Ib_AFE, Ic_AFE, Va_AFE, Vb_AFE, Vc_AFE, Vab_AFE, Vbc_AFE, Vca_AFE, Vno_AFE, Vpo_AFE;
double Ialpha_AFE, Ibeta_AFE, Id_AFE, Iq_AFE, Valpha_AFE, Vbeta_AFE, Vd_AFE, Vq_AFE, Valpha_AFE, Vbeta_AFE, VLL_d_AFE, VLL_q_AFE;
// variables for SVM
double Ma_1						= 0;
double Mk_1						= 0;

double Ma_2						= 0;
double Mk_2						= 0;

double duration1_1 				= 0;					//for intermediate converter 1 duty cycle calculation
double duration2_1 				= 0;					//for intermediate converter 1 duty cycle calculation
double duration0_1 				= 0;					//for intermediate converter 1 duty cycle calculation

double duration1_2 				= 0;					//for intermediate converter 2 duty cycle calculation
double duration2_2 				= 0;					//for intermediate converter 2 duty cycle calculation
double duration0_2 				= 0;					//for intermediate converter 2 duty cycle calculation

int region_SVPWM_1 				= 0;					// location region of reference voltage for SVPWM
double region_angle_SVPWM_1 	= 0;					// angle in location region of reference voltage for SVPWM

int region_SVPWM_2 				= 0;					// location region of reference voltage for SVPWM
double region_angle_SVPWM_2 	= 0;					// angle in location region of reference voltage for SVPWM


int region_DPWM_1 				= 0;					// location region of reference voltage for DPWM
double region_angle_DPWM_1 		= 0;					// angle in location region of reference voltage for DPWM

int region_DPWM_2 = 0;								// location region of reference voltage for DPWM
double region_angle_DPWM_2 		= 0;					// angle in location region of reference voltage for DPWM

double DutyA_1 					= 0;					//phase A duty cycle of converter 1
double DutyB_1 					= 0;					//phase B duty cycle of converter 1
double DutyC_1 					= 0;					//phase C duty cycle of converter 1

double DutyA_2 					= 0;					//phase A duty cycle of converter 2
double DutyB_2 					= 0;					//phase B duty cycle of converter 2
double DutyC_2 					= 0;					//phase C duty cycle of converter 2

void SVM(void);


// variable for on off and reset control
int Run[2] = {1,1};
int Run_flag = 0;
int Reset[2] = {1,1};
int Reset_flag = 0;
//---------------------------------------test function and variables for AD conversion---------------------------------------------------
void Sensor_scal(void);
// variables for measuring transfer function and loop gain
double Per_1 = 0;
double Per_2 = 0;

//=======================================================Handle PWM Interrupt==============================================================================
// y1 is phase A duty cycle, y2 is phase B duty cycle, y3 is phase C duty cycle
// u1 is 
void epwm1_isr(double *y1, double *y2, double *y3, double *y4, double *y5, double *y6, double *y7, double *y8, double *y9, double *y10, double u1, double u2, double u3, double u4, double u5, double u6, double u7, double u8, double u9, double u10, double u11, double u12, double u13, double u14, double u15, double u16, double u17, double u18, double u19, double u20, double u21, double u22, double u23, double u24)
{
    

// for VSI from sensor 1
	Ia_VSI		= -u1;					// 0.5 is because of two turns of the current sensors
	Ib_VSI		= -u2;
	Ic_VSI		= -Ia_VSI-Ib_VSI;
    Vab_VSI		= -u3;
	Vca_VSI		= -u4;                   // going from Per_2
	Vbc_VSI		= -Vab_VSI-Vca_VSI;
	Vdc_VSI		= u5;
        
    Ia_AFE		= -u6;					// 0.5 is because of two turns of the current sensors
	Ib_AFE		= -u7;                  // ?????????Module????????DSP?????????????
	Ic_AFE		= -Ia_AFE-Ib_AFE;
    Vab_AFE		= u8;                   // ??????????
	Vca_AFE		= u9;                   // going from Per_2
	Vbc_AFE		= -Vab_AFE-Vca_AFE;
	Vdc_AFE		= u10;
    
    
    Kpi_VSI		= u11;
	Kii_VSI		= u12;
    Kppq_VSI	= u13;
    Kipq_VSI	= u14;
    
    Kp_PLL_VSI  = u15;//1;//
    Ki_PLL_VSI  = u16;//2;//
    
    
    Kpi_AFE		= u17;
	Kii_AFE		= u18;
	Kpv_AFE		= u19;
	Kiv_AFE		= u20;
    
    Kp_PLL_AFE  = u21;//1;//
    Ki_PLL_AFE  = u22;//2;//
    
//     Pref_VSI    = u23;
//     Qref_VSI    = u24;
    idref_VSI    = u23;
    iqref_VSI    = u24;
    
    PLL_VSI();
    if (PLL_Voltage_OK_VSI == 1)
    {
//         idref_VSI=(P*Vd_VSI+Q*Vq_VSI)/(Vd_VSI*Vd_VSI+Vq_VSI*Vq_VSI);
//         iqref_VSI=(P*Vq_VSI-Q*Vd_VSI)/(Vd_VSI*Vd_VSI+Vq_VSI*Vq_VSI);
//         PQ_loop_VSI();
        I_loop_VSI();
        *y4 = 1;

    }
    else
    {
        *y4 = 0;
    }
    
    PLL_AFE();
    if (PLL_Voltage_OK_AFE == 1)
    {
// 			GpioDataRegs.GPCSET.bit.GPIO85 = 1;
//         V_loop_AFE();
		idref_AFE = 7.3485;
        I_loop_AFE();
        *y8 = 1;

    }
    else
    {
//             I_loop_AFE();
        *y8 = 0;
    }
    SVM();
    *y1 =   DutyA_2;
    *y2 =   DutyB_2;
    *y3 =   DutyC_2;
    *y5 =   DutyA_1;
    *y6 =   DutyB_1;
    *y7 =   FW_VP_VSI[1];//DutyC_1;
    *y9 =   Vd_out_VSI;//P_VSI;
    *y10 =  Vq_out_VSI;//Q_VSI;
}


void Sensor_scal(void)
{
//  Ia_AFE = u1;					// be careful about sensor direction and reference sign for AFE
// 	Ib_AFE = u2;					// be careful about sensor direction and reference sign for AFE
// 	Ic_AFE = -Ia_AFE-Ib_AFE;
// 	Vab_AFE = u1;
// 	Vca_AFE = u2;					// going from Per_1
// 	Vbc_AFE = -Vab_AFE-Vca_AFE;
// 	Vdc_AFE = u1;
// // for VSI from sensor 1
// 	Ia_VSI		= u1;				// 0.5 is because of two turns of the current sensors
// 	Ib_VSI		= u2;
// 	Ic_VSI		= -Ia_VSI-Ib_VSI;
// 	Vab_VSI		= u1;
// 	Vca_VSI		= u2;               // going from Per_2
// 	Vbc_VSI		= -Vab_VSI-Vca_VSI;
// 	Vdc_VSI		= u1;
// 	
// // for loop gain measurement 	
// 	Per_1	= u1;
// 	Per_2	= u2;
}

//========================================================inital parameters===========================================================

void inital_para(void)
{
// 	vdslop_VSI  			= 0;
	idslop_VSI				= 0;
	iqslop_VSI  			= 0;
	delta_id_antiwind_VSI 	= 0;
	delta_iq_antiwind_VSI 	= 0;
// 	delta_vd_antiwind_VSI 	= 0;
// 	delta_vq_antiwind_VSI 	= 0;
	//variables for Voltage Generation (Open loop or Close loop)
    SIN_VP_theta_VSI		= 0;
    COS_VP_theta_VSI		= 0;
    FW_VP_VSI[0]            = 0;
    FW_VP_VSI[1]            = 0;
    PLL_Voltage_OK_VSI      = 0;
    counter_PLL_VSI         = 0;	
	DQ_theta_VSI 			= 0;				//reference voltage angle in DQ axis
	Vd_out_VSI		  		= 0;				//reference voltage D value in DQ axis
	Vq_out_VSI		  		= 0;				//reference voltage Q value in DQ axis
	out_vector_alpha_VSI 	= 0;
	out_vector_beta_VSI		= 0;
	In_ud_VSI		  		= 0;
	In_uq_VSI 	 			= 0;
	In_idref_VSI			= 0;
	In_iqref_VSI			= 0;
	idref_VSI				= 0;
	iqref_VSI				= 0;	
	M_Index_VSI 			= 0;
	
	vdcslop_AFE  			= 129;
	idslop_AFE				= 0;
	iqslop_AFE  			= 0;
	delta_id_antiwind_AFE 	= 0;
	delta_iq_antiwind_AFE 	= 0;
	Vout_max_AFE			= 0;
	//variables for Voltage Generation (Open loop or Close loop)
	SIN_VP_theta_AFE		= 0;
    COS_VP_theta_AFE		= 0;
    FW_VP_AFE[0]            = 0;
    FW_VP_AFE[1]            = 0;
    PLL_Voltage_OK_AFE      = 0;
    counter_PLL_AFE         = 0;	
	DQ_theta_AFE 			= 0;				//reference voltage angle in DQ axis
	
	Vd_out_AFE		  		= 0;				//reference voltage D value in DQ axis
	Vq_out_AFE		  		= 0;				//reference voltage Q value in DQ axis
	Vd_out_AFE_final		= 0;
	Vq_out_AFE_final		= 0;
	out_vector_alpha_AFE 	= 0;
	out_vector_beta_AFE		= 0;
	Vout_max_AFE			= 0;
	Vout_vector_AFE			= 0;
	Vout_vector_AFE_final	= 0;
	In_ud_AFE		  		= 0;
	In_uq_AFE 	 			= 0;
	In_idref_AFE			= 0;
	idref_AFE				= 0;
	M_Index_AFE 			= 0;
	
//	output_vector_angle_1 	= 0;					//reference voltage angle
//	output_vector_angle_2 	= 0;					//reference voltage angle
	
	//variables for motor control
	f_line					= 60;					//line frequency (Hz)
	// variables for PWM duty cycle calculation
//	M_Index_1 				= 0;					//Modulation Index for converter 1;
//	M_Index_2				= 0;					//Modulation Index for convertet 2;
	
	Ma_1					= 0;
	Mk_1					= 0;

	Ma_2					= 0;
	Mk_2					= 0;
	
	DutyA_1 				= 0;					//phase A duty cycle of converter 1
	DutyB_1 				= 0;					//phase B duty cycle of converter 1
	DutyC_1 				= 0;					//phase C duty cycle of converter 1
	
	DutyA_2 				= 0;					//phase A duty cycle of converter 2
	DutyB_2 				= 0;					//phase B duty cycle of converter 2
	DutyC_2 				= 0;					//phase C duty cycle of converter 2

}

void protection_VSI(void)
{
	
	if ((abs(Ia_VSI) > I_peak_limit)||(abs(Ib_VSI) > I_peak_limit)||(abs(Ic_VSI) > I_peak_limit))
	{
		Run_flag = 0;
	}
	if ((abs(Vab_VSI) > V_peak_limit)||(abs(Vbc_VSI) > V_peak_limit)||(abs(Vca_VSI) > V_peak_limit))
	{
		Run_flag = 0;
	}
}

// void V_loop_VSI(void)
// {
// 	Va_VSI=One_3*(Vab_VSI-Vca_VSI);
// 	Vb_VSI=One_3*(Vbc_VSI-Vab_VSI);
// 	Vc_VSI=One_3*(Vca_VSI-Vbc_VSI);
// 	Valpha_VSI = (Va_VSI-Vb_VSI*0.5-Vc_VSI*0.5)*SQRT_2over3;
// 	Vbeta_VSI = (SQRT3_2*Vb_VSI-SQRT3_2*Vc_VSI)*SQRT_2over3;
// 	Vd_VSI = COS_VP_theta_VSI*Valpha_VSI+SIN_VP_theta_VSI*Vbeta_VSI;
// 	Vq_VSI = -SIN_VP_theta_VSI*Valpha_VSI+COS_VP_theta_VSI*Vbeta_VSI;
// 	if (vdslop_VSI < vdref_VSI)
// 	{	
// 		vdslop_VSI = vdslop_VSI+1;
// 	}
// 	else 
// 	{
// 		vdslop_VSI = vdref_VSI;
// 	}
// 	In_idref_VSI=In_idref_VSI+Kiv_VSI*(vdslop_VSI-Vd_VSI-delta_vd_antiwind_VSI)*Tsw_VSI;
// 	In_iqref_VSI=In_iqref_VSI+Kiv_VSI*(vqref_VSI-Vq_VSI-delta_vq_antiwind_VSI)*Tsw_VSI;
// 	idref_VSI=Kpv_VSI*(vdslop_VSI-Vd_VSI)+In_idref_VSI;		
// 	iqref_VSI=Kpv_VSI*(vqref_VSI-Vq_VSI)+In_iqref_VSI;
// }

void PLL_VSI(void)
{
	Va_VSI          = One_3*(Vab_VSI-Vca_VSI);
	Vb_VSI          = One_3*(Vbc_VSI-Vab_VSI);
	Vc_VSI          = One_3*(Vca_VSI-Vbc_VSI);
	Valpha_VSI      = (Va_VSI-Vb_VSI*0.5-Vc_VSI*0.5)*SQRT_2over3;
	Vbeta_VSI       = (SQRT3_2*Vb_VSI-SQRT3_2*Vc_VSI)*SQRT_2over3;
	Vd_VSI          = COS_VP_theta_VSI*Valpha_VSI+SIN_VP_theta_VSI*Vbeta_VSI;
	Vq_VSI          = -SIN_VP_theta_VSI*Valpha_VSI+COS_VP_theta_VSI*Vbeta_VSI;
    
    Ialpha_VSI					= (Ia_VSI-Ib_VSI*0.5-Ic_VSI*0.5)*SQRT_2over3;
	Ibeta_VSI					= (SQRT3_2*Ib_VSI-SQRT3_2*Ic_VSI)*SQRT_2over3;
	Id_VSI						= COS_VP_theta_VSI_ID*Ialpha_VSI+SIN_VP_theta_VSI_ID*Ibeta_VSI;
	Iq_VSI						= -SIN_VP_theta_VSI_ID*Ialpha_VSI+COS_VP_theta_VSI_ID*Ibeta_VSI;
    
    
	FW_VP_VSI[0]    = FW_VP_VSI[0]+Ki_PLL_VSI*Vq_VSI*Tsw_VSI;
    FW_VP_VSI[1]    = Kp_PLL_VSI*Vq_VSI+FW_VP_VSI[0];

	DQ_theta_VSI 	= DQ_theta_VSI+(FW_VP_VSI[1]+f_line*Two_PI)*Tsw_VSI;
	
	DQ_theta_VSI 	= DQ_theta_VSI>Two_PI?(DQ_theta_VSI-Two_PI):DQ_theta_VSI;
	DQ_theta_VSI 	= DQ_theta_VSI<0?(DQ_theta_VSI+Two_PI):DQ_theta_VSI;


	COS_VP_theta_VSI		= cos(DQ_theta_VSI);
	SIN_VP_theta_VSI		= sin(DQ_theta_VSI);
    
    COS_VP_theta_VSI_ID		= cos(DQ_theta_VSI+N*FW_VP_VSI[1]);
	SIN_VP_theta_VSI_ID		= sin(DQ_theta_VSI+N*FW_VP_VSI[1]);
		
	if (PLL_Voltage_OK_VSI == 0)
	{
		counter_PLL_VSI++;
		if (counter_PLL_VSI > counter_PLL_Max_VSI)	
		{
			PLL_Voltage_OK_VSI = 1;
			counter_PLL_VSI = 0;
		}		
	}
	else
	{
		PLL_Voltage_OK_VSI = 1;
		counter_PLL_VSI = 0;
		
	}
}

void PQ_loop_VSI(void)
{
	P_VSI = Vd_VSI*Id_VSI + Vq_VSI*Iq_VSI;
	Q_VSI = Vd_VSI*Iq_VSI - Vq_VSI*Id_VSI ;
// 	if (vdslop_VSI < vdref_VSI)
// 	{	
// 		vdslop_VSI = vdslop_VSI+1;
// 	}
// 	else 
// 	{
// 		vdslop_VSI = vdref_VSI;
// 	}
	In_idref_VSI=In_idref_VSI+Kipq_VSI*(Pref_VSI-P_VSI)*Tsw_VSI;
	In_iqref_VSI=In_iqref_VSI+Kipq_VSI*(Qref_VSI-Q_VSI)*Tsw_VSI;
	idref_VSI=Kppq_VSI*(Pref_VSI-P_VSI)+In_idref_VSI;		
	iqref_VSI=Kppq_VSI*(Qref_VSI-Q_VSI)+In_iqref_VSI;
}


void I_loop_VSI(void)
{
// 	Ialpha_VSI					= (Ia_VSI-Ib_VSI*0.5-Ic_VSI*0.5)*SQRT_2over3;
// 	Ibeta_VSI					= (SQRT3_2*Ib_VSI-SQRT3_2*Ic_VSI)*SQRT_2over3;
// 	Id_VSI						= COS_VP_theta_VSI_ID*Ialpha_VSI+SIN_VP_theta_VSI_ID*Ibeta_VSI;
// 	Iq_VSI						= -SIN_VP_theta_VSI_ID*Ialpha_VSI+COS_VP_theta_VSI_ID*Ibeta_VSI;
// 	if (idslop_VSI < idref_VSI)
// 	{	
// 		idslop_VSI = idslop_VSI+0.01;
// // 		iqslop_VSI = iqslop_VSI+0.05;
// 	}
// 	else 
// 	{
// 		idslop_VSI = idref_VSI;
// // 		iqslop_VSI = iqref_VSI;
// 	}
	idslop_VSI				= idref_VSI;
	iqslop_VSI				= iqref_VSI;
	In_ud_VSI				= In_ud_VSI+Kii_VSI*(idslop_VSI-Id_VSI-delta_id_antiwind_VSI)*Tsw_VSI;
	In_uq_VSI				= In_uq_VSI+Kii_VSI*(iqslop_VSI-Iq_VSI-delta_iq_antiwind_VSI)*Tsw_VSI;
	Vd_out_VSI				= K_decoup_VSI*Two_PI*f_line*Lout_VSI*Iq_VSI+Kpi_VSI*(idslop_VSI-Id_VSI)+In_ud_VSI;		
	Vq_out_VSI				= -K_decoup_VSI*Two_PI*f_line*Lout_VSI*Id_VSI+Kpi_VSI*(iqslop_VSI-Iq_VSI)+In_uq_VSI;
	out_vector_alpha_VSI	= COS_VP_theta_VSI_ID*Vd_out_VSI-SIN_VP_theta_VSI_ID*Vq_out_VSI;
	out_vector_beta_VSI		= SIN_VP_theta_VSI_ID*Vd_out_VSI+COS_VP_theta_VSI_ID*Vq_out_VSI;
	Ma_2					= out_vector_alpha_VSI*SQRT_3over2;///Vdc_VSI;							// Valpha*sqrt(3/2)/Vdc
	Mk_2					= out_vector_beta_VSI*SQRT_1over2;///Vdc_VSI;							// Vbeta*sqrt(3/2)/Vdc/sqrt(3)
//	output_vector_angle_2	= atan2(out_vector_beta_VSI, out_vector_alpha_VSI);

//	if(output_vector_angle_2<0)
//	{
//		output_vector_angle_2=output_vector_angle_2+2*PI;
//	}					

}

void O_loop_VSI(void)
{
	out_vector_alpha_VSI	= COS_VP_theta_VSI*Vd_out_VSI-SIN_VP_theta_VSI*Vq_out_VSI;
	out_vector_beta_VSI		= SIN_VP_theta_VSI*Vd_out_VSI+COS_VP_theta_VSI*Vq_out_VSI;
	Ma_2					= out_vector_alpha_VSI*SQRT_3over2/Vdc_VSI;							// Valpha*sqrt(3/2)/Vdc
	Mk_2					= out_vector_beta_VSI*SQRT_1over2/Vdc_VSI;							// Vbeta*sqrt(3/2)/Vdc/sqrt(3)					
}

void protection_AFE(void)
{
	
	if ((abs(Ia_AFE) > I_peak_limit)||(abs(Ib_AFE) > I_peak_limit)||(abs(Ic_AFE) > I_peak_limit))
	{
		Run_flag = 0;
	}
	if ((abs(Vab_AFE) > V_peak_limit)||(abs(Vbc_AFE) > V_peak_limit)||(abs(Vca_AFE) > V_peak_limit)||Vdc_AFE > V_dc_limit)
	{
		Run_flag = 0;
	}
}

void PLL_AFE(void)
{
	Va_AFE          = One_3*(Vab_AFE-Vca_AFE);
	Vb_AFE          = One_3*(Vbc_AFE-Vab_AFE);
	Vc_AFE          = One_3*(Vca_AFE-Vbc_AFE);
	Valpha_AFE      = (Va_AFE-Vb_AFE*0.5-Vc_AFE*0.5)*SQRT_2over3;
	Vbeta_AFE       = (SQRT3_2*Vb_AFE-SQRT3_2*Vc_AFE)*SQRT_2over3;
	Vd_AFE          = COS_VP_theta_AFE*Valpha_AFE+SIN_VP_theta_AFE*Vbeta_AFE;
	Vq_AFE          = -SIN_VP_theta_AFE*Valpha_AFE+COS_VP_theta_AFE*Vbeta_AFE;    
    
	FW_VP_AFE[0]    = FW_VP_AFE[0]+Ki_PLL_AFE*Vq_AFE*Tsw_AFE;
    FW_VP_AFE[1]    = Kp_PLL_AFE*Vq_AFE+FW_VP_AFE[0];

	DQ_theta_AFE 	= DQ_theta_AFE+(FW_VP_AFE[1]+60*Two_PI)*Tsw_AFE;
	
	DQ_theta_AFE 	= DQ_theta_AFE>Two_PI?(DQ_theta_AFE-Two_PI):DQ_theta_AFE;
	DQ_theta_AFE 	= DQ_theta_AFE<0?(DQ_theta_AFE+Two_PI):DQ_theta_AFE;


	COS_VP_theta_AFE		= cos(DQ_theta_AFE);
	SIN_VP_theta_AFE		= sin(DQ_theta_AFE);
		
	if (PLL_Voltage_OK_AFE == 0)
	{
		counter_PLL_AFE++;
		if (counter_PLL_AFE > counter_PLL_Max_AFE)	
		{
			PLL_Voltage_OK_AFE = 1;
			counter_PLL_AFE = 0;
		}		
	}
	else
	{
		PLL_Voltage_OK_AFE = 1;
		counter_PLL_AFE = 0;
		
	}
}

void I_loop_AFE(void)
{
	Ialpha_AFE					= (Ia_AFE-Ib_AFE*0.5-Ic_AFE*0.5)*SQRT_2over3;
	Ibeta_AFE					= (SQRT3_2*Ib_AFE-SQRT3_2*Ic_AFE)*SQRT_2over3;
	Id_AFE						= COS_VP_theta_AFE*Ialpha_AFE+SIN_VP_theta_AFE*Ibeta_AFE;
	Iq_AFE						= -SIN_VP_theta_AFE*Ialpha_AFE+COS_VP_theta_AFE*Ibeta_AFE;
	if (idslop_AFE < idref_AFE)
	{	
		idslop_AFE = idslop_AFE+0.001;
	}
	else 
	{
		idslop_AFE = idref_AFE;
		iqslop_AFE = iqref_AFE;
	}
// 	idslop_AFE				= idref_AFE;
// 	iqslop_AFE				= iqref_AFE;
	In_ud_AFE				= In_ud_AFE+Kii_AFE*(idslop_AFE-Id_AFE-delta_id_antiwind_AFE)*Tsw_AFE;
	In_uq_AFE				= In_uq_AFE+Kii_AFE*(iqslop_AFE-Iq_AFE-delta_iq_antiwind_AFE)*Tsw_AFE;
	Vd_out_AFE				= K_decoup_AFE*Two_PI*f_line*Lout_AFE*Iq_AFE+Kpi_AFE*(idslop_AFE-Id_AFE)+In_ud_AFE;		
	Vq_out_AFE				= -K_decoup_AFE*Two_PI*f_line*Lout_AFE*Id_AFE+Kpi_AFE*(iqslop_AFE-Iq_AFE)+In_uq_AFE;
    
// 	Vd_out_AFE				= -0.343*270;		
// 	Vq_out_AFE				= -0.032*270;			
	out_vector_alpha_AFE	= COS_VP_theta_AFE*Vd_out_AFE-SIN_VP_theta_AFE*Vq_out_AFE;
	out_vector_beta_AFE		= SIN_VP_theta_AFE*Vd_out_AFE+COS_VP_theta_AFE*Vq_out_AFE;
	Ma_1					= out_vector_alpha_AFE*SQRT_3over2;///Vdc_AFE;							// Valpha*sqrt(3/2)/Vdc
	Mk_1					= out_vector_beta_AFE*SQRT_1over2;///Vdc_AFE;							// Vbeta*sqrt(3/2)/Vdc/sqrt(3)
//	output_vector_angle_1	= atan2(out_vector_beta_AFE, out_vector_alpha_AFE)+f_line*Tsw_AFE*Two_PI;
//	if(output_vector_angle_1<0)
//	{
//		output_vector_angle_1=output_vector_angle_1+2*PI;
//	}
}
void V_loop_AFE(void)
{
	if (vdcslop_AFE< vdcref_AFE)
	{	
		vdcslop_AFE=vdcslop_AFE+0.01;
	}
	else 
	{
		vdcslop_AFE=vdcref_AFE;
		
	}
	In_idref_AFE=In_idref_AFE+Kiv_AFE*(vdcslop_AFE-Vdc_AFE-delta_vdc_antiwind_AFE)*Tsw_AFE;
	idref_AFE=((vdcslop_AFE-Vdc_AFE)*Kpv_AFE+In_idref_AFE)*Vdc_AFE*One_Vpd;
				
	if (idref_AFE >= 200 )
	{
		Y=200;
		delta_vdc_antiwind_AFE=(idref_AFE-Y)*R_load_AFE/C_dc_AFE;
		
	}
	else
	{
		if (idref_AFE < -200 )
		{
			Y=-200;
			delta_vdc_antiwind_AFE=(idref_AFE-Y)*R_load_AFE/C_dc_AFE;
			
		}
		else
		{
		    Y=idref_AFE;
			delta_vdc_antiwind_AFE=(idref_AFE-Y)*R_load_AFE/C_dc_AFE;
			
		}
		
	}
}

void SVM(void)
{
	//----------------------------------------------duty cycle calculation based on reference voltage-------------------------------------
	// converter 1 dutycycle calculation
	switch (PWM_Mode)
	{
		case 1:						// center aligned continious PWM
		{
			/*deterimine the region and the phase_angle in the region for SVPWM*/
			/*     010    110
					|_  II _|
				 III  |___|	 I
			  011_______|________100
					   _|_
				 IV  _|	  |_ VI
					|   V   |
				   001     101
			*/
			
			// find the section
			if (Ma_1 > 0)
			{
				if (Mk_1 > 0)
				{
					if (Ma_1 > Mk_1)
					{
						region_SVPWM_1	= 1;
						duration1_1		= Ma_1-Mk_1;
						duration2_1		= 2*Mk_1;
					}
					else
					{
						region_SVPWM_1	= 2;
						duration1_1		= Ma_1+Mk_1;
						duration2_1		= -Ma_1+Mk_1;
					}
				}
				else
					{
						if (Ma_1 > -Mk_1)
						{
							region_SVPWM_1	= 6;
							duration1_1		= -2*Mk_1;
							duration2_1		= Ma_1+Mk_1;
						}
						else
						{
							region_SVPWM_1	= 5;
							duration1_1		= -Ma_1-Mk_1;
							duration2_1		= Ma_1-Mk_1;
						}
					}
			}
			else
			{
				if (Mk_1 > 0)
				{
					if (-Ma_1 > Mk_1)
					{
						region_SVPWM_1	= 3;
						duration1_1		= 2*Mk_1;
						duration2_1		= -Ma_1-Mk_1;
					}
					else
					{
						region_SVPWM_1	= 2;
						duration1_1		= Ma_1+Mk_1;
						duration2_1		= -Ma_1+Mk_1;
					}
				}
				else
				{
					if (Ma_1 < Mk_1)
					{
						region_SVPWM_1	= 4;
						duration1_1		= -Ma_1+Mk_1;
						duration2_1		= -2*Mk_1;
					}
					else
					{
						region_SVPWM_1	= 5;
						duration1_1		= -Ma_1-Mk_1;
						duration2_1		= Ma_1-Mk_1;
					}
				}
			}
			
			duration0_1=1-duration1_1-duration2_1;
			
			//calculate duty cycle based on location of output vector
			switch (region_SVPWM_1)
			{
				case 1:	
				{
					DutyA_1=1-duration0_1/2;
					DutyB_1=duration2_1+duration0_1/2;
					DutyC_1=duration0_1/2;
					break;
				}
				case 2:
				{
					DutyA_1=duration1_1+duration0_1/2;
					DutyB_1=1-duration0_1/2;
					DutyC_1=duration0_1/2;	
					break;			
				}
				case 3:
				{
					DutyA_1=duration0_1/2;
					DutyB_1=1-duration0_1/2;
					DutyC_1=duration2_1+duration0_1/2;
					break;				
				}
				case  4:
				{
					DutyA_1=duration0_1/2;
					DutyB_1=duration1_1+duration0_1/2;
					DutyC_1=1-duration0_1/2;
					break;				
				}
				case  5:
				{
					DutyA_1=duration2_1+duration0_1/2;
					DutyB_1=duration0_1/2;
					DutyC_1=1-duration0_1/2;
					break;				
				}
				case 6:
				{
					DutyA_1=1-duration0_1/2;
					DutyB_1=duration0_1/2;
					DutyC_1=duration1_1+duration0_1/2;
					break;				
				}
			}
			break;
		}
		case 2:					// discontinious PWM
		{
			/*deterimine the region and the phase_angle in the region for DPWM*/
			/*     010    110
					|_  II _|
				 III  |___|	 I
			  011_______|________100
					   _|_
				 IV  _|	  |_ VI
					|   V   |
				   001     101
			*/
//			region_DPWM_1=(int)((output_vector_angle_1)/PIover6)+1;
//			region_SVPWM_1=(int)((output_vector_angle_1)/PIover3)+1;
//			region_angle_DPWM_1=output_vector_angle_1-(region_SVPWM_1-1)*PIover3;
//			duration1_1=M_Index_1*sin((PIover3)-region_angle_DPWM_1)/SIN_PIover3;
//			duration2_1=M_Index_1*sin(region_angle_DPWM_1)/SIN_PIover3;
//			duration0_1=1-duration1_1-duration2_1;
		
			switch (region_DPWM_1)
				{
					case 1:
					{
						DutyA_1=1;
						DutyB_1=duration2_1+duration0_1;
						DutyC_1=duration0_1;
//						PWM_pattern_1=0;
						break;
					}
					case 2:
					{
						DutyA_1=1-duration0_1;
						DutyB_1=duration2_1;
						DutyC_1=0;
//						PWM_pattern_1=1;
						break;	
					}
					case 3:
					{
						DutyA_1=duration1_1;
						DutyB_1=1-duration0_1;
						DutyC_1=0;
//						PWM_pattern_1=1;
						break;
					}
					case 4:
					{
						DutyA_1=duration1_1+duration0_1;
						DutyB_1=1;
						DutyC_1=duration0_1;
//						PWM_pattern_1=0;
						break;	
					}				
					case 5:
					{
						DutyA_1=duration0_1;
						DutyB_1=1;
						DutyC_1=duration2_1+duration0_1;
//						PWM_pattern_1=0;
						break;
					}
					case 6:
					{
						DutyA_1=0;
						DutyB_1=1-duration0_1;
						DutyC_1=duration2_1;
//						PWM_pattern_1=1;
						break;	
					}
					case 7:
					{
						DutyA_1=0;
						DutyB_1=duration1_1;
						DutyC_1=1-duration0_1;
//						PWM_pattern_1=1;
						break;
					}
					case 8:
					{
						DutyA_1=duration0_1;
						DutyB_1=duration1_1+duration0_1;
						DutyC_1=1;
//						PWM_pattern_1=0;
						break;
					}
					case 9:
					{
						DutyA_1=duration2_1+duration0_1;
						DutyB_1=duration0_1;
						DutyC_1=1;
//						PWM_pattern_1=0;
						break;
					}
					case 10:
					{
						DutyA_1=duration2_1;
						DutyB_1=0;
						DutyC_1=1-duration0_1;
//						PWM_pattern_1=1;
						break;
					}			
					case 11:
					{
						DutyA_1=1-duration0_1;
						DutyB_1=0;
						DutyC_1=duration1_1;
//						PWM_pattern_1=1;
						break;
					}
					case 12:
					{
						DutyA_1=1;
						DutyB_1=duration0_1;
						DutyC_1=duration1_1+duration0_1;
//						PWM_pattern_1=0;
						break;
					}				
				}
			break;
		}
		default :				// others, disable PWM GPIO 85(DSPIO6)=1, GPIO 29,28(DSPIO 7,8) =1,1
		{
//			GpioDataRegs.GPCSET.bit.GPIO85 =1;
			//GpioDataRegs.GPASET.bit.GPIO29 =1;
			//GpioDataRegs.GPASET.bit.GPIO28 =1;	
			break;
		}
	}
	
	// converter 2 dutycycle calculation
	switch (PWM_Mode)
	{
		case 1:						// center aligned continious PWM
		{
			/*deterimine the region and the phase_angle in the region for SVPWM*/
			/*     010    110
					|_  II _|
				 III  |___|	 I
			  011_______|________100
					   _|_
				 IV  _|	  |_ VI
					|   V   |
				   001     101
			*/
			// find the section
			if (Ma_2 > 0)
			{
				if (Mk_2 > 0)
				{
					if (Ma_2 > Mk_2)
					{
						region_SVPWM_2	= 1;
						duration1_2		= Ma_2-Mk_2;
						duration2_2		= 2*Mk_2;
					}
					else
					{
						region_SVPWM_2	= 2;
						duration1_2		= Ma_2+Mk_2;
						duration2_2		= -Ma_2+Mk_2;
					}
				}
				else
					{
						if (Ma_2 > -Mk_2)
						{
							region_SVPWM_2	= 6;
							duration1_2		= -2*Mk_2;
							duration2_2		= Ma_2+Mk_2;
						}
						else
						{
							region_SVPWM_1	= 5;
							duration1_2		= -Ma_2-Mk_2;
							duration2_2		= Ma_2-Mk_2;
						}
					}
			}
			else
			{
				if (Mk_2 > 0)
				{
					if (-Ma_2 > Mk_2)
					{
						region_SVPWM_2	= 3;
						duration1_2		= 2*Mk_2;
						duration2_2		= -Ma_2-Mk_2;
					}
					else
					{
						region_SVPWM_2	= 2;
						duration1_2		= Ma_2+Mk_2;
						duration2_2		= -Ma_2+Mk_2;
					}
				}
				else
				{
					if (Ma_2 < Mk_2)
					{
						region_SVPWM_2	= 4;
						duration1_2		= -Ma_2+Mk_2;
						duration2_2		= -2*Mk_2;
					}
					else
					{
						region_SVPWM_2	= 5;
						duration1_2		= -Ma_2-Mk_2;
						duration2_2		= Ma_2-Mk_2;
					}
				}
			}
			duration0_2=1-duration1_2-duration2_2;
			//calculate duty cycle based on location of output vector
			switch (region_SVPWM_2)
			{
				case 1:	
				{
					DutyA_2=1-duration0_2/2;
					DutyB_2=duration2_2+duration0_2/2;
					DutyC_2=duration0_2/2;
					break;
				}
				case 2:
				{
					DutyA_2=duration1_2+duration0_2/2;
					DutyB_2=1-duration0_2/2;
					DutyC_2=duration0_2/2;	
					break;			
				}
				case 3:
				{
					DutyA_2=duration0_2/2;
					DutyB_2=1-duration0_2/2;
					DutyC_2=duration2_2+duration0_2/2;
					break;				
				}
				case  4:
				{
					DutyA_2=duration0_2/2;
					DutyB_2=duration1_2+duration0_2/2;
					DutyC_2=1-duration0_2/2;
					break;				
				}
				case  5:
				{
					DutyA_2=duration2_2+duration0_2/2;
					DutyB_2=duration0_2/2;
					DutyC_2=1-duration0_2/2;
					break;				
				}
				case 6:
				{
					DutyA_2=1-duration0_2/2;
					DutyB_2=duration0_2/2;
					DutyC_2=duration1_2+duration0_2/2;
					break;				
				}
			}
			break;
		}
		case 2:					// discontinious PWM
		{
			/*deterimine the region and the phase_angle in the region for DPWM*/
			/*     010    110
					|_  II _|
				 III  |___|	 I
			  011_______|________100
					   _|_
				 IV  _|	  |_ VI
					|   V   |
				   001     101
			*/
//			region_DPWM_2=(int)((output_vector_angle_2)/PIover6)+1;
//			region_SVPWM_2=(int)((output_vector_angle_2)/PIover3)+1;
//			region_angle_DPWM_2=output_vector_angle_2-(region_SVPWM_2-1)*PIover3;
//			duration1_2=M_Index_2*sin((PIover3)-region_angle_DPWM_2)/SIN_PIover3;
//			duration2_2=M_Index_2*sin(region_angle_DPWM_2)/SIN_PIover3;
//			duration0_2=1-duration1_2-duration2_2;
		
			switch (region_DPWM_2)
				{
					case 1:
					{
						DutyA_2=1;
						DutyB_2=duration2_2+duration0_2;
						DutyC_2=duration0_2;
//						PWM_pattern_1=0;
						break;
					}
					case 2:
					{
						DutyA_2=1-duration0_2;
						DutyB_2=duration2_2;
						DutyC_2=0;
//						PWM_pattern_1=1;
						break;	
					}
					case 3:
					{
						DutyA_2=duration1_2;
						DutyB_2=1-duration0_2;
						DutyC_2=0;
//						PWM_pattern_1=1;
						break;
					}
					case 4:
					{
						DutyA_2=duration1_2+duration0_2;
						DutyB_2=1;
						DutyC_2=duration0_2;
//						PWM_pattern_1=0;
						break;	
					}				
					case 5:
					{
						DutyA_2=duration0_2;
						DutyB_2=1;
						DutyC_2=duration2_2+duration0_2;
//						PWM_pattern_1=0;
						break;
					}
					case 6:
					{
						DutyA_2=0;
						DutyB_2=1-duration0_2;
						DutyC_2=duration2_2;
//						PWM_pattern_1=1;
						break;	
					}
					case 7:
					{
						DutyA_2=0;
						DutyB_2=duration1_2;
						DutyC_2=1-duration0_2;
//						PWM_pattern_1=1;
						break;
					}
					case 8:
					{
						DutyA_2=duration0_2;
						DutyB_2=duration1_2+duration0_2;
						DutyC_2=1;
//						PWM_pattern_1=0;
						break;
					}
					case 9:
					{
						DutyA_2=duration2_2+duration0_2;
						DutyB_2=duration0_2;
						DutyC_2=1;
//						PWM_pattern_1=0;
						break;
					}
					case 10:
					{
						DutyA_2=duration2_2;
						DutyB_2=0;
						DutyC_2=1-duration0_2;
//						PWM_pattern_1=1;
						break;
					}			
					case 11:
					{
						DutyA_2=1-duration0_2;
						DutyB_2=0;
						DutyC_2=duration1_2;
//						PWM_pattern_1=1;
						break;
					}
					case 12:
					{
						DutyA_2=1;
						DutyB_2=duration0_2;
						DutyC_2=duration1_2+duration0_2;
//						PWM_pattern_1=0;
						break;
					}				
				}
			break;
		}
		default :				// others, disable PWM GPIO 85(DSPIO6)=1, GPIO 29,28(DSPIO 7,8) =1,1
		{
			//GpioDataRegs.GPCSET.bit.GPIO85 =1;
			//GpioDataRegs.GPASET.bit.GPIO29 =1;
			//GpioDataRegs.GPASET.bit.GPIO28 =1;	
			break;
		}
	}	
}
//===========================================================================
// No more.
//===========================================================================
