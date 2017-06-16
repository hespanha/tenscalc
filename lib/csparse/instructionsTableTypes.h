typedef enum instructionType_e {
    I_set=1,
    I_load=2,
    I_sum=3,
    I_sumprod=4,
    I_inv=5,
    I_minus_inv_sqr=6,
    I_2_inv_cube=7,
    I_div=8,
    I_minus_dot=9,
    I_minus_dot_div=10,
    I_plus_minus_dot=11,
    I_plus_sqr=12,
    I_plus_minus_dot_div=13,
    I_min=14,
    I_min0=15,
    I_max=16,
    I_max0=17,
    I_max_abs=18,
    I_clp=19,
    I_exp=20,
    I_log=21,
    I_cos=22,
    I_minus_cos=23,
    I_sin=24,
    I_minus_sin=25,
    I_abs=26,
    I_sqrt=27,
    I_Dsqrt=28,
    I_DDsqrt=29,
    I_sqr=30,
    I_2times=31,
    I_2=32,
    I_cube=33,
    I_3sqr=34,
    I_6times=35,
    I_atan=36,
    I_Datan=37,
    I_DDatan=38,
    I_luS2A=39,
    I_luS2Asym=40,
    I_mldivideA2F1=41,
    I_mldivideA2Fn=42,
    I_Mnorm2=43,
    I_Mnorminf=44,
    I_Mplus=45,
    I_Mmtimes=46,
    I_Mtimes=47,
    I_Msum=48,
    I_Mmin=49,
    I_Mones=50,
    I_Mzeros=51,
    I_Meye=52,
    I_Mdiag=53,
    I_Mclp=54,
    I_Mctranspose=55,
    I_Msubsref=56,
    I_Mtprod=57,
    I_Mtprod_matlab=58,
    I_Mfull=59,
    I_Mreshape=60,
    I_Mrepmat=61,
    I_Mcat=62,
    I_Mlu=63,
    I_Mldl=64,
    I_Mchol=65,
    I_Mmldivide_l1=66,
    I_Mmldivide_u=67,
    I_Mmldivide_u1=68,
    I_Mmldivide_d=69,
    I_Mrdivide=70,
    I_Mcompose=71,
} instructionType_t;
