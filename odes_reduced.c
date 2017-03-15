/*=======================================================================

 User-supplied f (ODEs) routine. Compute the function ydot = f(t,y)

-----------------------------------------------------------------------*/

/* Header files with descriptions of the contents used */

#include <math.h>                    /* Standard math functions */
#include <cvode/cvode.h>             /* CVODE functions and constants */
#include <cvode/cvode_dense.h>       /* Prototype for CVDense solver */
#include <nvector/nvector_serial.h>  /* Serial N_Vector types, functions, macros */
#include <sundials/sundials_dense.h> /* Definition of type DlsMat (dense matrix) */
#include <sundials/sundials_types.h> /* Definition of type realtype */

/*-----------------------------------------------------------------------*/

/* Type definition for user-supplied data passed to the solver functions */

typedef struct {
  realtype *rate, n_H, T_g, x_e;
} *User_Data;

/*-----------------------------------------------------------------------*/

int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *x, *ode, *rate;
  realtype n_H, x_e, loss, form;
  User_Data data;

  /* Obtain pointers to the y and ydot vector data arrays */
  x = NV_DATA_S(y);
  ode = NV_DATA_S(ydot);

  /* Retrieve the array of reaction rate coefficients and
   * the total number density from the user-supplied data */
  data = (User_Data) user_data;
  rate = data->rate;
  n_H = data->n_H;

  /* The electron abundance is a conserved quantity, given by the sum
   * of the abundances of all ionized species in the chemical network */
  x_e = x[0]+x[1]+x[3]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10]+x[12]+x[13]+x[14]+x[15]+x[18]+x[19]+x[22]+x[23]+x[26];

  /* The ODEs created by MakeRates begin here... */
  loss = -rate[67]*x[24]*n_H-rate[68]*x[20]*n_H-rate[69]*x[16]*n_H-rate[70]*x[21]*n_H-rate[71]*x[29]*n_H-rate[72]*x[11]*n_H-rate[73]*x[28]*n_H-rate[74]*x[17]*n_H-rate[75]*x[2]*n_H-rate[76]*x[27]*n_H-rate[216]*x_e*n_H-rate[217]*x_e*n_H-rate[266]-rate[267];
  form = +rate[48]*x[5]*x[30]*n_H;
  ode[0] = form+x[0]*loss;
  loss = -rate[49]*x[30]*n_H-rate[77]*x[20]*n_H-rate[78]*x[16]*n_H-rate[79]*x[16]*n_H-rate[80]*x[21]*n_H-rate[81]*x[11]*n_H-rate[82]*x[11]*n_H-rate[83]*x[11]*n_H-rate[84]*x[11]*n_H-rate[85]*x[28]*n_H-rate[86]*x[17]*n_H-rate[87]*x[17]*n_H-rate[88]*x[27]*n_H-rate[89]*x[27]*n_H-rate[90]*x[27]*n_H-rate[91]*x[4]*n_H-rate[155]*x[31]*n_H-rate[156]*x[31]*n_H-rate[169]*x[30]*n_H-rate[177]*x[24]*n_H-rate[178]*x[20]*n_H-rate[179]*x[11]*n_H-rate[180]*x[17]*n_H-rate[181]*x[4]*n_H-rate[242]*x_e*n_H;
  form = +rate[296]*x[25];
  ode[1] = form+x[1]*loss;
  loss = -rate[75]*x[0]*n_H-rate[147]*x[3]*n_H-rate[166]*x[18]*n_H-rate[184]*x[10]*n_H-rate[190]*x[26]*n_H-rate[198]*x[19]*n_H-rate[210]*x[12]*n_H-rate[213]*x[15]*n_H-rate[288]-rate[311];
  form = +rate[248]*x[9]*x_e*n_H;
  ode[2] = form+x[2]*loss;
  loss = -rate[47]*x[31]*n_H-rate[95]*x[24]*n_H-rate[107]*x[20]*n_H-rate[119]*x[16]*n_H-rate[130]*x[29]*n_H-rate[140]*x[28]*n_H-rate[146]*x[17]*n_H-rate[147]*x[2]*n_H-rate[148]*x[27]*n_H-rate[229]*x_e*n_H-rate[230]*x_e*n_H;
  form = +rate[58]*x[5]*x[11]*n_H+rate[59]*x[6]*x[30]*n_H+rate[72]*x[0]*x[11]*n_H+rate[132]*x[6]*x[11]*n_H+rate[133]*x[8]*x[11]*n_H+rate[258]*x[19]*x[30]*n_H;
  ode[3] = form+x[3]*loss;
  loss = -rate[7]*x[31]*n_H-rate[14]*x[30]*n_H-rate[18]*x[24]*n_H-rate[23]*x[20]*n_H-rate[31]*x[16]*n_H-rate[91]*x[1]*n_H-rate[100]*x[10]*n_H-rate[101]*x[10]*n_H-rate[114]*x[26]*n_H-rate[115]*x[26]*n_H-rate[124]*x[23]*n_H-rate[168]*x[18]*n_H-rate[176]*x[5]*n_H-rate[181]*x[1]*n_H-rate[204]*x[7]*n_H-rate[206]*x[6]*n_H-rate[209]*x[8]*n_H-rate[212]*x[12]*n_H-rate[214]*x[14]*n_H-rate[290]-rate[291]-rate[313]-rate[314]-rate[319]*x[31]*n_H-rate[325]*x[30]*n_H;
  form = +rate[37]*x[28]*x[29]*n_H+rate[186]*x[15]*x[24]*n_H+rate[192]*x[15]*x[20]*n_H+rate[197]*x[15]*x[16]*n_H+rate[213]*x[2]*x[15]*n_H+rate[264]*x[29]*x[29]*n_H;
  ode[4] = form+x[4]*loss;
  loss = -rate[48]*x[30]*n_H-rate[50]*x[24]*n_H-rate[52]*x[20]*n_H-rate[54]*x[16]*n_H-rate[56]*x[29]*n_H-rate[58]*x[11]*n_H-rate[60]*x[11]*n_H-rate[61]*x[28]*n_H-rate[63]*x[17]*n_H-rate[65]*x[27]*n_H-rate[154]*x[31]*n_H-rate[170]*x[20]*n_H-rate[171]*x[16]*n_H-rate[172]*x[11]*n_H-rate[173]*x[28]*n_H-rate[174]*x[17]*n_H-rate[175]*x[27]*n_H-rate[176]*x[4]*n_H-rate[215]*x_e*n_H-rate[265];
  form = +rate[153]*x[18]*x[30]*n_H+rate[169]*x[1]*x[30]*n_H+rate[249]*x[18]*x[31]*n_H+rate[250]*x[18]*x[31]*n_H+rate[266]*x[0]+rate[295]*x[30];
  ode[5] = form+x[5]*loss;
  loss = -rate[46]*x[31]*n_H-rate[59]*x[30]*n_H-rate[127]*x[29]*n_H-rate[132]*x[11]*n_H-rate[135]*x[17]*n_H-rate[137]*x[27]*n_H-rate[206]*x[4]*n_H-rate[226]*x_e*n_H-rate[227]*x_e*n_H;
  form = +rate[47]*x[3]*x[31]*n_H+rate[70]*x[0]*x[21]*n_H+rate[163]*x[11]*x[18]*n_H+rate[172]*x[5]*x[11]*n_H+rate[179]*x[1]*x[11]*n_H+rate[199]*x[7]*x[11]*n_H+rate[205]*x[11]*x[14]*n_H;
  ode[6] = form+x[6]*loss;
  loss = -rate[57]*x[30]*n_H-rate[104]*x[20]*n_H-rate[126]*x[11]*n_H-rate[128]*x[28]*n_H-rate[160]*x[31]*n_H-rate[187]*x[20]*n_H-rate[193]*x[16]*n_H-rate[199]*x[11]*n_H-rate[200]*x[28]*n_H-rate[201]*x[17]*n_H-rate[203]*x[27]*n_H-rate[204]*x[4]*n_H-rate[247]*x_e*n_H-rate[259]*x[24]*n_H;
  form = +rate[85]*x[1]*x[28]*n_H+rate[90]*x[1]*x[27]*n_H+rate[91]*x[1]*x[4]*n_H+rate[101]*x[4]*x[10]*n_H+rate[161]*x[18]*x[29]*n_H+rate[162]*x[18]*x[29]*n_H+rate[202]*x[14]*x[29]*n_H+rate[298]*x[29];
  ode[7] = form+x[7]*loss;
  loss = -rate[62]*x[30]*n_H-rate[92]*x[24]*n_H-rate[105]*x[20]*n_H-rate[118]*x[16]*n_H-rate[129]*x[29]*n_H-rate[133]*x[11]*n_H-rate[134]*x[11]*n_H-rate[139]*x[28]*n_H-rate[141]*x[17]*n_H-rate[143]*x[27]*n_H-rate[188]*x[20]*n_H-rate[194]*x[16]*n_H-rate[207]*x[17]*n_H-rate[209]*x[4]*n_H-rate[228]*x_e*n_H-rate[285];
  form = +rate[56]*x[5]*x[29]*n_H+rate[57]*x[7]*x[30]*n_H+rate[71]*x[0]*x[29]*n_H+rate[87]*x[1]*x[17]*n_H+rate[164]*x[18]*x[28]*n_H+rate[173]*x[5]*x[28]*n_H+rate[200]*x[7]*x[28]*n_H+rate[208]*x[14]*x[28]*n_H+rate[283]*x[28];
  ode[8] = form+x[8]*loss;
  loss = -rate[248]*x_e*n_H;
  form = +rate[75]*x[0]*x[2]*n_H+rate[147]*x[2]*x[3]*n_H+rate[166]*x[2]*x[18]*n_H+rate[184]*x[2]*x[10]*n_H+rate[190]*x[2]*x[26]*n_H+rate[198]*x[2]*x[19]*n_H+rate[210]*x[2]*x[12]*n_H+rate[213]*x[2]*x[15]*n_H+rate[288]*x[2]+rate[311]*x[2];
  ode[9] = form+x[9]*loss;
  loss = -rate[51]*x[30]*n_H-rate[93]*x[28]*n_H-rate[94]*x[28]*n_H-rate[97]*x[17]*n_H-rate[100]*x[4]*n_H-rate[101]*x[4]*n_H-rate[182]*x[20]*n_H-rate[183]*x[16]*n_H-rate[184]*x[2]*n_H-rate[243]*x_e*n_H-rate[244]*x_e*n_H-rate[245]*x_e*n_H-rate[252]*x[31]*n_H-rate[256]*x[30]*n_H-rate[262]*x[29]*n_H-rate[263]*x[29]*n_H;
  form = +rate[41]*x[26]*x[31]*n_H+rate[77]*x[1]*x[20]*n_H+rate[78]*x[1]*x[16]*n_H+rate[88]*x[1]*x[27]*n_H+rate[89]*x[1]*x[27]*n_H+rate[177]*x[1]*x[24]*n_H+rate[185]*x[14]*x[24]*n_H+rate[186]*x[15]*x[24]*n_H+rate[268]*x[24]+rate[297]*x[24]+rate[300]*x[24]+rate[302]*x[26];
  ode[10] = form+x[10]*loss;
  loss = -rate[3]*x[31]*n_H-rate[22]*x[20]*n_H-rate[28]*x[16]*n_H-rate[36]*x[29]*n_H-rate[39]*x[28]*n_H-rate[45]*x[18]*n_H-rate[58]*x[5]*n_H-rate[60]*x[5]*n_H-rate[72]*x[0]*n_H-rate[81]*x[1]*n_H-rate[82]*x[1]*n_H-rate[83]*x[1]*n_H-rate[84]*x[1]*n_H-rate[126]*x[7]*n_H-rate[132]*x[6]*n_H-rate[133]*x[8]*n_H-rate[134]*x[8]*n_H-rate[136]*x[12]*n_H-rate[138]*x[14]*n_H-rate[163]*x[18]*n_H-rate[172]*x[5]*n_H-rate[179]*x[1]*n_H-rate[199]*x[7]*n_H-rate[205]*x[14]*n_H-rate[280]-rate[281]-rate[282]-rate[308];
  form = +rate[11]*x[21]*x[30]*n_H+rate[32]*x[21]*x[21]*n_H+rate[33]*x[21]*x[28]*n_H+rate[35]*x[17]*x[21]*n_H+rate[95]*x[3]*x[24]*n_H+rate[107]*x[3]*x[20]*n_H+rate[119]*x[3]*x[16]*n_H+rate[140]*x[3]*x[28]*n_H+rate[146]*x[3]*x[17]*n_H+rate[147]*x[2]*x[3]*n_H+rate[148]*x[3]*x[27]*n_H+rate[206]*x[4]*x[6]*n_H+rate[230]*x[3]*x_e*n_H;
  ode[11] = form+x[11]*loss;
  loss = -rate[64]*x[30]*n_H-rate[96]*x[24]*n_H-rate[109]*x[20]*n_H-rate[120]*x[16]*n_H-rate[131]*x[29]*n_H-rate[136]*x[11]*n_H-rate[142]*x[28]*n_H-rate[149]*x[17]*n_H-rate[150]*x[27]*n_H-rate[189]*x[20]*n_H-rate[195]*x[16]*n_H-rate[210]*x[2]*n_H-rate[212]*x[4]*n_H-rate[231]*x_e*n_H-rate[232]*x_e*n_H-rate[233]*x_e*n_H;
  form = +rate[61]*x[5]*x[28]*n_H+rate[62]*x[8]*x[30]*n_H+rate[73]*x[0]*x[28]*n_H+rate[139]*x[8]*x[28]*n_H+rate[140]*x[3]*x[28]*n_H+rate[145]*x[22]*x[28]*n_H+rate[165]*x[17]*x[18]*n_H+rate[174]*x[5]*x[17]*n_H+rate[180]*x[1]*x[17]*n_H+rate[201]*x[7]*x[17]*n_H+rate[207]*x[8]*x[17]*n_H+rate[211]*x[14]*x[17]*n_H+rate[287]*x[17];
  ode[12] = form+x[12]*loss;
  loss = -rate[98]*x[24]*n_H-rate[111]*x[20]*n_H-rate[121]*x[16]*n_H-rate[234]*x_e*n_H-rate[235]*x_e*n_H-rate[236]*x_e*n_H-rate[237]*x_e*n_H;
  form = +rate[63]*x[5]*x[17]*n_H+rate[64]*x[12]*x[30]*n_H+rate[74]*x[0]*x[17]*n_H+rate[108]*x[17]*x[26]*n_H+rate[130]*x[3]*x[29]*n_H+rate[134]*x[8]*x[11]*n_H+rate[135]*x[6]*x[17]*n_H+rate[136]*x[11]*x[12]*n_H+rate[141]*x[8]*x[17]*n_H+rate[142]*x[12]*x[28]*n_H+rate[146]*x[3]*x[17]*n_H+rate[149]*x[12]*x[17]*n_H+rate[152]*x[17]*x[22]*n_H;
  ode[13] = form+x[13]*loss;
  loss = -rate[66]*x[30]*n_H-rate[112]*x[20]*n_H-rate[122]*x[16]*n_H-rate[138]*x[11]*n_H-rate[144]*x[28]*n_H-rate[151]*x[17]*n_H-rate[167]*x[31]*n_H-rate[185]*x[24]*n_H-rate[191]*x[20]*n_H-rate[196]*x[16]*n_H-rate[202]*x[29]*n_H-rate[205]*x[11]*n_H-rate[208]*x[28]*n_H-rate[211]*x[17]*n_H-rate[214]*x[4]*n_H-rate[238]*x_e*n_H;
  form = +rate[93]*x[10]*x[28]*n_H+rate[100]*x[4]*x[10]*n_H+rate[102]*x[15]*x[24]*n_H+rate[103]*x[26]*x[29]*n_H+rate[104]*x[7]*x[20]*n_H+rate[106]*x[26]*x[28]*n_H+rate[114]*x[4]*x[26]*n_H+rate[175]*x[5]*x[27]*n_H+rate[203]*x[7]*x[27]*n_H+rate[259]*x[7]*x[24]*n_H+rate[262]*x[10]*x[29]*n_H+rate[263]*x[10]*x[29]*n_H+rate[299]*x[27];
  ode[14] = form+x[14]*loss;
  loss = -rate[102]*x[24]*n_H-rate[116]*x[20]*n_H-rate[186]*x[24]*n_H-rate[192]*x[20]*n_H-rate[197]*x[16]*n_H-rate[213]*x[2]*n_H-rate[240]*x_e*n_H;
  form = +rate[128]*x[7]*x[28]*n_H+rate[129]*x[8]*x[29]*n_H+rate[131]*x[12]*x[29]*n_H+rate[168]*x[4]*x[18]*n_H+rate[176]*x[4]*x[5]*n_H+rate[181]*x[1]*x[4]*n_H+rate[204]*x[4]*x[7]*n_H+rate[206]*x[4]*x[6]*n_H+rate[209]*x[4]*x[8]*n_H+rate[212]*x[4]*x[12]*n_H+rate[214]*x[4]*x[14]*n_H+rate[291]*x[4]+rate[314]*x[4];
  ode[15] = form+x[15]*loss;
  loss = -rate[1]*x[31]*n_H-rate[10]*x[30]*n_H-rate[15]*x[24]*n_H-2*rate[24]*x[16]*n_H-rate[25]*x[29]*n_H-rate[26]*x[29]*n_H-rate[27]*x[29]*n_H-rate[28]*x[11]*n_H-rate[29]*x[28]*n_H-rate[30]*x[28]*n_H-rate[31]*x[4]*n_H-rate[42]*x[18]*n_H-rate[54]*x[5]*n_H-rate[69]*x[0]*n_H-rate[78]*x[1]*n_H-rate[79]*x[1]*n_H-rate[118]*x[8]*n_H-rate[119]*x[3]*n_H-rate[120]*x[12]*n_H-rate[121]*x[13]*n_H-rate[122]*x[14]*n_H-rate[123]*x[22]*n_H-rate[158]*x[18]*n_H-rate[171]*x[5]*n_H-rate[183]*x[10]*n_H-rate[193]*x[7]*n_H-rate[194]*x[8]*n_H-rate[195]*x[12]*n_H-rate[196]*x[14]*n_H-rate[197]*x[15]*n_H-rate[272]-rate[273]-rate[303]-rate[304];
  form = +rate[2]*x[21]*x[31]*n_H+rate[9]*x[20]*x[30]*n_H+rate[22]*x[11]*x[20]*n_H+rate[32]*x[21]*x[21]*n_H+rate[34]*x[21]*x[28]*n_H+rate[130]*x[3]*x[29]*n_H+rate[134]*x[8]*x[11]*n_H+rate[222]*x[19]*x_e*n_H+rate[227]*x[6]*x_e*n_H+rate[255]*x[24]*x[30]*n_H+rate[276]*x[21]+rate[281]*x[11]+rate[306]*x[21]+rate[308]*x[11];
  ode[16] = form+x[16]*loss;
  loss = -rate[5]*x[31]*n_H-rate[35]*x[21]*n_H-rate[38]*x[29]*n_H-rate[63]*x[5]*n_H-rate[74]*x[0]*n_H-rate[86]*x[1]*n_H-rate[87]*x[1]*n_H-rate[97]*x[10]*n_H-rate[108]*x[26]*n_H-rate[110]*x[26]*n_H-rate[135]*x[6]*n_H-rate[141]*x[8]*n_H-rate[146]*x[3]*n_H-rate[149]*x[12]*n_H-rate[151]*x[14]*n_H-rate[152]*x[22]*n_H-rate[165]*x[18]*n_H-rate[174]*x[5]*n_H-rate[180]*x[1]*n_H-rate[201]*x[7]*n_H-rate[207]*x[8]*n_H-rate[211]*x[14]*n_H-rate[286]-rate[287]-rate[310]-rate[318]*x[31]*n_H-rate[324]*x[30]*n_H;
  form = +rate[13]*x[28]*x[30]*n_H+rate[30]*x[16]*x[28]*n_H+rate[31]*x[4]*x[16]*n_H+rate[34]*x[21]*x[28]*n_H+rate[39]*x[11]*x[28]*n_H+rate[40]*x[28]*x[28]*n_H+rate[111]*x[13]*x[20]*n_H+rate[121]*x[13]*x[16]*n_H+rate[189]*x[12]*x[20]*n_H+rate[195]*x[12]*x[16]*n_H+rate[210]*x[2]*x[12]*n_H+rate[212]*x[4]*x[12]*n_H+rate[237]*x[13]*x_e*n_H+rate[254]*x[28]*x[31]*n_H;
  ode[17] = form+x[17]*loss;
  loss = -rate[42]*x[16]*n_H-rate[45]*x[11]*n_H-rate[153]*x[30]*n_H-rate[157]*x[20]*n_H-rate[158]*x[16]*n_H-rate[159]*x[21]*n_H-rate[161]*x[29]*n_H-rate[162]*x[29]*n_H-rate[163]*x[11]*n_H-rate[164]*x[28]*n_H-rate[165]*x[17]*n_H-rate[166]*x[2]*n_H-rate[168]*x[4]*n_H-rate[241]*x_e*n_H-rate[249]*x[31]*n_H-rate[250]*x[31]*n_H;
  form = +rate[49]*x[1]*x[30]*n_H+rate[83]*x[1]*x[11]*n_H+rate[86]*x[1]*x[17]*n_H+rate[94]*x[10]*x[28]*n_H+rate[154]*x[5]*x[31]*n_H+rate[155]*x[1]*x[31]*n_H+rate[156]*x[1]*x[31]*n_H+rate[160]*x[7]*x[31]*n_H+rate[167]*x[14]*x[31]*n_H+rate[265]*x[5]+rate[267]*x[0]+rate[271]*x[26]+rate[285]*x[8]+rate[292]*x[31]+rate[293]*x[30];
  ode[18] = form+x[18]*loss;
  loss = -rate[44]*x[31]*n_H-rate[125]*x[29]*n_H-rate[198]*x[2]*n_H-rate[222]*x_e*n_H-rate[223]*x_e*n_H-rate[224]*x_e*n_H-rate[225]*x_e*n_H-rate[246]*x_e*n_H-rate[258]*x[30]*n_H-rate[277]-rate[279];
  form = +rate[45]*x[11]*x[18]*n_H+rate[46]*x[6]*x[31]*n_H+rate[54]*x[5]*x[16]*n_H+rate[55]*x[23]*x[30]*n_H+rate[60]*x[5]*x[11]*n_H+rate[69]*x[0]*x[16]*n_H+rate[84]*x[1]*x[11]*n_H+rate[118]*x[8]*x[16]*n_H+rate[119]*x[3]*x[16]*n_H+rate[120]*x[12]*x[16]*n_H+rate[121]*x[13]*x[16]*n_H+rate[123]*x[16]*x[22]*n_H+rate[126]*x[7]*x[11]*n_H+rate[127]*x[6]*x[29]*n_H+rate[159]*x[18]*x[21]*n_H+rate[275]*x[21]+rate[305]*x[21];
  ode[19] = form+x[19]*loss;
  loss = -rate[0]*x[31]*n_H-rate[9]*x[30]*n_H-rate[19]*x[29]*n_H-rate[20]*x[29]*n_H-rate[21]*x[29]*n_H-rate[22]*x[11]*n_H-rate[23]*x[4]*n_H-rate[52]*x[5]*n_H-rate[68]*x[0]*n_H-rate[77]*x[1]*n_H-rate[104]*x[7]*n_H-rate[105]*x[8]*n_H-rate[107]*x[3]*n_H-rate[109]*x[12]*n_H-rate[111]*x[13]*n_H-rate[112]*x[14]*n_H-rate[113]*x[22]*n_H-rate[116]*x[15]*n_H-rate[157]*x[18]*n_H-rate[170]*x[5]*n_H-rate[178]*x[1]*n_H-rate[182]*x[10]*n_H-rate[187]*x[7]*n_H-rate[188]*x[8]*n_H-rate[189]*x[12]*n_H-rate[191]*x[14]*n_H-rate[192]*x[15]*n_H-rate[257]*x[30]*n_H-rate[269]-rate[270]-rate[301]-rate[316]*x[31]*n_H-rate[322]*x[30]*n_H-rate[326]*x[29]*n_H;
  form = +rate[1]*x[16]*x[31]*n_H+rate[8]*x[24]*x[30]*n_H+2*rate[15]*x[16]*x[24]*n_H+rate[16]*x[24]*x[28]*n_H+rate[24]*x[16]*x[16]*n_H+rate[25]*x[16]*x[29]*n_H+rate[30]*x[16]*x[28]*n_H+rate[122]*x[14]*x[16]*n_H+rate[190]*x[2]*x[26]*n_H+rate[219]*x[23]*x_e*n_H+rate[223]*x[19]*x_e*n_H+rate[224]*x[19]*x_e*n_H+rate[251]*x[24]*x[31]*n_H+rate[273]*x[16]+rate[278]*x[21]+rate[282]*x[11]+rate[304]*x[16]+rate[307]*x[21];
  ode[20] = form+x[20]*loss;
  loss = -rate[2]*x[31]*n_H-rate[11]*x[30]*n_H-2*rate[32]*x[21]*n_H-rate[33]*x[28]*n_H-rate[34]*x[28]*n_H-rate[35]*x[17]*n_H-rate[70]*x[0]*n_H-rate[80]*x[1]*n_H-rate[159]*x[18]*n_H-rate[275]-rate[276]-rate[278]-rate[305]-rate[306]-rate[307];
  form = +rate[3]*x[11]*x[31]*n_H+rate[10]*x[16]*x[30]*n_H+rate[22]*x[11]*x[20]*n_H+rate[24]*x[16]*x[16]*n_H+2*rate[28]*x[11]*x[16]*n_H+rate[29]*x[16]*x[28]*n_H+rate[36]*x[11]*x[29]*n_H+rate[39]*x[11]*x[28]*n_H+rate[83]*x[1]*x[11]*n_H+rate[132]*x[6]*x[11]*n_H+rate[135]*x[6]*x[17]*n_H+rate[136]*x[11]*x[12]*n_H+rate[137]*x[6]*x[27]*n_H+rate[138]*x[11]*x[14]*n_H+rate[198]*x[2]*x[19]*n_H+rate[226]*x[6]*x_e*n_H+rate[229]*x[3]*x_e*n_H+rate[246]*x[19]*x_e*n_H+rate[257]*x[20]*x[30]*n_H+rate[280]*x[11];
  ode[21] = form+x[21]*loss;
  loss = -rate[99]*x[24]*n_H-rate[113]*x[20]*n_H-rate[123]*x[16]*n_H-rate[145]*x[28]*n_H-rate[152]*x[17]*n_H-rate[239]*x_e*n_H;
  form = +rate[65]*x[5]*x[27]*n_H+rate[66]*x[14]*x[30]*n_H+rate[76]*x[0]*x[27]*n_H+rate[97]*x[10]*x[17]*n_H+rate[98]*x[13]*x[24]*n_H+rate[110]*x[17]*x[26]*n_H+rate[112]*x[14]*x[20]*n_H+rate[115]*x[4]*x[26]*n_H+rate[116]*x[15]*x[20]*n_H+rate[117]*x[23]*x[29]*n_H+rate[122]*x[14]*x[16]*n_H+rate[124]*x[4]*x[23]*n_H+rate[125]*x[19]*x[29]*n_H+rate[137]*x[6]*x[27]*n_H+rate[138]*x[11]*x[14]*n_H+rate[143]*x[8]*x[27]*n_H+rate[144]*x[14]*x[28]*n_H+rate[148]*x[3]*x[27]*n_H+rate[150]*x[12]*x[27]*n_H+rate[151]*x[14]*x[17]*n_H+rate[326]*x[20]*x[29]*n_H;
  ode[22] = form+x[22]*loss;
  loss = -rate[43]*x[31]*n_H-rate[55]*x[30]*n_H-rate[117]*x[29]*n_H-rate[124]*x[4]*n_H-rate[219]*x_e*n_H-rate[220]*x_e*n_H-rate[221]*x_e*n_H-rate[274];
  form = +rate[44]*x[19]*x[31]*n_H+rate[52]*x[5]*x[20]*n_H+rate[53]*x[26]*x[30]*n_H+rate[68]*x[0]*x[20]*n_H+rate[82]*x[1]*x[11]*n_H+rate[105]*x[8]*x[20]*n_H+rate[107]*x[3]*x[20]*n_H+rate[109]*x[12]*x[20]*n_H+rate[111]*x[13]*x[20]*n_H+rate[113]*x[20]*x[22]*n_H+rate[158]*x[16]*x[18]*n_H+rate[171]*x[5]*x[16]*n_H+rate[183]*x[10]*x[16]*n_H+rate[193]*x[7]*x[16]*n_H+rate[194]*x[8]*x[16]*n_H+rate[195]*x[12]*x[16]*n_H+rate[196]*x[14]*x[16]*n_H+rate[197]*x[15]*x[16]*n_H+rate[256]*x[10]*x[30]*n_H+rate[272]*x[16]+rate[277]*x[19]+rate[303]*x[16];
  ode[23] = form+x[23]*loss;
  loss = -rate[8]*x[30]*n_H-rate[15]*x[16]*n_H-rate[16]*x[28]*n_H-rate[17]*x[28]*n_H-rate[18]*x[4]*n_H-rate[50]*x[5]*n_H-rate[67]*x[0]*n_H-rate[92]*x[8]*n_H-rate[95]*x[3]*n_H-rate[96]*x[12]*n_H-rate[98]*x[13]*n_H-rate[99]*x[22]*n_H-rate[102]*x[15]*n_H-rate[177]*x[1]*n_H-rate[185]*x[14]*n_H-rate[186]*x[15]*n_H-rate[251]*x[31]*n_H-rate[255]*x[30]*n_H-rate[259]*x[7]*n_H-rate[260]*x[29]*n_H-rate[261]*x[29]*n_H-rate[268]-rate[297]-rate[300];
  form = +rate[0]*x[20]*x[31]*n_H+rate[6]*x[27]*x[31]*n_H+rate[19]*x[20]*x[29]*n_H+rate[90]*x[1]*x[27]*n_H+rate[108]*x[17]*x[26]*n_H+rate[112]*x[14]*x[20]*n_H+rate[182]*x[10]*x[20]*n_H+rate[183]*x[10]*x[16]*n_H+rate[184]*x[2]*x[10]*n_H+rate[218]*x[26]*x_e*n_H+rate[220]*x[23]*x_e*n_H+rate[221]*x[23]*x_e*n_H+rate[225]*x[19]*x_e*n_H+rate[238]*x[14]*x_e*n_H+rate[243]*x[10]*x_e*n_H+rate[244]*x[10]*x_e*n_H+rate[245]*x[10]*x_e*n_H+rate[270]*x[20]+rate[271]*x[26]+rate[289]*x[27]+rate[301]*x[20]+rate[312]*x[27]+rate[316]*x[20]*x[31]*n_H+rate[322]*x[20]*x[30]*n_H;
  ode[24] = form+x[24]*loss;
  loss = -rate[296];
  form = +rate[49]*x[1]*x[30]*n_H+rate[77]*x[1]*x[20]*n_H+rate[78]*x[1]*x[16]*n_H+rate[79]*x[1]*x[16]*n_H+rate[80]*x[1]*x[21]*n_H+rate[81]*x[1]*x[11]*n_H+rate[82]*x[1]*x[11]*n_H+rate[83]*x[1]*x[11]*n_H+rate[84]*x[1]*x[11]*n_H+rate[85]*x[1]*x[28]*n_H+rate[86]*x[1]*x[17]*n_H+rate[87]*x[1]*x[17]*n_H+rate[88]*x[1]*x[27]*n_H+rate[89]*x[1]*x[27]*n_H+rate[90]*x[1]*x[27]*n_H+rate[91]*x[1]*x[4]*n_H+rate[155]*x[1]*x[31]*n_H+rate[156]*x[1]*x[31]*n_H+rate[169]*x[1]*x[30]*n_H+rate[177]*x[1]*x[24]*n_H+rate[178]*x[1]*x[20]*n_H+rate[179]*x[1]*x[11]*n_H+rate[180]*x[1]*x[17]*n_H+rate[181]*x[1]*x[4]*n_H+rate[242]*x[1]*x_e*n_H;
  ode[25] = form+x[25]*loss;
  loss = -rate[41]*x[31]*n_H-rate[53]*x[30]*n_H-rate[103]*x[29]*n_H-rate[106]*x[28]*n_H-rate[108]*x[17]*n_H-rate[110]*x[17]*n_H-rate[114]*x[4]*n_H-rate[115]*x[4]*n_H-rate[190]*x[2]*n_H-rate[218]*x_e*n_H-rate[271]-rate[302];
  form = +rate[42]*x[16]*x[18]*n_H+rate[43]*x[23]*x[31]*n_H+rate[50]*x[5]*x[24]*n_H+rate[51]*x[10]*x[30]*n_H+rate[67]*x[0]*x[24]*n_H+rate[79]*x[1]*x[16]*n_H+rate[80]*x[1]*x[21]*n_H+rate[81]*x[1]*x[11]*n_H+rate[92]*x[8]*x[24]*n_H+rate[95]*x[3]*x[24]*n_H+rate[96]*x[12]*x[24]*n_H+rate[99]*x[22]*x[24]*n_H+rate[157]*x[18]*x[20]*n_H+rate[170]*x[5]*x[20]*n_H+rate[178]*x[1]*x[20]*n_H+rate[182]*x[10]*x[20]*n_H+rate[187]*x[7]*x[20]*n_H+rate[188]*x[8]*x[20]*n_H+rate[189]*x[12]*x[20]*n_H+rate[191]*x[14]*x[20]*n_H+rate[192]*x[15]*x[20]*n_H+rate[252]*x[10]*x[31]*n_H+rate[269]*x[20]+rate[274]*x[23]+rate[279]*x[19];
  ode[26] = form+x[26]*loss;
  loss = -rate[6]*x[31]*n_H-rate[65]*x[5]*n_H-rate[76]*x[0]*n_H-rate[88]*x[1]*n_H-rate[89]*x[1]*n_H-rate[90]*x[1]*n_H-rate[137]*x[6]*n_H-rate[143]*x[8]*n_H-rate[148]*x[3]*n_H-rate[150]*x[12]*n_H-rate[175]*x[5]*n_H-rate[203]*x[7]*n_H-rate[289]-rate[299]-rate[312];
  form = +rate[17]*x[24]*x[28]*n_H+rate[18]*x[4]*x[24]*n_H+rate[20]*x[20]*x[29]*n_H+rate[21]*x[20]*x[29]*n_H+rate[23]*x[4]*x[20]*n_H+rate[26]*x[16]*x[29]*n_H+rate[27]*x[16]*x[29]*n_H+rate[31]*x[4]*x[16]*n_H+rate[94]*x[10]*x[28]*n_H+rate[99]*x[22]*x[24]*n_H+rate[101]*x[4]*x[10]*n_H+rate[113]*x[20]*x[22]*n_H+rate[123]*x[16]*x[22]*n_H+rate[145]*x[22]*x[28]*n_H+rate[152]*x[17]*x[22]*n_H+rate[167]*x[14]*x[31]*n_H+rate[185]*x[14]*x[24]*n_H+rate[191]*x[14]*x[20]*n_H+rate[196]*x[14]*x[16]*n_H+rate[202]*x[14]*x[29]*n_H+rate[205]*x[11]*x[14]*n_H+rate[208]*x[14]*x[28]*n_H+rate[211]*x[14]*x[17]*n_H+rate[214]*x[4]*x[14]*n_H+rate[239]*x[22]*x_e*n_H+rate[260]*x[24]*x[29]*n_H+rate[261]*x[24]*x[29]*n_H;
  ode[27] = form+x[27]*loss;
  loss = -rate[4]*x[31]*n_H-rate[13]*x[30]*n_H-rate[16]*x[24]*n_H-rate[17]*x[24]*n_H-rate[29]*x[16]*n_H-rate[30]*x[16]*n_H-rate[33]*x[21]*n_H-rate[34]*x[21]*n_H-rate[37]*x[29]*n_H-rate[39]*x[11]*n_H-2*rate[40]*x[28]*n_H-rate[61]*x[5]*n_H-rate[73]*x[0]*n_H-rate[85]*x[1]*n_H-rate[93]*x[10]*n_H-rate[94]*x[10]*n_H-rate[106]*x[26]*n_H-rate[128]*x[7]*n_H-rate[139]*x[8]*n_H-rate[140]*x[3]*n_H-rate[142]*x[12]*n_H-rate[144]*x[14]*n_H-rate[145]*x[22]*n_H-rate[164]*x[18]*n_H-rate[173]*x[5]*n_H-rate[200]*x[7]*n_H-rate[208]*x[14]*n_H-rate[254]*x[31]*n_H-rate[283]-rate[284]-rate[309]-rate[317]*x[31]*n_H-rate[323]*x[30]*n_H;
  form = +rate[5]*x[17]*x[31]*n_H+rate[6]*x[27]*x[31]*n_H+rate[7]*x[4]*x[31]*n_H+rate[12]*x[29]*x[30]*n_H+2*rate[14]*x[4]*x[30]*n_H+rate[19]*x[20]*x[29]*n_H+rate[23]*x[4]*x[20]*n_H+rate[25]*x[16]*x[29]*n_H+rate[35]*x[17]*x[21]*n_H+rate[36]*x[11]*x[29]*n_H+2*rate[38]*x[17]*x[29]*n_H+rate[86]*x[1]*x[17]*n_H+rate[96]*x[12]*x[24]*n_H+rate[109]*x[12]*x[20]*n_H+rate[114]*x[4]*x[26]*n_H+rate[120]*x[12]*x[16]*n_H+rate[124]*x[4]*x[23]*n_H+rate[126]*x[7]*x[11]*n_H+rate[127]*x[6]*x[29]*n_H+rate[149]*x[12]*x[17]*n_H+rate[150]*x[12]*x[27]*n_H+rate[151]*x[14]*x[17]*n_H+rate[188]*x[8]*x[20]*n_H+rate[194]*x[8]*x[16]*n_H+rate[207]*x[8]*x[17]*n_H+rate[209]*x[4]*x[8]*n_H+rate[233]*x[12]*x_e*n_H+rate[235]*x[13]*x_e*n_H+rate[236]*x[13]*x_e*n_H+rate[253]*x[29]*x[31]*n_H+rate[286]*x[17]+rate[310]*x[17]+rate[318]*x[17]*x[31]*n_H+rate[324]*x[17]*x[30]*n_H;
  ode[28] = form+x[28]*loss;
  loss = -rate[12]*x[30]*n_H-rate[19]*x[20]*n_H-rate[20]*x[20]*n_H-rate[21]*x[20]*n_H-rate[25]*x[16]*n_H-rate[26]*x[16]*n_H-rate[27]*x[16]*n_H-rate[36]*x[11]*n_H-rate[37]*x[28]*n_H-rate[38]*x[17]*n_H-rate[56]*x[5]*n_H-rate[71]*x[0]*n_H-rate[103]*x[26]*n_H-rate[117]*x[23]*n_H-rate[125]*x[19]*n_H-rate[127]*x[6]*n_H-rate[129]*x[8]*n_H-rate[130]*x[3]*n_H-rate[131]*x[12]*n_H-rate[161]*x[18]*n_H-rate[162]*x[18]*n_H-rate[202]*x[14]*n_H-rate[253]*x[31]*n_H-rate[260]*x[24]*n_H-rate[261]*x[24]*n_H-rate[262]*x[10]*n_H-rate[263]*x[10]*n_H-2*rate[264]*x[29]*n_H-rate[298]-rate[326]*x[20]*n_H;
  form = +rate[4]*x[28]*x[31]*n_H+rate[7]*x[4]*x[31]*n_H+rate[16]*x[24]*x[28]*n_H+rate[18]*x[4]*x[24]*n_H+rate[29]*x[16]*x[28]*n_H+rate[33]*x[21]*x[28]*n_H+rate[40]*x[28]*x[28]*n_H+rate[88]*x[1]*x[27]*n_H+rate[89]*x[1]*x[27]*n_H+rate[91]*x[1]*x[4]*n_H+rate[92]*x[8]*x[24]*n_H+rate[100]*x[4]*x[10]*n_H+rate[102]*x[15]*x[24]*n_H+rate[105]*x[8]*x[20]*n_H+rate[115]*x[4]*x[26]*n_H+rate[116]*x[15]*x[20]*n_H+rate[118]*x[8]*x[16]*n_H+rate[133]*x[8]*x[11]*n_H+rate[139]*x[8]*x[28]*n_H+rate[141]*x[8]*x[17]*n_H+rate[142]*x[12]*x[28]*n_H+rate[143]*x[8]*x[27]*n_H+rate[144]*x[14]*x[28]*n_H+rate[160]*x[7]*x[31]*n_H+rate[187]*x[7]*x[20]*n_H+rate[193]*x[7]*x[16]*n_H+rate[199]*x[7]*x[11]*n_H+rate[200]*x[7]*x[28]*n_H+rate[201]*x[7]*x[17]*n_H+rate[203]*x[7]*x[27]*n_H+rate[204]*x[4]*x[7]*n_H+rate[228]*x[8]*x_e*n_H+rate[231]*x[12]*x_e*n_H+rate[232]*x[12]*x_e*n_H+rate[234]*x[13]*x_e*n_H+rate[238]*x[14]*x_e*n_H+2*rate[240]*x[15]*x_e*n_H+rate[247]*x[7]*x_e*n_H+rate[284]*x[28]+rate[285]*x[8]+rate[289]*x[27]+2*rate[290]*x[4]+rate[309]*x[28]+rate[312]*x[27]+2*rate[313]*x[4]+rate[317]*x[28]*x[31]*n_H+2*rate[319]*x[4]*x[31]*n_H+rate[323]*x[28]*x[30]*n_H+2*rate[325]*x[4]*x[30]*n_H;
  ode[29] = form+x[29]*loss;
  loss = -rate[8]*x[24]*n_H-rate[9]*x[20]*n_H-rate[10]*x[16]*n_H-rate[11]*x[21]*n_H-rate[12]*x[29]*n_H-rate[13]*x[28]*n_H-rate[14]*x[4]*n_H-rate[48]*x[5]*n_H-rate[49]*x[1]*n_H-rate[51]*x[10]*n_H-rate[53]*x[26]*n_H-rate[55]*x[23]*n_H-rate[57]*x[7]*n_H-rate[59]*x[6]*n_H-rate[62]*x[8]*n_H-rate[64]*x[12]*n_H-rate[66]*x[14]*n_H-rate[153]*x[18]*n_H-rate[169]*x[1]*n_H-rate[255]*x[24]*n_H-rate[256]*x[10]*n_H-rate[257]*x[20]*n_H-rate[258]*x[19]*n_H-rate[293]-rate[294]-rate[295]-rate[315]*x[31]*n_H-rate[320]*x_e*n_H-2*rate[321]*x[30]*n_H-rate[322]*x[20]*n_H-rate[323]*x[28]*n_H-rate[324]*x[17]*n_H-rate[325]*x[4]*n_H-rate[328];
  form = +rate[0]*x[20]*x[31]*n_H+rate[1]*x[16]*x[31]*n_H+rate[2]*x[21]*x[31]*n_H+rate[3]*x[11]*x[31]*n_H+rate[4]*x[28]*x[31]*n_H+rate[5]*x[17]*x[31]*n_H+rate[27]*x[16]*x[29]*n_H+rate[41]*x[26]*x[31]*n_H+rate[42]*x[16]*x[18]*n_H+rate[43]*x[23]*x[31]*n_H+rate[44]*x[19]*x[31]*n_H+rate[45]*x[11]*x[18]*n_H+rate[46]*x[6]*x[31]*n_H+rate[47]*x[3]*x[31]*n_H+rate[60]*x[5]*x[11]*n_H+rate[67]*x[0]*x[24]*n_H+rate[68]*x[0]*x[20]*n_H+rate[69]*x[0]*x[16]*n_H+rate[70]*x[0]*x[21]*n_H+rate[71]*x[0]*x[29]*n_H+rate[72]*x[0]*x[11]*n_H+rate[73]*x[0]*x[28]*n_H+rate[74]*x[0]*x[17]*n_H+rate[75]*x[0]*x[2]*n_H+rate[76]*x[0]*x[27]*n_H+rate[78]*x[1]*x[16]*n_H+rate[80]*x[1]*x[21]*n_H+rate[81]*x[1]*x[11]*n_H+rate[82]*x[1]*x[11]*n_H+rate[98]*x[13]*x[24]*n_H+rate[106]*x[26]*x[28]*n_H+rate[110]*x[17]*x[26]*n_H+rate[125]*x[19]*x[29]*n_H+rate[131]*x[12]*x[29]*n_H+rate[154]*x[5]*x[31]*n_H+rate[170]*x[5]*x[20]*n_H+rate[171]*x[5]*x[16]*n_H+rate[172]*x[5]*x[11]*n_H+rate[173]*x[5]*x[28]*n_H+rate[174]*x[5]*x[17]*n_H+rate[175]*x[5]*x[27]*n_H+rate[176]*x[4]*x[5]*n_H+rate[217]*x[0]*x_e*n_H+rate[221]*x[23]*x_e*n_H+rate[223]*x[19]*x_e*n_H+rate[225]*x[19]*x_e*n_H+rate[229]*x[3]*x_e*n_H+rate[232]*x[12]*x_e*n_H+rate[234]*x[13]*x_e*n_H+rate[236]*x[13]*x_e*n_H+rate[267]*x[0]+rate[278]*x[21]+rate[279]*x[19]+rate[281]*x[11]+rate[282]*x[11]+rate[307]*x[21]+rate[308]*x[11]+rate[321]*x[30]*x[30]*n_H+rate[322]*x[20]*x[30]*n_H+rate[323]*x[28]*x[30]*n_H+rate[324]*x[17]*x[30]*n_H+rate[325]*x[4]*x[30]*n_H+rate[327]*x[31]*n_H;
  ode[30] = form+x[30]*loss;
  loss = -rate[0]*x[20]*n_H-rate[1]*x[16]*n_H-rate[2]*x[21]*n_H-rate[3]*x[11]*n_H-rate[4]*x[28]*n_H-rate[5]*x[17]*n_H-rate[6]*x[27]*n_H-rate[7]*x[4]*n_H-rate[41]*x[26]*n_H-rate[43]*x[23]*n_H-rate[44]*x[19]*n_H-rate[46]*x[6]*n_H-rate[47]*x[3]*n_H-rate[154]*x[5]*n_H-rate[155]*x[1]*n_H-rate[156]*x[1]*n_H-rate[160]*x[7]*n_H-rate[167]*x[14]*n_H-rate[249]*x[18]*n_H-rate[250]*x[18]*n_H-rate[251]*x[24]*n_H-rate[252]*x[10]*n_H-rate[253]*x[29]*n_H-rate[254]*x[28]*n_H-rate[292]-rate[315]*x[30]*n_H-rate[316]*x[20]*n_H-rate[317]*x[28]*n_H-rate[318]*x[17]*n_H-rate[319]*x[4]*n_H-2*rate[327]*n_H;
  form = +rate[8]*x[24]*x[30]*n_H+rate[9]*x[20]*x[30]*n_H+rate[10]*x[16]*x[30]*n_H+rate[11]*x[21]*x[30]*n_H+rate[12]*x[29]*x[30]*n_H+rate[13]*x[28]*x[30]*n_H+rate[17]*x[24]*x[28]*n_H+rate[20]*x[20]*x[29]*n_H+rate[21]*x[20]*x[29]*n_H+2*rate[26]*x[16]*x[29]*n_H+rate[37]*x[28]*x[29]*n_H+rate[48]*x[5]*x[30]*n_H+rate[49]*x[1]*x[30]*n_H+rate[50]*x[5]*x[24]*n_H+rate[51]*x[10]*x[30]*n_H+rate[52]*x[5]*x[20]*n_H+rate[53]*x[26]*x[30]*n_H+rate[54]*x[5]*x[16]*n_H+rate[55]*x[23]*x[30]*n_H+rate[56]*x[5]*x[29]*n_H+rate[57]*x[7]*x[30]*n_H+rate[58]*x[5]*x[11]*n_H+rate[59]*x[6]*x[30]*n_H+rate[60]*x[5]*x[11]*n_H+rate[61]*x[5]*x[28]*n_H+rate[62]*x[8]*x[30]*n_H+rate[63]*x[5]*x[17]*n_H+rate[64]*x[12]*x[30]*n_H+rate[65]*x[5]*x[27]*n_H+rate[66]*x[14]*x[30]*n_H+rate[75]*x[0]*x[2]*n_H+rate[77]*x[1]*x[20]*n_H+rate[79]*x[1]*x[16]*n_H+rate[81]*x[1]*x[11]*n_H+rate[84]*x[1]*x[11]*n_H+rate[85]*x[1]*x[28]*n_H+rate[87]*x[1]*x[17]*n_H+rate[93]*x[10]*x[28]*n_H+rate[97]*x[10]*x[17]*n_H+rate[103]*x[26]*x[29]*n_H+rate[104]*x[7]*x[20]*n_H+rate[117]*x[23]*x[29]*n_H+rate[128]*x[7]*x[28]*n_H+rate[129]*x[8]*x[29]*n_H+rate[147]*x[2]*x[3]*n_H+rate[153]*x[18]*x[30]*n_H+rate[157]*x[18]*x[20]*n_H+rate[158]*x[16]*x[18]*n_H+rate[159]*x[18]*x[21]*n_H+rate[161]*x[18]*x[29]*n_H+rate[162]*x[18]*x[29]*n_H+rate[163]*x[11]*x[18]*n_H+rate[164]*x[18]*x[28]*n_H+rate[165]*x[17]*x[18]*n_H+rate[166]*x[2]*x[18]*n_H+rate[168]*x[4]*x[18]*n_H+2*rate[215]*x[5]*x_e*n_H+3*rate[216]*x[0]*x_e*n_H+rate[217]*x[0]*x_e*n_H+rate[218]*x[26]*x_e*n_H+rate[219]*x[23]*x_e*n_H+2*rate[220]*x[23]*x_e*n_H+rate[222]*x[19]*x_e*n_H+2*rate[224]*x[19]*x_e*n_H+rate[225]*x[19]*x_e*n_H+rate[226]*x[6]*x_e*n_H+2*rate[227]*x[6]*x_e*n_H+rate[228]*x[8]*x_e*n_H+rate[230]*x[3]*x_e*n_H+2*rate[231]*x[12]*x_e*n_H+rate[233]*x[12]*x_e*n_H+rate[234]*x[13]*x_e*n_H+2*rate[235]*x[13]*x_e*n_H+rate[237]*x[13]*x_e*n_H+rate[239]*x[22]*x_e*n_H+rate[241]*x[18]*x_e*n_H+rate[265]*x[5]+rate[266]*x[0]+rate[270]*x[20]+rate[273]*x[16]+rate[274]*x[23]+rate[276]*x[21]+rate[277]*x[19]+rate[280]*x[11]+rate[282]*x[11]+rate[284]*x[28]+rate[286]*x[17]+rate[293]*x[30]+2*rate[294]*x[30]+rate[301]*x[20]+rate[302]*x[26]+rate[304]*x[16]+rate[306]*x[21]+rate[309]*x[28]+rate[310]*x[17]+3*rate[315]*x[30]*x[31]*n_H+2*rate[316]*x[20]*x[31]*n_H+2*rate[317]*x[28]*x[31]*n_H+2*rate[318]*x[17]*x[31]*n_H+rate[319]*x[4]*x[31]*n_H+2*rate[320]*x[30]*x_e*n_H+2*rate[321]*x[30]*x[30]*n_H+rate[322]*x[20]*x[30]*n_H+rate[323]*x[28]*x[30]*n_H+rate[324]*x[17]*x[30]*n_H+2*rate[328]*x[30];
  ode[31] = form+x[31]*loss;

  /* Store the electron abundance in the user data */
  data->x_e = x_e;

  return(0);
}
/*=======================================================================*/
