#include <iostream> 
#include "z_update.h"

float z_update(float* sv, int node){
// fast sub-variables = protein level regulations:

float s[5] = {0.0};
/*
 * Each biological factor is associated with an index number as follows
 * Wnt Dsh IGFI  R-Smad Ihh Gli2 Bcat Lef/Tcf Runx2 Sox9 PTHrP PPR Col-X PKA
 * 1   2    3     4     5   6    7     8       9    10    11   12  13    14
 * MEF2C FGF FGFR3 STAT1 Smadcomplex Col-II Nkx3.2 ERK1/2 TGFbeta 
 *  15   16   17    18       19       20     21     22      23
 * MMP13 Smad7 Smad3 FGFR1 ATF2 NFkb HDAC4 CCND1 Dlx5  BMP  p38
 *  24   25    26    27    28   29   30    31    32   33    34
 *  GSK DC PP2A Akt PI3K Elk-1 
 *  35  36  37  38   39   40   
 * Ras IGF-IR Msx2 EF-1 ATF4 HIF-2alpha  GREM1 FRZB DKK1 Frizzled-LRP5/7
 * 41   42    43   44   45     46         47   48   49      50     
 * Cytokines  ALK1 ALK5 Rinflam TAK1 JNK Proteoglycans IkB-a Socs FOXO1
 *  51        52   53    54     55   56     57          58   59    60
*/
  switch (node){
    case 1:
      return 1.0 - sv[49] - sv[48];
    case 2:
      return sv[50] - 0.3*sv[4];
    case 3:
      return 1.0;
    case 4:
      return sv[52] - 0.5*sv[25] + sv[52];
    case 5:
      return 1.0; 
    case 6:
      return  sv[5]/1.5;
    case 7:
      return 1.0 - sv[36];
    case 8:
      return sv[7]- sv[10];
    case 9:
      return sv[38] + sv[22] + sv[32] - sv[10] - sv[14]- sv[26]*sv[30] - sv[31] - sv[43];
    case 10:
      return sv[14] + sv[26] - sv[7] - sv[9];
    case 11:
      return 1.0;
    case 12:
      return sv[11];
    case 13:
      return 1.0;
    case 14:
      return (sv[12] + sv[3])*(s[2]);
    case 15:
      return sv[34] + sv[19] - sv[30] - sv[26];
    case 16:
      return 1.0;
    case 17:
      return sv[16];
    case 18:
      return (sv[17] + sv[27] + 0.75*sv[42])*s[3] - sv[19] + sv[54]*(1-sv[17]);
    case 19:
      return sv[4] - 0.25 * sv[22];
    case 20:
      return 1.0;
    case 21:
      return 1.0;
    case 22:
      return sv[41] - sv[37] + sv[55];
    case 23:
      return 1.0;
    case 24:
      return 1.0; 
    case 25:
      return 1.0;
    case 26:
      return sv[53] - 0.5*(s[3]*(sv[22]+ sv[31]+ sv[25]));
    case 27:
      return sv[16];
    case 28:
      return (sv[34] + sv[34]* sv[26])*s[2];
    case 29:
      return (sv[27] + sv[38] + sv[55] - sv[58])*s[2];
    case 30:
      return sv[37] - sv[35];
    case 31:
      return 1.0 - 0.5*(sv[35] + sv[18]);
    case 32:
      return sv[34] + sv[19] - sv[43]; 
    case 33:
      return 1.0-sv[47];
    case 34:
      return (sv[23] + sv[33] + sv[55])*s[3];
    case 35: // free GSK3B:
      return 1.0 - sv[14] - sv[38];
    case 36: // GSK3B insulated in DC:
      return (1 + sv[37] - sv[2] - 0.5*sv[22]);
    case 37:
      return sv[14];
    case 38:
      return sv[39] - 0.5*sv[37];
    case 39:
      return (sv[41] + sv[42] + sv[54])*s[3];
    case 40:
      return s[2] *(sv[22] + sv[56]);
    case 41:
      return (sv[1] + sv[33] +  sv[17] + sv[27])*s[4];
    case 42:
      return sv[3];
    case 43:
      return 1.0 - sv[32];
    case 44:
      return 1.0;
    case 45:
      return (sv[14] + sv[22])*s[2];
    case 46:
      return 1.0;
    case 47:
      return 1.0;
    case 48:
      return 1.0;
    case 49:
      return 1.0;
    case 50:
      return sv[1] - sv[49];
    case 51:
      return 1.0;
    case 52:
      return sv[23] + sv[33];
    case 53:
      return sv[23];
    case 54:
      return sv[51] - sv[59];
    case 55:
      return sv[54] - 0.5*sv[34];
    case 56:
      return sv[55];
    case 57:
      return 1.0 - sv[24];
    case 58:
      return 1.0 - sv[55];
    case 59:
      return 1.0;
    case 60:
      return 0.75 + sv[56] - sv[38] - sv[22];

    default:
      return 0.0;
      std::cout<< "Unknown state variable. \n" ;
    }
}