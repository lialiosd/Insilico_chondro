import numpy as np

def z_update_fast(state, node, s):
  if node == 0:
    y = 1.0 - state[48] - state[47]
  elif node == 1:
    y = state[49] - 0.3 * state[3]
  elif node == 2:
    y = 1.0
  elif node == 3:
    y = state[51] - 0.5 * state[24] + state[51]
  elif node == 4:
    y = 1.0
  elif node == 5:
    y = (2.0/3.0) * state[4]
  elif node == 6:
    y = 1.0 - state[35]
  elif node == 7:
    y = state[6]- state[9]
  elif node == 8:
    y = state[37] + state[21] + state[31] - state[9] - state[13] -\
          state[25] * state[29] - state[30] - state[42]
  elif node == 9:
    y = state[13] + state[25] - state[6] - state[8] 
  elif node == 10:
    y = 1.0
  elif node == 11:
    y = state[10]
  elif node == 12:
    y = 1.0 
  elif node == 13:
    y = (state[11] + state[2]) * s[1]
  elif node == 14:
    y = state[33] + state[18] - state[29] - state[25]
  elif node == 15:
    y = 1.0 
  elif node == 16:
    y = state[15]
  elif node == 17:
    y = (state[16] + state[26] + 0.75*state[41]) * s[2] -\
          state[18] + state[53]*(1-state[16])
  elif node == 18:
    y = state[3] - 0.25 * state[21]  
  elif node == 19:
    y = 1.0
  elif node == 20:
    y = 1.0
  elif node == 21:
    y = state[40] - state[36] + state[54]
  elif node == 22:
    y = 1.0
  elif node == 23:
    y = 1.0
  elif node == 24:
    y = 1.0
  elif node == 25:
    y = state[52] - 0.5 * s[2] * (state[21] + state[30] + state[24])
  elif node == 26:
    y = state[15]   
  elif node == 27:
    y = (state[33] + state[33] * state[25]) * s[1]
  elif node == 28:
    y = (state[26] + state[37] + state[54] - state[57]) * s[1];
  elif node == 29:
    y = state[36] - state[34]
  elif node == 30:
    y = 1.0 - 0.5 * (state[34] + state[17])
  elif node == 31:
    y = state[33] + state[18] - state[42] 
  elif node == 32:
    y = 1.0 - state[46]
  elif node == 33:
    y = (state[22] + state[32] + state[54]) * s[2]
  elif node == 34:
    y = 1.0 - state[13] - state[37]
  elif node == 35:
    y = (1.0 + state[36] - state[1] - 0.5 * state[21])
  elif node == 36:
    y = state[13]
  elif node == 37:
    y = state[38] - 0.5 * state[36]
  elif node == 38:
    y = (state[40] + state[41] + state[53])*s[2]
  elif node == 39:
    y = s[1] *(state[21] + state[55])
  elif node == 40:
    y = (state[0] + state[32] + state[16] + state[26]) * s[3];
  elif node == 41:
    y = state[2]
  elif node == 42:
    y = 1.0 - state[31]
  elif node == 43:
    y = 1.0
  elif node == 44:
    y = (state[13] + state[21]) * s[1]
  elif node == 45:
    y = 1.0
  elif node == 46:    
    y = 1.0
  elif node == 47:        
    y = 1.0
  elif node == 48:        
    y = 1.0
  elif node == 49:    
    y = state[0] - state[48]       
  elif node == 50:
    y = 1.0
  elif node == 51:    
    y = state[22] + state[32]
  elif node == 52:
    y = state[22]
  elif node == 53:
    y = state[50] - state[58]
  elif node == 54:
    y = state[53] - 0.5 * state[33]
  elif node == 55:
    y = state[54]     
  elif node == 56:       
    y = 1.0 - state[23]
  elif node == 57:
    y = 1.0 - state[54]
  elif node == 58:
    y = 1.0
  elif node == 59:
    y = 0.75 + state[55] - state[37] - state[21]
  else:
    print("Node out of bounds!")
  
  return np.clip(y, 0.0, 1.0)

def z_update_slow(state, node, s):
  if node == 0:
    y = 2.0 * state[5] - state[42]
  elif node == 1:
    y = 1.0;
  elif node == 2:
    y = state[13] 
  elif node == 3:
    y = 1.0 - state[44] + state[43]
  elif node == 4:
    y = (state[8] + state[18] + state[28] + state[44]  - state[43] - state[16]) * s[1]   
  elif node == 5:
    y = 1.0 - state[16] 
  elif node == 6:
    y = 1.0
  elif node == 7:
    y = (1.0 + state[8]+ state[6]) * s[2]
  elif node == 8:
    y = (state[45] + state[7] + state[8] + state[14] + state[31] -\
         state[20] * state[18] - state[42] + state[5] - state[13] + state[39]) * s[3] 
  elif node == 9:
    y = (state[13] - state[28] + state[33] + state[28] + state[9] + state[20]) * s[3]
  elif node == 10:
    y = (state[5] * state[9] + state[9] + state[25]) * s[2] 
  elif node == 11:
    y = (state[5] + state[9] + state[18] - state[17])*s[1]
  elif node == 12:
    y = (state[8] + state[3] + state[14] - state[13] -\
         state[10] + state[45] - state[17] + state[6]) * s[1] 
  elif node == 13:
    y = 1.0
  elif node == 14:
    y = (state[8] + state[18]) * s[1]
  elif node == 15:
    y = (state[6] + state[8]) * s[1] 
  elif node == 16:
    y = state[9] - state[39] * state[21] 
  elif node == 17:
    y = 1.0
  elif node == 18:
    y = 1.0
  elif node == 19:
    y = (state[9] + state[26] + state[20] - state[6] + state[59]) * s[2] 
  elif node == 20:
    y = (state[9] + state[13]) * s[1]
  elif node == 21:
    y = state[14]
  elif node == 22:
    y = state[5] - state[39] + state[25]
  elif node == 23:
    y = (state[8] + state[28] +state[45] + state[5] + state[45] * state[8] * state[33] +\
         state[18] + state[39] - state[25] + state[39] * state[27] * state[55] - state[59]) * 0.2 
  elif node == 24:
    y = (state[17] + state[18] + state[25] + 2.0 * state[28]) * s[4] 
  elif node == 25:
    y = 1.0
  elif node == 26:
    y = state[8] - 0.5 * state[18]
  elif node == 27:
    y = 1.0
  elif node == 28:
    y = 1.0
  elif node == 29:
    y = 1.0
  elif node == 30:
    y = (state[27] + state[5] + 1.5 * state[13] + state[59] * state[25]) * s[3] 
  elif node == 31:
    y = (state[14] + state[18] - state[25])
  elif node == 32:
    y = (state[5] + state[28]) * s[1] 
  elif node == 33:
    y = 1.0 - 0.4 * state[43]
  elif node == 34:
    y = 1.0
  elif node == 35:
    y = 1.0
  elif node == 36:
    y = (1.0 - state[6]) * s[1]
  elif node == 37:
    y = state[8]
  elif node == 38:
    y = state[8]
  elif node == 39:
    y = 1.0
  elif node == 40:
    y = 1.0
  elif node == 41:
    y = (state[25] + state[28] + state[17]) * s[2]
  elif node == 42:
    y = (state[18] + state[25] + state[42] - state[31]) * s[1] 
  elif node == 43:
    y = (state[25] + state[28] + state[39] - state[18] + state[3]) * s[2]
  elif node == 44:
    y = 1.0
  elif node == 45:
    y = state[28]
  elif node == 46:
    y = state[9] - state[6]
  elif node == 47:
    y = 0.5 * state[18] - state[6] + state[9] -state[28] -state[27]
  elif node == 48:
    y = 0.5 * state[18] - state[6] + state[9] -state[28] -state[27]
  elif node == 49:
    y = 1.0
  elif node == 50:
    y = state[39] * state[27] * state[55] + state[28]
  elif node == 51:
    y = s[2] * (0.4 + state[7] - state[25] + state[8] + state[27] * state[55])
  elif node == 52:
    y = s[1] * (0.4 + state[25] - state[28] + state[9])
  elif node == 53:
    y = 1.0 - state[59]
  elif node == 54:
    y = 1.0
  elif node == 55:
    y = 1.0
  elif node == 56:
    y = s[2] * (state[9] + state[59] - state[28])
  elif node == 57:
    y = state[28]
  elif node == 58:
    y = s[1] * (state[41] + state[16])
  elif node == 59:
    y = state[25] - state[39] * state[27] - state[28]
  else:
    print("Node out of bounds!")
    
  return np.clip(y, 0.0, 1.0)