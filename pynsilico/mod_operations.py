import numpy as np
from mod_updates import z_update_fast, z_update_slow

def check_state(new_state, attractors, atol=0.01):
  '''Returns True if the attractor candidate's state already exists
  in the known states. Incrimates the frequency of the encountered attractor'''
  for attractor in attractors:
    if np.allclose(attractor.state, new_state, atol=atol):
      attractor.freq += 1
      return True
  return False

def generate_new_state(n, slowid, fastid):
  new_state = np.zeros([3, n])
  new_state[0, :] = np.random.rand(n)
  new_state[1, :] = np.sqrt(new_state[0,:]) # Slow nodes
  new_state[2, :] = np.sqrt(new_state[0,:]) # Fast nodes
  
  #Full control over fast and slow nodes
  
  new_state[1, slowid] = 1.0
  new_state[2, slowid] = new_state[0, slowid]
  
  new_state[2, fastid] = 1.0
  new_state[1, fastid] = new_state[0, fastid]
  return new_state
      
def find_attractor(state, s, tol=1e-5):
  attr = 0
  cnt  = 0
  new_state = state
  n = np.size(new_state, 1)
  while attr == 0:
    attr = 1
    for i in np.random.permutation(n):
      attrf = 0
      while attrf == 0:
        attrf = 1
        for j in np.random.permutation(n):
          changed_node = z_update_fast(new_state[0, :], j, s)
          if (abs(new_state[1, j] - changed_node) > tol):
            new_state[1, j] = changed_node
            new_state[0, j] = min([1, new_state[1, j] * new_state[2, j]])
            attrf = 0
            break
      changed_node = z_update_slow(new_state[0, :], i, s)
      if (abs(new_state[2, i] - changed_node) > tol):
        new_state[2, i] = changed_node
        new_state[0, i] = min([1, new_state[1, i] * new_state[2, i]])
        attr = 0
        cnt += 1
        break
  return new_state