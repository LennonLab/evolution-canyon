import numpy as np
import DiversityTools as dt
import bide


""" code for running an individual-based simulation model for the Evolution
    Canyon project. """

def microbide(nPatches, lgp, northVal, southVal, env_optima, dorm_decay, btn):

    """
    imRate  :  number individuals immigrating from the regional pool
               each time step
    
    lgp  :  log-series parameter, typically close to 0.99
    
    state  :  determines whether local communities will have different
              have different environmental values. Value can be:
              'heterogeneous' or 'homogeneous'
            
    time  :  number of time steps to run simulation. Needs to be large
             enough to allow the community to reach a relatively stable state
    """
    
    ceil = 50 #     number of individuals per patch
    optima_dict = {}
    decay_dict = {}
    
    nCOM = [list([]) for _ in range(nPatches)]
    sCOM = [list([]) for _ in range(nPatches)]
    nCOM, sCOM = bide.immigration([nCOM, sCOM], optima_dict,
                                    env_optima, decay_dict, dorm_decay, ceil)
    
    Nx = np.random.beta(12.0, 1.2, nPatches).tolist()
    Sx = np.random.beta(1.2, 12.0, nPatches).tolist()

    t, im  = 1, 1 # number of time steps & individuals immigrating per time step
    while t <= 1000:
        
        """ Immigration from regional pool. At time = 1, the community patches
        are populated with the same large number of individuals. Then a zero-sum
        dynamic is enforced where death, immigration, emigration, etc. occur but
        never change the abundance in a patch. The zero-sum dynamic is our way
        of keeping patch abundance from exploding and/or drastically shrinking.
        
        Species IDs of immigrants are chosen at random from a log-series
        distribution, i.e., a monotonically decreasing frequency distribution.
        """
        
        """ Immigration """
        nCOM, sCOM = bide.immigration([nCOM, sCOM], optima_dict, env_optima, 
                                        decay_dict, dorm_decay, im)
    
        """ Reproduction """
        nCOM, sCOM = bide.reproduction([nCOM, sCOM])
            
        """ Dormancy """
        nCOM, sCOM = bide.dormancy(nCOM, northVal, optima_dict)
        
        """ Between patch dispersal """
        nCOM, sCOM = bide.disperse(nCOM, sCOM, southVal, btn)
            
        """ Emigration & Death """
        nCOM, sCOM = bide.death_emigration(nCOM, northVal, optima_dict, ceil)
            
        t += 1
        
        N, S, D, A, BP, Evar, nABV = dt.getNSDA(nCOM)
        print 'time:',t, 'north N:',N,'S:',S,'D:',D,'A:',A,'BP:',BP,'Evar:',Evar
        
        # ABV stands for 'abundance vector'; a list of species identities and
        # their numerical abundances
        
        N, S, D, A, BP, Evar, sABV = dt.getNSDA(sCOM)
        print 'time:',t, 'south N:',N,'S:',S,'D:',D,'A:',A,'BP:',BP,'Evar:',Evar
        
        BC = dt.Bray_Curtis(nABV, sABV)
        print 'time:',t, 'Bray-Curtis as percent similarity',BC, '\n'
    
    #print 'number of sites in north and south:',len(nCOM), len(sCOM)    
    return [nCOM, sCOM]