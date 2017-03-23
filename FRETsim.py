#!/usr/bin/env python

#########################################
# Copyright 2017, Mario Rosasco
#
# Special thanks to Eric Senning 
# for inspiration and helpful discussion
#########################################

import matplotlib.pyplot as plt
import numpy as np

class FretPair:
    def __init__(self, dist, p1mobility, p2mobility):
        self.P1 = 0
        self.P2 = dist # in Angstrom
        self.dist = dist # in Angstrom
        self.excited = True
        self.emitted = False
        self.quenched = False
        self.p1mobility = p1mobility
        self.p2mobility = p2mobility
        
    def relaxPair(self):
        self.P1 = 0 + np.random.triangular(-self.p1mobility/2.0, 0, self.p1mobility/2.0)
        self.P2 = self.dist + np.random.triangular(-self.p2mobility/2.0, 0, self.p2mobility/2.0)
        
    def getDist(self):
        return float(self.P2 - self.P1)
        
    def checkExcited(self, dT, tau, R_0):
        k_T = ((1.0/tau)*((R_0/self.getDist())**6))
        k_Donor = 1.0/tau
        
        P_quenched = 1-(np.exp(-k_T * dT))
        P_emitted = 1-(np.exp(-k_Donor * dT))
        
        randTest = np.random.rand()
        # Using the same random variable to test for both
        # quenching and emission, so need to account for multiple
        # testing. Is this the right way to account for the TIErr ...?
        P_emitted_adj = P_emitted*(1-P_quenched)
        
        if randTest < P_quenched:
            self.quenched = True
            self.excited = False
            
        elif randTest > (1-P_emitted_adj):
            self.emitted = True
            self.excited = False

def runSim(N, avgDist, sigma, p1mob, p2mob, dT, finalT, tauDonor, R_0, verbose=True):
    '''
    Runs a simulation of time-dependency of FRET for two particles that are mobile on the timescale of dT
    
    Note that future simulations may make the mobility itself a time-dependent function; for now users should
    ensure that dT is large enough to reasonably allow for any movements of P1 and P2 to occur, but short enough
    that each timepoint provides a reasonable approximation of the probability of energy transfer
    relative to the probability of emission.
    
    Parameters
    ----------
    N: int
        Number of pairs to start the simulation with
    
    avgDist: numerical
        The mean distance between the pairs of points
    
    sigma: numerical
        The stdev of the distances between the pairs of points
        
    p1mob, p2mob: numerical
        The mobility of each point. Currently interpreted as the width of a triangular distribution centered at 0.
        
    dT: numerical
        time step (in ns)
        
    finalT: numberical
        length of simulation (in ns)
        
    tauDonor: numerical
        time constant for the donor fluorescence emission (in ns)
        
    R_0: numerical
        distance between donor and acceptor at which energy transfer occurs with 50% efficiency (in Angstrom)
        
    verbose: boolean
        if True, will report the time at each cycle of the simulation
        
    Returns
    -------
    (startingDists, remainingDists, emittedDists, quenchedDists): 4 numerical lists
        Contain the distances for the starting pairs, the remaining pairs, the emitted pairs, and the quenched pairs, respectively
    '''
    # Initialize N pairs of points with separation following a normal distribution centered around some position
    startingDists = np.random.normal(loc = avgDist, scale = sigma, size = N)

    # Build the list of pairs
    pairList = []
    for dist in startingDists:
        pairList.append(FretPair(dist, p1mob, p2mob))

    # 2) Iterate through the number of steps needed to complete the simulation length
    for currT in np.arange(0, finalT, dT):
        if verbose:
            print "t =", currT, "ns"
        
        for pair in pairList:
            if pair.excited:
                # a) Allow each site to move, based on the range of movements predicted
                pair.relaxPair()

                # b) Compute the FRET and emission probability for the pair
                pair.checkExcited(dT, tauDonor, R_0)

    # Collect all the pairs that are still excited
    remainingDists = []
    emittedDists = []
    quenchedDists = []
    for pair in pairList:
        if pair.excited:
            remainingDists.append(pair.getDist())
        elif pair.quenched:
            quenchedDists.append(pair.getDist())
        elif pair.emitted:
            emittedDists.append(pair.getDist())
        else:
            print "Error: Pair found that was not categorized appropriately. Please check the simulation code."
            quit()

    n_final = len(remainingDists)
    n_emitted = len(emittedDists)
    n_transferred = N - n_final - n_emitted

    k_D = 1.0/tauDonor
    F_D = N*(1 - np.exp(-k_D*finalT))
    F_DA = float(n_emitted)

    print "Of", N, "starting pairs", n_final, "remained."
    print n_emitted, "emitted during the simulation,"
    print n_transferred, "underwent energy transfer."
    print "Calculated transfer efficiency is:", 1.0-F_DA/F_D
    
    return (startingDists, remainingDists, emittedDists, quenchedDists)
    
if __name__ == '__main__':    
    # Set up the simulation environment
    N = 100000
    avgDist = 6 # in Angstrom
    sigma = 5 # in Angstrom
    p1mob = 1 # currently interpreted as width of triangular distribution. 
    p2mob = 0.01 # currently interpreted as width of triangular distribution. 

    dT = 0.1 # in ns
    finalT = 10.0 # in ns
    tauDonor = 2.0 # in ns

    R_0 = 10.0 # in Angstrom
    
    # Run the simulation
    startingDists, remainingDists, emittedDists, quenchedDists = runSim(N, avgDist, sigma, p1mob, p2mob, dT, finalT, tauDonor, R_0, verbose=True)
    
    # Show distance histograms
    ax1 = plt.subplot(311)
    plt.hist(startingDists,100)    
    plt.ylabel("# At t=0")
    
    ax2 = plt.subplot(312, sharex=ax1, sharey=ax1)
    plt.hist(emittedDists,100)
    plt.ylabel("# Emitted")
    
    ax3 = plt.subplot(313, sharex=ax1, sharey=ax1)
    plt.hist(quenchedDists,100)
    plt.ylabel("# Quenched")
    
    plt.xlabel("Distances (Angstrom)")
    plt.locator_params(axis='y', numticks=3) #Why doesn't this work here? Hmm...
    plt.show()
