#!/usr/bin/env python

#########################################
# Copyright 2017, Mario Rosasco
#
# Special thanks to Eric Senning 
# for inspiration and helpful discussion
#########################################

import matplotlib.pyplot as plt
import numpy as np
import sys

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
        if True, will report the time at each cycle of the simulation and some final summary information
        
    Returns
    -------
    E: float
        The calculated energy transfer efficiency. Note that for this to be meaningful, the simulation must have been run long enough
        for most of the particles to have undergone energy transfer or emission.
    startingDists, remainingDists, emittedDists, quenchedDists: 4 numerical lists
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
    E = 1.0-F_DA/F_D

    if verbose:
        print "Of", N, "starting pairs", n_final, "remained."
        print n_emitted, "emitted during the simulation,"
        print n_transferred, "underwent energy transfer."
        print "Calculated transfer efficiency is:", E
    
    return (E, startingDists, remainingDists, emittedDists, quenchedDists)
    
if __name__ == '__main__':    
    # Run the simulation
    while(True):
        print '''
        Please select from the following options:
        1 - Run simulation for single average distance (user defined values)
        2 - Generate FRET curves (multiple simulations, user defined values)
        3 - Run simulation for single average distance (demonstration values)
        4 - Quit
        '''
        try:
            selection=int(raw_input('Please choose an option above: '))
        except ValueError:
            print "Invalid selection."
            continue
            
        # quit
        if selection == 4:
            quit()
        
        elif selection in [1,2,3]:
            # run single sim     
            if selection == 1:
                try:
                    # Set up the simulation environment
                    N = int(raw_input('Please enter number of FRET pairs: '))
                    avgDist = float(raw_input('Please enter average distance between FRET pairs (in Angstrom): '))
                    sigma = float(raw_input('Please enter std deviation of distance between FRET pairs (in Angstrom): '))
                    p1mob = float(raw_input('Please enter mobility of FRET donor (width of distribution, in Angstrom): '))
                    p2mob = float(raw_input('Please enter mobility of FRET acceptor (width of distribution, in Angstrom): '))

                    dT = float(raw_input('Please enter time step for simulation (in ns): '))
                    finalT = float(raw_input('Please enter total length of simulation (in ns): '))
                    tauDonor = float(raw_input('Please enter the time constant of the FRET donor (in ns): '))

                    R_0 = float(raw_input('Please enter the R_0 of the FRET pair(in Angstrom): '))
                except ValueError:
                    print "Error: could not parse input. Please confirm input data type and try again."
                    quit()
            
                # Run the simulation
                E, startingDists, remainingDists, emittedDists, quenchedDists = runSim(N, avgDist, sigma, p1mob, p2mob, dT, finalT, tauDonor, R_0, verbose=True)
                
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
                plt.show()
            
            # run mutiple sims to make FRET curve   
            elif selection == 2:
                try:
                    # Set up the simulation environment
                    N = int(raw_input('Please enter number of FRET pairs: '))
                    sigma = float(raw_input('Please enter std deviation of distance between FRET pairs (in Angstrom): '))
                    p1mob = float(raw_input('Please enter mobility of FRET donor (width of distribution, in Angstrom): '))
                    p2mob = float(raw_input('Please enter mobility of FRET acceptor (width of distribution, in Angstrom): '))

                    dT = float(raw_input('Please enter time step for simulation (in ns): '))
                    finalT = float(raw_input('Please enter total length of simulation (in ns): '))
                    tauDonor = float(raw_input('Please enter the time constant of the FRET donor (in ns): '))

                    R_0 = float(raw_input('Please enter the R_0 of the FRET pair(in Angstrom): '))
                except ValueError:
                    print "Error: could not parse input. Please confirm input data type and try again."
                    quit()
                    
                Evals_immobile = []
                Evals_mobile = []
                rvals = []
                
                for avgDist in range(1,20,1):
                    print "Running simulation with average distance =", avgDist
                    rvals.append(avgDist)
                    p1mob = 5
                    p2mob = 5
                    Evals_mobile.append(runSim(N, avgDist, sigma, p1mob, p2mob, dT, finalT, tauDonor, R_0, verbose=False)[0])
                    p1mob = 0.00001
                    p2mob = 0.00001
                    Evals_immobile.append(runSim(N, avgDist, sigma, p1mob, p2mob, dT, finalT, tauDonor, R_0, verbose=False)[0])
                    
                plt.plot(rvals, Evals_mobile, label="mobile donor/acceptor")
                plt.plot(rvals, Evals_immobile, label="immobile donor/acceptor")
                plt.ylabel("FRET Efficiency")
                plt.xlabel("Distance (Angstrom)")
                plt.legend(['mobile donor/acceptor', 'immobile donor/acceptor'])
                plt.show()
            
            # run single sim with predefined, demonstration values
            elif selection == 3:
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
                E, startingDists, remainingDists, emittedDists, quenchedDists = runSim(N, avgDist, sigma, p1mob, p2mob, dT, finalT, tauDonor, R_0, verbose=True)
                
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
                plt.show()
            
        else:
            print "Invalid selection."
            continue
    
    
        
    
