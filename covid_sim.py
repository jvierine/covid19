#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt

# 
# Death-rate: 1-2%, estimate cruise-ship 7/700, population bias towards old people
# Flock-immunity: 
# 
# 
class population():
    
    def __init__(self,
                 N,
                 pop_density=1e4,
                 p_infection=0.05,  # infection probability
                 L_infection=0.04,  # infection length
                 p_death=0.01,      # probability of death (Cruise ship)
                 t_restrictions=7.0, # time when restrictions on human interactions are placedst_
                 t_return_to_normal=100.0, # return to business as usual
                 T_infected=7.0):

        n.random.seed(0)
        
        self.N_pop = long(N)
        self.p_infection = p_infection
        self.t_infected = n.zeros(self.N_pop)
        self.t_infected[:] = 0.0

        self.t_return_to_normal = t_return_to_normal
        self.p_infection_orig = p_infection
        self.L_infection_orig = L_infection

        # infection length
        self.L_infection=L_infection
        
        # time after which patient is cured
        self.T_infected=T_infected
        
        self.infected=n.zeros(self.N_pop)
        
        self.infected[0:3]=1.0
        self.t_infected[0:3]=0.0
        
        self.immune=n.zeros(self.N_pop)
        self.immune[:]=0.0
        
        self.dead=n.zeros(self.N_pop)
        self.dead[:]=0.0
        self.p_death=p_death
        
        self.L_x=n.sqrt(N/pop_density)
        self.L_y=n.sqrt(N/pop_density)        
        
        self.pos_x=n.random.rand(self.N_pop)*self.L_x
        self.pos_y=n.random.rand(self.N_pop)*self.L_y
        self.t_now=0.0
        self.t_list=[]
        self.n_inf=[]
        self.n_dead=[]
        self.n_cured=[]
        self.n_cumulative_infected=[]
        self.flights=True
        self.n_total_infected=n.sum(self.infected)
        self.t_restrictions=t_restrictions


    def timestep(self):
        """
        use rules
        """
        print("t_now %1.2f infectious %d immune %d dead %d"%(self.t_now,n.sum(self.infected),n.sum(self.immune),n.sum(self.dead)))
        
        # cured after time sick, becomes immune
        i_idx=n.where( (self.infected>0) & ( (self.t_now - self.t_infected) > self.T_infected ) )[0]
#        print(i_idx)
        self.infected[i_idx]=0.0
        self.immune[i_idx]=1.0
        self.t_infected[i_idx]=0
        
        for i in i_idx:
            # probability of death
            if n.random.rand(1) < self.p_death:
                self.dead[i] = 1.0

        # model infections
        i_idx=n.where( (self.infected>0) )[0]

        
        for i in i_idx:
            p0_x=self.pos_x[i]
            p0_y=self.pos_y[i]
            # distance is less than infection range. can't get reinfected
            p_idx=n.where( (n.sqrt((self.pos_x - p0_x)**2.0 + (self.pos_y-p0_y)**2.0) < self.L_infection) & (self.infected < 1.0) & (self.immune < 1.0) )[0]
            for pi in p_idx:
                # if gets infected
                if n.random.rand(1) < self.p_infection:
                    self.infected[pi] = 1.0
                    self.t_infected[pi] = self.t_now
                    self.n_total_infected+=1.0
        self.t_list.append(self.t_now)
        self.n_inf.append(n.sum(self.infected))
        self.n_cured.append(n.sum(self.immune))
        self.n_dead.append(n.sum(self.dead))
        self.n_cumulative_infected.append(self.n_total_infected)

        if self.t_now > self.t_restrictions and self.t_now < self.t_return_to_normal:
            print("implementing quarantine policy")
            self.p_infection=0.0003
            self.L_infection=0.04
            self.flights=False
        elif self.t_now > self.t_return_to_normal:
            self.p_infection = self.p_infection_orig
            self.L_infection = self.L_infection_orig
            self.flights=True
            
            

        if self.flights:
            self.fly()
        
        self.t_now += 1.0
            
    def fly(self,n_fly=0.01):
        n_flights=int(self.N_pop*n_fly)
        for i in range(n_flights):
            pi=int(n.random.rand(1)*self.N_pop)

            self.pos_x[pi]=n.random.rand(1)*self.L_x
            self.pos_y[pi]=n.random.rand(1)*self.L_y
        
    def plot(self,n_plot_max=10000,show=False):

        if n_plot_max == None:
            n_plot_max=len(self.infected)
        pidx=n.arange(0,n_plot_max)
            
        i_idx=pidx[n.where(self.infected[pidx])[0]]
        plt.clf()
        plt.subplot(211)
        plt.plot(self.pos_x[pidx],self.pos_y[pidx],".",color="blue")
        plt.plot(self.pos_x[i_idx],self.pos_y[i_idx],".",color="red")
        if self.t_now <= self.t_restrictions:
            plt.title("No restrictions in place")
        else:
            plt.title("Interaction and travel restrictions in place")
            
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.subplot(212)
        plt.plot(self.t_list,100.0*n.array(self.n_inf)/self.N_pop,color="red",label="infected")
        plt.plot(self.t_list,100.0*n.array(self.n_cumulative_infected)/self.N_pop,color="brown",label="total infected")        
        plt.plot(self.t_list,100.0*n.array(self.n_cured)/self.N_pop,color="green",label="cured")
        plt.plot(self.t_list,100.0*n.array(self.n_dead)/self.N_pop,color="black",label="dead")
        plt.ylabel("Percent of population")
        plt.axvline(self.t_restrictions,color="red")
        plt.axvline(self.t_return_to_normal,color="green")        
        plt.legend()
        plt.xlabel("Time (days)")
        if show:
            plt.show()
        else:
            plt.pause(0.01)

p=population(4e4)
for i in range(150):
    p.plot()
    p.timestep()
p.plot(show=True)

