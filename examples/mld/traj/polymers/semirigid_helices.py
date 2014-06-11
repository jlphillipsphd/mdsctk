#!/usr/bin/python
# 
#  
#                 This source code is part of
#  
#                         M D S C T K
#  
#        Molecular Dynamics Spectral Clustering ToolKit
#  
#                         VERSION 1.2.0
#  Written by Joshua L. Phillips.
#  Copyright (c) 2013, Joshua L. Phillips.
#  check out http://github.com/douradopalmares/mdsctk/ for more information.
# 
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#  
#  If you want to redistribute modifications, please consider that
#  derived work must not be called official MDSCTK. Details are found
#  in the README & LICENSE files - if they are missing, get the
#  official version at <website>.
#  
#  To help us fund MDSCTK development, we humbly ask that you cite
#  the papers on the package - you can find them in the top README file.
#  
#  For more info, check our website at http://github.com/douradopalmares/mdsctk/
#  
# 

# Bring in useful modules
from __future__ import division
import random
import math
import sys

random.seed(0)
rnoise=random.Random(0)
rshift=random.Random(0)

def cross_pi(val1,val2):
    minval=min(val1,val2)
    maxval=max(val1,val2)
    while (minval<-math.pi):
        minval+=math.pi
        maxval+=math.pi
    while (maxval>math.pi):
        minval-=math.pi
        maxval-=math.pi
    if (minval*maxval<0.0): return 1
    return 0

if len(sys.argv)!=6:
    print "Usage: %s <number_of_structures> <length> <sd_theta> <sd_phi> <num_of_domains>"%sys.argv[0]
    exit()
############################################################################
ntrials=int(sys.argv[1])         # Number of configurations to generate
nlinks=int(sys.argv[2])         # Number of links in chain
theta_sd=float(sys.argv[3])
phi_sd=float(sys.argv[4])
ndomains=int(sys.argv[5])
theta_fixed=1.537648     # For an alpha-helix
phi_fixed=0.8290859      # For an alpha-helix
bondlength=3.88     # Typical C-alpha to C-alpha
theta_buffer=0.08725  # phi angles are indiscernable when theta gets too close to zero or pi
theta_init=[]
phi_init=[]
domain=[]
theta_last=[]
phi_last=[]
temp_theta_last=[]
temp_phi_last=[]
riemann_sheet=[]
riemann_last=[]

## Get domain rates
theta_domain_rate=[]
theta_domain_rate_sd=[]
theta_domain_sd=[]
theta_domain_noise=[]
theta_domain_flip=[]
phi_domain_rate=[]
phi_domain_rate_sd=[]
phi_domain_sd=[]
phi_domain_noise=[]
phi_domain_flip=[]
for i in range(ndomains):
    theta_domain_rate.append(-abs(float(sys.stdin.readline())))
    theta_domain_rate_sd.append(abs(float(sys.stdin.readline())))    
    theta_domain_sd.append(abs(float(sys.stdin.readline())))        
    print "Domain %d theta rate:  %f"%(i,theta_domain_rate[i])
    print "Domain %d theta rate noise: %f"%(i,theta_domain_rate_sd[i])
    print "Domain %d theta noise: %f"%(i,theta_domain_sd[i])
    theta_domain_noise.append(0.0)
    theta_domain_flip.append(False)
    phi_domain_rate.append(-abs(float(sys.stdin.readline())))
    phi_domain_rate_sd.append(abs(float(sys.stdin.readline())))
    phi_domain_sd.append(abs(float(sys.stdin.readline())))
    print "Domain %d phi   rate:  %f"%(i,phi_domain_rate[i])
    print "Domain %d phi   rate noise: %f"%(i,phi_domain_rate_sd[i])
    print "Domain %d phi   noise: %f"%(i,phi_domain_sd[i])
    phi_domain_noise.append(0.0)   
    phi_domain_flip.append(False)

mydomain=-1
for i in range(nlinks-1):
    if (i % int((nlinks-1)/ndomains) == 0 and mydomain+1 < ndomains):
        mydomain += 1

    domain.append(mydomain)
    theta_init.append(theta_fixed)  # Initial theta value
    phi_init.append(phi_fixed)  # Initial phi value
    
    theta_last.append(theta_init[i])
    phi_last.append(phi_init[i])

    temp_theta_last.append(theta_init[i])
    temp_phi_last.append(phi_init[i])

    riemann_sheet.append(0)  # Counter of number of times theta passes 0 or pi
    riemann_last.append(0)

#add one more place in the riemann_sheet variable
riemann_sheet.append(0)
riemann_last.append(0)

print "Link domain assignment:",domain

############################################################################

################ Define Python classes for simulation ########################
# This is a utility class for managing points
class point:
    def __init__(self):
        self.x=0.0
        self.y=0.0
        self.z=0.0
    def set_value(self,x,y,z):
        self.x=x
        self.y=y
        self.z=z
    def add(self, apoint):
        self.set_value(self.x+apoint.x,self.y+apoint.y,self.z+apoint.z)
    def subtract(self,apoint):
        self.set_value(self.x-apoint.x,self.y-apoint.y,self.z-apoint.z)
    def distance(self, apoint):
        return math.sqrt((self.x-apoint.x)**2+(self.y-apoint.y)**2+(self.z-apoint.z)**2)
    def mag(self):
        return math.sqrt(self.x**2+self.y**2+self.z**2)
    def randomize(self, length):
        self.x=2.*random.random()-1.
        self.y=2.*random.random()-1.
        self.z=2.*random.random()-1.
        mag=self.mag()
        self.x*=length/mag
        self.y*=length/mag
        self.z*=length/mag
        #print 'Magnitude of new link =%lf'%(self.mag())
        #self.printpoint()
    def printpoint(self):
        print '(%6.4lf, %6.4lf, %6.4lf)' % (self.x,self.y,self.z)
    def theta(self):
        theta=math.atan2(math.sqrt(self.x*self.x+self.y*self.y),self.z)
        return theta
    def phi(self):
        phi=math.atan2(self.y,self.x)
        return phi
    def rotate_y(self,theta):
        c=math.cos(theta)
        s=math.sin(theta)
        bpoint=point()
        bpoint.x=s*self.z+c*self.x
        bpoint.y=self.y
        bpoint.z=c*self.z-s*self.x
        self.x=bpoint.x
        self.y=bpoint.y
        self.z=bpoint.z
    def rotate_z(self,theta):
        c=math.cos(theta)
        s=math.sin(theta)
        bpoint=point()
        bpoint.x=c*self.x-s*self.y
        bpoint.y=c*self.y+s*self.x
        bpoint.z=self.z
        self.x=bpoint.x
        self.y=bpoint.y
        self.z=bpoint.z
    def genlink(self,base,vector,penult,theta,phi,flip_axis):
        #print '*****************************'
        apoint=point()       
        #Create a new link pointing in the direction of the last link
        apoint.set_value(vector.x,vector.y,vector.z)
        #Rotate this by theta around the normal to penult x vector
        axis=penult.cross_product(vector)
	#print 'Axis before',
	#axis.printpoint()
        if (flip_axis):
            #print 'Flipping axis'
            axis.reverse()
	#print 'Axis after',
	#axis.printpoint()
        bpoint=apoint.rotate(axis,-theta)
        #print 'Point after theta rotation',
        #bpoint.printpoint()
        #Rotate the result by phi around the final vector
        apoint=bpoint.rotate(vector,phi)
        #print 'Point after phi rotation',
        #apoint.printpoint()
        #print 'Base:',
        #base.printpoint()
        apoint.add(base)
        return apoint
    def cross_product(self,vecb):
        apoint=point()
        apoint.x=self.y*vecb.z-self.z*vecb.y
        apoint.y=self.z*vecb.x-self.x*vecb.z
        apoint.z=self.x*vecb.y-self.y*vecb.x
        return apoint
    def dot(self,apoint):
	return self.x*apoint.x+self.y*apoint.y+self.z*apoint.z
    def reverse(self):
        self.set_value(-self.x,-self.y,-self.z)
    def rotate(self,axis,theta):
        c=math.cos(theta)
        s=math.sin(theta)
        #print 'Axis:',
        #axis.printpoint()
        mag=axis.mag()
        if mag == 0:
            return self
        x=axis.x/mag
        y=axis.y/mag
        z=axis.z/mag
        omc=1.-c
        #Create rotate matrix around axis
        q11=x*x*omc+c
        q12=x*y*omc-z*s
        q13=x*z*omc+y*s
        q21=x*y*omc+z*s
        q22=y*y*omc+c
        q23=y*z*omc-x*s
        q31=x*z*omc-y*s
        q32=y*z*omc+x*s
        q33=z*z*omc+c
        #rotate self and store in new point
        apoint=point()
        apoint.x=q11*self.x+q12*self.y+q13*self.z
        apoint.y=q21*self.x+q22*self.y+q23*self.z
        apoint.z=q31*self.x+q32*self.y+q33*self.z
        return apoint

# This is a utility class for managing chains
class chain:
    def __init__(self):
        self.points=[]
        self.n_points=0
        apoint=point()
        self.append(apoint)
    def append(self,apoint):
        (self.points).append(apoint)
        self.n_points+=1
        #print 'Appending point (%lf,%lf,%lf), npoints=%lf'\
        #      %(apoint.x,apoint.y,apoint.z,self.n_points)
    def grow(self,theta,phi,flip_axis):
        apoint=point()
        #print '**********************************'
        #print 'Entering grow, n_points=%d'%(self.n_points)
        if self.n_points>2:
            #self.printchain()
            #print 'Final_vector',
            bpoint=self.final_vector()
            #bpoint.printpoint()
            #print 'Penult_vector',
            cpoint=self.penult_vector()
            #cpoint.printpoint()
            bpoint=apoint.genlink(self.end(),self.final_vector(),self.penult_vector(),theta,phi,flip_axis)
            #print 'Generated point'
            #bpoint.printpoint()
            self.append(bpoint)
        else:
            if self.n_points==2:
                bpoint=point()
                bpoint.x=bondlength*math.cos(theta)+bondlength
                bpoint.y=bondlength*math.sin(theta)
                self.append(bpoint)
            else:
                bpoint=point()
                bpoint.set_value(bondlength,0.,0.)
                self.append(bpoint)
    def distance(self,point1,point2):
        return self.points[point1].distance(self.points[point2])
    def printchain(self):
        print 'n_points=%d' % self.n_points
        for i in self.points:
            i.printpoint()
    def printchainpdb(self,filenumber):
        pdbfilename='chain.%05d.pdb'%(filenumber)
        pdbfile=open(pdbfilename,'w')
        print>>pdbfile,'COMMENT SYSTEM OF CA ATOMS'
        icount=1
        for iatom in self.points:
            print>>pdbfile,'ATOM  % 5d  CA       % 4d    %8.3lf%8.3lf%8.3lf'%\
                                 (icount,icount,iatom.x,iatom.y,iatom.z)
            icount+=1
        pdbfile.close()
    def printchainpsf(self):
        pdbfilename='chain.psf'
        psffile=open(pdbfilename,'w')
        print>>psffile,'PSF\n'
        print>>psffile,'       1 !NTITLE'
        print>>psffile,'REMARKS Dummy PSF file to make VMD movies of trajectories look nice\n'
        print>>psffile,'%7d !NATOM'%(len(self.points))
        for iatom in range(len(self.points)):
            print>>psffile,'%8d%7d          C          0.000000       12.0110           0'%\
                                 (iatom+1,iatom+1)
        print>>psffile,''
        print>>psffile,'%8d !NBOND: bonds'%(len(self.points)-1)
        outstring=''
        for iatom in range(len(self.points)-1):
            newstring='%8d%8d'%(iatom+1,iatom+2)
            outstring=outstring+newstring
            if (iatom+1) % 4 == 0:
                print>>psffile,'%s'%outstring
                outstring=''
        print>>psffile,'%s'%outstring
        psffile.close()
    def end(self):
        return self.points[self.n_points-1]
    def final_vector(self):
        if self.n_points<1:
            print 'Cannot calculate final_vector for one link'
            exit(1)
        #print 'In final_vector, n_points=%d'%(self.n_points)
        apoint=point()
        apoint.x=self.points[self.n_points-1].x-self.points[self.n_points-2].x
        apoint.y=self.points[self.n_points-1].y-self.points[self.n_points-2].y
        apoint.z=self.points[self.n_points-1].z-self.points[self.n_points-2].z
        return apoint
    def penult_vector(self):
        if self.n_points<2:
            print 'Cannot calculate penult_vector for 2 links'
            exit(1)
        #print 'In final_vector, n_points=%d'%(self.n_points)
        apoint=point()
        apoint.x=self.points[self.n_points-2].x-self.points[self.n_points-3].x
        apoint.y=self.points[self.n_points-2].y-self.points[self.n_points-3].y
        apoint.z=self.points[self.n_points-2].z-self.points[self.n_points-3].z
        apoint.reverse()
        return apoint
        

# This is a utility class for managing distances
class dists:
    def __init__(self,ndists):
        self.ndists=ndists
        self.dists=[]
        for i in range(ndists):
            self.dists.append([])
    def add_value(self,distn,dist):
        self.dists[distn].append(dist)
    def printdists(self):
        f=open("dists.txt","w")
        for i in range(len(self.dists[0])):
            for j in range(self.ndists):
                f.write('% 10.5lf ' % self.dists[j][i])
            f.write('\n')

################ Main Program ##############
end_to_end=0
end_to_end_sqr=0
theta=theta_init
phi=phi_init

for itrial in range(ntrials):
    
    #Create polymer
    achain=chain()

    # Get noise for random walk in domain movement
    for i in range(ndomains):
        theta_domain_noise[i]=rshift.gauss(0.0,theta_domain_rate_sd[i])
        phi_domain_noise[i]=rshift.gauss(0.0,phi_domain_rate_sd[i])
        theta_domain_flip[i]=False
        phi_domain_flip[i]=False

    for i in range(nlinks-1):

        ## print "Trial %d - link %d - domain %d - noise %f"%(itrial,i,domain[i],theta_domain_noise[domain[i]])
        theta[i]=theta_last[i]+theta_domain_rate[domain[i]]+theta_domain_noise[domain[i]]
        phi[i]=phi_last[i]+phi_domain_rate[domain[i]]+phi_domain_noise[domain[i]]

        # Don't get too extended or we lose clarity for phi
        if (theta[i] < theta_buffer):
            theta[i] = theta_buffer
            theta_domain_flip[domain[i]]=True
        elif (theta[i] > theta_fixed):
            theta[i] = theta_fixed
            theta_domain_flip[domain[i]]=True

        # if (phi[i] < 0):
        #     phi[i] = 0.0
        #     phi_domain_flip[domain[i]]=True
        # elif (phi[i] > phi_fixed):
        #     phi[i] = phi_fixed
        #     phi_domain_flip[domain[i]]=True

        # Don't need to flip past pi because it changes the
        # handedness of the helix...
        if (phi[i] < phi_fixed):
            phi[i] = phi_fixed
            phi_domain_flip[domain[i]]=True
        elif (phi[i] > math.pi):
            phi[i] = math.pi
            phi_domain_flip[domain[i]]=True
            
        #Check to see if we've crossed zero or n*pi
        ##if (cross_pi(theta[i],theta_last[i])==1): 
            ##riemann_sheet[i+1]+=1
            #print 'Crossed pi theta=%lf theta_last=%lf'%(theta[i],theta_last[i])
	#print 'Chain after adding all links'
        #achain.printchain()
        #print '#########################################'

        flip_axis=0
        #print 'riemann_sheet[i]mod2=%d'%(riemann_sheet[i]%2)
        ##if (riemann_last[i]%2==0): flip_axis=0
	#if (flip_axis==1):
        #    print 'Flipping axis in timestep=%d, riemann_sheet=%d'%(itrial,riemann_sheet[i])
        temp_theta=abs(theta[i]+rnoise.gauss(0.0,theta_domain_sd[domain[i]])+rnoise.gauss(0.0,theta_sd))
        temp_phi=phi[i]+rnoise.gauss(0.0,phi_domain_sd[domain[i]])+rnoise.gauss(0.0,phi_sd) 

        # Temp constraints... the only real problem is theta getting too small...
        # or too close to pi...
        rel_theta = (temp_theta - (math.pi/2))
        temp_theta = rel_theta - (round(rel_theta / math.pi) * math.pi) + (math.pi/2)
        if (temp_theta < theta_buffer):
            temp_theta += theta_buffer
        elif (temp_theta > math.pi-theta_buffer):
            temp_theta -= theta_buffer

        achain.grow(temp_theta,temp_phi,flip_axis)
        theta_last[i]=theta[i]
        phi_last[i]=phi[i]
        temp_theta_last[i]=temp_theta
        temp_phi_last[i]=temp_phi
    #
    for i in range(nlinks):
        riemann_last[i]=riemann_sheet[i]

    for i in range(ndomains):
        if (theta_domain_flip[i]):
            theta_domain_rate[i] = -theta_domain_rate[i]
        if (phi_domain_flip[i]):
            phi_domain_rate[i] = -phi_domain_rate[i]

    #Calculate end-to-end distance
    length=achain.distance(0,nlinks-1)        
    end_to_end+=length
    end_to_end_sqr+=length*length
    if itrial==0:
        achain.printchainpsf()
    achain.printchainpdb(itrial)

    ## Output angles...
    # for i in range(1,len(theta)):
    #     print "%9.6f"%(theta[i]),
    # for i in range(2,len(phi)):
    #     print "%9.6f"%(phi[i]),
    # print


#print 'Ave end-to-end length=%8.3lf' % (end_to_end/ntrials)
#print 'Ave end-to-end^2 length=%8.3lf' % (end_to_end_sqr/ntrials)
