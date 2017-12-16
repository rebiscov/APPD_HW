#!/usr/bin/env python2.7
import sys
import os
import math

# Link parameters
link_latency = "10us"
link_bandwidth = 10
link_bandwidth_unit = "Gbps"

# XML generation functions
def issueHead():
        head = ("<?xml version='1.0'?>\n"
                "<!DOCTYPE platform SYSTEM \"http://simgrid.gforge.inria.fr/simgrid/simgrid.dtd\">\n"
                "<platform version=\"4.1\">\n\n")

        config_clause = ("<!--  WARNING:  This <config></config> clause below\n"
                       "makes it so that NO COMPUTATION TIME is simulated. This is because\n"
                       "in this module, for pedagogic purposes, we don't want to muddy the\n"
                       "(simulation) waters with computational times. As a results, this\n"
                       "XML platform file may not be suitable for running other\n"
                       "simulations, unless you remove the <config></config> clause.\n"
                       "-->\n"
                       "<config>\n"
                       "<prop id=\"smpi/simulate-computation\" value=\"0\"></prop>\n"
                       "<prop id=\"smpi/host-speed\" value=\""+str(real_compute_power)+"\"></prop>\n"
                       "</config>\n\n")

        AS_head = "<zone id=\"AS0\" routing=\"Full\">\n"

        return head + config_clause + AS_head


def issueTail():
	return "</zone>\n</platform>\n"

def issueLink1(x):
	return "  <link id=\"link-"+str(x)+"\" latency=\""+str(link_latency)+"\" bandwidth=\""+link_bandwidth+"\"/>\n"

def issueLink2(x,y):
	return "  <link id=\"link-"+str(x)+"-"+str(y)+"\" latency=\""+str(link_latency)+"\" bandwidth=\""+link_bandwidth+"\"/>\n"

def issueLink3(x,y,bw):
	return "  <link id=\"link-"+str(x)+"-"+str(y)+"\" latency=\""+str(link_latency)+"\" bandwidth=\""+str(bw)+link_bandwidth_unit+"\"/>\n"

def issueHost(index):
	return "  <host id=\"host-"+str(index)+"."+hostname+"\" speed=\""+sim_compute_power+"\"/>\n"

def issueRouteHead(index1, index2):
	return "  <route src=\"host-"+str(index1)+"."+hostname+"\" dst=\"host-"+str(index2)+"."+hostname+"\">\n"
def issueRouteTail():
	return "  </route>\n"

def issueRouteLink1(x):
	return "\t<link_ctn id=\"link-"+str(x)+"\"/>\n"

def issueRouteLink2(x,y):
	return "\t<link_ctn id=\"link-"+str(x)+"-"+str(y)+"\"/>\n"

######################################################################
# Parse command-line arguments
if (len(sys.argv) != 6):
	print >> sys.stderr, "Usage: smpi-generate-clique.py <num hosts> <real-machine-compute-power> <simulation-node-compute-power> <simulation-link-bandwidth> <simulation-link-latency> \n"
	print >> sys.stderr, "Example: smpi-generate-clique.py 32 1000Gf 100Gf 10Gbps 10us \n"
	print >> sys.stderr, "  Will generate a clique_<num hosts>.xml and hostfile_<num hosts>.txt file\n"
	exit(1)

num_hosts = int(sys.argv[1])
sim_compute_power = sys.argv[2]+"Gf"
real_compute_power = int(sys.argv[3])*1000000000
link_bandwidth = sys.argv[4]
link_latency = sys.argv[5]
hostname = "nimportequoi.fr"

###############################################################
# Generate CLIQUE XML file
filename = "./clique-"+str(num_hosts)+"-platform.xml"
fh = open(filename, 'w')
fh.write(issueHead())

# Create hosts
for i in range(0,num_hosts):
	fh.write(issueHost(i))

# Create links
for i in range(0,num_hosts*num_hosts):
	fh.write(issueLink1(i))

# Create routes
for i in range (0,num_hosts):
	for j in range(i+1,num_hosts):
		fh.write(issueRouteHead(i,j))
  		fh.write(issueRouteLink1(i*num_hosts + j))
    		fh.write(issueRouteTail())

fh.write(issueTail())
fh.close()
print >> sys.stderr, "Clique XML platform file created: "+filename

###############################################################
## Generate host file
filename = "./clique-"+str(num_hosts)+"-hostfile.txt"
fh = open(filename, 'w')

for i in range(0,num_hosts):
	fh.write("host-"+str(i)+"."+hostname+"\n")

fh.close()
print >> sys.stderr, "Hostfile created: "+filename
