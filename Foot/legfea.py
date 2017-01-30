#PFEA imports
import pfea.solver
from pfea.geom import octet

import numpy as np
from math import *

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

subdiv = 1

mat_matrix = [[[0,0,0,0,0,0,0],
			   [0,0,0,0,0,0,0],
			   [0,0,0,0,0,0,0],
			   [0,0,0,0,0,0,0],
			   [0,0,0,0,0,0,0]],
			  [[0,0,0,0,0,0,0],
			   [0,1,1,1,1,1,0],
			   [0,1,1,1,1,1,0],
			   [0,1,1,1,1,1,0],
			   [0,0,0,0,0,0,0]],
			  [[0,0,0,0,0,0,0],
			   [0,0,0,0,0,0,0],
			   [0,0,0,0,0,0,0],
			   [0,0,0,0,0,0,0],
			   [0,0,0,0,0,0,0]]]

if subdiv > 1:
	dims = np.shape(mat_matrix)
	new_mat_matrix = np.zeros(((dims[0]-2)*subdiv+2,(dims[1]-2)*subdiv+2,(dims[2]-2)*subdiv+2))

	for i in range(1,dims[0]-1):
		for j in range(1,dims[1]-1):
			for k in range(1,dims[2]-1):
				if mat_matrix[i][j][k] == 1:
					dex = ((i-1)*subdiv+1,(j-1)*subdiv+1,(k-1)*subdiv+1)
					for l in range(0,subdiv):
						for m in range(0,subdiv):
							for n in range(0,subdiv):
								new_mat_matrix[l+dex[0]][m+dex[1]][n+dex[2]] = 1

	#print(new_mat_matrix)

	mat_matrix = new_mat_matrix

# Material Properties

#Physical Voxel Properties
vox_pitch = 0.61/subdiv #m

frame_props = {"nu"  : 0.33, #poisson's ratio
			   "d1"	 : 0.02, #m
			   "d2"	 : 0.02, #m
			   "th"  : 0,
			   "E"   :  8e9, #Pa
			   "G"   :   42000000000,  #Pa
			   "rho" :  500, #kg/m^3
			   "beam_divisions" : 10,
			   "cross_section"  : 'rectangular',
			   "roll": 0}

nodes,frames,node_frame_map,dims = pfea.geom.octet.from_material(mat_matrix,vox_pitch)
#print(node_frame_map.shape)
frame_props["Le"] = pfea.geom.octet.frame_length(vox_pitch)

constraints = []
loadids = []
loads = [] 

nodes = np.array(nodes)
nodemin = np.min(nodes.T[2])
nodemax = np.max(nodes.T[2])

for nidnum,node in enumerate(nodes):
	if abs(nodemin-nodes[nidnum][2])<0.1*vox_pitch:
		constraints.append({'node':nidnum,'DOF':0, 'value':0})
		constraints.append({'node':nidnum,'DOF':1, 'value':0})
		constraints.append({'node':nidnum,'DOF':2, 'value':0})
		constraints.append({'node':nidnum,'DOF':3, 'value':0})
		constraints.append({'node':nidnum,'DOF':4, 'value':0})
		constraints.append({'node':nidnum,'DOF':5, 'value':0})
	if abs(nodemax-nodes[nidnum][2])<0.1*vox_pitch:
		loadids.append(nidnum)

totalload = 2500*9.8 # N
pernodeload = totalload/len(loadids)

for loadid in loadids:
	loads.append({'node':loadid,'DOF':2, 'value':-pernodeload})

#Group frames with their characteristic properties.
out_frames = [(np.array(frames),[],{'E'   : frame_props["E"],
								 	  'rho' : frame_props["rho"],
								 	  'nu'  : frame_props["nu"],
								 	  'd1'  : frame_props["d1"],
								 	  'd2'  : frame_props["d2"],
								 	  'th'  : frame_props["th"],
								 	  'beam_divisions' : frame_props["beam_divisions"],
								 	  'cross_section'  : frame_props["cross_section"],
								 	  'roll': frame_props["roll"],
								 	  'loads':{'element':0},
								 	  'prestresses':{'element':0},
								 	  'Le': frame_props["Le"]})]
#Format node positions
out_nodes = np.array(nodes)



global_args = {"gravity":[0,0,0]}

res_displace,C,Q = pfea.solver.analyze_System(out_nodes, global_args, out_frames, constraints,loads)

global_args["debug_plot"] = True

#print(Q)

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

if global_args["debug_plot"]:
	### Right now the debug plot only does x-y-z displacements, no twisting
	xs = []
	ys = []
	zs = []

	rxs = []
	rys = []
	rzs = []

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.set_aspect('equal')
	frame_coords = []
	factor = 10

	#print(matplotlib.projections.get_projection_names())
	for i,node in enumerate(nodes):
		xs.append(node[0])
		ys.append(node[1])
		zs.append(node[2])
		rxs.append(node[0]+res_displace[i][0]*factor)
		rys.append(node[1]+res_displace[i][1]*factor)
		rzs.append(node[2]+res_displace[i][2]*factor)

	frame_args = out_frames[0][2]
	st_nrg = Q.T[0,:]#0.5*frame_args["Le"]/frame_args["E"]*(Q.T[0,:]**2/frame_args["Ax"]+Q.T[4,:]**2/frame_args["Iy"]+Q.T[5,:]**2/frame_args["Iz"])
	st_nrg = st_nrg.T
	qmax = np.max(abs(st_nrg))
	#print st_nrg
	print qmax
	print np.mean(abs(st_nrg))

	factors = np.zeros((len(st_nrg),2))

	for i, strain in enumerate(st_nrg):
		if strain > 0:
			if abs(strain)/(pi**3*frame_args["E"]*frame_args["d1"]**4/(4*frame_args["Le"]**2)) > 1.0: 	
				print("**DANGER AT FRAME ID {0}".format(i))
			#print "element {0} safety factor".format(abs(strain)/(pi**3*frame_args["E"]*frame_args["d1"]**4/(4*frame_args["Le"]**2)))
			factors[i][0] = abs(strain)/(pi**3*frame_args["E"]*frame_args["d1"]**4/(4*frame_args["Le"]**2))
		if strain < 0:
			if abs(strain)/(37.5e6*frame_args["d1"]**2) > 1.0:
				print("**DANGER AT FRAME ID {0}".format(i))
			#print "element {0} safety factor".format(abs(strain)/(37.5e6*frame_args["d1"]**2))
			factors[i][1] = abs(strain)/(37.5e6*frame_args["d1"]**2)
	
	max_ten = np.max(factors.T[1])
	max_com = np.max(factors.T[0])
	for i,frame in enumerate(frames):
		nid1 = int(frame[0])
		nid2 = int(frame[1])
		start = [xs[nid1],ys[nid1],zs[nid1]]
		end   = [xs[nid2],ys[nid2],zs[nid2]]
		rstart = [rxs[nid1],rys[nid1],rzs[nid1]]
		rend   = [rxs[nid2],rys[nid2],rzs[nid2]]

		#ax.plot([start[0],end[0]],[start[1],end[1]],[start[2],end[2]],color='r', alpha=0.5)
		if(st_nrg[i] < 0 ):
			ax.plot([rstart[0],rend[0]],[rstart[1],rend[1]],[rstart[2],rend[2]],color='r', alpha=(1.0*factors[i][1]/max_ten)**3)
		else:
			ax.plot([rstart[0],rend[0]],[rstart[1],rend[1]],[rstart[2],rend[2]],color='b', alpha=(1.0*factors[i][0]/max_com)**3)

	'''
	for dframe in dframes:
		nid1 = int(dframe[0])
		nid2 = int(dframe[1])
		dstart = [dxs[nid1],rys[nid1],rzs[nid1]]
		dend   = [dxs[nid2],rys[nid2],rzs[nid2]]
		ax.plot([dstart[0],dend[0]],[dstart[1],dend[1]],[dstart[2],dend[2]],color='b', alpha=0.1)
	'''	

	axisEqual3D(ax)
	#ax.scatter(xs,ys,zs, color='r',alpha=0)
	#ax.scatter(rxs,rys,rzs, color='b',alpha=0.3)
	plt.show()
	#print(frames)