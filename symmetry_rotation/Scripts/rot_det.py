# Backend Code for Rotation Detection of Segmented Vertebrae
# best if used in 3d Slicer Jupyter env

'''
Author: Dillon Haughton
Version: 1
'''

# Necessary Modules. Make sure to install in Slicer env
import slicer
import pandas as pd
import numpy as np
import SegmentStatistics
import vtk
import vtk.util.numpy_support
#import k3d
import trimesh
import scipy
import vg
import matplotlib
import matplotlib.pyplot as plt
#import JupyterNotebooksLib as slicernb
import os



#=========================================================================================
# Following equation derived from 
# "Finding Mirror Symmetry via Registration and Optimal Symmetric Pairwise Assignment of Curves"
# Cicconet et al.
def Sv(normal):
	# These must be unit vectors to work.
	normal = normal / np.linalg.norm(normal) # normal = normal of mirror plane
	return np.identity(3) - 2*np.outer(normal, np.transpose(normal))

def distance(point, normal):
	normal = normal / np.linalg.norm(normal)
	return np.dot(np.subtract(point, np.array([0,0,0])), normal)

def reflect_points(x,p,v):
	v = v / np.linalg.norm(v)
	# Reflects all points of x across plane with point (p) and normal (v)
	# This function may be slow depending on the size of x. 
	# Points of mesh are reflected and triangles / faces still coorelate.
	reflect_coord = []
	for i in x:
		reflect_point = Sv(v).dot(i) + 2*distance(p,v)*v
		reflect_coord += [reflect_point]
	return reflect_coord     

def reflect_vector(V,N):
	R = V -2*(np.dot(V,N))*N
	return R


#=========================================================================================

def get_projection(axis, normal):
	# Get element of an axis vector projected onto surface using normal of surface.
	v = (axis.dot(normal)) * normal
	proj = axis - v
	return proj


class vertebrae:
	def __init__(self, volume, segmentation_file, segments=['all']):
		self.volume = volume
		self.segmentation_file = segmentation_file
		self.segments = segments
		#self.plot = k3d.plot(grid_visible=False)
		self.large_data = pd.DataFrame({'Volume':[], 'SegmentID':[], 'trimesh':[], 'Princ_x':[], 'Princ_y':[], 'Princ_z':[], 'Center':[]})
		self.data = pd.DataFrame({'Volume':[], 'SegmentID':[]})

		scan = self.volume
		#print(active_read_log)

		# Read log of current files in Data Folder
		# Limit analysis to exclude Cervicals
		#to_analyze = ['Segment_{}'.format(i) for i in np.arange(8,26,1)]

		

		'''
		try:
			active_read_log.remove('.DS_Store')
		except:
			pass

		# Limit number you want to include in current run
		for j, ww in enumerate(active_read_log):
			break_up = ww.split('/')
			print(break_up)
			if scan in break_up:
				number = j
			else:
				pass
				
		active_read_log = [active_read_log[number]]	# Error here is no number produced.
		'''
		# Initiate loop over values you want to look at
		'''
		for i in active_read_log:
			
			# Will find appropriate files since each has their own format
			get_File    = os.listdir(os.path.join('Data',i))
			try:
				get_File.remove('.DS_Store')
			except:
				pass
			

			# This part has issues across folder formats. Should just 
			# input the files naturally. 
			for ll in get_File:
				if ll.split('.')[-2] == 'nii':

					if ll.split('_')[-1] == 'seg.nii.gz':
						Segments = os.path.join(os.getcwd(), 'Data/{}/{}'.format(i,ll))

					elif ll.split('_')[-1] == 'image.nii':
						Volume = os.path.join(os.getcwd(), 'Data/{}/{}'.format(i,ll))	
					
					#else:
					#	Volume = os.path.join(os.getcwd(), 'Data/{}/{}'.format(i,ll))
				
			#Volume      = os.path.join(os.getcwd(), 'Data/{}/{}_CT-iso.nii.gz'.format(i,i))
			#Segments    = os.path.join(os.getcwd(), 'Data/{}/{}_CT-iso_seg.nii.gz'.format(i,i))

			# Load Volumes and Segmentations
			#print('Loading Volumes: ', Volume.split('/')[-1], ' with Segment: ', Segments.split('/')[-1])
			#print('------------------------------------------------------------------------------------')
			'''


		#slicer.util.loadVolume(self.volume)
		#seg = slicer.util.loadSegmentation(self.segmentation_file)

		displayNode = self.segmentation_file.GetDisplayNode()
		displayNode.SetAllSegmentsVisibility(True)

		# Load Proper Statistics to Calculate Principal Inertial Axes
		#print('Initializing Parameter acquisition')
		segStatLogic = SegmentStatistics.SegmentStatisticsLogic()
		segStatLogic.getParameterNode().SetParameter("Segmentation", self.segmentation_file.GetID())
		segStatLogic.getParameterNode().SetParameter("LabelmapSegmentStatisticsPlugin.centroid_ras.enabled",str(True))
		segStatLogic.getParameterNode().SetParameter("LabelmapSegmentStatisticsPlugin.principal_axis_x.enabled",str(True))
		segStatLogic.getParameterNode().SetParameter("LabelmapSegmentStatisticsPlugin.principal_axis_y.enabled",str(True))
		segStatLogic.getParameterNode().SetParameter("LabelmapSegmentStatisticsPlugin.principal_axis_z.enabled",str(True))
		segStatLogic.computeStatistics()
		stats = segStatLogic.getStatistics()

		if self.segments != ['all']:
			to_analyze = self.segments
		else:
			to_analyze = list(stats["SegmentIDs"])
			self.segments = to_analyze
		#print(type(to_analyze))	

		# Store center mass, and 3 principle inertial axes for each segment
		# and save them into pandas dataframe
		
		#print('Initializing Loop')
		#print()
		#frame = pd.DataFrame({'Volumes':[], 'SegmentID':[], 'Center':[], 'Princ_x':[], 'Princ_y':[], 'Princ_z':[]})


		for segmentID in stats["SegmentIDs"]:


			if segmentID in to_analyze:
			#if segmentID == 'Segment_20' or segmentID == 'Segment_21' or segmentID == 'Segment_22' or segmentID == 'Segment_23' or segmentID == 'Segment_24' or segmentID == 'Segment_25':
			#if segmentID == 'Segment_22':
			#try:
				center = np.array(stats[segmentID,"LabelmapSegmentStatisticsPlugin.centroid_ras"])
				princ_x = np.array(stats[segmentID,"LabelmapSegmentStatisticsPlugin.principal_axis_x"])
				princ_y = np.array(stats[segmentID,"LabelmapSegmentStatisticsPlugin.principal_axis_y"])
				princ_z = np.array(stats[segmentID,"LabelmapSegmentStatisticsPlugin.principal_axis_z"])

				all_vecs = np.array([princ_x, princ_y, princ_z])

				# Need to confirm each principal axis. ** want them standardized pointing one direction **
				# Need to get away from principal axis. First one should be easy enough and can mostly stick with princ_y.
				# Would like one that points from centroid to tip of spinous process but its unlikely to find a good other vector

				# z symetry should be max min from bounds centered at centroid.
				# x symetry should be from the orthogonal of that one. two much symetry in other planes to miss these.

				princ_x = all_vecs[np.argmax(np.abs(all_vecs[:,1]))]
				princ_y = all_vecs[np.argmax(np.abs(all_vecs[:,0]))]
				princ_z = all_vecs[np.argmax(np.abs(all_vecs[:,2]))]


				# Create 3D representation of Segment and convert into x,y,z coordinates
				self.segmentation_file.CreateClosedSurfaceRepresentation()
				segmentPolyData = vtk.vtkPolyData()
				self.segmentation_file.GetClosedSurfaceRepresentation(segmentID, segmentPolyData)
				polyData = segmentPolyData.GetPolys().GetData()

				# Retrieve triangles of mesh
				values = [int(polyData.GetTuple1(i)) for i in range(polyData.GetNumberOfTuples())]
				triangles = []
				while values:
					n = values[0] # number of points in the polygon
					triangles.append(values[1:n+1])
					del values[0:n+1]

				# Retrieve x,y,z point data
				pointData = segmentPolyData.GetPoints().GetData()
				pointCoordinates = vtk.util.numpy_support.vtk_to_numpy(pointData)

				# Place pointCoordinates' center at the origin
				pointCoordinates = pointCoordinates - center

				# Store as trimesh python object
				tri_mesh = trimesh.Trimesh(pointCoordinates, np.array(triangles))

				# For each segment. Storage of Volume path, Segment path, Segment ID, Center, Principle Axis, and trimesh 
				to_frame = pd.DataFrame({'Volume':[scan], 'SegmentID':[segmentID], 'trimesh':[tri_mesh], 'Princ_x':[princ_x], 
										 'Princ_y':[princ_y], 'Princ_z':[princ_z], 'Center':[center]})

				self.large_data = pd.concat([self.large_data, to_frame], ignore_index=True)
				self.data = pd.concat([self.data, to_frame[['Volume', 'SegmentID']]])
				#print('Successful for ', segmentID)
			#except:
			#    print('Unable to generate model for: ', segmentID)

		#print('Complete!')
		

		# Data acquisition moving forward

		#self.data = pd.DataFrame({'SegmentID':[],
		#                     	  'FE':[], 'Rot':[], 'Side':[], 'Warning':[],
		#                     	  'norm_y':[], 'y':[],
		#                     	  'norm_x':[], 'x':[],
		#                     	  'norm_z':[], 'z':[],
		#                   		  'max_zb':[], 'min_zb':[], 'center':[], 'box1':[], 'cylinder':[],
		#                     	  'line1':[], 'line2':[], 'line3':[], 'v_body':[],
		#                      	  'proj_x':[], 'proj_y':[], 'proj_z':[]})



		#segmentIDs = frame['SegmentID'].values
		#segmentIndex= frame['SegmentID'].index

	def calculate_symmetry_rot(self):	    
		# Limit to omit Cervicals defined in the begining
		frame = self.large_data
		#print(type(self.segments))
		slice_view = self.large_data[self.large_data['SegmentID'].isin(self.segments)]
		#slice_view = frame[frame['SegmentID'].isin(['Segment_18'])]

		segmentIDs   = slice_view['SegmentID'].values
		segmentIndex = slice_view.index
		rots = []
		bodys= []
		spinos=[]
		sym_vs=[]
		sym_ps=[]
		for i,j in zip(segmentIndex, segmentIDs):

			warning = ''
			#print('-------------- ',j,' --------------')

			looking_at = frame.loc[i, 'trimesh'] 
			trimesh.repair.fix_inversion(looking_at)

			#print('1.1 Reflecting ')
			reflect_y = reflect_points(looking_at.vertices, np.array([0,0,0]), frame.loc[i,'Princ_y'])
			tri_reflect = trimesh.Trimesh(reflect_y, looking_at.faces)
			trimesh.repair.fix_inversion(tri_reflect)

			#print('2.1 Registration by ICP ')
			registration_y = trimesh.registration.icp(tri_reflect.vertices, looking_at.vertices)
			registration_roty = registration_y[0][:3,:3]
			registration_trany = registration_y[0][:3,-1]

			svroty = np.dot(Sv(frame.loc[i,'Princ_y']), np.transpose(registration_roty))
			eigvaly, eigvecy = np.linalg.eig(svroty)

			#print('3.1 Determing Y-Symmetry Plane on entire vertebrae')
			# Syms did have
			sym_normaly = eigvecy[:,list(np.where(np.min(eigvaly)))[0][0]].real
			sym_pointy  = 0.5 * (np.dot(registration_roty, 2* distance(np.array([0,0,0]),frame.loc[i,'Princ_y'])*frame.loc[i,'Princ_y']) + registration_trany).real

			sym_vs += [sym_normaly]
			sym_ps += [sym_pointy]

			#print('3.2 Projecting Y axis on symmetry plane for Rotation ')
			projection_y = get_projection(np.array([0,1,0]), sym_normaly)
			#print('3.3 Projecting Z axis on symmetry plane for Side Bending ')

			rots += [vg.signed_angle(np.array([0,1,0]), projection_y.real, look=np.array([0,0,-1]))]


			#print('3.4 Rotation calculated succesfully ')
			if np.abs(rots[i]) > 20:
				print('WARNING: Obscene Rotation detected!!')
				warning += ' Obscene Rotation '

			slicedup = trimesh.intersections.mesh_plane(looking_at, sym_normaly, sym_pointy)
			slice_k3d = trimesh.load_path(slicedup)
			center_path = slice_k3d.centroid
			both = slice_k3d.discrete
			
			# Will yield two slices one for the spinous process and one for vertebral body
			both_sum = []
			for y in both:
				
				both_sum += [np.max(y[:,1])]

			both_sum = np.array(both_sum)    
			# Determine which is vertebral body by seeing which is more positive
			body = both[np.argsort(both_sum)[-1]]
			spinous = both[np.argsort(both_sum)[-2]]  

			bodys += [body]
			spinos+= [spinous]  

			#frame = pd.DataFrame({'Rotation':[rots]})
		self.data['Rotation'] = rots      
		self.large_data['Rotation'] = rots
		self.large_data['Sym_V'] = sym_vs
		self.large_data['Sym_P'] = sym_ps
		self.large_data['Body'] = bodys
		self.large_data['Spinous'] = spinos

		

		#self.plot += k3d.line(body + center, name='body', width=0.3)
		#self.plot += k3d.line(spinous + center, name='spinous',width=0.3)

	def show_plot(self):
		# Will need to update as features continue to develop.
		values = self.large_data
		# Can skip this since the segmentation are loaded already.
		#if self.large_data.columns.isin(['trimesh']).any():
		#	for index, row in self.large_data.iterrows():
		#		#print(row)
		#		center = row['Center']
		#		mesh = row['trimesh']
		#		self.plot += k3d.mesh(mesh.vertices + center, mesh.faces,color=0xcaeaff)


		if self.large_data.columns.isin(['Rotation']).any():
			# Displaying this could be feasibly done bt loading the 
			for index, row in self.large_data.iterrows():
				#print(row)
				center = row['Center']
				sym_point = row['Sym_P']
				body = row['Body']
				spines=row['Spinous']

				combined = np.concatenate([body , spines])
				combined = combined + center


				vtk_points1 = vtk.vtkPoints()
				vtk_points1.SetData(vtk.util.numpy_support.numpy_to_vtk(combined))


				polydata1 = vtk.vtkPolyData()
				polydata1.SetPoints(vtk_points1)


				# Create the vtkSphereSource object.
				sphere1 = vtk.vtkSphereSource()
				sphere1.SetRadius(1.0)

				# Create the vtkGlyph3D object.
				glyph1 = vtk.vtkGlyph3D()
				glyph1.SetInputData(polydata1)
				glyph1.SetSourceConnection(sphere1.GetOutputPort())


				pointCloudModelNode1 = slicer.modules.models.logic().AddModel(glyph1.GetOutputPort())

		if self.large_data.columns.isin(['FE']).any():
			# This is spitting out the inner ring of the vertebrae
			for index, row in self.large_data.iterrows():
				#print(row)
				center = row['Center']
				body_z = row['body_z']
				body_x = row['body_x']
				for i in range(len(row['body_z'])):
					try:
						body_z = body_z[i] + center

					#combined = np.concatenate([body_z , body_x])
					#combined = combined + center


						vtk_points2 = vtk.vtkPoints()
						vtk_points2.SetData(vtk.util.numpy_support.numpy_to_vtk(body_z))


						polydata2 = vtk.vtkPolyData()
						polydata2.SetPoints(vtk_points2)


						# Create the vtkSphereSource object.
						sphere2 = vtk.vtkSphereSource()
						sphere2.SetRadius(1.0)

						# Create the vtkGlyph3D object.
						glyph2 = vtk.vtkGlyph3D()
						glyph2.SetInputData(polydata2)
						glyph2.SetSourceConnection(sphere2.GetOutputPort())	
					except:
						pass		

					pointCloudModelNode2 = slicer.modules.models.logic().AddModel(glyph2.GetOutputPort())

				#============================
				for i in range(len(row['body_x'])):
					try:
						body_x = body_x[i] + center

						vtk_points3 = vtk.vtkPoints()
						vtk_points3.SetData(vtk.util.numpy_support.numpy_to_vtk(body_x))


						polydata3 = vtk.vtkPolyData()
						polydata3.SetPoints(vtk_points3)


						# Create the vtkSphereSource object.
						sphere3 = vtk.vtkSphereSource()
						sphere3.SetRadius(1.0)

						# Create the vtkGlyph3D object.
						glyph3 = vtk.vtkGlyph3D()
						glyph3.SetInputData(polydata3)
						glyph3.SetSourceConnection(sphere3.GetOutputPort())	

						pointCloudModelNode3 = slicer.modules.models.logic().AddModel(glyph3.GetOutputPort())
					except:
						pass	

				#self.plot += k3d.line(body + center) 		
				#self.plot += k3d.line(spines + center)

				#origins = np.array(list(sym_point + center))
				#axises  = np.array([0,1,0]) * 40

				#self.plot +=  k3d.vectors(origins, axises, name='Axis')

		# Can implement later is necessary		
		#if self.large_data.columns.isin(['V_Body']).any():
		#	for index, row in self.large_data.iterrows():
		#		#print(row)
		#		center = row['Center']
		#		mesh = row['V_Body']
		#		self.plot += k3d.mesh(mesh.vertices + center, mesh.faces,color=0x42EBFF)		
		print(self.data)

		#self.plot.display()		


	def seperate_symmetry_body(self):
		frame = self.large_data
		slice_view = self.large_data[self.large_data['SegmentID'].isin(self.segments)]
		#slice_view = frame[frame['SegmentID'].isin(['Segment_18'])]

		segmentIDs   = slice_view['SegmentID'].values
		segmentIndex = slice_view.index
		bodys = []
		bodys_c = []

		for i,j in zip(segmentIndex, segmentIDs):
			looking_at = frame.loc[i, 'trimesh'] 
			trimesh.repair.fix_inversion(looking_at)

			sym_normaly = frame.loc[i, 'Sym_V']
			sym_pointy  = frame.loc[i, 'Sym_P']
			

		
			#print('4.1 Seperating Vertebral body for x and z based symmetry plane ')
			# Yield slice from symmetry plane
			slicedup = trimesh.intersections.mesh_plane(looking_at, sym_normaly, sym_pointy)
			slice_k3d = trimesh.load_path(slicedup)
			center_path = slice_k3d.centroid
			both = slice_k3d.discrete
				
				# Will yield two slices one for the spinous process and one for vertebral body
			both_sum = []
			for y in both:
					
				both_sum += [np.max(y[:,1])]

			both_sum = np.array(both_sum)    
				# Determine which is vertebral body by seeing which is more positive
			body = both[np.argsort(both_sum)[-1]]
			spinous = both[np.argsort(both_sum)[-2]]

				# Get box of vertebral body giving points to estimate initial values
				# for vertebral body seperation
			body_path = trimesh.load_path(body)
			body_path_2d = body_path.to_planar()
			bounds_2d = trimesh.bounds.oriented_bounds_2D(body_path_2d[0].vertices)
			box_2d = trimesh.path.creation.rectangle(bounds_2d[1])
			box    = box_2d.apply_transform(np.linalg.inv(bounds_2d[0]))
			box_3d = box.to_3D(transform = body_path_2d[1])
			bodys_center = box_3d.centroid
			box_3d = box_3d.discrete
			self.mob = (box_3d[0][0] + box_3d[0][1] + box_3d[0][2] + box_3d[0][3])/4

			max_ys = spinous[np.argmax(spinous[:,1])]
			min_yb = body[np.argmin(body[:,1])]
			mid_point = ((max_ys + min_yb)/2)

			box_3d = np.delete(box_3d, -1, axis=1)
			yb = box_3d[0][np.argsort(box_3d[0][:,1])[:2]]
			min_zb = yb[np.argsort(yb[:,2])[0]]
			max_zb = yb[np.argsort(yb[:,2])[1]]

			orthogs_r = np.cross(sym_normaly, max_zb - min_zb)
			orthogs_r = orthogs_r / np.linalg.norm(orthogs_r)

			# Initial planes to reflect for x and z
			self.for_x = max_zb - min_zb
			self.for_x = self.for_x / np.linalg.norm(self.for_x)
				
			start = (max_zb + min_zb)/2
			startv = mid_point - start
			startv = startv/np.linalg.norm(startv)
				
			self.for_z = np.cross(sym_normaly, self.for_x)
			self.for_z = self.for_z / np.linalg.norm(self.for_z)

			right = trimesh.intersections.slice_mesh_plane(looking_at,  sym_normaly, [mid_point[0], mid_point[1], mid_point[2]], cap=True)
			left  = trimesh.intersections.slice_mesh_plane(looking_at, -sym_normaly, [mid_point[0], mid_point[1], mid_point[2]], cap=True)

				
			def reduction_function_end():
				# This optimization is terrible need to adjust
				# Function to seperate vertebral body using two planes
				# correct position and orientation is determined by highest percent of 
				# volume / volume of bounding box
				index= 0 
				ratio = []
				values = []
				
				# These values may need to be tinkered with in order to optimize results
				# may try substuting reflection instead of orthogs.
				backwards = np.arange(0,20,1)
				inwards   = np.arange(-0.3,0.3,0.1)
				for lu in backwards:
					for wu in inwards:
						#plot = k3d.plot()
						plane1  = np.array([orthogs_r[0]-wu, orthogs_r[1], orthogs_r[2]])
						plane1 = plane1/np.linalg.norm(plane1)

						plane2  = -reflect_vector(plane1, -self.for_z)
						plane2  = plane2/np.linalg.norm(plane2)
						
						center = lu*(startv) + start
						initial1 = trimesh.intersections.slice_mesh_plane(right, -plane1, center, cap=True)
						initial2 = trimesh.intersections.slice_mesh_plane(left, -plane2, center, cap=True)

						final = trimesh.util.concatenate(initial1, initial2)
						#plot += k3d.points(center)
						#plot += k3d.mesh(final.vertices, final.faces)
						#plot.display()
						
						volume  = final.volume
						bounds  = trimesh.bounds.oriented_bounds(final)[1]

						V       = np.prod(bounds)
						surface = volume/V
						values += [[lu, wu]]
						ratio += [surface]
					
					
				find = np.argmax(ratio)
				return values[find]

			lu, wu = reduction_function_end()   

			orthogs_r = np.array([orthogs_r[0]-wu, orthogs_r[1], orthogs_r[2]])
			orthogs_r = orthogs_r / np.linalg.norm(orthogs_r)
			orthogs_l = -reflect_vector(orthogs_r, -self.for_z)
			#orthogs_l = 2 * np.dot(sym_normaly, orthogs_r) * sym_normaly - orthogs_r
			#orthogs_l = np.array([orthogs_r[0]+wu, orthogs_r[1], orthogs_r[2]])
			orthogs_l = orthogs_l/np.linalg.norm(orthogs_l)

			true_value = lu*(startv) + start

			first = trimesh.intersections.slice_mesh_plane(right, -orthogs_r, true_value, cap=True)
			second = trimesh.intersections.slice_mesh_plane(left, -orthogs_l, true_value, cap=True)

			v_body = trimesh.util.concatenate(first, second)
			
			v_body_center = v_body.center_mass
			
			bodys += [v_body]
			bodys_c += [v_body_center]

		self.large_data['V_Body'] = bodys
		self.large_data['V_Body_Center'] = bodys_c

	
	def export_data(self, append_to):
		# Only want part of the volume as input
		try:
			appendage = pd.read_csv(append_to, index_col=0)
		except:
			print('No File: Generating Data Log')
			appendage = pd.DataFrame({'Volume':[], 'SegmentID':[], 'Rotation':[], 'L/R':[]})	
			appendage.to_csv(append_to)

		gather = pd.concat([appendage, self.data], ignore_index=True)
		print(gather)
		gather.to_csv(append_to)



	def calculate_symmetry_side_flex(self):
		# Requires previously calculated vertebral body
		# Limit to omit Cervicals defined in the begining
		frame = self.large_data
		slice_view = self.large_data[self.large_data['SegmentID'].isin(self.segments)]

		segmentIDs   = slice_view['SegmentID'].values
		segmentIndex = slice_view.index
		sides = []
		flexes = []
		bodys_z= []
		bodys_x= []
		
		for i,j in zip(segmentIndex, segmentIDs):

			warning = ''
			#print('-------------- ',j,' --------------')

			looking_at = frame.loc[i, 'trimesh'] 
			v_body = frame.loc[i, 'V_Body']
			v_body_center = frame.loc[i, 'V_Body_Center']
			trimesh.repair.fix_inversion(v_body)

			#print('5.1 Reflecting Vertebral Body')
			reflect_x = reflect_points(v_body.vertices, np.array(self.mob), np.array(self.for_x))
			reflect_z = reflect_points(v_body.vertices, np.array(self.mob), np.array(self.for_z))
			tri_reflectx = trimesh.Trimesh(reflect_x, v_body.faces)
			tri_reflectz = trimesh.Trimesh(reflect_z, v_body.faces)
			trimesh.repair.fix_inversion(tri_reflectx)
			trimesh.repair.fix_inversion(tri_reflectz)

			#print('6.1 Registering Vertebral Body by ICP ')
			# This section seems to have the most trouble may need to loop it with different 
			# initial ICP values in order to get better results. 
			registration_x = trimesh.registration.icp(tri_reflectx.vertices, v_body.vertices)
			registration_z = trimesh.registration.icp(tri_reflectz.vertices, v_body.vertices)

			registration_rotx = registration_x[0][:3,:3]
			registration_tranx = registration_x[0][:3,-1]

			registration_rotz = registration_z[0][:3,:3]
			registration_tranz = registration_z[0][:3,-1]

			svrotx = np.dot(Sv(self.for_x), np.transpose(registration_rotx))
			eigvalx, eigvecx = np.linalg.eig(svrotx)

			svrotz = np.dot(Sv(self.for_z), np.transpose(registration_rotz))
			eigvalz, eigvecz = np.linalg.eig(svrotz)

			sym_normalx = eigvecx[:,list(np.where(np.min(eigvalx)))[0][0]].real
			sym_pointx  = 0.5 * (np.dot(registration_rotx, 2* distance(v_body_center,self.for_x)*self.for_x) + registration_tranx).real

			sym_normalz = eigvecz[:,list(np.where(np.min(eigvalz)))[0][0]].real
			sym_pointz  = 0.5 * (np.dot(registration_rotz, 2* distance(v_body_center,self.for_z)*self.for_z) + registration_tranz).real

			slicedup2 = trimesh.intersections.mesh_plane(v_body, sym_normalz, sym_pointz)
			slice_k3d2 = trimesh.load_path(slicedup2)
			path_discrete = slice_k3d2.discrete
			print(path_discrete)
			slicedup3 = trimesh.intersections.mesh_plane(v_body, sym_normalx, sym_pointx)
			slice_k3d3 = trimesh.load_path(slicedup3)
			path_discrete3 = slice_k3d3.discrete
			
			# Sometimes can get multiple lines in this cut
			if len(path_discrete) > 1:
				print('WARNING: Multiple cuts seen in z symetry plane')
				warning += ' Multiple cuts in z '

			#checking_x1 = vg.signed_angle(sym_normalx, sym_normaly, look=np.array([0,1,0]))
			#checking_x2 = vg.signed_angle(sym_normalx, sym_normalz, look=np.array([-1,0,0]))

			# This is the hardest warning to catch
			#if np.abs(checking_x1) <= 45 or np.abs(checking_x2) <= 45:
			#	print('WARNING: Close estimation between x,y,z')
			#	warning += ' Close estimation in x '

			#print('7.1 Symetry Planes of Vertebral Body determined ')

			#print('7.2 Projecting x axis onto x symmetry plane for Side Bending')
			projection_x = get_projection(np.array([1,0,0]), sym_normalx)
			#print('7.3 Projecting z axis onto z symmetry plane for Flexion/ Extension')
			projection_z = get_projection(np.array([0,0,1]), sym_normalz)

			flex = vg.signed_angle(np.array([0,0,1]), projection_z.real, look=np.array([-1,0,0]))
			side= vg.signed_angle(np.array([1,0,0]), projection_x.real, look=np.array([0,1,0]))	
			

			sides += [side]
			flexes += [flex]
			bodys_z+= [path_discrete]
			bodys_x+= [path_discrete3]

		self.data['Side'] = sides
		self.data['FE'] = flexes      
		self.large_data['Side'] = sides
		self.large_data['FE'] = flexes
		#self.large_data['Sym_V'] = sym_vs
		#self.large_data['Sym_P'] = sym_ps
		self.large_data['body_z'] = bodys_z
		self.large_data['body_x'] = bodys_x	
		

	def calculate_symmetry_flex(self, body):
		print('Upcoming update')	
'''
				
				print('5.1 Reflecting Vertebral Body')
				reflect_x = reflect_points(v_body.vertices, np.array(mob), for_x)
				reflect_z = reflect_points(v_body.vertices, np.array(mob), np.array(for_z))
				tri_reflectx = trimesh.Trimesh(reflect_x, v_body.faces)
				tri_reflectz = trimesh.Trimesh(reflect_z, v_body.faces)
				trimesh.repair.fix_inversion(tri_reflectx)
				trimesh.repair.fix_inversion(tri_reflectz)

				print('6.1 Registering Vertebral Body by ICP ')
				# This section seems to have the most trouble may need to loop it with different 
				# initial ICP values in order to get better results. 
				registration_x = trimesh.registration.icp(tri_reflectx.vertices, v_body.vertices)
				registration_z = trimesh.registration.icp(tri_reflectz.vertices, v_body.vertices)

				registration_rotx = registration_x[0][:3,:3]
				registration_tranx = registration_x[0][:3,-1]

				registration_rotz = registration_z[0][:3,:3]
				registration_tranz = registration_z[0][:3,-1]

				svrotx = np.dot(Sv(for_x), np.transpose(registration_rotx))
				eigvalx, eigvecx = np.linalg.eig(svrotx)

				svrotz = np.dot(Sv(for_z), np.transpose(registration_rotz))
				eigvalz, eigvecz = np.linalg.eig(svrotz)

				sym_normalx = eigvecx[:,list(np.where(np.min(eigvalx)))[0][0]].real
				sym_pointx  = 0.5 * (np.dot(registration_rotx, 2* distance(v_body_center,for_x)*for_x) + registration_tranx).real

				sym_normalz = eigvecz[:,list(np.where(np.min(eigvalz)))[0][0]].real
				sym_pointz  = 0.5 * (np.dot(registration_rotz, 2* distance(v_body_center,for_z)*for_z) + registration_tranz).real

				slicedup2 = trimesh.intersections.mesh_plane(looking_at, sym_normalz, sym_pointz)
				slice_k3d2 = trimesh.load_path(slicedup2)
				path_discrete = slice_k3d2.discrete
				
				slicedup3 = trimesh.intersections.mesh_plane(looking_at, sym_normalx, sym_pointx)
				slice_k3d3 = trimesh.load_path(slicedup3)
				path_discrete3 = slice_k3d3.discrete
				
				# Sometimes can get multiple lines in this cut
				if len(path_discrete) > 1:
					print('WARNING: Multiple cuts seen in z symetry plane')
					warning += ' Multiple cuts in z '

				checking_x1 = vg.signed_angle(sym_normalx, sym_normaly, look=np.array([0,1,0]))
				checking_x2 = vg.signed_angle(sym_normalx, sym_normalz, look=np.array([-1,0,0]))

				# This is the hardest warning to catch
				if np.abs(checking_x1) <= 45 or np.abs(checking_x2) <= 45:
					print('WARNING: Close estimation between x,y,z')
					warning += ' Close estimation in x '

				print('7.1 Symetry Planes of Vertebral Body determined ')

				print('7.2 Projecting x axis onto x symmetry plane for Side Bending')
				projection_x = get_projection(np.array([1,0,0]), sym_normalx)
				print('7.3 Projecting z axis onto z symmetry plane for Flexion/ Extension')
				projection_z = get_projection(np.array([0,0,1]), sym_normalz)

				flex = vg.signed_angle(np.array([0,0,1]), projection_z.real, look=np.array([-1,0,0]))
				side= vg.signed_angle(np.array([1,0,0]), projection_x.real, look=np.array([0,1,0]))


				append = pd.DataFrame({'SegmentID':[j],
									   'FE':[flex], 'Rot':[rot], 'Side':[side], 'Warning':[warning],
									   'norm_y':[sym_normaly], 'y':[sym_pointy],
									   'norm_x':[sym_normalx], 'x':[sym_pointx],
									   'norm_z':[sym_normalz], 'z':[sym_pointz],
									   'max_zb':[max_zb], 'min_zb':[min_zb], 'center':[mob], 'box1':[box_3d],
									   'line1':[body], 'line2':[spinous], 'line3':[path_discrete], 'v_body': [v_body],
									   'proj_x':[projection_x], 'proj_y':[projection_y], 'proj_z':[projection_z]})


				print('8.1 Updating Data')
				data = data.append(append, ignore_index=True)
				print()
				
				complete = pd.merge(frame, data, on='SegmentID')

				#except:
					#print('----------  Aborted due to error ----------')
					##print()
					## Needs to be dropped or everything will be off
					#warning += 'Abortion Detected'
					##frame = frame.drop(i)
					##print(frame)
			
				# zero offset from L5...

				lookat = j
				view = complete[complete['SegmentID'] == lookat]


				mesh = view['trimesh'].values[0]
				#center = view['center'].values[0]
				v_body = view['v_body'].values[0]

				line1 = view['line1'].values[0]
				line2 = view['line2'].values[0]
				line3 = view['line3'].values[0]
				sym_norm = view['norm_x'].values[0]
				sym_point= view['x'].values[0]

				sym_normz = view['norm_z'].values[0]
				sym_pointz = view['z'].values[0]

				sym_pointy= view['y'].values[0]
				projection_x = view['proj_x'].values[0]
				projection_y = view['proj_y'].values[0]
				projection_z = view['proj_z'].values[0]

				max_ys = view['max_zb'].values[0]
				min_yb = view['min_zb'].values[0]

				center = view['Center'].values[0]
				princ_x = view['Princ_x'].values[0]
				princ_y = view['Princ_y'].values[0]
				princ_z = view['Princ_z'].values[0]

				box1 = view['box1'].values[0]
				cylinder = view['cylinder'].values[0]

				#center = view['center'].values[0]

				origins = np.array(list(sym_pointz + center) + list(sym_pointy + center) + list(sym_point + center))
				axises  = np.array([0,0,1,0,1,0,1,0,0]) * 40
				projs   = np.array(list(projection_z) + list(projection_y) + list(projection_x)) * 40
				colors  = np.array([0xFF0000] * 6)
				

				p_ori   = np.array(list([0,0,0]) + list([0,0,0]) + list([0,0,0]))
				princ   = np.array(list(princ_x) + list(princ_y) + list(princ_z)) * 40
				p_color = np.array([0xFFFFFF, 0xFFFFFF, 0x00FF00, 0x00FF00, 0x00FFFF, 0x00FFFF])

				#plot = k3d.plot()
				plot += k3d.mesh(mesh.vertices + center, mesh.faces, opacity=0.3, name='Main Mesh', color=0xcaeaff)
				plot += k3d.mesh(v_body.vertices + center, v_body.faces, opacity=0.5, color=0x6db0ba, name='Vertebral Body')
				#plot += k3d.mesh(cylinder.vertices, cylinder.faces)
				#plot += k3d.points(reflect_z)
				plot += k3d.line(line1 + center, name='body', width=0.3)
				plot += k3d.line(line2 + center, name='spinous',width=0.3)
				plot += k3d.line(line3[0] + center, color=0x0000ff, name='z_plane', width=0.3)
				for lll in path_discrete3:
					plot += k3d.line(lll + center, color=0x0000ff, width=0.3)

				#plot += add_plane_k3d(origin = sym_point, normal = sym_norm, mesh = [-50,50,-50,50])
				#plot += add_plane_k3d(origin = true_value, normal = orthogs_r, mesh = [-50,50,-50,50])
				#plot += add_plane_k3d(origin = true_value, normal = orthogs_l, mesh = [-50,50,-50,50])
				#plot += k3d.points(box1)
				#plot += k3d.points(center)

				plot += k3d.vectors(origins, axises, name='Axis')
				plot += k3d.vectors(origins, projs, colors= colors, name='Projections')
				#plot += k3d.vectors(p_ori, princ, colors=p_color, name='Principals')

				

				
				print('FE: ',view['FE'].values[0])
				print('Rot: ',view['Rot'].values[0])
				print('Side: ',view['Side'].values[0])
				print('Warning: ',view['Warning'].values[0])
				print()
		   
			
			print('COMPLETE')
			plot.display()

			# Find a way to numerically compare the two halves to see how symmetric they are

			#masterRepresentation.Initialized()
			#masterRepresentation.Modified()


'''
