from dolfin import *
from mshr import *
from tube import *


'''
FunctionSpaces for the mesh domain. Also,local and global Node numbers and global Element numbers of the mesh. 
'''

#nodewise functionspace:
def FS_node(bmesh):
	
	FS=FunctionSpace(bmesh,'CG',1)
	return FS


#elementwise functionspace:
def FS_elem(bmesh):
	FS=FunctionSpace(bmesh,'DG',0)
	return FS

#vector Function space
def FS_vec(bmesh):
	FS=VectorFunctionSpace(bmesh, 'Lagrange', 2)
	return FS



# element number
def elem_num(bmesh):
	element_numbers = bmesh.num_cells()
	return element_numbers



#global node number (better for 1 processor)
def glob_node(bmesh):
	node_numbers = bmesh.num_vertices()
	return node_numbers


#local node number (better for more than 1 processor)
def loc_node(bmesh):
	local_range = FS_node(bmesh).dofmap().ownership_range()
	local_dim = local_range[1] - local_range[0]
	return local_dim

# coordinates of each node. 
def coord(bmesh):
	cor=FS_node(bmesh).tabulate_dof_coordinates()
	return cor







