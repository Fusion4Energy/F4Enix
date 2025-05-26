from mesh2vtk_2 import MeshtalParser
from mesh2vtk_2 import FMesh

meshtal = r'meshes\meshtam'
meshtam = r'meshes\meshtan'

mesh1 = MeshtalParser(meshtal)
mesh2 = MeshtalParser(meshtam)

data1 = mesh1.get_FMesh(64)
data2 = mesh2.get_FMesh(64)

# operate with meshdata
sum_data = 10 * data1 + 20 * data2


# convert MeshData to Fmesh object
sum_mesh = FMesh(sum_data, 'sum_mesh')

# write vtk mesh
sum_mesh.write_vtk('sum_mesh')
