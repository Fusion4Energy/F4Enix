from mesh2vtk_2 import MeshtalParser
from mesh2vtk_2 import FMesh

filename = r'meshes\meshtam'
meshtal = MeshtalParser(filename)

mesh = meshtal.get_FMesh(54)

#get mesh corresponding to a specific energy and/or time bin

# return MeshData object with energy bin 0 and last time bin
etbin_data = mesh.get_etbin_data(0,mesh.nt-1)

# return MeshData object with energy last bin and all time bins
ebin_data = mesh.get_etbin_data(ebin=-1)

# return MeshData object with time last bin and all energy bins
tbin_data = mesh.get_etbin_data(tbin=-1)

# convert MeshData to Fmesh object
etbin_mesh = FMesh(etbin_data, 'e0_t-1_mesh')
ebin_mesh = FMesh(ebin_data, 'e_last_all_time_mesh')
tbin_mesh = FMesh(tbin_data, 't_last_all_time_mesh')

# write vtk mesh
mesh.write_vtk('mesh',binlabels=('time_energy','error'))
etbin_mesh.write_vtk('et_mesh',binlabels=('e0_t4','error'))
ebin_mesh.write_vtk('e_mesh',binlabels=('elast_time','error'))
tbin_mesh.write_vtk('t_mesh',binlabels=('tlast_energy','error'))
