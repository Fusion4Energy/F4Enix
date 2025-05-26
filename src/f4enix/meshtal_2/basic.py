from mesh2vtk_2 import CDGSMeshParser, CUVMeshParser, MeshtalParser

cdgsfile = r"meshes\cdgs"
meshfile = r"meshes\meshtal"
cuvfile = r"meshes\cuvmsk2"

# open mesh files
cdgs = CDGSMeshParser(cdgsfile)
mesh = MeshtalParser(meshfile)
cuv = CUVMeshParser(cuvfile)

# get mesh tally from file
cdgs_1 = cdgs.get_FMesh(1)
mesh_24 = mesh.get_FMesh(24)
cuv_4104 = cuv.get_FMesh(4104)

# write mesh in vtk file (vtk extension is set by the writer method)
cdgs_1.write_vtk("vtk_cdgs")
mesh_24.write_vtk("vtk_mesh")
cuv_4104.write_vtk("vtk_cuv")
