import numpy
import vtk
import math
from vtk.util import numpy_support

def makeVTKarray(array):
    return numpy_support.numpy_to_vtk(numpy.array(array).ravel(),array_type=vtk.VTK_FLOAT)

def get_labels(nt,ne,labels):
    valstr,errstr = labels
      
    it_label = []
    if nt > 1 and ne> 1:
       for it in range(nt-1):
          ie_label = []
          for ie in range(ne-1):
             val_label = f'{valstr}_{it+1:03d}_{ie+1:03d}'
             err_label = f'{errstr}_{it+1:03d}_{ie+1:03d}'
             ie_label.append((val_label, err_label))
          val_label = f'{valstr}_{it+1:03d}_Tot'
          err_label = f'{errstr}_{it+1:03d}_Tot'  
          ie_label.append((val_label, err_label))
          it_label.append(ie_label)
       
       ie_label = []
       for ie in range(ne-1):
           val_label = f'{valstr}_Tot_{ie+1:03d}'
           err_label = f'{errstr}_Tot_{ie+1:03d}'
           ie_label.append((val_label, err_label))
       val_label = f'{valstr}_Tot_Tot'
       err_label = f'{errstr}_Tot_Tot'  
       ie_label.append((val_label, err_label))
       it_label.append(ie_label)
                
    elif nt > 1 :
       for it in range(nt-1):
          val_label = f'{valstr}_{it+1:03d}'
          err_label = f'{errstr}_{it+1:03d}'  
          ie_label = [(val_label, err_label)]
          it_label.append(ie_label)
       
       val_label = f'{valstr}_Tot'
       err_label = f'{errstr}_Tot'
       ie_label  = [(val_label, err_label)]
       it_label.append(ie_label)                       
        
    elif ne > 1 :
       ie_label = []
       for ie in range(ne-1):
          val_label = f'{valstr}_{ie+1:03d}'
          err_label = f'{errstr}_{ie+1:03d}'  
          ie_label.append((val_label, err_label))
         
       
       val_label = f'{valstr}_Tot'
       err_label = f'{errstr}_Tot'
       ie_label.append((val_label, err_label))
       it_label.append(ie_label)                       

    else:
       val_label = f'{valstr}'
       err_label = f'{errstr}'  
       ie_label = [(val_label, err_label)]
       it_label.append(ie_label)

    return it_label   

def rectilinear_grid(mesh, labels=('value','error')):
    
    rgrid = vtk.vtkRectilinearGrid()
    rgrid.SetDimensions(mesh.nx1+1,mesh.nx2+1,mesh.nx3+1) 

    bin1 = makeVTKarray(mesh.x1bin)
    bin2 = makeVTKarray(mesh.x2bin)
    bin3 = makeVTKarray(mesh.x3bin)      

    rgrid.SetXCoordinates(bin1)
    rgrid.SetYCoordinates(bin2)
    rgrid.SetZCoordinates(bin3)     

    pdata = rgrid.GetCellData()  
    
    bin_labels = get_labels(mesh.nt,mesh.ne,labels)

    for it in range(mesh.nt):
       for ie in range(mesh.ne):
          vtkVal = makeVTKarray(mesh.data[it,ie,:,:,:,0]) 
          vtkErr = makeVTKarray(mesh.data[it,ie,:,:,:,1]) 
          vtkVal.SetName(bin_labels[it][ie][0])
          vtkErr.SetName(bin_labels[it][ie][1])
          pdata.AddArray(vtkVal)
          pdata.AddArray(vtkErr)
    return rgrid  

def structured_grid(mesh, trsf=None, labels=('value','error')):
    
    if trsf is None:
        origin = (0.,0.,0.)
    else:
        origin = trsf[0]    
    
    
    sgrid = vtk.vtkStructuredGrid()
    sgrid.SetDimensions(mesh.nx1+1,mesh.nx2+1,mesh.nx3+1) 

    pts = vtk.vtkPoints()
    pts.SetDataTypeToFloat()
    
    if mesh.geom == 'cyl' :
       if mesh.x3bin[1] > 0.5:
           mesh.x3bin = numpy.insert(mesh.x3bin,1,0.5,axis=0)
           mesh.nx3 += 1
           mesh.data = numpy.insert(mesh.data,0,mesh.data[:,:,0,:,:,:], axis=2)

    npts = (mesh.nx1+1)*(mesh.nx2+1)*(mesh.nx3+1) 
    pts.SetNumberOfPoints(npts)
    
    n=0
    if mesh.geom == 'cyl':
       for t in mesh.x3bin:  
           st = math.sin(t*2*math.pi)
           ct = math.cos(t*2*math.pi)
           for z in mesh.x2bin:
               z0 = z + origin[2]
               for r in mesh.x1bin:
                   x0 = r*ct + origin[0]
                   y0 = r*st + origin[1]
                   pts.InsertPoint(n,x0, y0, z0)
                   n += 1     
    else:
       for z in mesh.x3bin:
           z0 = z + origin[2]
           for y in mesh.x2bin:
               y0 = y + origin[1]
               for x in mesh.x1bin:
                   pts.InsertPoint(n,x+origin[0], y0, z0)
                   n += 1  
        
    sgrid.SetPoints(pts) 
 
    pdata = sgrid.GetCellData()  
    bin_labels = get_labels(mesh.nt,mesh.ne,labels)

    for it in range(mesh.nt):
       for ie in range(mesh.ne):
          vtkVal = makeVTKarray(mesh.data[it,ie,:,:,:,0]) 
          vtkErr = makeVTKarray(mesh.data[it,ie,:,:,:,1]) 
          vtkVal.SetName(bin_labels[it][ie][0])
          vtkErr.SetName(bin_labels[it][ie][1])
          pdata.AddArray(vtkVal)
          pdata.AddArray(vtkErr)
    return sgrid  