import numpy
from .utils import ExtraBin
from .read_data import get_column_data, get_block_data
from .read_data import get_cuv_element, get_cdgs_element

def get_rotation_matrix(axs,vec):
    if vec is None:
        one = numpy.argmin(abs(axs))
        dummy = numpy.array((0,0,0))
        dummy[one] = 1.
        vecX = axs.cross(dummy)
        vecX /= numpy.linalg.norm(vecX)
    else:
        vecX = vec    
    
    vecY = numpy.cross(axs,vecX)
    rotmat = numpy.array((vecX,vecY,axs))
    return rotmat       

def get_transformation(line, cyl=True):
    
    if cyl :
        if "VEC direction" in line :
            fmcnp6 = True
            org0 = line.index('at') + 2
            org1 = line.index('axis')

            axs0 = line.index('axis in') + 7
            axs1 = line.index('direction')

            vec0 = line.index('direction',axs1 + 1) + 9
        else:
            fmcnp6 = False     
            org0 = line.index('at') + 2
            org1 = line.index(',')

            axs0 = line.index('axis in') + 7
            axs1 = line.index('direction')

        origin = numpy.array(line[org0:org1].split(), dtype=numpy.double)   
        axis   = numpy.array(line[axs0:axs1].split(), dtype=numpy.double) 
        if fmcnp6:
            vec = numpy.array(line[vec0:].split(), dtype=numpy.double)
        else:
            vec = None
    
        axisZ = abs(axis.dot(numpy.array((0,0,1.)))-1.0) < 1e-12
        if vec is None:
            axisX =None
        else:    
            axisX = abs(axis.dot(numpy.array((1.,0,0)))-1.0) < 1e-12 
        
        if axisZ and (axisX is None or axisX) :
            rotmat = None
        else: 
            rotmat = get_rotation_matrix(axis,vec)    
    else:
        org0 = 0
        org1 = 44
        vecx0 = 45
        vecx1 = 87
        vecy0 = 88
        origin = numpy.array(line[org0:org1].split(), dtype=numpy.double)
        vecX = numpy.array(line[vecx0:vecx1].split(), dtype=numpy.double)
        vecY = numpy.array(line[vecy0:].split(), dtype=numpy.double)
        axisX = abs(vecX.dot(numpy.array((1.,0,0.)))-1.0) < 1e-12    
        axisY = abs(vecY.dot(numpy.array((0,1.,0)))-1.0) < 1e-12 
        if axisX and axisY :
            rotmat = None
        else: 
            rotmat = get_rotation_matrix(vecX,vecY)    

    if numpy.linalg.norm(origin) < 1e-12 and rotmat is None :
           return None
    return (origin, rotmat)

def get_etbin_tag(line):
    if "Energy" in line:
        return "erg"
    elif "Decay time" in line:
        return "dtme"
    elif "Time" in line:
        return 'tme'
    elif "Nuclide" in line:
        return 'nuc'
    elif "Cell" in line:
        return 'cel'
    

def scan_meshfile(fic):
    tally=dict()
    while True:
        pos = fic.tell()
        line = fic.readline()[0:30]
        if " Mesh Tally Number" == line[0:18] :
            tnum = line[18:].strip()
            boundaries, data = get_mesh_block(fic)
            block_pos = (pos, boundaries, data)
            tally[tnum] = block_pos
        elif line == '': 
            break    
    return tally

def scan_cuvfile(fic):
    cuvmesh=dict()
    while True:
        pos = fic.tell()
        line = fic.readline()[0:30]
        if " Mesh Tally Number" == line[0:18] :
            tnum = line[18:].strip()
            boundaries, data = skip_cuvdata(fic)
            cuvmesh[tnum] = (pos, boundaries, data)
        elif line == '': 
            break     
    return cuvmesh

def scan_cdgsfile(fic):
    srcmesh=dict()
    line = fic.readline().split()
    nmesh= int(float(line[1]))
    line = fic.readline().split()
    total_strength = float(line[1])
     
    for i in range(nmesh) :
        pos = fic.tell()
        line = fic.readline()
        if "mesh_id" == line[0:7]:
            mshnum = line.split()[1]
            srcmesh[mshnum] = [pos]
            while True:
                pos = fic.tell()
                line = fic.readline()[0:11]
                if line == "source_data":
                    srcmesh[mshnum].append(pos)
                    break
            while True:  
                line = fic.readline()[0:15]
                if line == "end_source_data": break
    return srcmesh

def skip_cuvdata(fic):
    
    line = fic.readline()[0:30]
    while  line[0:22] != " Tally bin boundaries:":
        bound_pos = fic.tell()
        line = fic.readline()[0:30]
        if line == '':
            break
    
    geom, trsf, meshbins= get_mesh_boundaries(fic)
    x1bin, x2bin, x3bin, ebin, tbin = meshbins
    n1 = len(x1bin) - 1
    n2 = len(x2bin) - 1
    n3 = len(x3bin) - 1
    nelemts = n1*n2*n3

    data_pos = fic.tell()
    fic.skipline(1)
    for i in range(nelemts):
        headline = fic.readline().split()
        index = int(float(headline[0]))
        ncel = int(float(headline[-1]))
        fic.skipline(3*ncel)
    return (bound_pos, data_pos)    

def get_mesh_block(fic):
    blk_line = 0
    while blk_line < 2 :      
        line = fic.readline()
        pos = fic.tell()
        if line.strip(' ') == '\n': 
           if blk_line == 0 : 
               bound_pos = pos
           elif blk_line == 1: 
               data_pos = pos  
           blk_line += 1
        elif line == '':
            break
    return (bound_pos, data_pos)    

def get_header(fic,position):
    particle_list  = ('neutron', 'photon','electron','proton','deuteron')
    fic.seek(position)
    fic.skipline()
    comments = ''
    while True:
        line = fic.readline()
        if "mesh tally." in line : break
        comments += line
    if comments == '' : comments = None

    part = None
    for p in particle_list:
        if p in line:
            part = p
            break
    return part, comments    

def get_cdgsheader(fic,position):
    fic.seek(position)
    fic.skipline()
    comments = fic.readline()
    line = fic.readline().split()
    cooling = float(line[1])
    line = fic.readline().split()
    strength = float(line[1])

    return cooling, strength, comments
        

def get_mesh_type(fic,position,geom='rec', explicit_time=False):
    fic.seek(position)
    line = fic.readline()       

    if "Result" in line :
        return 'col'
    elif "GEOMETRY" in line :
        return 'oldCUV'
    
    if explicit_time :
        fic.skipline(4)
    else:
        fic.skipline()     
    line = fic.readline()
    if geom == 'rec' :
        if "Z bin" in line :
            return 'ij'
        elif "Y bin" in line:
            return 'ik'
        elif "X bin" in line:
            return 'jk'
        else:
            return 'bad'
    else:    
        if "Theta bin" in line :
            return 'ij'
        elif "Z bin" in line:
            return 'ik'
        elif "R bin" in line:
            return 'jk'
        else:
            return 'bad'        
        
def get_mesh_boundaries(fic, position=None, cuv=None):
    if position is not None:
        fic.seek(position)
        fic.skipline()
    line = fic.readline()       

    if "origin at" in line :
        trsf = get_transformation(line)
        geom = 'cyl'
        line = fic.readline()
    else:
        trsf = None
        geom = 'rec'    

    init = line.index(":") + 1
    x1bin = numpy.array(line[init:].split(), dtype=numpy.double)    

    line = fic.readline()
    init = line.index(":") + 1
    x2bin = numpy.array(line[init:].split(), dtype=numpy.double)    

    line = fic.readline()
    init = line.index(":") + 1
    x3bin = numpy.array(line[init:].split(), dtype=numpy.double)    

    line = fic.readline()
    if "number:" in line:
        line = fic.readline()
    
    init = line.index(":") + 1

    bindata = numpy.array(line[init:].split(), dtype=numpy.double)
    bintag  = get_etbin_tag(line[0:init])
    tmpbin = ExtraBin(bindata, bintag)

    line = fic.readline()
    if line.strip(' ') == '\n':
        ebin = tmpbin
        bindata = numpy.array((0,1e20), dtype=numpy.double)
        bintag  = 'tme'
        tbin = ExtraBin(bindata, bintag, explicitbin = False)
    else:
        tbin = tmpbin
        init = line.index(":") + 1
        bindata = numpy.array(line[init:].split(), dtype=numpy.double)
        bintag  = get_etbin_tag(line[0:init])
        ebin = ExtraBin(bindata, bintag)

    return geom, trsf, (x1bin, x2bin, x3bin, ebin, tbin)

def get_cdgsmesh_boundaries(fic,position):
    fic.seek(position)
    fic.skipline(4)
    line = fic.readline().split()
    ergbin = line[1] == 'bins'
    line = fic.readline().split()
    nerg = int(float(line[1]))
  
    etag = 'erg' if ergbin else 'line'
    ebin = ExtraBin(fic.get_multilines(nerg), etag)
    tbin = ExtraBin(numpy.array((0,1e20)), 'tme', explicitbin = False)

    line = fic.readline().split()
    geom = line[1]

    line = fic.readline().split()
    nx1 = float(int(line[1]))  
    nx2 = float(int(line[2]))
    nx3 = float(int(line[3]))

    line = fic.readline()
    trsf = get_transformation(line, cyl = geom=='cyl')

    x1bin = fic.get_multilines(nx1)
    x2bin = fic.get_multilines(nx2)
    x3bin = fic.get_multilines(nx3)

    return geom, trsf, (x1bin, x2bin, x3bin, ebin, tbin)


def get_mesh_data(fic, position, geom, explicit_time, shape):
    
    type = get_mesh_type(fic,position,geom, explicit_time)
    if type == 'col' or type == 'cf':
        data = get_column_data(fic,position,type,shape)
    elif type == 'oldCUV':
        # in mcnp5 version of d1suned CUV was written in meshtal files with CDGS format
        data = get_CDGSmesh_data(fic,position,shape,oldCUV=True)
    else :   
        data = get_block_data(fic,position,type,explicit_time,shape)
    return data

def get_CUVmesh_data(fic, position, shape, norm=None, filter=None):
    fic.seek(position)
    fic.skipline()

    nt,ne,nx3,nx2,nx1,_ = shape
    nelement = nx1*nx2*nx3
    data = numpy.ndarray(shape)
    val_array = numpy.ndarray((nelement,ne))
    err_array = numpy.ndarray((nelement,ne))
       
    for it in range(nt):
        for i in range(nelement):
            val,err = get_cuv_element(fic, norm, filter)
            val_array[i,:] = val[:]
            err_array[i,:] = err[:]
            
        val_array = numpy.transpose(val_array.reshape((nx1,nx2,nx3,ne)))
        err_array = numpy.transpose(err_array.reshape((nx1,nx2,nx3,ne)))
        data[it,:,:,:,:,0] = val_array[:,:,:,:]
        data[it,:,:,:,:,1] = err_array[:,:,:,:]        
    return data

def get_CDGSmesh_data(fic, position, shape, oldCUV=False):
    fic.seek(position)
    if oldCUV:
        fic.skipline(4)
    else:    
        fic.skipline()

    nt,ne,nx3,nx2,nx1,_ = shape
    nelement = nx1*nx2*nx3
    data = numpy.ndarray(shape)
    val_array = numpy.ndarray((nelement,ne))
    err_array = numpy.ndarray((nelement,ne))
    zero = numpy.zeros(ne)   
    nextcell = True   
    end_data = False
    ne_read = ne if ne == 1 else ne-1
    for it in range(nt):
        for i in range(nelement):
            if nextcell:
                index,val,err = get_cdgs_element(fic,ne_read)
                if index is None : 
                    end_data = True
                else:    
                    index -= 1

                nextcell = False
            if i == index:
                nextcell = not end_data
                val_array[i,:] = val[:]
                err_array[i,:] = err[:]
            else:
                val_array[i,:] = zero[:]
                err_array[i,:] = zero[:]

        val_array = numpy.transpose(val_array.reshape((nx1,nx2,nx3,ne)))
        err_array = numpy.transpose(err_array.reshape((nx1,nx2,nx3,ne)))
        data[it,:,:,:,:,0] = val_array[:,:,:,:]
        data[it,:,:,:,:,1] = err_array[:,:,:,:]        
    return data