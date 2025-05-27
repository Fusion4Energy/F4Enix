import numpy

def read_block(fic,nl,array):
    # read value and error blocks
    fic.skipline(2)
    for j in range(2):
      fic.skipline(2)
      for i in range(nl) :
         line = fic.readline()
         array[i,:,j] = [ float(x) for x in line.split()[1:] ]
      fic.skipline()  

def get_block_data(fic,position,type,explicit_time,shape, dtme=False):
    fic.seek(position)   
    if type == 'ij':
        nc = shape[4]  # X
        nl = shape[3]  # Y
        nb = shape[2]  # Z
    elif type == 'ik':
        nc = shape[4]  # X
        nl = shape[2]  # Z
        nb = shape[3]  # Y
    elif type == 'jk':
        nc = shape[3]  # Y
        nl = shape[2]  # Z
        nb = shape[4]  # X

    ne = shape[1]
    nt = shape[0]
    blockdata  = numpy.ndarray((nl,nc,2))
    wrkshape = (nt,ne,nb,nl,nc,2)
    data  = numpy.ndarray(wrkshape)

# for ttag != dtme outer loop over energies
    if not explicit_time:
        for ie in range(ne):
            fic.skipline(2)
            for ib in range(nb):
                read_block(fic,nl,blockdata)
                data[0,ie,ib,:,:,:] = blockdata[:,:,:]
                fic.skipline() 
            fic.skipline()    
               
    else :
        for ie in range(ne):
            fic.skipline(3)
            for it in range(nt):
                fic.skipline(2)
                for ib in range(nb):
                    read_block(fic,nl,blockdata)
                    data[it,ie,ib,:,:,:] = blockdata[:,:,:]
                    fic.skipline()
                fic.skipline()     
   

    if type == 'jk':
        data = numpy.moveaxis(data,2,4)
    elif type == 'ik':               
        data = numpy.swapaxes(data,2,3)

    return data

def get_column_data(fic,position,type,shape):
    fic.seek(position) 
       
    nt, ne, nx3, nx2, nx1, _ = shape
    if type == 'col':
       ival = -2
       ierr = -1
    else:
       ival = -4
       ierr = -3
  
    fic.skipline(1)
    col_shape = (nt,ne,nx1,nx2,nx3,2)
    data = numpy.ndarray(col_shape)
    for ie in range(ne):
       for it in range(nt):
          for i1 in range(nx1):
             for i2 in range(nx2):
                for i3 in range(nx3):
                     strvalues = fic.readline().split() 
                     val = float(strvalues[ival])
                     err = float(strvalues[ierr])
                     data[it,ie,i1,i2,i3,:] = (val,err)
    data = numpy.swapaxes(data,2,4)                 
    return data  

def read_cuv_cell_info(fic,ncell):
    cell_info = []
    tot_value = []
    tot_error = []
    for ic in range(ncell):
        line = fic.readline().split()
        cell = int(float(line[0]))
        volf = float(line[1])
        cell_info.append((cell,volf))
        if len(line) == 4:
            tot_value.append(float(line[2]))
            tot_error.append(float(line[3])) 
   
    addtot = len(tot_value) > 0
    
    val_tab = []
    err_tab = []
    for ic in range(ncell):
        line = fic.readline().split()
        values = numpy.array(line,dtype=numpy.double)
        line = fic.readline().split()
        errors = numpy.array(line,dtype=numpy.double)

        if addtot: 
            numpy.append(values,tot_value[ic])
            numpy.append(errors,tot_error[ic])
        val_tab.append(values)
        err_tab.append(errors)
         
    return cell_info, val_tab, err_tab    

def read_cdgs_cell_info(fic,ncell,ne):
    cell_info = []
    val_tab = []
    err_tab = []

    for ic in range(ncell):
        headline = fic.readline().split()
        cell = int(float(headline[0]))
        volf = float(headline[1])          
        cell_info.append((cell,volf))
        
        values = fic.get_multilines(ne)
        errors = fic.get_multilines(ne)
        if ne != 1 : 
            nh = len(headline)
            if nh == 2 :
                totval = numpy.sum(values)
                errval = 0.
            elif nh==3:   
                totval = float(headline[2])
                errval = 0.
            else:
                totval = float(headline[2])
                errval = float(headline[3])
     
            values = numpy.append(values, totval)
            errors = numpy.append(errors, errval)
         
        val_tab.append(values)
        err_tab.append(errors)

    return cell_info, val_tab, err_tab    

def select_cells(cell_info, values, errors, filter):
    new_cell_info = []
    new_values_tab = []
    new_errors_tab = []
    sumf = 0.

    for cinfo,val,err in zip(cell_info, values, errors):
        if cinfo[0] in filter:
            new_cell_info.append(cinfo)
            new_values_tab.append(val)
            new_errors_tab.append(err)
            sumf += cinfo[1]
    return sumf, new_cell_info, new_values_tab, new_errors_tab         

def get_mean_value(cell_info,val,err):
    val_data = numpy.array(val,dtype=numpy.double)
    err_data = numpy.array(err,dtype=numpy.double)
    err_data = val_data * err_data
    tmp = numpy.array(cell_info)
    frac = tmp[:,1]
    val_mean = numpy.matmul(frac,val_data)
    err_mean = numpy.where(val_mean != 0., numpy.matmul(frac,err_data) / val_mean, 0.)
    return val_mean, err_mean

def get_cuv_element(fic, norm, filter):
    line = fic.readline().split()
    volume = float(line[4])
    ncell = int(float(line[5]))
    cell_info,values, errors = read_cuv_cell_info(fic,ncell)

    ne = len(values[0]) 

    if filter is not None:
        sumf, cell_info, values, errors = select_cells(cell_info, values, errors, filter)
    else:
        sumf = 1.

    if sumf == 0. :
        values = numpy.zero(ne)
        errors = numpy.zero(ne)
        return values, errors
    
    if ncell > 1:
        values, errors = get_mean_value(cell_info, values,errors)
    else:    
        values = values[0]
        errors = errors[0]
    if norm == 'total':
        values = values*volume
    elif norm == 'fraction':    
        values = values/sumf

    return  values, errors  

def get_cdgs_element(fic, ne, norm = None, filter= None):
    line = fic.readline().split()
    if "end_source_data" == line[0] or "" == line[0]:
        return None,None,None
    
    index = int(float(line[0]))
    volume = float(line[-2])
    ncell = int(float(line[-1]))
    
    cell_info,values, errors = read_cdgs_cell_info(fic,ncell,ne)

    if filter is not None:
        sumf, cell_info, values, errors = select_cells(cell_info, values, errors, filter)
    else:
        sumf = 0.
        for cf in cell_info:
            sumf += cf[1]

    if sumf == 0. :
        values = numpy.zero(ne)
        errors = numpy.zero(ne)
        return values, errors

    if ncell > 1:
        values, errors = get_mean_value(cell_info, values, errors)
    else:
        values = values[0]
        errors = errors[0]

    if norm == 'total':
        values = values*volume
    elif norm == 'fraction':    
        values = values/sumf

    return  index,values, errors  
        