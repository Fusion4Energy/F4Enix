import numpy
from ..Modules.functions.utils import ExtraBin

def format_XYZ_Dim(etbin:ExtraBin, nval:int = 6) -> str:
    
    def format_XYZ_Dim_inter(vec,nval=6):
        interfound = False
        xtab = []
        dx0  = 0
        x1   = vec[0]
        i    = 0
        for x2 in vec[1:]:
            dx = round(x2-x1,2)
            if dx == dx0 :
                i += 1
                interfound = True
            else:
                xtab.append([i,x1])
                dx0 = dx
                i = 0
            x1 = x2
        xtab.append([i,x2])

        if not interfound :
            return format_XYZ_Dim_long(vec)
        line=''
        newline = f'  {xtab[0][1]:10.3e}'

        for i,v in enumerate(xtab[1:]):
            newline += f' {v[0]:4d}I {v[1]:10.3e}'
            if (i+1)%nval == 0:
                line += newline+'\n'
                newline = '                 '

        if (i+1)%nval != 0:
            line += newline
        else:
            line = line[:-1]
        return line

    def format_XYZ_Dim_long(vec,nval=8):
        line=''
        newline = ''
        for i,v in enumerate(vec):
            newline += f' {v:10.3e}'
            if (i+1)%nval == 0:
                line += newline+'\n'
                newline = '                 '

        if (i+1)%nval != 0:
            line += newline
        else:
            line = line[:-1]
        return line

    if len(etbin) <= nval :
        return format_XYZ_Dim_long(etbin)
    else:
        return format_XYZ_Dim_inter(etbin) 

def format_ET_bin(etbin:ExtraBin) -> str:
    line  = f'         flag : {etbin.type}\n'
    line += f'    bin index : bin range  \n'

    if etbin.type == 'erg':
        if etbin.binbound:
            for i in range(etbin.nvalue) :
                line += f'       {i:4d}   :      {etbin[i]} - {etbin[i+1]}  MeV\n'
        else:
            for i in range(etbin.nvalue) :
                line += f'       {i:4d}   :      {etbin[i]}  MeV\n'  
        if etbin.totalbin:
                 line += f'       {etbin.nvalue:4d}   :      Total\n'    

    elif etbin.type == 'tme':
        if etbin.binbound:
            for i in range(etbin.nvalue) :
                line += f'       {i:4d}   :      {etbin[i]} - {etbin[i+1]}  shakes\n'
        else:
            for i in range(etbin.nvalue) :
                line += f'       {i:4d}   :      {etbin[i]}  shakes\n'  
        if etbin.totalbin:
                 line += f'       {etbin.nvalue:4d}   :      Total\n'            
       
    elif etbin.type == 'dtme':
        for i in range(etbin.nvalue) :
            line += f'       {i:4d}   :      {etbin[i]}  shakes\n'  

    else:
        binlist=[]
        ibin = 0
        for b in etbin:
            decpart,intpart = numpy.modf(b)
            binnum = int(numpy.rint(abs(decpart)*1000))-1
            if binnum == ibin :
                binlist.append([int(intpart)])
                ibin += 1
            else:
                binlist[binnum].append(int(intpart))

        ib = 0
        line += f'       {ib:4d}   :  Other \n'
        for b in binlist[1:] :
            ib += 1
            tmp = f'       {ib:4d}   :'
            rang = 0
            for i,v in enumerate(b):
                if i < len(b)-1 :
                    if b[i+1] < 0:
                        rang = 1
                    if rang == 0:
                        tmp += f' {v},'
                elif rang == 1:
                    tmp += f' {v}-'
                    rang = 2
                else :
                    tmp += f'{abs(v)},'
                    rang = 0
            line += f'{tmp[:-1]}\n' 

        ib += 1
        line += f'       {ib:4d}   :  Total \n'

    return line
