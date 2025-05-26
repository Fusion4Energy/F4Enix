from .binformat_functions import format_XYZ_Dim, format_ET_bin

class format_info:

    @staticmethod
    def print_mesh_info(info:dict) -> str:
       
        geom = 'rectangular' if info['geom'] == 'rec' else 'cylindrical'
        meshinfo = f" Tally          : {info['tally']}\n"

        if info['comments'] is not None:
            meshinfo += f" Comments       : \n {info['comments']}\n"

        meshinfo +=  f" Particle       : {info['particle']}\n"    
        meshinfo +=  f" Mesh geometry  : {geom}\n"
        
        print(meshinfo)
  
        xyzinfo = ''
        if info['geom'] == 'rec':
            xyzinfo += f" X dimensions   :{format_XYZ_Dim(info['x1bin'])}\n"
            xyzinfo += f" Y dimensions   :{format_XYZ_Dim(info['x2bin'])}\n"
            xyzinfo += f" Z dimensions   :{format_XYZ_Dim(info['x3bin'])}\n"
        else:
            xyzinfo += f" R dimensions   :{format_XYZ_Dim(info['x1bin'])}\n"
            xyzinfo += f" Z dimensions   :{format_XYZ_Dim(info['x2bin'])}\n"
            xyzinfo += f" T dimensions   :{format_XYZ_Dim(info['x3bin'])}\n"
     
        print(xyzinfo)

        if info['ebin'].explicit :
            line = f" Energy bins    :\n{format_ET_bin(info['ebin'])}"
            print(line)

        if info['tbin'].explicit :
            line = f" Energy bins    :\n{format_ET_bin(info['tbin'])}"
            print(line)

    @staticmethod
    def display_meshinfo(meshparser):
        line =  f" Meshtally file : {meshparser.filename}\n"     
        tallylist = meshparser.get_meshlist()
        for t in tallylist:
            particle,comments = meshparser.get_header(t)   
            if comments is not None:
                firstlinecom = comments.split('\n')[0].strip()
                line += f"   Tally {t} : {particle}  '{firstlinecom}' "
            else:
                line += f"   Tally {t} : {particle}  "
        print(line)

    @staticmethod 
    def display_tallyinfo(meshparser,tally):
        particle,comments = meshparser.get_header(tally)
        geom, trsf, meshbins = meshparser.get_boundaries(tally)
        info ={
            'geom':geom,
            'tally':tally,
            'comments':comments,
            'particle':particle,
            'x1bin':meshbins[0],
            'x2bin':meshbins[1],
            'x3bin':meshbins[2],
            'ebin':meshbins[3],
            'tbin':meshbins[4],
        }
        format_info.print_mesh_info(info)