import os

from .menu import Mesh2vtkMenu as menu
from ..Modules.mesh_parser  import CDGSMeshParser, MeshtalParser, CUVMeshParser
from ..Modules.FMesh import FMesh
from .showinfo import format_info
   
class Interface:   
    def __init__(self):
        self.meshfiles = dict()   # contains meshfile parser objects
        self.FMesh = dict()       # contains mesh extract from list
        self.meshList = dict()    # contains list of available mesh in opened meshfiles
        
    def run(self):
        menu.clear_screen()
        menu.display_menu('principal')
        ans,optname=menu.answer_loop('principal')
       
        while True:
            if ans[0:4] == 'open':
               self.load_meshfile(optname,ans[4:])
                       
  # Print meshtally information
            elif ans == 'info' :
                if len(self.meshfiles) == 0:
                    print(' No meshtally file')
                else:
                    self.display_info()

  # Write mesh to VTK
            elif ans == 'write':
                if len(self.meshfiles) == 0:
                    print(' No meshtally file')
                else:
                    self.write_vtk()

  # operate on mesh
            elif ans == 'operate':
                if len(self.meshfiles) == 0:
                    print(' No meshtally file')
                else:
                    self.operate()
            else:
                break
            ans,optname=menu.answer_loop('principal')
            menu.clear_screen()
            menu.display_menu('principal')

    def update_meshlist(self):
        for key in self.meshfiles.keys():
            if key not in self.meshList.keys():
                tally_list = self.meshfiles[key].get_meshlist()
                self.meshList[key] = tally_list
                self.FMesh[key] = dict()

    def append_meshresults(self,meshname):
        if '_results_' in self.meshList.keys(): 
            self.meshList['_results_'].append(meshname)
            self.meshList['_results_'].sort() 
        else:
             self.meshList['_results_'] = [meshname] 
             self.FMesh['_results_'] = dict()              

    def check_filename(self,name):
        if name is None:
            fname = input(' enter Meshtally file name:')
        else:
            fname = name
        if fname != '' and fname not in self.meshfiles.keys() :
            if os.path.isfile(fname):
                self.meshfilenames.append(fname)
            else:
                print (f' {fname} not found')
            print('\n Input files :')
            for f in self.meshfilenames:
                print(f" - {f}")

    def load_meshfile(self,filename,type):
        if filename is None:
            fname = input(' enter Meshtally file name:')
        else:
            fname = filename

        if fname != '' and fname not in self.meshfiles.keys() :
            if not os.path.isfile(fname):
                print (f' {fname} not found')
                return   
        else:
           return        
        
        if type == 'CUV' :
            self.meshfiles[filename] = CUVMeshParser(filename)
        elif type == 'CDGS':
            self.meshfiles[filename] = CDGSMeshParser(filename)
        elif type == '':
            self.meshfiles[filename] = MeshtalParser(filename)
        
        self.update_meshlist()         
        print('\n Input files :')
        for f in self.meshList.keys():
            print(f" - {f}") 

    def get_mesh(self,mshfile,tally):
        if mshfile not in self.FMesh.keys():
            print(f'meshfile {mshfile} in not loaded')
            return None
        
        if tally not in self.FMesh[mshfile].keys():
            msh = self.meshfiles[mshfile].get_FMesh(tally)
            self.FMesh[mshfile][tally] = msh
        else:    
            msh = self.FMesh[mshfile][tally]
        return msh     
       
    def operate(self):
        menu.display_menu('operate')       
        ans,meshname = menu.answer_loop('operate')
        onefile = menu.display_meshlist(self.meshList)
        if ans == 'scale' : 
            mlist = menu.answer_meshlist('4',ans,self.meshList,onefile)
            if len(mlist) > 0 :
                mshfile, tally, sfact =  mlist[0]
                msh = self.get_mesh(mshfile, tally) 
                smesh = msh.scale(sfact)
                self.append_meshresults(meshname)
                self.FMesh['_results_'][meshname] = FMesh(smesh,meshname,msh.trsf)
                
        
        elif ans == 'sum' : 
            mlist = menu.answer_meshlist('4',ans,self.meshList,onefile)
            if len(mlist) > 1 :
                mshfile,tally,sfact = mlist[0]
                smesh = self.get_mesh(mshfile, tally)
                trsf = smesh.trsf
                if sfact != 1.0: smesh = smesh.scale(sfact)
                
                for mshfile,tally,sfact in mlist[1:]:
                    cmesh = self.get_mesh(mshfile, tally)
                    if sfact != 1.0: cmesh = cmesh.scale(sfact)
                    smesh = smesh + cmesh   
                self.append_meshresults(meshname)     
                self.FMesh['_results_'][meshname] = FMesh(smesh,meshname,trsf) 
                  

    def write_vtk(self):
        onefile = menu.display_meshlist(self.meshList)
        mlist = menu.answer_meshlist('3',None,self.meshList,onefile)
        print()
        
        for mshfile,tally,_ in mlist:
           mshname = mshfile.strip('_')
           outname = f'{mshname}_{tally}'
           print(f'write vtk file : {outname}')
           msh = self.get_mesh(mshfile, tally)
           msh.write_vtk(outname)
   
    def display_info(self):
        menu.display_menu('info')       
        ans0,ans1 = menu.answer_loop('info')
    
        if len(self.meshfiles) > 1 :
            print(' Input files :')
            for f in self.meshfiles.keys():
                print(" - {}".format(f))
            while True:
                sfile = input(' Select file #:')
                if sfile in self.meshfiles.keys():
                    meshparser = self.meshfiles[sfile]
                    break
                else:
                    print (' bad filename')
        else:
            meshparser = next(self.meshfiles.values().__iter__())

        if   ans0 == 'meshfile_info' :
            format_info.display_meshinfo(meshparser)

        # print tally information
        elif ans0 == 'tally_info' :
            tallylist = meshparser.get_meshlist()
            tally = ans1
            if tally not in tallylist:
                print(' Tallies :')
                for f in tallylist:
                    print(" - tally {}".format(f))

            tally = input(' Select tally #:')        
            while tally not in tallylist :
                if tally.lower in ('end','stop','exit','cancel'):
                    return
                else:    
                    print (' bad tally number')
                    tally = input(' Select tally #:')

            format_info.display_tallyinfo(meshparser, tally)

