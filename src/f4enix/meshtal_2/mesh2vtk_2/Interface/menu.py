import os

principal_menu="""
 ***********************************************
   Process meshtally to VTK
 ***********************************************

 * Append mesh file         (open, openCUV, openCDGS)
 * Display mesh information (info)
 * Write VTK file           (write)
 * Mesh operation           (operate) 
 * Exit                     (end)
"""

info_menu="""
 * Mesh information         (meshfile_info)
 * Tally information        (tally_info)
 * Exit                     (end)
"""

operate_menu="""
 * Mesh scale               (scale)
 * Mesh sum                 (sum)
"""

menu_dict = {'principal' : principal_menu,
             'info' : info_menu,
             'operate' : operate_menu}

principal_keys = ['open','openCUV','openCDGS','info','write','operate','end']
info_keys = ['meshfile_info','tally_info','end']
operate_keys = ['scale','sum']

def check_tally(*args):
    if len(args) == 2:
        tally = args[0]
        meshList = args[1]
        msh = next(meshList.keys().__iter__())
    elif len(args) == 3:    
        msh = args[0]
        tally = args[1]
        meshList = args[2]
    
    if msh not in meshList.keys():
        print(f' {msh} meshfile not loaded.')
        return True
    
    if tally not in meshList[msh]:
        print(f' {tally} tally not in {msh} file.')
        return True
    return False

def get_float(val):
    try :
        val = float(val) 
    except:
        print (' scale factor should be float')
        return True    
    return val,False

class Mesh2vtkMenu:
    
    @staticmethod
    def display_menu(menu0):  
        if menu0 in menu_dict.keys():
            print(menu_dict[menu0])  
    
    @staticmethod
    def clear_screen():
        if os.name == 'nt':
            os.system('cls')
        else:
            os.system('clear')           
    
    @staticmethod  
    def answer_loop(menu):
        while True:
            ans = input(" enter action :")
            ans = ans.split()
            ans0 = None
            ans1 = None
            if len(ans) > 0:
                ans0 = ans[0]
                if len(ans) > 1 :
                    ans1 = ans[1]

            if menu == 'operate' and ans1 == None:
                if ans0 not in operate_keys:
                    print(' bad operation keyword')
                else:
                    print(' name of the result mesh is missing')
            elif ans0 in menu_dict[menu] : 
                break
            else:
                print(' not expected keyword')
        return ans0,ans1  
   
    @staticmethod
    def display_meshlist(meshList):
        # display mesh list
        if len(meshList) > 1 :
            onefile = False
            for msh in meshList.keys():
                if msh == '_results_': continue
                print(f" - {msh}")
                for tally in meshList[msh]:
                    print(f"     - {tally}")
        else:
            onefile = True
            msh = next(meshList.keys().__iter__())
            for tally in meshList[msh]:
                    print(f"     - {tally}")

        if "_results_" in meshList.keys():
            print(" - Result mesh (_results_)")
            for tally in meshList['_results_']: 
                print(f"     - {tally}")

        return onefile

    @staticmethod
    def answer_meshlist(menu,answer,meshList,onefile=True):
        label_00 = ' select mesh, finish "end":'
        label_01 = ' tally and scaling factor "tally factor":'
        label_02 = ' mesh and scaling factor "tally factor", finish "end":'
        label_03 = ' first mesh  "tally":'
        label_04 = ' second mesh "tally":'
        label_05 = ' select mesh :'
        label_10 = ' select mesh "file tally", finish "end":'
        label_11 = ' meshtally file, tally and scaling factor "file tally factor":'
        label_12 = ' meshtally file, tally and scaling factor "file tally factor", finish "end":'
        label_13 = ' first mesh  "file tally":'
        label_14 = ' second mesh "file tally":'
        label_15 = ' select mesh "file tally":'

        label_one = [label_00,label_01,label_02,label_03,label_04,label_05]
        label_mlt = [label_10,label_11,label_12,label_13,label_14,label_05]
  

        if menu == '3' :
            labelindex =  0 
            nans       =  2
            nvalues    = -1
        elif menu == '4' :
            if answer == 'scale' :
                labelindex =  1 
                nans       =  3
                nvalues    =  1
            elif answer == 'sum' :
                labelindex =  2 
                nans       =  3
                nvalues    = -1

        finish = False
        mlist = []
        f1 = 1. 

        if onefile :
            nans = nans - 1
        label = label_one if onefile else label_mlt
     
        while not finish:
            while True:
                error = False 
                sfile = input(label[labelindex])
                sfile = sfile.split()
                if len(sfile) == 1  and sfile[0] == 'end': 
                    finish = True
                    break
                if len(sfile) == nans :
                    if sfile[0] == 'all': 
                        for msh in meshList.keys():
                            for tally in meshList[msh].keys():
                                mlist.append([msh,tally,f1])
                        finish = True
                        break
                    
                    if onefile:
                        msh = next(meshList.keys().__iter__())
                        tally = sfile[0]
                        error1 = check_tally(tally,meshList)
                        error2 = False
                        if nans == 2 : 
                            f1,error2 = get_float(sfile[1])
                    else:  
                        msh = sfile[0]   
                        tally = sfile[1]
                        error1 = check_tally(msh,tally,meshList)
                        error2 = False
                        if nans == 3 : 
                            f1,error2 = get_float(sfile[2])
                    error = error1 or error2    
                else:
                    error = True
                    print(f' {len(sfile)} values are entered,{nans} are expected')

                if not error :
                    mlist.append([msh,tally,f1])
                    if nvalues == len(mlist) : finish = True  
                    break
        return mlist