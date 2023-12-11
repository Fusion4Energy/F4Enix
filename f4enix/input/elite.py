from f4enix.input.MCNPinput import Input
import os
import pandas as pd
import math
import copy
from numjuggler import parser
import re

class ELite_Input(Input):
    
    def initialize_elite(self, inpfile: os.PathLike, 
                         check_Elite: bool = True):
        # checks if the input is actually e-lite and initializes some variables
        # check is optional
        self.block_structure = pd.read_excel(inpfile)
        self.elite_cells_list = list(self.cells.keys())

        # check the envelopes
        # if set(self.elite_cells_list) != set(self.block_structure['Cell'
        #                                                             ].tolist()):
        #     msg = 'MCNP input is not an E-Lite file, or the Excel Block' + \
        #             ' Structure file is not compatible'
        #     print(msg)
        #     exit(0)
        
        self.sectors_env_cells_names = {}
        self.sectors_universes = {}
        self.sectors_names = [1, '2 & 3', 4, 5, 6, 7, 8, 9]
        self.sectors_cells = {}
        self.sectors_cells_names = {}
        self.sectors_L0_cells = {}
        self.sectors_L0_cells_names = {}
        self.tally_cards_types = ['F', '+F', '*F', 'FC', 'E', 'T', 'C', 'FQ',
                                  'FM', 'DE', 'DF', 'EM', 'TM', 'CM', 'CF',
                                  'SF', 'FS', 'SD', 'FU', 'FT', 'TF', 'NOTRN']
        
        for sec in self.sectors_names:
            # you will have to modify outer cell to put them in the sectors
            self.sectors_env_cells_names[sec] = self.block_structure[self.block_structure['Sector'] == sec]['Cell'].tolist()
            self.sectors_universes[sec] = self.block_structure[self.block_structure['Sector'] == sec]['Universe ID'].tolist()
            self.sectors_cells[sec] = []
            self.sectors_cells_names[sec] = []
            self.sectors_L0_cells[sec] = []
            self.sectors_L0_cells_names[sec] = []     

        for cell in self.cells.values():
            try:
                uni = cell._get_value_by_type('u')
                for sec in self.sectors_names:
                    if uni in self.sectors_universes[sec]:
                        self.sectors_cells[sec].append(cell)
                        self.sectors_cells_names[sec].append(cell.name)
            except:
                for sec in self.sectors_names:
                    if cell.name in self.sectors_env_cells_names[sec]:
                        self.sectors_cells[sec].append(cell)
                        self.sectors_cells_names[sec].append(cell.name)
                        self.sectors_L0_cells[sec].append(cell)
                        self.sectors_L0_cells_names[sec].append(cell.name)

        self.sector_boundaries = {1:       [427016, 427024],
                                '2 & 3': [437543, -437544],
                                4:       [451016, 451024],
                                5:       [459016, 459024],
                                6:       [467016, 467024],
                                7:       [475016, 475024],
                                8:       [483016, 483024],
                                9:       [491016, 491024],}

        self.sector_boundaries_angles = {1:      [50, 10],
                                        '2 & 3': [130, 50],
                                        4:      [170, 130], 
                                        5:      [-150, 170],
                                        6:      [-110, -150],
                                        7:      [-70, -110],
                                        8:      [-30, -70], 
                                        9:      [10, -30],}
        # this can be refined
        self.plasma_cells = {1:       427001,      
                            '2 & 3': 435001,
                            4:       451001,
                            5:       459001,
                            6:       467001,
                            7:       475001,
                            8:       483001,
                            9:       491001,}
        return

    def extract_sector(self, sectors, outfile: os.PathLike = 'sector'):
        
        if not isinstance(sectors, list):
            sectors = [sectors]

        boundaries_angles = self._get_boundaries_angles(sectors)

        # logging.info('Collecting the cells, surfaces, materials and transf.')
        # make sure these are str
        cells = []
        for sector in sectors:
            cells += self.sectors_L0_cells_names[sector]

        cset = set(cells)

        # first, get all surfaces needed to represent the cn cell.
        sset = set()  # surfaces
        mset = set()  # material
        # tset = set()  # transformations
        
        # duplicate the final set and work on a dynamic set that contains only 
        # new cells at each loop
        cell_set = copy.deepcopy(cset)

        # next runs: find all other cells:
        again = True
        while again:
            again = False
            new_set = set()
            uni_set = set()
            # loop over cells to extract
            for cell_num in cell_set:
                c = self.cells[str(cell_num)]
                # get hash cells in the cells that have to be extracted
                cref = c.get_refcells()
                # add the hash cells to extraction list
                new_set |= cref
                # collect universes in the definition of cells
                fill = c.get_f()
                if fill is not None:
                    uni_set.add(fill)
            # if one wants to extract also lower levels, loop over universes 
            # and collect their cells
            for key, c in self.cells.items():
                if c.get_u() in uni_set:
                    new_set.add(c.values[0][0])
            # get the new set with the cells to be checked
            cell_set = new_set - cell_set
            # check if loop is to be repeated
            if cell_set:
                again = True
                cset |= cell_set
        
        # sort the set
        cset = list(cset)
        cset.sort()


        # Get all surfaces and materials
        cells_cards = self.get_cells_by_id(cset)
        for i, (_, cell) in enumerate(cells_cards.items()):
            for v, t in cell.values:
                if t == 'sur':
                    sset.add(v)
                elif t == 'mat':
                    if int(v) != 0:  # void material is not defined in a card
                        mset.add('M'+str(v))

        # order surfaces
        sset = list(sset)
        sset.sort()

        self.modified_data_cards = copy.deepcopy(self.other_data)
        pattern_str = '|'.join(self.tally_cards_types)
        pattern = re.compile(f'^({pattern_str})\d+$')
        # for now remove all tally cards
        for key in self.other_data.keys():  
            
            # Explanation of the pattern:
            # ^           - Start of the string
            # (pattern)   - A group containing possible patterns
            # \d+         - One or more digits
            # $           - End of the string
            if pattern.match(key):
                self.modified_data_cards.pop(key)

        self.modified_surfaces = copy.deepcopy(self.surfs)
        # Future implementation
        self._set_sdef(sectors)

        self.modified_boundaries = self._set_boundaries(self._get_planes_angles(
                                                        boundaries_angles))
        self._set_periodic_boundaries()
        # extract tallies based on comments
        # self._extract_tallies(sectors)
        # logging.info('write MCNP reduced input')
        self._modify_graveyard()
        # logging.info('write MCNP reduced input')
        with open(outfile, 'w') as outfile:
            # Add the header lines
            for line in self.header:
                outfile.write(line)
            # Add the cells
            outfile.writelines(self._print_cards(cells_cards))
            outfile.writelines(self._print_cards(self.new_gy_outer))
            # Add a break
            outfile.write('\n')
            # Add the surfaces
            surfs = self.get_surfs_by_id(sset)
            outfile.writelines(self._print_cards(surfs))
            # Add a break
            outfile.write('\n')
            # Add materials
            materials = self.get_materials_subset(mset)
            outfile.write(materials.to_text()+'\n')
            outfile.writelines(self._print_cards(self.transformations))
            # Add the rest of the datacards
            outfile.writelines(self._print_cards(self.modified_data_cards, wrap=True))
            # Add a break
            outfile.write('\n')
        # logging.info('input written correctly')
        return

    def _set_sdef(self, sectors):
        new_si = 'SI70 L '
        new_sp = 'SP70 '
        for k, sector in enumerate(sectors):
            new_si = new_si + str(self.plasma_cells[sector]) + ' '
            if sector != '2 & 3':
                new_sp += '1 '
            else:
                new_sp += '2 '
        self.modified_data_cards['SI70'].input = new_si[:-1]
        self.modified_data_cards['SP70'].input = new_sp[:-1]
        return
    
    def _set_boundaries(self, angle_dic):
        # define the angles at which there are the planes cutting the sectors\
        bounds = list(angle_dic.keys())
        self.boundary_angles = angle_dic[bounds[0]] + angle_dic[bounds[1]]
        # derive the coefficients in plane equation of such planes
        self.coeff_x = []
        self.coeff_y = []
        ref_x = 0
        ref_y = 1
        for l, angle in enumerate(self.boundary_angles):
            self.coeff_x.append(ref_x * math.cos(math.radians(angle)) - ref_y * math.sin(math.radians(angle)))
            self.coeff_y.append(ref_x * math.sin(math.radians(angle)) + ref_y * math.cos(math.radians(angle)))
        # set a tolerance to coefficients' values to check if the planes are 
        self.tol = 1e-5
        # these dicts will group the planes parallel to each reference cutting plane
        self.boundary_surfs = {num: [] for num in self.boundary_angles}
        self.boundary_surfs_names = {num: [] for num in self.boundary_angles}
        
        # For each angle, fulfill the list of parallel planes
        for l, angle in enumerate(self.boundary_angles):
            # Check all surfaces for each angle
            for surf in self.surfs.values():
                # check if the surface is a plane
                if surf.stype == 'p':
                    # check if the surface is a plane parallel to z
                    if surf.scoefs[2] == 0 and surf.scoefs[3] == 0:
                        # check if the plane has a transformation
                        try:
                            # if so, apply it and check if it's a plane parallel to
                            # the boundary
                            tr = surf._get_value_by_type('tr')
                            tr_name = 'TR' + str(tr)
                            trans = self.transformations[tr_name]
                            # check if it's only a translation
                            if len(trans.values) < 12:
                                continue
                            if trans.unit == '*':
                                p_x_coeff = (surf.scoefs[0]*math.cos(math.radians(trans.values[4][0]))-surf.scoefs[1]*math.cos(math.radians(trans.values[5][0])))/math.sqrt(surf.scoefs[0]**2 + surf.scoefs[1]**2)
                                p_y_coeff = (-surf.scoefs[0]*math.cos(math.radians(trans.values[7][0]))+surf.scoefs[1]*math.cos(math.radians(trans.values[8][0])))/math.sqrt(surf.scoefs[0]**2 + surf.scoefs[1]**2)
                            else:
                                p_x_coeff = (surf.scoefs[0]*trans.values[4][0]-surf.scoefs[1]*trans.values[5][0])/math.sqrt(surf.scoefs[0]**2 + surf.scoefs[1]**2)
                                p_y_coeff = (-surf.scoefs[0]*trans.values[7][0]+surf.scoefs[1]*trans.values[8][0])/math.sqrt(surf.scoefs[0]**2 + surf.scoefs[1]**2)
                            if self.coeff_x[l]-self.tol <= (p_x_coeff)/math.sqrt(surf.scoefs[0]**2 + surf.scoefs[1]**2) <= self.coeff_x[l]+self.tol and self.coeff_y[l]-self.tol <= (p_y_coeff)/math.sqrt(surf.scoefs[0]**2 + surf.scoefs[1]**2)<= self.coeff_y[l] + self.tol:
                                    self.boundary_surfs[angle].append(surf)
                                    self.boundary_surfs_names[angle].append(surf.name)
                        except UnboundLocalError:
                            # if not, just check if it's parallel to the boundary
                            if self.coeff_x[l]-self.tol <= (surf.scoefs[0])/math.sqrt(surf.scoefs[0]**2 + surf.scoefs[1]**2) <= self.coeff_x[l]+self.tol and self.coeff_y[l]-self.tol <= (surf.scoefs[1])/math.sqrt(surf.scoefs[0]**2 + surf.scoefs[1]**2) <= self.coeff_y[l]+self.tol:
                                self.boundary_surfs[angle].append(surf)
                                self.boundary_surfs_names[angle].append(surf.name)
        surfs_lis = []                     
        for ang_bounds in angle_dic.values():
            for ang_bound in ang_bounds:
                surfs_lis += self.boundary_surfs_names[ang_bound]
        surfs_lis = list(set(surfs_lis))
        return surfs_lis
    
    def _get_planes_angles(self, boundaries_angles):
        mirroring = [360, 180, -180, -360]
        angle_dic = {}
        for elem in boundaries_angles:
            if elem > 0:
                mirror = mirroring[1:]
            else:
                mirror = mirroring[:-1]
            ang_lis = [elem]
            for ang in mirror:
                ang_lis.append(elem + ang)
            angle_dic[elem] = ang_lis
        return angle_dic
    
    def _get_boundaries_angles(self, sectors):

        boundaries_angles = []

        for k, sec in enumerate(sectors):
            bounds = self.sector_boundaries_angles[sec]
            if k == 0:
                boundaries_angles.append(bounds[1])
            if k == (len(sectors)-1):
                boundaries_angles.append(bounds[0])

        return boundaries_angles

    def _set_periodic_boundaries(self):
        
        for surf_n in self.modified_boundaries:
            surf = self.modified_surfaces[str(surf_n)]
            surf.input[0] = '*' + surf.input[0]
        return
    
    # def _extract_tallies(self, sectors):
        
    #     sector = 
    #     for key, card in self.other_data:
    #         tally_commment = self.other_data['FC'+ str(tally_num)]
    #         pattern = rf'^({"|".join(self.tally_cards_types)})\d+$'
    #     return

    def _modify_graveyard(self, sectors):

        for k, sector in enumerate(sectors):
            if k == 0:
                self.new_gy = union_cell(self.cells['801'], -self.sector_boundaries[sector][0], None)
                self.new_outercell = cut_cell(self.cells['800'], self.sector_boundaries[sector][0], None)
            elif k == len(sectors) - 1:
                self.new_gy = union_cell(self.new_gy, -self.sector_boundaries[sector][1], None)
                self.new_outercell = cut_cell(self.new_outercell, self.sector_boundaries[sector][1], None)
            else:
                continue
        self.new_gy_outer = {'800': self.new_outercell,
                             '801': self.new_gy}
        return
        
            


# initialize cut_cell function
def union_cell(cell:parser.Card, union_surface:int, new_cell_num:int = None):

    if new_cell_num is None:
        new_cell_num = cell.name

    new_cell = copy.deepcopy(cell)
    # Introduce parentheses before the third word in the first row
    first_row = new_cell.input[0].split()

    if new_cell._get_value_by_type('mat') == 0:
        first_row[2] = '(' + first_row[2]
    else:
        first_row[3] = '(' + first_row[3]

    new_cell.input[0] = ' '.join(first_row)

    # Check all rows if there are letters in the row
    for i in range(len(new_cell.input)):

        row = new_cell.input[i].split()
        param_cards_idx = -1

        for m, words in enumerate(row):
            if any(c.isalpha() for c in words):
                param_cards_idx = m
                break

        if param_cards_idx != -1:
            if union_surface >= 0:
                row.insert(param_cards_idx, 
                           ') :{:<' + str(len(str(union_surface))) + '} ' )
            else:
                row.insert(param_cards_idx, 
                           ') :-{:<' + str(len(str(union_surface))) + '} ' )
                
            new_cell.input[i] = ' '.join(row)

            if new_cell.input[i][:5] != '     ' and i != 0:
                new_cell.input[i] = '     ' + new_cell.input[i]
            break

    for k in range(len(new_cell.values)-1, -1, -1):
        if new_cell.values[k][1] == 'sur' or new_cell.values[k][1] == 'cel':
            last_sur_idx = k
            break
        
    new_cell.values.insert(k+1, (abs(union_surface),'sur'))

    new_cell.name = new_cell_num
    new_cell._set_value_by_type('cel', new_cell_num)
    
    return new_cell

# initialize cut_cell function
def cut_cell(cell:parser.Card, split_surface:int, new_cell_num:int = None):

    if new_cell_num is None:
        new_cell_num = cell.name

    new_cell = copy.deepcopy(cell)
    # Introduce parentheses before the third word in the first row
    first_row = new_cell.input[0].split()

    if new_cell._get_value_by_type('mat') == 0:
        first_row[2] = '(' + first_row[2]
    else:
        first_row[3] = '(' + first_row[3]

    new_cell.input[0] = ' '.join(first_row)

    # Check all rows if there are letters in the row
    for i in range(len(new_cell.input)):

        row = new_cell.input[i].split()
        param_cards_idx = -1

        for m, words in enumerate(row):
            if any(c.isalpha() for c in words):
                param_cards_idx = m
                break

        if param_cards_idx != -1:
            if split_surface >= 0:
                row.insert(param_cards_idx, 
                           ') {:<' + str(len(str(split_surface))) + '} ' )
            else:
                row.insert(param_cards_idx, 
                           ') -{:<' + str(len(str(split_surface))) + '} ' )
                
            new_cell.input[i] = ' '.join(row)

            if new_cell.input[i][:5] != '     ' and i != 0:
                new_cell.input[i] = '     ' + new_cell.input[i]
            break

    for k in range(len(new_cell.values)-1, -1, -1):
        if new_cell.values[k][1] == 'sur' or new_cell.values[k][1] == 'cel':
            last_sur_idx = k
            break
        
    new_cell.values.insert(k+1, (abs(split_surface),'sur'))

    new_cell.name = new_cell_num
    new_cell._set_value_by_type('cel', new_cell_num)
    
    return new_cell