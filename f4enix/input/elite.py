from f4enix.input.MCNPinput import Input

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
        return

    def extract_sector(self, sectors, outfile: os.PathLike = 'sector',
                    only_envelopes: bool = False):
        
        if not isinstance(sectors, list):
            sectors = [sectors]

        boundaries_angles = self._get_boundaries_angles(sectors)

        # logging.info('Collecting the cells, surfaces, materials and transf.')
        # make sure these are str
        cells = []
        for sector in sectors:
            cells.append(self.sectors_cells_names[sector])

        cset = []
        for cell in cells:
            cset.append(str(cell))
        cset = set(cset)

        # first, get all surfaces needed to represent the cn cell.
        sset = set()  # surfaces
        mset = set()  # material
        # tset = set()  # transformations
        
        again = True
        cell_set = cset
        # next runs: find all other cells:
        while again:
            again = False
            new_set = set()
            for cell_num in cell_set:
                c = self.cells[cell_num]
                cref = c.get_refcells()
                new_set = new_set.union(cref)
            cell_set = new_set.difference(cell_set)
            if cell_set:
                again = True
                cset.union(cell_set)

        # Get all surfaces and materials
        cells_cards = self.get_cells_by_id(cset)
        for i, (_, cell) in enumerate(cells_cards.items()):
            for v, t in cell.values:
                if t == 'sur':
                    sset.add(v)
                elif t == 'mat':
                    if int(v) != 0:  # void material is not defined in a card
                        mset.add('M'+str(v))

        self.modified_data_cards = {**self.other_data}
        # Future implementation
        # self._set_sdef(sectors)

        modified_boundaries = self._set_boundaries(sector, _get_planes_angles(self, boundaries_angles))
    
        # logging.info('write MCNP reduced input')
        with open(outfile, 'w') as outfile:
            # Add the header lines
            for line in self.header:
                outfile.write(line)
            # Add the cells
            outfile.writelines(self._print_cards(cells_cards))
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

    def _set_sdef(self, sector):
        if sector == 1:
            sec = 1
        elif sector == '2 & 3':
            sec = 2
        else:
            sec = sector - 1
        cell_dist = self.modified_data_cards['SDEF'].split()
        for n, opt in enumerate(cell_dist):
            if opt == 'CEL':
                break
        distr_card = ''.join(filter(str.isnumeric, cell_dist[n+1]))
        
        inp_si = self.modified_data_cards['SI'+distr_card].input
        inp_sp = self.modified_data_cards['SP'+distr_card].input

        new_inp_si = []
        for l, word in enumerate(inp_si.split()):
            if l  in [0, 1, int(sec) + 1]:
                new_inp_si.append(word)
        new_inp_si = ' '.join(new_inp_si)

        new_inp_sp = []
        for l, word in enumerate(inp_sp.split()):
            if l  in [0, int(sec)]:
                new_inp_sp.append(word)
        new_inp_sp = ' '.join(new_inp_si)

        self.modified_data_cards['SI'+distr_card].input = new_inp_si
        self.modified_data_cards['SP'+distr_card].input = new_inp_sp

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
                                angle_rad = math.radians(trans.values[4][0])
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
                if k == 0:
                    bounds = sector_boundaries[sec]
                    boundaries_angles.append(sector_boundaries_angles[sec][1])
                if k == (len(sectors_input)-1):
                    bounds = sector_boundaries[sec]
                    boundaries_angles.append(sector_boundaries_angles[sec][0])

            return boundaries_angles