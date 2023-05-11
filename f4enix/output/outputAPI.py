"""
Parsing of MCNP output file (.o)

Only a few features are implemented for the moment being.
"""
import os
import re
import logging
import pandas as pd
import pyvista as pv

# -- Identifiers --
SURFACE_ID = 'currently being tracked has reached surface'
CELL_ID = 'other side of the surface from cell'
POINT_ID = 'x,y,z coordinates:'
# -- Patterns --
NUM_PAT = re.compile(r'\d+')
SCIENTIFIC_PAT = re.compile(r'-*\d.\d+E[+|-]\d+')
# patComments = re.compile(r'(?i)C\s+')
# patUniverse = re.compile(r'(?i)u=\d+')
# patNPS = re.compile(r'(?i)nps')
# patNPS_value = re.compile(r'\s+[eE.0-9]+\s*')
# patXYZ = re.compile(r'\s+[eE.0-9]+\s[eE.0-9]+\s[eE.0-9]+\s*')
# patSDEF = re.compile(r'(?i)sdef')
# patSDEFsur = re.compile(r'(?i)sur=\d+')


class Output:
    def __init__(self, filepath: os.PathLike) -> None:
        """Object representing and MCNP output file

        Parameters
        ----------
        filepath : os.PathLike
            path to the output file to be parsed

        Attributes
        ----------
        filepath: os.PathLike
            path to the original output file
        name: str
            name of the original file
        lines: list[str]
            list of all the lines of the file

        Examples
        --------
        >>> from f4enix.output.outputAPI import Output
        ... # parse the output
        ... outp = Output('test.o')
        ... # print excel and .vtk file containing info on lost particles
        ... outp.print_lp_debug('outfile_name')
        ... # Get the results of the statistical checks
        ... outp.get_statistical_checks()
        {11: 'All zeros',
         12: 'All zeros',
         21: 'Passed',
         22: 'Missed',
         ...
         66: 'Passed',
         76: 'Missed',
         86: 'Passed',
         96: 'Passed',
         106: 'Missed'}

        """
        self.filepath = filepath
        self.name = os.path.basename(filepath)
        logging.info('reading {} ...'.format(self.name))
        self.lines = self._read_lines()
        logging.info('reading completed')

    def _read_lines(self) -> list[str]:
        lines = []
        with open(self.filepath, 'r', errors="surrogateescape") as infile:
            for line in infile:
                lines.append(line)
        return lines

    def print_lp_debug(self, outpath: os.PathLike,
                       input_model: os.PathLike = None) -> None:
        """prints both an excel ['LPdebug_{}.vtp'] and a vtk cloud point file
        ['LPdebug_{}.vtp'] containing information about the lost particles
        registered in an MCNP run.

        Parameters
        ----------
        outpath : os.PathLike
            path to the folder where outputs will be dumped.
        """

        # -- Variables --
        surfaces = []
        cells = []
        x = []
        y = []
        z = []

        logging.info(
            'Recovering lost particles surfaces and cells in '+self.name)

        for i, line in enumerate(self.lines):

            if line.find(SURFACE_ID) != -1:  # LP in surface
                surfaces.append(NUM_PAT.search(line).group())

            if line.find(CELL_ID) != -1:  # LP in cell
                cells.append(NUM_PAT.search(line).group())

            if line.find(POINT_ID) != -1:  # LP in cell
                point = SCIENTIFIC_PAT.findall(line)  # [0:3]
                x.append(float(point[0]))
                y.append(float(point[1]))
                z.append(float(point[2]))
            #     if '***' in self.lines[i-10]:
            #         gp = SCIENTIFIC_PAT.findall(self.lines[i-9])[0:3]
            #     else:
            #         gp = SCIENTIFIC_PAT.findall(self.lines[i-10])[0:3]
            #     try:
            #         gp[2]
            #         for i in range(len(gp)):
            #             gp[i] = gp[i][0:-3]+'E'+gp[i][-3:]
            #         gp = '   '.join(gp)
            #         globalpointList.append(gp)
            #     except:
            #         globalpointList.append('NO')

        # Building the df
        df = pd.DataFrame()
        df['Surface'] = surfaces
        df['Cell'] = cells
        df['x'] = x
        df['y'] = y
        df['z'] = z

        # Get a complete set of locations
        loc = df.drop_duplicates().set_index(['Surface', 'Cell']).sort_index()

        # Count multiple occasions
        df['count'] = 1
        global_df = df[['Cell', 'Surface', 'count']]
        global_df = global_df.groupby(['Cell', 'Surface']).sum()
        global_df = global_df.sort_values(by='count', ascending=False)
        logging.debug('global df was built')

        # --- Printing to excel ---
        logging.info('printing to excel')

        # dump in the excel output
        outfile = os.path.join(outpath, 'LPdebug_{}.xlsx'.format(self.name))
        with pd.ExcelWriter(outfile) as writer:
            # global df
            global_df.to_excel(writer, sheet_name='global')
            # loc.to_excel(writer, sheet_name='Locations')

        # visualize the cloud point
        logging.info('building the cloud point')
        point_cloud = pv.PolyData(loc[['x', 'y', 'z']].values)
        point_cloud.plot()
        outfile = os.path.join(outpath, 'LPdebug_{}.vtp'.format(self.name))
        point_cloud.save(outfile)

        logging.info('dump completed')

    def get_statistical_checks(self) -> dict[int, str]:
        """
        Retrieve the result of the 10 statistical checks for all tallies.
        They are registered as either 'Missed', 'Passed' or 'All zeros' in a
        dictionary indicized using the tallies numbers.

        Returns
        -------
        stat_checks : dict[int, str]
            keys are the tally numbers, values the result of the statistical
            checks.

        """
        # Some global key words and patterns
        start_stat_check = 'result of statistical checks'
        miss = 'missed'
        passed = 'passed'
        allzero = 'no nonzero'
        pat_tnumber = re.compile(r'\s*\t*\s*\d+')
        end = 'the 10 statistical checks are only'

        # Recover statistical checks
        statcheck_flag = False
        stat_checks = {}
        for line in self.lines:
            if line.find(start_stat_check) != -1:
                statcheck_flag = True

            elif statcheck_flag:
                # Control if is a tally line
                tallycheck = pat_tnumber.match(line)
                if tallycheck is not None:
                    tnumber = int(tallycheck.group())
                    if line.find(miss) != -1:
                        result = 'Missed'
                    elif line.find(passed) != -1:
                        result = 'Passed'
                    elif line.find(allzero) != -1:
                        result = 'All zeros'
                    else:
                        print('Warning: tally n.'+str(tnumber) +
                              ' not retrieved')

                    stat_checks[tnumber] = result

                elif line.find(end) != -1:
                    # Exit from loop when all checks are read
                    break

        return stat_checks
