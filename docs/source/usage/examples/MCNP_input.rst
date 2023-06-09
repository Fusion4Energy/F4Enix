Reading and setting values on numjuggler cards
**********************************************

The F4Enix MCNP input parser is based on numjuggler, hence the cells,
surfaces and all data cards are stored as numjuggler.Card.
The recommended way to read and set values on this cards is to use the
:code:`._get_value_by_type()` and :code:`._set_value_by_type()`.

.. todo::
    complete with the possible types inferring them from numjuggler


Change densities of cells based on condition
********************************************
Changing the cell densities of an MCNP input is a common operation.
F4Enix allow to do so with different kinds of logics.
In this example only the cells
belonging to a specific universe will be changed by a factor that depends
on the material itself. Void cells are not considered and all other cells
in the universe for which a factor has not be specified are scaled by
a constant factor.

.. code-block:: python

    from f4enix.output.MCNPinput import Input

    # Load the input file
    inp = Input.from_input('test_file.i')

    # Set some density correction factors for specific materials
    density_factors = {400: 1.1, 300: 10}

    # Select the universe of interest
    universe = 100

    # Cycle on the cells dictionary
    for idx_cell, cell in self.cells.items():
        # check for correct universe
        if cell._get_value_by_type('u') == universe:
            # Change density based on material
            mat = cell._get_value_by_type('mat')
            # get the density of the cell
            rho = cell.get_d()
            if mat in [400, 300]:
                cell.set_d(rho*density_factors[mat])
            elif mat == 0:
                pass  # do not change density of void cells
            else:
                # scale all other densitites by a constant factor
                cell.set_d(rho*0.5)
    
    # Print the modified file
    inp.write('new_input.i')

Here is another example in case the cell and density data was directly
provided in a .csv file

.. code-block:: python

    from f4enix.output.MCNPinput import Input
    import pandas as pd

    # Read a density .csv of two columns: 'cell' and 'factor'
    densities = pd.read_csv('densities.csv')

    # Load the input file
    inp = Input.from_input('test_file.i')

    # --- Modify all cells listed in the dataframe ---
    for _, row in densities.items():
        # get the cell
        cell = inp.cells[row['cell']]
        # get the current density
        rho = cell.get_d()
        # correct the density
        cell.set_d(rho*row['factor'])
    
    # Print the modified file
    inp.write('new_input.i')
    