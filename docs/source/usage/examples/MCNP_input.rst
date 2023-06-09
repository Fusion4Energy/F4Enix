==========
MCNP input
==========

The following are a series of usage examples of the MCNP input API.

Change densities of cells based on condition
********************************************
For different reasons, it may be useful to change the densities of cells
in an MCNP input based on some kind of logic. In this example only the cells
belonging to a specific universe will be changed, by a factor that depends
on the material itself. Many other logics can of course be implemented in
a similar way from users, including reading the cells affected and relative
density correction factors from a file.

.. code:: 
    from f4enix.output.outputAPI