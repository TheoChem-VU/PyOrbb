PyOrbb Command-Line Tools
=========================

PyOrbb offers two command-line interface programs that allow the user to quickly get started in their bonding analysis.
The first, ``pyorbb open``, starts the interactive user-interface that shows the most important orbital mixing situations found by PyOrbb. The second, ``pyorbb excel`` writes orbital information to a usefull Excel spreadsheet that shows the complete information extracted by PyOrbb.


``pyorbb``
----------

**usage:** ``pyorbb [-h] {excel,open} ...``

**options:**
    ``-h, --help``: show this help message and exit

PyOrbb command-line scripts:
    ``{excel,open}``

    ``excel``: Read orbital information from an ADF calculation and write them to an Excel file.

    ``open``: Start the PyOrbb analysis program


``pyorbb open``
---------------

Start the PyOrbb analysis program

**usage:** ``pyorbb open [-h] [rkf]``

**positional arguments:**
    ``rkf``: The path to the `adf.rkf` file to analyse. If omitted, starts an empty PyOrbb GUI.

**options:**
    ``-h, --help``: show this help message and exit


``pyorbb excel``
----------------
Read orbital information from an ADF calculation and write them to an Excel file.

**usage:** ``pyorbb excel [-h] [-o OUTPUT] rkf``


**positional arguments:**
    ``rkf``: The path to the ``adf.rkf`` file to summarize in an Excel file.

**options:**
    ``-h, --help``: show this help message and exit
    ``-o OUTPUT, --output OUTPUT``: Set the output Excel file to write to.