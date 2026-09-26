Installation Manual PyOrbb |ProjectVersion|
##################################################

PyOrbb Graphical User Interface
===============================

Download and install the latest PyOrbb version here:

`<https://github.com/TheoChem-VU/PyOrbb/releases/latest>`_

Select the appropriate installer for your operating system and follow the installation instructions.



PyOrbb Python Library
=====================

Installation of PyOrbb for use with Python is straightforwardly done via ``pip``:

.. code-block::

	python -m "pip install pyorbb"


Expert users who would like to work on the code may also install the source code:

.. code-block::
	
	git clone git@github.com:TheoChem-VU/PyOrbb.git
	cd PyOrbb
	python -m "pip install -e ."


To check if PyOrbb is correctly installed try to run the following:

.. code-block::

	pyorbb --help

which should output:

.. code-block::

	usage: pyorbb [-h] {excel,open} ...

	options:
	  -h, --help    show this help message and exit

	PyOrbb command-line scripts:
	  {excel,open}
	    excel       Read orbital information from an ADF calculation and write them to an Excel file.
	    open        Start the PyOrbb analysis program
