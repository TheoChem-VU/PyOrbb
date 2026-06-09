Installation Guide
==================

Installation of PyOrbb is straightforwardly done via ``pip``:

.. code-block::

	python -m "pip install pyorbb"


For expert users who would like to work on the code we also offer installation of the source code:

.. code-block::
	
	git clone git@github.com:TheoChem-VU/PyOrbb.git
	cd PyOrbb
	python -m "pip install -e ."


To check if PyOrbb is correctly installed try to run the following:

.. code-block::

	pyorbb --help

which should output:

.. code-block::

	usage: pyorbb [-h] {excel,analyse} ...

	options:
	  -h, --help       show this help message and exit

	PyOrbb command-line scripts:
	  {excel,analyse}
	    excel          Read orbital information from an ADF calculation and write them to an Excel file.
	    analyse        Start an interactive PyOrbb orbital diagram.
