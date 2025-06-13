rm -r _build
rm -r api
sphinx-apidoc -f -e -o ./api ../src/pyfmo -t _templates
make html