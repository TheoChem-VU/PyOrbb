rm -rf build
rm -rf dist

sed -i '' 's/%REPLACE_VERSION/$(git describe --tags --abbrev=0)/g' pyproject.toml

briefcase build -r

cp new_vtk.py build/pyfmo/macos/app/PyOrbb.app/Contents/Resources/app_packages/vtk.py
cat removed_files.txt | while read line 
do
   rm $line
done

briefcase package --adhoc-sign
