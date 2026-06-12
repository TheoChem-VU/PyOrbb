rm -rf build
rm -rf dist

sed -i "s/%REPLACE_VERSION/$(git describe --tags --abbrev=0 | cut -c 2-)/g" pyproject.toml
git describe --tags --abbrev=0 | cut -c 2- > src/pyorbb/VERSION

briefcase create -a PyOrbb --no-input
briefcase build -r

briefcase package --adhoc-sign
mv dist/PyOrbb*.deb dist/PyOrbb.deb