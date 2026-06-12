rm -rf build
rm -rf dist

apt install python3.11

sed -i "s/%REPLACE_VERSION/$(git describe --tags --abbrev=0 | cut -c 2-)/g" pyproject.toml
git describe --tags --abbrev=0 | cut -c 2- > src/pyorbb/VERSION

# briefcase create --no-input
briefcase build -r

briefcase package --adhoc-sign
mv dist/PyOrbb*.deb dist/PyOrbb.deb