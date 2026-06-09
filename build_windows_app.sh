rm -rf build;
rm -rf dist;

sed -i "s/%REPLACE_VERSION/$(git describe --tags --abbrev=0 | cut -c 2-)/g" pyproject.toml;
git describe --tags --abbrev=0 | cut -c 2- > src/pyorbb/VERSION

briefcase build -r;
briefcase package --adhoc-sign;

mv dist/PyOrbb*.msi dist/PyOrbb.msi;

sed -i "s/version= \"$(git describe --tags --abbrev=0 | cut -c 2-)\"/version= \"%REPLACE_VERSION\"/g" pyproject.toml;