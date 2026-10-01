rm -rf build
rm -rf dist

sudo apt install python3.11
sudo apt -y install socat


sed -i "s/%REPLACE_VERSION/$(git describe --tags --abbrev=0 | cut -c 2-)/g" pyproject.toml
git describe --tags --abbrev=0 | cut -c 2- > src/pyorbb/VERSION


# briefcase create --no-input
briefcase build -r --target ubuntu:24.04

briefcase package --adhoc-sign --target ubuntu:24.04
mv dist/PyOrbb*.deb dist/PyOrbb.deb