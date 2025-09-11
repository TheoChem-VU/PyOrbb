rm -rf build
rm -rf dist

briefcase build -r
briefcase package --adhoc-sign
