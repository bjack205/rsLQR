#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
CUR_DIR=$PWD
CUR_BRANCH=$(git rev-parse --abbrev-ref HEAD)

echo "Deploying to gh-pages branch"
git config user.email "bjack205@gmail.com"
git config user.name "Brian Jackson"
echo $PWD
git checkout origin/gh-pages
cp -r build/docs/html/* html
git add html/*
git commit -m "Update documentation"
git push -f origin gh-pages
git checkout $CUR_BRANCH

cd $CUR_DIR

