#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
CUR_DIR=$PWD
CUR_BRANCH=$(git rev-parse --abbrev-ref HEAD)

echo "Deploying to gh-pages branch"
git config user.email "bjack205@gmail.com"
git config user.name "Brian Jackson"
echo $PWD
git fetch
git checkout origin/gh-pages
echo "Copying html files"
cp -r build/docs/html/* html
sleep 0.5
echo "Adding html folder"
git add html/*
echo "Committing changes"
git commit -m "Update documentation"
echo "pushing to origin"
git push -f origin gh-pages
git checkout $CUR_BRANCH

cd $CUR_DIR

