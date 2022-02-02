#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
CUR_DIR=$PWD
CUR_BRANCH=$(git rev-parse --abbrev-ref HEAD)

echo "Deploying to gh-pages branch"
git config user.email "bjack205@gmail.com"
git config user.name "Brian Jackson"
echo $PWD
git fetch
git checkout origin/gh-pages && exit 1
echo "Copying html files"
cp -r build/docs/html/* html && exit 1
sleep 0.5
echo "Adding html folder"
git add html/* && exit 1
echo "Committing changes"
git commit -m "Update documentation" && exit 1
echo "pushing to origin"
git push -f origin gh-pages && exit 1
git checkout $CUR_BRANCH && exit 1

cd $CUR_DIR

