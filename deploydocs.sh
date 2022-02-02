#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
CUR_DIR=$PWD
CUR_BRANCH=$(git rev-parse --abbrev-ref HEAD)

echo "Deploying to gh-pages branch"
git config user.email "bjack205@gmail.com"
git config user.name "Brian Jackson"
git fetch
git checkout gh-pages && exit 1
cp -r build/docs/html/* html && exit 1
sleep 0.5
git add html/* && exit 1
git commit -m "Update documentation" && exit 1
git push -f origin gh-pages && exit 1
git checkout $CUR_BRANCH && exit 1

cd $CUR_DIR

