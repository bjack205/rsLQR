#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
CUR_DIR=$PWD
CUR_BRANCH=$(git rev-parse --abbrev-ref HEAD)

echo "Deploying to gh-pages branch"
git config user.email "bjack205@gmail.com"
git config user.name "Brian Jackson"
echo "Fetch"
git fetch
echo "Checkout"
git checkout gh-pages || exit 1
echo "Copy"
cp -r build/docs/html/* html || exit 1
sleep 0.5
echo "Add"
git add html/* || exit 1
echo "Commit"
git commit -m "Update documentation" || exit 1
echo "Push"
git push -f origin gh-pages || exit 1
echo "Checkout"
git checkout $CUR_BRANCH || exit 1

cd $CUR_DIR

