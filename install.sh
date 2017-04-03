#!/bin/bash
mvn package
cd src/main/webapp
npm install
npm run build
cd ../../../

