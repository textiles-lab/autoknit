rm -f misc-cactus.st misc-cactus.st misc-cactus.js misc-cactus.k misc-cactus.dat
node Maekfile.js
dist/interface obj:../autoknit-tests/models/misc-cactus.obj constraints:../autoknit-tests/constraints/misc-cactus.cons obj-scale:10.0 stitch-width:3.66 stitch-height:1.73 save-traced:misc-cactus.st peel-test:-1
dist/schedule st:misc-cactus.st js:misc-cactus.js
node misc-cactus.js
node ../knitout-backend-swg/knitout-to-dat.js misc-cactus.k misc-cactus.dat
