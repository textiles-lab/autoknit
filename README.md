# Autoknit

A re-implementation of ["Automatic Machine Knitting of 3D Meshes"](https://textiles-lab.github.io/publications/2018-autoknit/).
Does not match the code used in the paper exactly, but is close in most regards.

The latest version of this code is available at https://github.com/textiles-lab/autoknit .

## License

This code is placed in the public domain.

## Pre-Built Releases

This code is automatically built using Github Actions; check for releases at https://github.com/textiles-lab/autoknit/releases .
However, it is also reasonably straightforward to build on your own using the steps below.

## <a name="building"></a>Building

Building autoknit is handled with the single-file build tool [maek](https://github.com/ixchow/maek), and uses the [nest-libs](https://github.com/15-466/nest-libs/releases) pre-built library package (for SDL2 and glm), along with the Eigen linear algebra library, which you will need to fetch separately.

### <a name="mac"></a>MacOS setup
```
#install nodejs (used for build script and to post-process scheduled output):
brew install node

#extract nest-libs package as a sibling of autoknit folder:
curl 'https://github.com/15-466/nest-libs/releases/download/v0.13/nest-libs-macos-v0.13.tar.gz' -L -O
tar xvfz nest-libs-macos-v0.13.tar.gz

#clone repository:
git clone git@github.com:textiles-lab/autoknit
cd autoknit
git submodule init
git submodule update
```

### <a name="linux"></a>Linux setup

```
#install nodejs (used for build script and to post-process scheduled output):
sudo apt-get install nodejs

#extract nest-libs package as a sibling of autoknit folder:
curl 'https://github.com/15-466/nest-libs/releases/download/v0.13/nest-libs-linux-v0.13.tar.gz' -L -O
tar xvfz nest-libs-macos-v0.13.tar.gz

#clone repository:
git clone git@github.com:textiles-lab/autoknit
cd autoknit
git submodule init
git submodule update
```

### <a name="windows"></a>Windows setup

First, make sure that git is installed in such a way that it can be run from a command prompt; also install [nodejs](https://nodejs.org/), which is used by the build script and to post-process scheduled output.

Then, from a `Visual Studio 2022 > x64 Native Tools Command Prompt for VS 2022` command prompt do:
```
#clone repository:
git clone git@github.com:textiles-lab/autoknit
cd autoknit
git submodule init
git submodule update
```

Finally, download the nest-libs package (https://github.com/15-466/nest-libs/releases/download/v0.13/nest-libs-windows-v0.13.zip) and place the contained `nest-libs/` folder next to the `autoknit` folder where you checked out the code.

### Linux/Windows/MacOS build

NOTE: on windows, be sure to use a `Visual Studio 2022 > x64 Native Tools Command Prompt for VS 2022`.

```
cd autoknit
node Maekfile.js
```

Command line flags for maek include:
- `-j8` run with 8 compilation jobs in parallel. Adjust the number to suit your machine. Defaults to number of cores + 1.
- `-q` quit on first error
- Specifically, `node Maekfile.js -j1 -q` is useful when you want to work on one error at a time during development!

Note that maek is a small build system, entirely contained in [Maekfile.js](Maekfile.js). You can read the file to see more about command line options or to see how the build is structured.

## <a name="usage"></a>Usage

Step-by-step instructions for creating knitting machine instructions for the [misc-cactus.obj](https://github.com/textiles-lab/autoknit-tests/raw/master/models/misc-cactus.obj) model from the [autoknit-tests](https://github.com/textiles-lab/autoknit-tests) repository.

### <a name="constraints"></a>Step 1: Constraints

Launch the interface, telling it to load from ```misc-cactus.obj``` and to save constraints to ```misc-cactus.cons```:

```
cd dist
./interface obj:misc-cactus.obj constraints:misc-cactus.cons
```
(Note that the `interface` executable is built in the `dist/` subdirectory, so you will need to change to that directory before running it.)

You will see a 3D view of the loaded model:

![](usage/constraints-01.png)

You can <a name="rotate"></a>rotate this model with the right mouse button, zoom with the mouse wheel, and <a name="pan"></a>pan with <kbd>shift</kbd> + right mouse button.

The point on the surface your mouse is over will be highlighted with a grey sphere (the red, green, and blue spheres show the location of the corners of the triangle:

![](usage/constraints-02.png)

Pressing the <a name="constraint-keybinding"></a><kbd>c</kbd> key will add a constraint:

![](usage/constraints-03.png)

Pressing the <kbd>c</kbd> key while hovering over a constraint will add a connected constraint point, which you can press left mouse button to place:

![](usage/constraints-04.png)

You can click and drag constraint points to move them:

![](usage/constraints-05.png)

Pressing the <a name="red-cons-keybinding"></a><kbd>+</kbd> key while hovering over a constraint will move it later in time (redder) while pressing the <a name="blue-cons-keybinding"></a><kbd>-</kbd> key will move it earlier in time (bluer):

![](usage/constraints-06.png)

Pressing the <a name="delete-keybinding"></a><kbd>X</kbd> key while hovering over a constraint will delete it:

![](usage/constraints-07.png)

Before you can proceed to the next step, you will need to create constraints for (at least) all of the boundaries of the cactus:

![](usage/constraints-08.png)

You can also place constraints elsewhere on the model to control the knitting direction, and use the <a name="cut-keybinding"></a><kbd>R</kbd> key to cut out a region around a constraint (useful for starting/ending on meshes without boundaries).

### Step 2: <a name="peeling"></a>Peeling/Linking

Now that constraints are specified, the rest of the steps proceed automatically. However, the interface can provide visualization tools to show you what is happening.

#### Manual method:

Load the cactus object and the constraints into the interface. The <a name="obj-scale"></a>```obj-scale``` parameter tells the interface how much to scale the object, while the <a name="stitch-width"></a>```stitch-width``` and <a name="stitch-height"></a>```stitch-height``` parameters give the stitch size relative to the scaled object. The <a name="save-traced"></a>```save-traced:``` parameter tells the interface where to save its traced stitches.

```
./interface obj:misc-cactus.obj load-constraints:misc-cactus.cons obj-scale:10.0 stitch-width:3.66 stitch-height:1.73 save-traced:misc-cactus.st
```

Press the <a name="peel-keybinding"></a><kbd>p</kbd> key to step through peeling:

![](usage/peel-anim.gif)

During peeling, you can use the <a name="graph-keybinding"></a><kbd>g</kbd> key to show or hide the portions of the row-column graph created so far, and the <a name="toggleview-keybinding"></a><kbd>s</kbd> key to toggle whether the original model, interpolated value, or current slice model are being shown:

![](usage/peel-model.png) ![](usage/peel-times.png) ![](usage/peel-slice.png)
![](usage/peel-model-graph.png) ![](usage/peel-times-graph.png) ![](usage/peel-slice-graph.png)

Once the peeling has finished, you can press <a name="traced-keybinding"></a><kbd>t</kbd> to create and save the traced path:

![](usage/peel-done.png) ![](usage/traced.png)

#### <a name="peel-step"></a>Automatic method:

If you don't want to press <kbd>p</kbd> a whole lot, you can just pass the ```peel-step:N``` option to do ```N``` steps of peeling. ```peel-step:-1``` will peel until the mesh is finished.

```
./interface obj:misc-cactus.obj load-constraints:misc-cactus.cons obj-scale:10.0 stitch-width:3.66 stitch-height:1.73 save-traced:misc-cactus.st peel-step:-1
```

### <a name="scheduling"></a>Step 3: Scheduling

Now that the traced stitches have been created, they need to be assigned knitting machine needles. We call this step scheduling, and it has its own executable, called ```schedule```.
The only parameters used by schedule are ```st:``` which gives an input stitches file and ```js:``` which gives an output javascript file:

```
./schedule st:misc-cactus.st js:misc-cactus.js
```

Schedule doesn't have an UI; it just does a relatively large combinatorial search and then dumps its output into a javascript file.

### <a name="knitout"></a>Step 4: Knitout

Running the javascript file created by ```schedule``` will create knitout instructions, ready for use on your machine:
```
NODE_PATH=.. node misc-cactus.js out:misc-cactus.k
```

Note that the javascript file created by ```schedule``` uses some helper functions defined in the [```node_modules/autoknit.js```](node_modules/autoknit.js) file to do things like cast on tubes, bring in/out yarns, and perform transfers. You may want to customize ```autoknit.js``` for your machine.

## Status By Pipeline Step

This implementation is mostly complete, but is not fully working.

- Interface/wrapper - working.
- Model (obj) loading - working.
- Constraint specification - working.
- Peeling - mostly working.
    - Could be improved to handle ending with short rows.
    - Could be improved to deal with orphaned chains.
- Linking (including split/merge cases) - working.
- Tracing - mostly working.
    - Could be improved with more extensive ancestor traversal when tucking at the ends of short rows.
    - Sometimes generates short yarns; next-stitch-picking heuristic could be improved.
    - Might want to add a lazy vs eager switch for moving to the next row after splits. (Currently, the behavior is eager).
- Scheduling - working for small cases only.
    - Need to add a greedy version (currently only optimal is used).
- Knitting instructions -- mostly working.
    - Need a better yarn-in function for split tubes that tucks on front/back and then drops later.
    - Should add the option to use separate yarn for starting tubes
	- Should add the option to tuck yarn in from the edge of the beds
