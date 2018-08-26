# Autoknit

A re-implementation of ["Automatic Machine Knitting of 3D Meshes"](https://textiles-lab.github.io/publications/2018-autoknit/).
Does not match the code used in the paper exactly, but is close in most regards.

## License

This code is placed in the public domain.

## Building

You will need Perforce's Jam/MR tool to build, along with the SDL2 library (opengl + mouse handling), the glm math library, and the Eigen linear algebra library.

### MacOS setup
```
#clone repository:
git clone git@github.com:textiles-lab/autoknit
cd autoknit
git submodule init
git submodule update

#install prerequisite libraries and build tool:
brew install ftjam sdl2 glm eigen
```

### Linux setup
TBD

### Window setup

First, install (ftjam)[https://www.freetype.org/jam/] somewhere in your ```%path%``` so you can run it from a command prompt. (Also make sure that git is installed in such a way that it can be run from a command prompt.)

Then, from a ```Visual Studio 2017 > x64 Native Tools Command Prompt for VS 2017``` do:
```
#clone repository:
git clone git@github.com:textiles-lab/autoknit
cd autoknit
git submodule init
git submodule update

#install pre-built versions of sdl2 and glm libraries:
git clone git@github.com:ixchow/kit-libs-win
#get a copy of the Eigen headers:
git clone git@github.com:eigenteam/eigen-git-mirror eigen
```

### Linux/Windows/MacOS build
```
cd autoknit
jam -j8
```

The ```-j8``` parameter means to run up to 8 compilation jobs in parallel. You may wish to adjust this for your particular machine.

## Usage

TODO

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
