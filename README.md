# Autoknit 2

A re-implementation of autoknit; very much in-progress.

## Status

Unless otherwise marked, items are incomplete.

- Mesh peeling
 - 3D Model Loading
 - Constraint specification (probably involves [smoothed?] shortest-path computations)
 - Interpolation
 - Geodesic extraction + trimming
- Mesh linking
- Tracing
- Scheduling
 - Partially re-implemented, but may have gone down the wrong path (in terms of scheduling tubes instead of connections).
 - Traced stitch format
  - `Stitch.*pp`
 - DAG embedding
  - `embed_DAG.*pp` -- untested
 - Transfer Planning
  - `plan_transfers*pp` -- mostly complete, pending testing (and potential modification for ordered decreases)

## Data

Models are vertices + triangles, and must be manifold.
Probably stored in .obj format.

Constraints are (connected) points on the surface of the model. Connections are (smoothed?) surface geodesics.


## Thoughts

It would be cool to figure out how to do the mesh extraction step on a signed distance field because that gives better feature size control (e.g. could pre-filter to avoid aliasing during stitch extraction).

Lots of alternate mesh extraction strategies probably exist. Anything that can chop a 3D model into connected oriented tubes provides a starting point.
