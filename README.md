# Fold Molding

This project is to design decomposition algorithm for free form geometry with the constraint of height field.


## Height Field
Height field is the constraint that decomposes the target object into several parts. Every part should be seperated in one direction without conflict.


## Algorithm

This algorithm aims to detect height field, avoid overlaping and inter-locking. It maps an object modeled by triangle mesh into a graph. Every triangle is mapped to a vertex in the graph. Edges are connected to the vetices whose corresponding triangles are adjacent. We assign weight to every edge according to the length of its corresponding triangle edge. Alpha expansion algorithm is applied to achieve the graph cut problem. We then adjust the segmentations by the graph cut results


## Run

This project depends on libigl. You should install libigl first, and then run the codd

    This project reads .obj or .off files only.

## Dependencies

We recommend you to install libigl using git via:

    git clone --recursive https://github.com/libigl/libigl.git

If you have installed libigl at `/path/to/libigl/` then a good place to clone
this library is `/path/to/libigl-example-project/`.
