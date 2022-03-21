# Parallel k-d tree

Parallel implementation of k-d tree using MPI, OpenMP and C++17.

## Summary

In computer science, a k-d dimensional tree is a space-partitioning data structure for
organizing points in a k-dimensional space. k-d trees are a useful data structure for several
applications, such as searches involving a multidimensional search key (e.g. range searches
and nearest neighbor searches) and creating point clouds. k-d trees are a special case of
binary space partitioning trees. In this brief report is presented a parallel implementations
of k-d tree which uses two standard frameworks for parallel programming, namely OpenMP
(shared memory) and MPI (distributed memory).

A k-d tree may be represented with a binary tree, whose nodes are defined approximately in the following way:

```cpp
struct knode {
    point_type point
    int axis
    knode *left;
    knode *right;
};
```

Our data can be represented using double and integer values

## Compilation

Clone the repository, a `Makefile` is available along with the following recipes:

| Recipe         | Result                                                                                                                                                              |
| ---------------| ------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `int`          | Produce a binary that accept as input a file compose with integers numbers                                                                                          |
| `double`       | Produce a binary that accept as input a file compose with double numbers                                                                                            |
| `debug_int`    | Print a k-d tree, is reccomended to use with small dataset (< 1000)                                                                                                 |
| `debug_double` | Print a k-d tree, is reccomended to use with small dataset (< 1000)                                                                                                 |

All the executable all are in the folder `bin`.

There are also an optional parameter which are used to make different executable

| Parameter | Values   | Explanation
|-----------|----------|----------------------------------|
| `src`     | `mpi`    | Compile using MPI                |
| `src`     | `omp`    | Compile using OMP                |
| `src`     | `serial` | Compile using a single processor |
|-----------|----------|----------------------------------|
