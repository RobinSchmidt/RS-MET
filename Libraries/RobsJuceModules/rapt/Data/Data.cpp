#include "Data.h"

namespace RAPT
{

#include "ArrayTools.cpp"
#include "MatrixTools.cpp"
#include "StandardContainerTools.cpp"
#include "Buffers.cpp"
#include "MultiArray.cpp"
#include "Trees.cpp"

}


/*

For more ideas for data structures and algorithms, see:
https://github.com/gibsjose/cpp-cheat-sheet/blob/master/Data%20Structures%20and%20Algorithms.md

https://www.youtube.com/watch?v=vVL6NFzr0Rg  Every Data Structure Explained in 20 Minutes!
-Explains: array, linked list, stack, queue, deqeue, hash map, hash set, tree, binary search tree,
 heap, graph, trie, disjoint set, bloom filter, LRU cache


ToDo:
-Maybe make a class that wraps around an existing C-style array and gives it the interface of 
 std::vector such that the std::algorithms can be used on it. Maybe rsVectorView or 
 rsPseudoContainer or rsPseudoStdVector or rsRawVector

*/