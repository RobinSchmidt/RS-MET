#include "CellularExperiments.h"


// ToDo: write classes like CellularAutomaton2N3S ("3I3S": each node has 3 inputs and 3 states)




void testCellularAutomaton1()
{


  int dummy = 0;
}



/*
For a cellular automaton where each node has N inputs (and therefore N-1 neighbours) and S 
states, the number of possible rules is: S^(S^N) because there are S^N possible different input 
patterns and for each such input pattern, there are S possible outcomes. With N=3 inputs (2 
neighbours) for each node and S=2 states, we have 2^3 = 8 different input patterns and therefore
2^8 = 256 possible rules, which can be enumerated according to Stephen Wolfram's scheme: ...
*/

/*
2D regular layout with 3 connections for each node:

  |         |         |         |
--O----O----O----O----O----O----O----O--  1st row
       |         |         |         |
       |         |         |         |
--O----O----O----O----O----O----O----O--  2nd row
  |         |         |         |
  |         |         |         |
--O----O----O----O----O----O----O----O--  3rd row
       |         |         |         |
       |         |         |         |
--O----O----O----O----O----O----O----O--  4th row
  |         |         |         |

...this "brickwall" pattern can be stretched out into a totally regular hexagonal pattern where the
corners of each hexagon are the nodes, like this:

  O   O   O   O       \
 / \ / \ / \ / \       |-- 1st row
    O   O   O   O     /
    |   |   |   |
    O   O   O   O     \
 \ / \ / \ / \ /       |-- 2nd row
  O   O   O   O       /
  |   |   |   |    
  O   O   O   O       \
 / \ / \ / \ / \       |-- 3rd row
    O   O   O   O     /
    |   |   |   |
    O   O   O   O     \
 \ / \ / \ / \ /       |-- 4th row
  O   O   O   O       /
  |   |   |   |    

That means, a regular pattern in 2D where each node has 3 connections, can conveniently be 
represented by a 2D array in code, where in each row with even index, the even column-indexed nodes
have vertical connections only to the row above and odd column-indexed nodes have vertical 
connections only to the row below. In odd indexed rows, it's the other way around. With periodic
boundaries (wraparound), thje number of rows and columns should be even. In order to have unique
neighbours in all directions, it must be >=3, the 1st even number >= is 4, so this is the minimum 
number for rows and columns




  tbc.
*/



/*
under construction:
3D, 3 connections: O node connects downward (into screen), X node connects upward (out of screen)

       |         |         |         |
  O----X    X----O    O----X    X----O   
  |         |         |         |
  |         |         |         |
--X    X----O    O----X    X----O    O--
  |    |
  |    |
  O    O    



       |         |         |         |     
--O    O----X    X----O    O----X    X--
  |         |         |         |
  |         |         |         | 
  X----O    O----X    X----O    O----X
       |         |         |         |
       |         |         |         |
  X    X----O    O----X    X----O    O--
  |         |         |         |
  |         |         |         |
  O----X    X----O    O----X    X----O
       |         |         |         |


         

4 basic elements? O----X, X----O

*/