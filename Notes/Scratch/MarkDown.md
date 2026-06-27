In this document, I give some examples of important MarkDown constructs suitable for copy-and-paste.

LaTeX Equations
---------------

$`E = mc^2`$,  $E = mc^2$, $`\boxed{  E = mc^2}`$, $\boxed{E = mc^2}$

Boxed and centered equations:

$$\boxed{
  \left( \sum_{k=1}^n a_k b_k \right)^2 
  \leq 
  \left( \sum_{k=1}^n a_k^2 \right) \left( \sum_{k=1}^n b_k^2 \right)
}$$


Comments
--------

Sometimes we may want to add some preliminary text or comments about what to do into out markdown
files which should not be rendered into the output. There are various ways to achieve this. To see
these, you need to look at the mardown source code. In a rendered output, you should see nothing
here:

<!---
your comment goes here
and here
-->

[comment]: <> (This is a comment, it will not be included)
[comment]: <> (in  the output file unless you use it in)
[comment]: <> (a reference style link.)

[//]: <> (This is also a comment.)

[//]: # (This may be the most platform independent comment)


Code Snippets
-------------


Resources
---------

https://www.markdownguide.org/basic-syntax/  
https://www.markdownguide.org/extended-syntax/