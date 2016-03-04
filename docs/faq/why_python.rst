.. _why_python:

=============
 Why Python?
=============

From Nipy's documentation::

The choice of programming language has many scientific and practical
consequences. Matlab is an example of a high-level language. Languages
are considered high level if they are able to express a large amount
of functionality per line of code; other examples of high level
languages are Python, Ruby, Octave, and R. In contrast, C is a
low-level language. Low level languages can achieve higher execution
speed, but at the cost of code that is considerably more difficult to
read. C++ and Java occupy the middle ground sharing the advantages and
the disadvantages of both levels.

Low level languages are a particularly ill-suited for exploratory
scientific computing, because they present a high barrier to access by
scientists that are not specialist programmers. Low-level code is
difficult to read and write, which slows development
([Prechelt2000ECS]_, [boehm1981]_, [Walston1977MPM]_) and makes it more
difficult to understand the implementation of analysis
algorithms. Ultimately this makes it less likely that scientists will
use these languages for development, as their time for learning a new
language or code base is at a premium. Low level languages do not
usually offer an interactive command line, making data exploration
much more rigid. Finally, applications written in low level languages
tend to have more bugs, as bugs per line of code is approximately
constant across many languages [brooks78].

In contrast, interpreted, high-level languages tend to have
easy-to-read syntax and the native ability to interact with data
structures and objects with a wide range of built-in
functionality. High level code is designed to be closer to the level
of the ideas we are trying to implement, so the developer spends more
time thinking about what the code does rather than how to write
it. This is particularly important as it is researchers and scientists
who will serve as the main developers of scientific analysis
software. The fast development time of high-level programs makes it
much easier to test new ideas with prototypes. Their interactive
nature allows researchers flexible ways to explore their data.

SPM is written in Matlab, which is a high-level language specialized
for matrix algebra. Matlab code can be quick to develop and is
relatively easy to read. However, Matlab is not suitable as a basis
for a large-scale common development environment. The language is
proprietary and the source code is not available, so researchers do
not have access to core algorithms making bugs in the core very
difficult to find and fix. Many scientific developers prefer to write
code that can be freely used on any computer and avoid proprietary
languages. Matlab has structural deficiencies for large projects: it
lacks scalability and is poor at managing complex data structures
needed for neuroimaging research. While it has the ability to
integrate with other languages (e.g., C/C++ and FORTRAN) this feature
is quite impoverished. Furthermore, its memory handling is weak and it
lacks pointers - a major problem for dealing with the very large data
structures that are often needed in neuroimaging. Matlab is also a
poor choice for many applications such as system tasks, database
programming, web interaction, and parallel computing. Finally, Matlab
has weak GUI tools, which are crucial to researchers for productive
interactions with their data.


.. [boehm1981]
   Boehm, Barry W. (1981) *Software Engineering Economics*. Englewood
   Cliffs, NJ: Prentice-Hall.

.. [Prechelt2000ECS]
   Prechelt, Lutz. 2000. An Empirical Comparison of Seven Programming
   Languages. *IEEE Computer* 33, 23--29.

.. [Walston1977MPM]
   Walston, C E, and C P Felix. 1977. A Method of Programming
   Measurement and Estimation. *IBM Syst J* 16, 54-73.
