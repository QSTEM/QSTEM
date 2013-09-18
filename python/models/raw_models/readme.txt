Files in this folder are simply data structures.  They don't know how to plot themselves.

If you want to implement a new way of plotting things (say, using Matplotlib), then subclass these classes and
add in your desired plotting methods, ideally to match the methods already present in the Chaco-based plot classes one folder up.