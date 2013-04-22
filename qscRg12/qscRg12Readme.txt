qscRg12, written by Ray Hu (GXH115@bham.ac.uk)

Purpose: a simple C program to help users to get the edge of scanning windows X and Y after rotation, so that we won't need the graphical user interface when we conduct the simulation under Linux OS remotely.

Example usage:

qscRg12 ABC.cfg X.tilting Y.tilting Z.tilting Spaced
for example.
qscRg12 ABC.cfg 10 20 0 0.5
would mean that outputting the X_start/end Y_start/end for X rotation 10 deg, Y rotation 20 deg. And also push the edges out for 0.5 Angstrom, because we don't want the scanning window cut right on the atoms on the rim.

The program would return information on screen like this.
Xmax-> 14.524338 Xmin-> 5.093791
Ymax-> 13.863408 Ymin-> 3.394938
Spaced Xmax-> 15.024338 Xmin-> 4.593791
Spaced Ymax-> 14.363408 Ymin-> 2.894938
15.024338 4.593791 14.363408 2.894938

I use it with another shell script qsc_AutoScriptMS. like
qsc_AutoScriptMS X.tilting Y.tilting Z.tilting Spaced
no need to specified which .cfg to use anymore. but make sure there is only one .cfg available at current folder.