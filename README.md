# ConfocalMicroscopy
Part of the Origami BioBandage project. The code is supposed to read the confocal microscopy pictures and recognize alive and dead cells.

Files lodepng.h and lodepng.cpp are required extra libraries. File stemscan.cpp contains a code that I have written myself. The code is supposed to read the pictures from the folder input-pictures (one picture at a time) and then it recognises alive and dead cells (alive cells are green, dead are red and black colour represents the background, I’ve calculated that cell of this kind (stem cell) is 16 pixels  wide on average) and marks them with white and blue dots, and returns modified picture.
