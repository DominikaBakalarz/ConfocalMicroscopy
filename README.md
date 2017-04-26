# ConfocalMicroscopy
Part of the Origami BioBandage project. The code is supposed to read the confocal microscopy pictures and recognize alive and dead cells.

Files lodepng.h and lodepng.cpp are required extra libraries. File stemscan.cpp contains a code that I have written myself. The code is supposed to read the pictures from the folder input-pictures (one picture at a time) and then it recognises alive and dead cells (alive cells are green, dead are red and black colour represents the background, I’ve calculated that cell of this kind (stem cell) is 16 pixels  wide on average) and marks them with white and blue dots, and returns modified picture.

File stemscan requires 4 arguments, given in this order matchPercent, interPercent, cellSize, numberOfTrials. First argument describes what percentage of pixels in given cirle needs to green/red in order to marks cell as alive/dead respectively, second describes the percentage of cells overlapping that we allow - since we assume that cells are circles, small overlapping gives more accurate results, third argument describes cell size - most accurate is radius of 8 pixels, forth describes how many times do we want to search for the cells.
