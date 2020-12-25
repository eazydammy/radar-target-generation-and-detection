# Radar Target Generation and Detection

## Implementation steps for the 2D CFAR process

* Determine the number of Training cells for each dimension. Similarly, pick the number of guard cells.
* Slide the cell under test across the complete matrix. Make sure the CUT has margin for Training and Guard cells from the edges.
* For every iteration sum the signal level within all the training cells. To sum convert the value from logarithmic to linear using db2pow function.
* Average the summed values for all of the training cells used. After averaging convert it back to logarithmic using pow2db.
* Further add the offset to it to determine the threshold.
* Next, compare the signal under CUT against this threshold.
* If the CUT level > threshold assign it a value of 1, else equate it to 0.

## Selection of Training, Guard cells and offset

Tr = 12;
Td = 6;

Gr = 6;
Gd = 3;

offset = 1.4; (by SNR value in dB)

The selection of values were obtained by tweaking for accuracy. Also, the higher the value of cells, the longer it took to run the program.

## Steps taken to suppress the non-thresholded cells at the edges

The few non-thresholded cells are set to zero by creating a logical index to cells that have not been thresholded to zero or one.
This is to keep the map size same as it was before CFAR

> RDM(RDM~=0 & RDM~=1) = 0
