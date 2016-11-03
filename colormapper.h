/*
 *
 *   2014
 *   Author:       Giulia Picciau - DIBRIS, Università degli studi di Genova
 *   Supervisors:  Leila De Floriani - DIBRIS, Università degli studi di Genova
 *                 Patricio Simari - Department of Electrical Engineering and Computer Science, The Catholic University of America
 *
 *   Title:          Fast and scalable mesh superfacets
 *   Submission to Pacific Graphics 2014
 *
 *
 **/
#ifndef COLORMAPPER_H
#define COLORMAPPER_H


#include <math.h>
#include <iostream>

/**
 *
 *   @author Giulia Picciau - Università degli studi di Genova
 *
 *   @date 2014
 *
 *   @brief ColorMapper class generates a colormap to assign a different color to faces of a mesh depending on
 *   which cluster they belong to
 *
 *
 */
class ColorMapper
{
public:
    //Constructor
    ColorMapper();

    /**
     * @brief ColorMapper Constructor
     * @param num The number of regions of the segmentation
     */
    ColorMapper(int num);

    /*
     * Functions to assign the color to a triangle depending on the region to which it belongs
     */
    void setRed(int);
    void setGreen(int);
    void setBlue(int);

    /*
     *Color intensities
     */
    float red;
    float green;
    float blue;

    //Number of different colors required (=number of regions)
    int numColors;

    //sampling interval on the hue circumference
    float stepHue;
    //Value and saturation (default=1.0)
    float value;
    float saturation;

    //hue of the triangle (depends on the region)
    float hue;
};

#endif // COLORMAPPER_H
