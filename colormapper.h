/*
 *
 *   2014
 *   Author:       Giulia Picciau - DIBRIS, Università degli studi di Genova
 *   Supervisors:  Leila De Floriani - DIBRIS, Università degli studi di Genova
 *                 Patricio Simari - Department of Electrical Engineering and Computer Science, The Catholic University of America
 *
 *   Title:          Fast and scalable mesh superfacets
 *   Submission to Eurographics Symposium on Geometry Processing 2014
 *
 *
 **/
#ifndef COLORMAPPER_H
#define COLORMAPPER_H


#include <math.h>
#include <iostream>

/**
 *
 *   @author Giulia Picciau
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
    ColorMapper();
    ColorMapper(int);

    void setRed(int);
    void setGreen(int);
    void setBlue(int);

    float red;
    float green;
    float blue;

    int numColors;

    float stepHue;
    float value;
    float saturation;

    float hue;
};

#endif // COLORMAPPER_H
