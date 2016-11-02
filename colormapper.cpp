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

#include "colormapper.h"


ColorMapper::ColorMapper()
{
}

/**
 * @brief ColorMapper::ColorMapper creates a subdivision of the hsv cone depending on the number of segments
 * @param totalC: number of segments in which the mesh is divided
 */
ColorMapper::ColorMapper(int num){
    numColors=num;

    saturation=1.0;
    value=1.0;

    stepHue=360.0/numColors;
}


/**
 * @brief ColorMapper::setRed assigns a value for the red channel
 * @param region: the index of the segment to which the triangle belongs
 */
void ColorMapper::setRed(int region){

    hue=region*stepHue;
    float hueSixty=hue/60.0;
    float X =(1.0-fabs(fmod(hueSixty,2.0)-1.0));

    if((hueSixty >=0.0 && hueSixty < 1.0) || (hueSixty >= 5.0 && hueSixty <6.0))
        red=saturation*value;
    else if((hueSixty >= 1.0 && hueSixty < 2.0) || (hueSixty >= 4.0 && hueSixty <5.0)){
        red=X;
    }
    else{
        red=0.0;
    }
}

/**
 * @brief ColorMapper::setGreen assigns a value for the green channel
 * @param region: the index of the segment to which the triangle belongs
 */
void ColorMapper::setGreen(int region){

    hue=region*stepHue;
    float hueSixty=hue/60.0;
    float X=(saturation*value) * (1.0-fabs(fmod(hueSixty,2.0)-1.0));

    if((hueSixty>=0.0 && hueSixty < 1.0) || (hueSixty >=3.0 && hueSixty < 4.0))
        green=X;
    else if(hueSixty >= 1.0 && hueSixty < 3.0)
        green=saturation*value;
    else
        green=0.0;
}


/**
 * @brief ColorMapper::setBlue assigns a value for the blue channel
 * @param region: the index of the segment to which the triangle belongs
 */
void ColorMapper::setBlue(int region){

    hue=region*stepHue;
    float hueSixty=hue/60.0;
    float X=(saturation*value) * (1.0-fabs(fmod(hueSixty,2.0)-1.0));

    if(hueSixty >= 0.0 && hueSixty < 2.0)
        blue=0.0;
    else if(hueSixty >= 3.0 && hueSixty < 5.0)
        blue=saturation*value;
    else
        blue=X;
}
