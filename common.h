#ifndef COMMON_H
#define COMMON_H

#include "normals.h"
#include "Timer.h"
#include "Reader.h"
#include <tr1/unordered_map>
#include <queue>
#include <tr1/unordered_set>
#include <numeric>
#include <QString>
#include <QStringList>
#include <cfloat>
#include <fstream>



/// types for indices of triangles and for adjacencies between them
typedef unsigned long int faceind;
typedef unsigned long long int edgekey;

/// types for grid initialization
typedef signed char coordind;
typedef int gridkey;

/*!
 * \brief The pointDist struct: used for the priority queue in the Dijkstra-based expansion
 */
struct pointDist{
    int indexP;
    double distanceP;
};

/*!
 * \brief The compare struct: redefinition of the > operator to be able to use a priority queue
 */
struct compare{

    bool operator ()(pointDist p1, pointDist p2){

        return p1.distanceP > p2.distanceP;
    }
};

#endif // COMMON_H
