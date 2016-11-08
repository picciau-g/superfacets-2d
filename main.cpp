#include "segmenter.h"
#include "meshvisualizer.h"
#include "vertexbasedsegmenter.h"
#include <stdio.h>
#include <QApplication>

using namespace std;

/**
 * @author Giulia Picciau
 * Submission to Eurographics Symposium on Geometric Processing 2014
 */


/**
 * @brief print_help prints instructions on screen
 */
void print_help(){
    printf("\n NAME: \n");
    printf("\t Superfacets Segmentation \n\n");
    printf("\t Author: Giulia Picciau");
    printf("\t Date: 04-15-2014");
    printf("\nUSAGE:\n");


    printf("  -m [input mesh]\n");
    printf("\t the mesh to be read as input\n\n");

    printf("  -r [radius]\n");
    printf("\t approximated radius of a single cluster (mutually exclusive with -nseg)\n\n");

    printf("  -nseg [segments]\n");
    printf("\t number of desired regions (mutually exclusive with option -r\n\n");

    printf("  -a [alpha]\n");
    printf("\t parameter that states how much the angular distance will be important in the total one\n\n");

    printf("  -e [etaconvex]\n");
    printf("\t value of the weight for convex angles\n\n");

    printf("  -vis\n");
    printf("\t if you want to visualize the output subdivision after the segmentation process\n\n");

    printf("  -out [fileout]\n");
    printf("\t file on which is written the output segmentation\n\n");

    printf("  -ov [segmfile]\n");
    printf("\t only calls the visualizer on an already computed segmentation (stored in segmfile)\n\n");

    printf("  -h:\n");
    printf("\t to put the header in the output segmentation file\n\n");

    printf("  -flood [initmode] :\n");
    printf("\t specify the initialization mode we want to use (flooding (initmode=1) or squares (initmode=0))\n\n");

    printf("  -mult [factor]: \n");
    printf("\t specify the number which multiplies the distance threshold in the expansion step (default=2)\n\n");

    printf("  -nIter [iters]\n");
    printf("\t how many iterative steps we want to perform at most (default=50)\n\n");

    printf("  -debug\n");
    printf("\t if this flag is specified, a lot of debug messages will be shown\n\n\n");
}

/**
 * @brief callVis opens the visualizer window
 * @param MN name of the mesh file
 * @param SN name of the segmentation file
 */
void callVis(string MN, string SN){

    MeshVisualizer *MV = new MeshVisualizer(MN, SN);
    MV->resize(800, 600);
    MV->show();
}

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    bool debugM=false;
    int timesR = 2;
    int nIters = 50;
    bool triBased = true;

    if(argc==1){
        print_help();
        exit(0);
    }


    //Parse arguments
    string meshfilename;
    string segFile;
    string fieldFile;
    double alpha;
    double radius;
    int number;
    double etaConvex=0.2; //Default value
    bool putHeader=false;
    bool visualize=false;
    bool justvisualize=false;
    bool floodI = true;

    for(int n_opt=1;n_opt<argc;n_opt+=2){

        char* option=argv[n_opt];

        if(!strcmp(option, "-m")){ //Mesh name
            QString input = argv[n_opt+1];
            meshfilename = input.toStdString();
        }

        else if(!strcmp(option, "-r")){  //Radius (implies radius-based init)
            QString input = argv[n_opt+1];
            radius = input.toDouble();
            number = -1;
        }

        else if(!strcmp(option, "-nseg")){  //Number of regions
            QString input = argv[n_opt+1];
            number = input.toInt();
            radius = -1.0;
        }

        else if(!strcmp(option, "-a")){ //Alpha
            QString input = argv[n_opt+1];
            alpha = input.toDouble();
        }

        else if(!strcmp(option, "-eta")){ //Value of eta if the angle is convex
            QString input = argv[n_opt+1];
            etaConvex = input.toDouble();
        }

        else if(!strcmp(option, "-h")){  //put the header in the file
            putHeader = true;
            n_opt--;
        }

        else if(!strcmp(option, "-vis")){  //Launch the visualizer when segmentation is done
            visualize = true;
            n_opt--;
        }

        else if(!strcmp(option, "-out")){  //The output file in which we write the segmentation
            QString input = argv[n_opt+1];
            segFile = input.toStdString();
        }

        else if(!strcmp(option, "-flood")){  //Radius-based, decide the strategy
            QString input = argv[n_opt+1];
            floodI = (bool)input.toInt();
        }

        else if(!strcmp(option, "-mult")){
            QString input = argv[n_opt+1];
            timesR = input.toInt();
        }
        else if(!strcmp(option, "-nIter")){
            QString input = argv[n_opt+1];
            nIters = input.toInt();
        }

        else if(!strcmp(option, "-ov")){  //Calls only the visualizer
            justvisualize = true;
            QString input = argv[n_opt+1];
            segFile = input.toStdString();
        }

        else if(!strcmp(option, "-debug")){
            debugM=true;
            cout<<"Debug mode selected"<<endl;
        }

        else if(!strcmp(option, "-VB")){
            triBased = false;
            justvisualize=false;
            QString input = argv[n_opt+1];
            fieldFile = input.toStdString();
        }
        else{
            cout<<"Uncorrect usage"<<endl;
            print_help();
            exit(0);
        }
    }

    if(justvisualize){  /// if we only want to call the visualizer on a segmentation

        callVis(meshfilename, segFile);
        return a.exec();
    }

    if(!triBased){

        cout<<"Loading vertices"<<endl;
        VertexBasedSegmenter *VBSuperSeg = new VertexBasedSegmenter;
        cout<<"Created"<<endl;
        VBSuperSeg->filename = meshfilename;
        VBSuperSeg->fieldfilename = fieldFile;
        VBSuperSeg->setNCluster(number);
        VBSuperSeg->setAlpha(alpha);
        VBSuperSeg->callLoad();

        VBSuperSeg->startSeg();
        cout<<"Converged, now writing"<<endl;
        VBSuperSeg->callWriter(segFile);

        return 0;
    }

    /// Create segmenter class
    Segmenter *SuperSeg = new Segmenter;

    /// set parameters
    SuperSeg->filename=meshfilename;  /// file to segment
    SuperSeg->setAlpha(alpha);   /// weight of the angular distance
    /// Spatial threshold / number of regions
    if(radius > 0)
        SuperSeg->setMaxD(radius);
    if(number > 0)
        SuperSeg->setNCluster(number);
    SuperSeg->setEtaConvex(etaConvex); /// weight of convex angle
    SuperSeg->putHeader=putHeader;
    SuperSeg->setTimesR(timesR);  /// factor which multiplies the threshold in the expansion step
    SuperSeg->setMaxIters(nIters); /// Number of maximum iterative steps to make
    SuperSeg->floodInit = floodI;  /// If we know the radius, to decide between flood and grid initialization
    SuperSeg->debugMode = debugM;  /// If we are running in debug mode (more messages will be displayed)

    /// Load the model
    SuperSeg->callLoad();
    bool haveMoved=true;
    int countIter=0;
    Timer TMR;
    TMR.start();

    /// Iterate until convergence
    while(haveMoved && countIter<SuperSeg->getMaxIters()){
        cout<<"Iteration "<<countIter<<"..."<<endl;
        if(haveMoved)
            SuperSeg->expansionStep();
        if(debugM)
            cout<<"Expanded"<<endl;
        haveMoved = SuperSeg->updateCenters();
        if(debugM)
            cout<<"Updated"<<endl;
        countIter++;
    }
    TMR.stop();
    double timeForSeg = TMR.getElapsedTimeInMilliSec();
    timeForSeg += SuperSeg->initTime;
    SuperSeg->setElapsedTime(timeForSeg);
    SuperSeg->setIters(countIter);

    /// Write to file
    if(strcmp(segFile.c_str(), ""))
        SuperSeg->writeSegmOnFile(segFile);
    else
        cout<<"No output file"<<endl;

    /// calls the visualizer after the segmentation
    if(visualize){

        callVis(meshfilename, segFile);
        return a.exec();
    }

    return 0;

}
