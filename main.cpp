#include "segmenter.h"
#include <stdio.h>
#include <QApplication>

#define MAX_ITERATIONS 10

/**
 * @author Giulia Picciau
 * Submission to Eurographics Symposium on Geometric Processing 2014
 */


/**
 * @brief print_help prints instructions on screen
 */
void print_help()
{
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
void callVis(string MN, string SN)
{

    // MeshVisualizer *MV = new MeshVisualizer(MN, SN);
    // MV->resize(800, 600);
    // MV->show();
}

Mesh<Triangle> OpenMeshFile(const QString& pMeshName)
{
    Mesh<Triangle> mesh;
    QStringList parsed = pMeshName.split(".");

    if(!parsed.back().compare("tri"))
    {
        Reader::readMeshFile(mesh, pMeshName.toStdString());
    }
    else if(!parsed.back().compare("off"))
    {
        Reader::readOFFMesh(mesh, pMeshName.toStdString());
    }
    else
    {
        std::cout << "Not a valid file format (it must be .off or .tri), instead it is " << parsed.back().toStdString() << std::endl;
        exit(0);
    }

    mesh.Build();

    return mesh;
}


void WriteOutputSegmentation(const std::vector<int>& pSegm, std::string pFileOut, bool pPutHeader = true)
{
    //open file in write mode
    ofstream ofs;
    ofs.open(pFileOut.c_str());

    if(ofs.is_open())
    {

        //if we want an header (visualization purposes)
        if(pPutHeader)
        {
            ofs << "#";
            // ofs << asctime(localT) << std::endl;
            // ofs << "#mesh input=" << m_MeshName << std::endl;
            // ofs << "#radius=" << m_RegionRadius << std::endl;
            // ofs << "#alpha=" << m_Alpha << std::endl;
        }
        for(size_t a=0; a<pSegm.size(); a++)
        {
            ofs << pSegm[a] << endl;
        }

        ofs.close();
    }
    else
    {
        std::cout << "Unable to write on file" << std::endl;
    }
}


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    bool debugM=false;
    int timesR = 2;
    int nIters = 50;
    bool triBased = true;

    if(argc==1)
    {
        print_help();
        exit(0);
    }


    //Parse arguments
    std::string meshfilename;
    std::string segFile = "Output.dat";
    double alpha;
    double radius;
    int number;
    bool floodI = true;

    for(int n_opt=1;n_opt<argc;n_opt+=2)
    {

        char* option=argv[n_opt];

        if(!strcmp(option, "-m"))
        { //Mesh name
            QString input = argv[n_opt+1];
            meshfilename = input.toStdString();
        }

        else if(!strcmp(option, "-r"))
        {  //Radius (implies radius-based init)
            QString input = argv[n_opt+1];
            radius = input.toDouble();
            number = -1;
        }

        else if(!strcmp(option, "-nseg"))
        {  //Number of regions
            QString input = argv[n_opt+1];
            number = input.toInt();
            radius = -1.0;
        }

        else if(!strcmp(option, "-a"))
        { //Alpha
            QString input = argv[n_opt+1];
            alpha = input.toDouble();
        }

        else if(!strcmp(option, "-eta"))
        { //Value of eta if the angle is convex
            QString input = argv[n_opt+1];
            //etaConvex = input.toDouble();
        }

        else if(!strcmp(option, "-h"))
        {  //put the header in the file
            //putHeader = true;
            n_opt--;
        }

        else if(!strcmp(option, "-vis"))
        {  //Launch the visualizer when segmentation is done
            //visualize = true;
            n_opt--;
        }

        else if(!strcmp(option, "-out"))
        {  //The output file in which we write the segmentation
            QString input = argv[n_opt+1];
            segFile = input.toStdString();
        }

        else if(!strcmp(option, "-flood"))
        {  //Radius-based, decide the strategy
            QString input = argv[n_opt+1];
            floodI = (bool)input.toInt();
        }

        else if(!strcmp(option, "-mult"))
        {
            QString input = argv[n_opt+1];
            timesR = input.toInt();
        }
        else if(!strcmp(option, "-nIter"))
        {
            QString input = argv[n_opt+1];
            nIters = input.toInt();
        }

        else if(!strcmp(option, "-debug"))
        {
            debugM=true;
            cout<<"Debug mode selected"<<endl;
        }
        else
        {
            cout<<"Uncorrect usage"<<endl;
            print_help();
            exit(0);
        }
    }


    auto procMesh = OpenMeshFile(QString::fromStdString(meshfilename));

    Segmenter segmenter = (number <= 0) ? Segmenter(procMesh, radius, alpha, floodI) : Segmenter(procMesh, number, alpha);
    segmenter.StartSegmentation();

    for(unsigned short idx = 0; idx < MAX_ITERATIONS; ++idx)
    {
        std::cout << "STEP " << idx << std::endl;
        segmenter.ClassificationStep();
        if(!segmenter.UpdateCenters())
        {
            std::cout << "Converged after " << idx << " Steps" << std::endl;
            break;
        }
    }

    WriteOutputSegmentation(segmenter.GetSegmentation(), segFile);


    return 0;

}
