/*
B. Naegel
5/11/2014

Input images :
-double nrrd original image
-int nrrd label image
*/

#include <sstream>
#include <cmath>
#include <algorithm>
#include "include/Image.h"
#include "include/ComponentTree.h"
#include "include/FlatSE.h"
#include "nifti/nifti1_io.h"


using namespace LibTIM;

const int magicConstant=1e8;
int classe;
int imageId;


// Main algorithm
// Input:
// - tree: initialized and attributed component-tree
// - alpha: parameter (floating number strictly comprised between 0 and 1)
// Output:
// -selectedNodes: list of nodes selected by the algorithm
vector <Node *> computeSolution(ComponentTree<U8> &tree, double alpha)
{
    vector <Node *> orderedNodes;
    vector <Node *> selectedNodes;

    // Breadth-first scan of the tree
    // The nodes are stored in the vector orderedNodes
    std::queue<Node *> fifo;
    fifo.push(tree.m_root);

    while(!fifo.empty())
    {
        Node *tmp=fifo.front();
        fifo.pop();

        orderedNodes.push_back(tmp);

        for(unsigned int i=0; i<tmp->childs.size(); i++)
        {
            fifo.push(tmp->childs[i]);
        }
    }

    double exprl,exprr;

    // Scan all the nodes from the leafs in reverse order
    // It ensures that all the nodes are processed before their father
    for(int i=orderedNodes.size()-1; i>=0; i--)
    {
        // Take the following node in the list
        Node *tmp=orderedNodes[i];

        // Compute left and right expressions (costs to keep or skip the node)
        exprl=alpha*tmp->n;
        exprr=(1-alpha)*tmp->ps;

        // if tmp is not the root
        if(tmp->father!=tmp)
        {
            // if tmp is a leaf
            if(tmp->childs.size()==0)
            {
                // The node is kept
                if(exprl<exprr)
                {
                    tmp->calpha=exprl;
                    selectedNodes.push_back(tmp);
                }
                // The node is skipped
                else
                {
                    tmp->calpha=exprr;
                }
            }
            // if tmp is not a leaf
            else
            {
                // Compute the sum of costs of all child nodes
                double sum=0.0;
                for(unsigned int i=0;i<tmp->childs.size();i++)
                {
                    sum+=tmp->childs[i]->calpha;
                }

                // Add this sum to exprr (cost to skip the node)
                exprr+=sum;

                if(exprl<exprr)
                {
                    // The node is kept
                    tmp->calpha=exprl;
                    selectedNodes.push_back(tmp);
                }
                else
                {
                    // The node is skipped
                    tmp->calpha=exprr;
                }
            }
        }
    }

    return selectedNodes;
}




int getMaxGrey(Node *n, int h)
{
    if(n->h > h) {
        h=n->h;
    }
    for(int i=0; i<n->childs.size(); i++) {
        h=getMaxGrey(n->childs[i],h);
    }
    return h;
}


void sampleTree(Node *node, vector<Node*> &nodes, int nmax)
{
    if((rand()%100)<10) {
            if(nodes.size()<nmax) {
                nodes.push_back(node);
            }
            else return;
        }

    for(int i=0; i<node->childs.size(); i++)
        sampleTree(node->childs[i],nodes,nmax);
}

template <typename T>
Image<T> readNifti(char *filename)
{

    nifti_image *nim=NULL;
    nim=nifti_image_read(filename,1);

    Image<T> inputImage(nim->nx,nim->ny,nim->nz);


    for(int i=0; i<inputImage.getBufSize(); i++) {
        inputImage(i)=((T *)nim->data)[i];
    }

    return inputImage;
}


void writeNifti(Image<unsigned char> &im, const char *filename)
{

    nifti_image *nim=NULL;
    int t[8]={3,im.getSizeX(), im.getSizeY(), im.getSizeZ(), 0,0,0,0};
    nim=nifti_make_new_nim(t,2,1);


    for(int i=0; i<im.getBufSize(); i++) {
        ((unsigned char *)(nim->data))[i]=(unsigned char)im(i);
    }

    int test=nifti_set_filenames(nim,filename,0,nim->byteorder);
    nifti_image_write(nim);
}


int main(int argc, char **argv)
{

    if (argc!=6) {
        fprintf(stderr, "usage: %s <input image (double)> <label image (int)> <output path> <alpha> <threshold (int)> \n",argv[0]);
        exit(1);
    }

    Image<unsigned char> inputImage=readNifti<unsigned char>(argv[1]);
    Image<unsigned char> labelImage=readNifti<unsigned char>(argv[2]);


    // Threshold the input image
    Image<unsigned char>::iterator it;
    Image<unsigned char>::iterator end=inputImage.end();
    int threshold = atoi(argv[5]);
    for(it=inputImage.begin(); it!=end; ++it){
        if (*it > threshold){
            *it = 255;
        }
        else{
            *it = 0;

        }
    }


    // alpha parameter (alpha is a floating number which must be comprised between 0 and 1)
    double alpha=atof(argv[4]);

    if(alpha<0 || alpha>1)
    {
        std::cout<<"Error: alpha must be comprised between 0 and 1\n";
        exit(1);
    }

    FlatSE connexity;
    //26-connexity
    connexity.make3DN26();

    ComponentTree<unsigned char> tree(inputImage, labelImage, connexity);


    // Compute the selected nodes
    vector<Node *> selectedNodes;
    selectedNodes=computeSolution(tree,alpha);

    std::cerr << "Number of nodes: " << selectedNodes.size() <<"\n";

    // Computation of the result image from the set of selected nodes
    Image<unsigned char> imRes;
    imRes.setSize(inputImage.getSize());
    imRes.fill(0);
    for(unsigned int i=0; i<selectedNodes.size(); i++)
    {
        // Draw the node in imRes
        // i.e. set all pixels belonging to the node to the grey-level of node
        tree.constructNode(imRes, selectedNodes[i]);
    }


    // Threshold the result to obtain a binary segmentation
    Image<unsigned char>::iterator it2;
    Image<unsigned char>::iterator end2 = imRes.end();

    for(it2=imRes.begin(); it2!=end2; ++it2){
        if (*it2 > 10){
            *it2 = 255;
        }

        else{
            *it2 = 0;

        }
    }

    writeNifti(imRes, argv[3]);
    return 0;
}

