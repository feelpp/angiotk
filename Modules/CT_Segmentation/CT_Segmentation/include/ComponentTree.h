#ifndef ComponentTree_h
#define ComponentTree_h

#include "Morphology.h"

namespace LibTIM {

using std::vector;

/** \defgroup ccTree Component-Tree Based Algorithms
	\ingroup Morpho
**/

/*@{*/

/** @brief Component tree representation of an image
 *	The component tree (or max-tree) is a particular image structure. Nodes represents flat zones.
 *	Leafs represents regional maxima (max-tree) or minima (min-tree). Father-son relation represents
 *	an inclusion relation (with respect to level sets representation of image).
 *	This structure is efficient for computing attribute openings (area, contrast, ...)
**/
const int localMax=std::numeric_limits<int>::max();
const int localMin=std::numeric_limits<int>::min();

//TODO: a generic structure that contains Node attributes
struct Node {
		Node()
		: label(-1),xmin(localMax),xmax(localMin),ymin(localMax),
		ymax(localMin),area(0),
		contrast(0), volume(0),  contourLength(0),
		complexity(0), subNodes(0),status(true),
		rSignature(0),m01(0),m10(0),m20(0),m02(0),n(0), ps(0), calpha(0),
		active(true), father(0)
			{
			pixels.reserve(7); childs.reserve(5);
			}
        //Attributes
		int label;
		// BN 2/5/09
		int uniqueId;
		int ori_h;
		int h;
		int xmin;
		int xmax;
		int ymin;
		int ymax;
		int area;
		int contrast;
		int volume;
		int contourLength;
		int complexity;
		double compacity;

		// BN 31/12/08
		// add attributes related to the grey-level mean of the component
		double mean;
		// mean of dilated region Id
		double meanId;
		double var;

		double meanNghb;
		double varNghb;
		double varId;
		int maxNghb;
		int minNghb;

		int subNodes;
		bool status;

        // BN 02/05/09
        std::vector<int> nghbNdg;

		//BN 19/02/08: shape descriptor based on R-signature (or F2d signature with r=0)
		double *rSignature;



		//experimental... moments
		long double m01;
		long double m10;
		long double m20;
		long double m02;
		//moment d'inertie (first Hu's invariant moment)
		double I;

		// BN 8/09/2009
		// attributes for nodes selection (see paper)
		//int p;
		int n;
		int ps;
		//int ns;
		// BN 9/09/2009
		double calpha;


		bool active;
		//Distance measure between the node and a reference vector
		//(and eventually a covariance matrix for Mahalanobis distance)
		double dist;

		//Common to all type of nodes:
		Node *father;
		typedef std::vector<TOffset> ContainerPixels;
		ContainerPixels pixels;
                ContainerPixels pixels_tmp;
		typedef  std::vector<Node *> ContainerChilds;
		ContainerChilds childs;
		};


typedef std::vector<std::vector<Node *> > IndexType;


template <class T>
class ComponentTreeStrategy;

template <class T>
class SalembierRecursiveImplementation;

template <class T>
class ComponentTree {
	public:
		ComponentTree() {};
		ComponentTree(Image <T> &img);
		ComponentTree(Image <T> &img, FlatSE &connexity);
		// BN 08/09/09
		// constructor based on binary ground-truth used for nodes selection (see paper)
		ComponentTree(Image <T> &img, Image <U8> &gt);
		ComponentTree(Image <T> &img, Image <U8> &gt,FlatSE &connexity);
		//Copy constructor
		ComponentTree(ComponentTree <T> &tree);
		~ComponentTree();

		enum ConstructionDecision {MIN,MAX,DIRECT};
		Image <T> constructImage(ConstructionDecision decision=MIN);
		Image <T> &constructImageOptimized();
		/**
		  * @brief Print tree on standard output
		**/

		void print();

		/**
		  * @brief Print size of tree (in kilobytes)
		**/
		void printSize();

		/**
		  * @brief Print tree in a .dot file to generate clean drawing
		**/

		int writeDot(const char *filename);

		/**
		  * @brief Restore original tree (i.e. clear all filtering)
		**/

		int restore();

		/**
		  * @brief Area filtering
		**/

		int areaFiltering(int tMin, int tMax=std::numeric_limits<int>::max());

		// BN 04/01/09
		int meanFiltering(int tMin, int tMax=std::numeric_limits<int>::max());

		/**
		  * @brief Contrast filtering
		**/

		//Keep only cc having certain contrast
		int contrastFiltering(int tMin, int tMax=std::numeric_limits<int>::max());


		/**
		  * @brief Volumic filtering
		**/

		int volumicFiltering(int tMin, int tMax=std::numeric_limits<int>::max());

		/**
		  * @brief Connected opening  (experimental))
		**/

		int connectedOpening(FlatSE &se);
		int connectedOpening2(FlatSE &se);

		/**
		  * @brief Shape filtering
		**/


		int complexityFiltering(int tMin, int tMax);
		int compacityFiltering(int tMin, int tMax);
		int intensityFiltering(int tMin, int tMax);

		int boundingBoxFiltering(int min, int max);

		//Experimental 16/07/07
		//N most brightest leafs
		vector<Node *> intensityPruning(int N);


		void setFalse();

		/**
		  * @brief Return the signature of a branch at point (x,y)
		  * (i.e. the value of all attributes along the path from root to the node containing (x,y))
		**/

		typedef std::map<int, Node *> SignatureType;
		SignatureType branchSignature(TCoord x, TCoord y);
		SignatureType branchSignature(Node *node);
		void activeBranch(Node *leaf);
		void activeBranch(SignatureType &signature);
		int writeSignature(SignatureType &signature, const char *file);

	//private:

		void erase_tree();

		//Helper functions for filtering
		std::vector <TOffset > merge_pixels(Node *tree);
		std::vector <TOffset > merge_pixelsFalseNodes(Node *tree);
		void merge_pixels(Node *tree, std::vector <TOffset > &res);

		bool isInclude(FlatSE &se, Node::ContainerPixels &pixels);

		Node * coordToNode(TCoord x, TCoord y);
		Node * offsetToNode(TOffset offset);

        Image <U8> coordToComponent(TCoord x, TCoord y);
		Image <U8> offsetToComponent(TOffset offset);

		void constructImageMin(Image <T> &res);
		void constructImageMax(Image<T> &res);
		void constructImageDirect(Image<T> &res);
		void constructImageDirectExpe(Image<T> &res);

		void constructNode(Image <T> &res, Node *node);
		void constructNodeDirect(Image <T> &res, Node *node);

		//Experimental 16/07/07
		void constructBranch(Image <T> &res,Node *leaf);

		// Internal structure
		// root node

		Node *m_root;
		//TSize *m_size;

		//original data
		Image<T> m_img;

		//STATUS image
		Image <int> STATUS;

		//max-tree index
		IndexType index;

        Image<Node *> indexNodes;

		int totalNodes;

};

/** @brief Abstract class for strategy to compute component tree
  *	Abstract class encapsulating the various strategies to compute Max-tree
  *
  *
**/
template <class T>
class ComponentTreeStrategy {
	public:

	ComponentTreeStrategy() {};
	virtual ~ComponentTreeStrategy() {};

	virtual Node *computeTree()=0;
	virtual int computeAttributes(Node *tree)=0;



};

/** @brief Salembier recursive implementation
 **/

template <class T>
class SalembierRecursiveImplementation: public ComponentTreeStrategy <T> {
	public:
	SalembierRecursiveImplementation(ComponentTree <T> *parent, FlatSE &connexity)
	:m_parent(parent)
		{
		this->totalNodes=0;
		this->init(m_parent->m_img, connexity);
		}
    //BN 08/09/09
    SalembierRecursiveImplementation(ComponentTree <T> *parent, Image<U8> &gt, FlatSE &connexity)
	:m_parent(parent)
		{
		this->totalNodes=0;
		this->init(m_parent->m_img, gt, connexity);
		}

	~SalembierRecursiveImplementation()
		{
		delete[] hq;
		}

	Node *computeTree();
	Node *computeTree(Image<U8> &gt);

	int computeAttributes(Node *tree);

	int computeArea(Node *tree);
	int computeContrast(Node *tree);
	int computeVolume(Node *tree);

	//BN 08/09/09

	int computeN(Node *tree);

	// BN 31/12/08
	double computeMean(Node *tree);

	int computeSubNodes(Node *tree);

	long double computeM01(Node *tree);
	long double computeM10(Node *tree);
	long double computeM02(Node *tree);
	long double computeM20(Node *tree);
	void computeInertiaMoment(Node *tree);

	//Shape-based attributes
	int computeContourLength();
	int computeComplexityAndCompacity(Node *tree);

	int computeBoundingBox(Node *tree);


	private:
		//Helper functions
		inline int update_attributes(Node *n, TOffset &imBorderOffset);
		inline int update_attributes(Node *n, TOffset &imBorderOffset, Image<U8> &gt);
		inline int flood(int m) ;
		inline int flood(int m, Image<U8> &gt) ;
		void link_node(Node *tree, Node *child) ;
		Node *new_node(int h, int n)  ;

		int init(Image <T> &img, FlatSE &connexity) ;
		//BN 08/09/09
		int init(Image <T> &img, Image<U8> &gt, FlatSE &connexity) ;

		int hToIndex(int h)  {return h-hMin;}
		int indexToH(int h)  {return h+hMin;}

		//members
		Image <T> imBorder;
		FlatSE se;
		TSize oriSize[3];

		static const T BORDER=T(0);
		static const int NOT_ACTIVE=-1;
		static const int ACTIVE=-2;
		static const int BORDER_STATUS=-3;
		TCoord back[3];
		TCoord front[3];
		/** @brief Hierarchical queue
		**/
		//typedef std::map <int, std::queue<TOffset> > HierarchicalQueue;
		typedef std::queue<TOffset> * HierarchicalQueue;

		 HierarchicalQueue hq;

		 T hMin;
		 T hMax;
		// size of different containers
		 int numberOfLevels;

		 int totalNodes;

		 Image <int> STATUS;
		 vector <int> number_nodes;
		 vector <bool> node_at_level;
		// For now, container for accessing nodes by level and cc number
		//typedef std::map <T, std::map<TLabel,  Node *> > IndexType;
		//typedef Node *** IndexType;

		IndexType index;
        Image<Node *> indexNodes;

		ComponentTree<T> *m_parent;
};

/** @brief Union-Find implementation (Najman-Couprie)
 **/


template <class T>
class UnionFindImplementation: public ComponentTreeStrategy <T> {




};




/*@}*/

}//end namespace

#include "ComponentTree.hxx"
#endif
